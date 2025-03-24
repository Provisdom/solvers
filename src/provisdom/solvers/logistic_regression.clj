(ns provisdom.solvers.logistic-regression
  (:require
    [clojure.spec.alpha :as s]
    [clojure.spec.gen.alpha :as gen]
    [provisdom.apache-math.apache-matrix :as apache-mx]
    [provisdom.math.core :as m]
    [provisdom.math.matrix :as mx]
    [provisdom.math.special-functions :as special-fns]
    [provisdom.math.tensor :as tensor]
    [provisdom.math.vector :as vector]
    [provisdom.utility-belt.anomalies :as anomalies]))

(s/def ::max-iter
  (s/with-gen (s/nilable ::m/int+)
              #(gen/one-of [(s/gen (s/int-in 100 1000))
                            (gen/return nil)])))

(s/def ::abs-accu
  (s/with-gen ::m/finite+
              #(s/gen (s/double-in :min m/tiny-dbl
                                   :max 1.0
                                   :NaN? false))))

(s/def ::ridge-lambda
  (s/with-gen (s/nilable ::m/open-prob)
              #(gen/one-of [(s/gen (s/double-in 1e-14 1e-3))
                            (gen/return nil)])))

(defn iterative-reweighted-least-squares
  "Can optionally set `ridge-lambda` parameter to do Ridge regression."
  ([x-mx y] (iterative-reweighted-least-squares x-mx y {}))
  ([x-mx y {::keys [max-iter abs-accu ridge-lambda] :or {abs-accu 1e-15}}]
   (let [y-sum (vector/kahan-sum y)
         data-count (count y)
         dependent-var-count (mx/columns x-mx)
         y-ave (/ y-sum data-count)
         weight (m/log (/ y-ave (m/one- y-ave)))
         xt (mx/transpose x-mx)
         max-iter (or max-iter 1000)]
     (let [w* (loop [w (vec (repeat dependent-var-count 0.0))
                     i 0]
                (if (>= i max-iter)
                  {::anomalies/category ::anomalies/no-solve
                   ::anomalies/fn       (var iterative-reweighted-least-squares)
                   ::anomalies/message  (str "Solution didn't converge within"
                                             " `abs-accu` in `max-iter`.")}
                  (let [sz (map (fn [x yi]
                                  (let [ni (+ weight (vector/dot-product w x))]
                                    (if-not (m/finite? ni)
                                      [m/nan m/nan]
                                      (let [mi (special-fns/logistic ni)
                                            si (* mi (m/one- mi))
                                            diff (- yi mi)]
                                        [si
                                         (+ ni
                                            (if (zero? diff)
                                              0.0
                                              (m/div diff
                                                     si
                                                     (* (m/sgn diff)
                                                        m/inf+))))]))))
                                x-mx
                                y)
                        s (map first sz)]
                    (if (some m/nan? s)
                      {::anomalies/category ::anomalies/no-solve
                       ::anomalies/fn       (var iterative-reweighted-least-squares)
                       ::anomalies/message  "No Solution: NaN."}
                      (let [z (mapv second sz)
                            s-mx (mx/diagonal-matrix s)
                            xsx (if ridge-lambda
                                  (let [lambda-diagonal (mx/diagonal-matrix
                                                          (vec
                                                            (repeat dependent-var-count
                                                                    ridge-lambda)))]
                                        (apache-mx/->apache-matrix
                                          (tensor/add (mx/mx* xt s-mx x-mx) lambda-diagonal)))
                                  (apache-mx/->apache-matrix (mx/mx* xt s-mx x-mx)))
                            xsx-1 (apache-mx/inverse xsx)]
                        (if-not xsx-1
                          {::anomalies/category ::anomalies/no-solve
                           ::anomalies/fn       (var iterative-reweighted-least-squares)
                           ::anomalies/message  "No Solution: no inverse."}
                          (let [xsx-1-mx (apache-mx/apache-matrix->matrix xsx-1)
                                w' (vector/to-vector
                                     (mx/mx* xsx-1-mx xt s-mx (mx/column-matrix z)))
                                w' (if ridge-lambda
                                     (tensor/add w'
                                                 (tensor/multiply (* 2.0 ridge-lambda)
                                                                  w))
                                     w')]
                            (if (tensor/roughly? w' w abs-accu)
                              w'
                              (recur w' (inc i))))))))))]
       w*))))

(s/def ::parameters
  (s/and ::vector/vector-finite
         (fn [params] (pos? (count params)))))

(s/def ::y
  (s/and (vector/vector-of-spec {:pred ::m/prob})
         (fn [y] (pos? (count y)))))

(s/def ::x-mx
  (s/and ::mx/matrix-finite
         (fn [x-mx]
           (pos? (mx/columns x-mx)))))

(s/fdef iterative-reweighted-least-squares
        :args (s/and (s/cat :x-mx ::x-mx
                            :y ::y
                            :opts (s/? (s/keys :opt [::max-iter
                                                     ::abs-accu
                                                     ::ridge-lambda])))
                     (fn [{:keys [x-mx y]}]
                       (and (= (mx/rows x-mx) (count y))
                            (>= (count y) 2)
                            (not (every? m/one? y))
                            (not (every? zero? y)))))
        :ret (s/or :parameters ::parameters
                   :anomaly ::anomalies/anomaly))

(defn log-likelihood-of-data-given-parameters
  ""
  [{::keys [x-mx y parameters]}]
  (reduce
    (fn [tot [yi x]]
      (let [xb (vector/dot-product x parameters)]
        (+ tot
           (- (* yi xb)
              (m/log (inc (m/exp xb)))))))
    0.0
    (partition 2 (interleave y x-mx))))

(s/fdef log-likelihood-of-data-given-parameters
        :args (s/with-gen
                (s/and (s/cat :args (s/keys :req [::x-mx ::y ::parameters]))
                       (fn [{:keys [args]}]
                         (let [{::keys [x-mx y parameters]} args]
                           (and (= (mx/rows x-mx) (count y))
                                (= (mx/columns x-mx) (count parameters))))))
                #(gen/return [{::x-mx       [[0 30] [1 2] [3 5] [7 23] [10 10]]
                               ::y          [0.0 0.56 0.8 0.88 1.0]
                               ::parameters [0.7607052740510958 -0.15751338511807894]}]))
        :ret ::m/number)
