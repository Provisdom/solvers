(ns provisdom.solvers.generalized-regression-neural-network
  (:require
    [clojure.spec.alpha :as s]
    [clojure.spec.gen.alpha :as gen]
    [provisdom.math.core :as m]
    [provisdom.math.intervals :as intervals]
    [provisdom.math.matrix :as mx]
    [provisdom.math.tensor :as tensor]
    [provisdom.math.vector :as vector]
    [provisdom.solvers.optimize-univariate :as opt-uni]
    [provisdom.utility-belt.anomalies :as anomalies]))

(s/def ::max-iter
  (s/with-gen (s/nilable ::m/int+)
    #(gen/one-of [(s/gen (s/int-in 10 100))
                  (gen/return nil)])))

(s/def ::abs-accu
  (s/with-gen ::m/finite+
    #(s/gen (s/double-in :min m/tiny-dbl
              :max 1.0
              :NaN? false))))

(s/def ::rel-accu
  (s/with-gen ::m/finite+
    #(s/gen (s/double-in :min m/tiny-dbl
              :max 1.0
              :NaN? false))))

(defn- negative-half-distances
  ""
  [x-mx x-new]
  (mapv (fn [x]
          (let [x-adj (tensor/subtract x x-new)]
            (* -0.5 (vector/dot-product x-adj x-adj))))
    x-mx))

(defn- solve*
  ""
  [neg-half-distances spread y]
  (let [weight (m/pow spread -2.0)
        kernels (mapv (fn [d]
                        (m/exp (* d weight)))
                  neg-half-distances)
        top (vector/dot-product kernels y)
        divisor (vector/kahan-sum kernels)]
    (m/div top divisor)))

(defn solve-for-spread
  ""
  ([args] (solve-for-spread args {}))
  ([{::keys [x-training-mx
             y-training
             x-test-mx
             y-test
             spread-guess
             spread-interval]}
    {::keys [max-iter abs-accu rel-accu] :or {abs-accu 1e-6 rel-accu 1e-14}}]
   (let [max-iter (or max-iter 100)
         negative-half-distances-mx (mapv
                                      (fn [x-test]
                                        (negative-half-distances
                                          x-training-mx x-test))
                                      x-test-mx)
         spread-f (fn [spread]
                    (reduce +
                      (map
                        (fn [neg-half-distances y*]
                          (m/sq (- y*
                                  (solve*
                                    neg-half-distances spread y-training))))
                        negative-half-distances-mx
                        y-test)))
         sol (opt-uni/optimize-univariate
               {::opt-uni/univariate-f           spread-f
                ::opt-uni/strict-finite-interval spread-interval
                ::opt-uni/guess                  spread-guess}
               {::opt-uni/max-iter max-iter
                ::opt-uni/rel-accu rel-accu
                ::opt-uni/abs-accu abs-accu})]
     (if (anomalies/anomaly? sol)
       sol
       (::opt-uni/point sol)))))

(s/def ::x-training-mx
  (s/and ::mx/matrix-finite
    (fn [mx]
      (pos? (mx/columns mx)))))

(s/def ::y-training ::vector/vector-finite)

(s/def ::x-test-mx
  (s/and ::mx/matrix-finite
    (fn [mx]
      (pos? (mx/columns mx)))))

(s/def ::y-test ::vector/vector-finite)
(s/def ::spread-guess ::m/finite+)
(s/def ::spread-interval ::intervals/finite+-interval)

(s/fdef solve-for-spread
  :args (s/with-gen
          (s/and (s/cat :args (s/keys :req [::x-training-mx
                                            ::y-training
                                            ::x-test-mx
                                            ::y-test
                                            ::spread-guess
                                            ::spread-interval])
                   :opts (s/? (s/keys :opt [::max-iter ::abs-accu ::rel-accu])))
            (fn [{:keys [args]}]
              (let [{::keys [x-training-mx y-training x-test-mx
                             y-test spread-guess spread-interval]} args]
                (and (= (mx/rows x-training-mx) (count y-training))
                  (= (mx/rows x-test-mx) (count y-test))
                  (= (mx/columns x-training-mx) (mx/columns x-test-mx))
                  (not (== (first spread-interval) (second spread-interval)))
                  (intervals/in-bounds?
                    (intervals/bounds spread-interval true false)
                    spread-guess)))))
          #(gen/bind (gen/tuple (s/gen ::x-training-mx)
                       (s/gen ::y-test)
                       (s/gen ::spread-interval))
             (fn [[x-training-mx y-test spread-interval]]
               (gen/fmap (fn [[y-training x-test-mx spread-guess]]
                           [{::x-training-mx   x-training-mx
                             ::y-training      y-training
                             ::x-test-mx       x-test-mx
                             ::y-test          y-test
                             ::spread-guess    spread-guess
                             ::spread-interval spread-interval}])
                 (gen/tuple
                   (gen/vector (s/gen ::m/finite) (mx/rows x-training-mx))
                   (gen/vector
                     (gen/vector (s/gen ::m/finite) (mx/columns x-training-mx))
                     (count y-test))
                   (gen/double* {:min (first spread-interval)
                                 :max (second spread-interval)}))))))
  :ret (s/or :solution ::m/number
         :anomaly ::anomalies/anomaly))

(defn solve
  ""
  [{::keys [x-mx y x-new spread]}]
  (solve* (negative-half-distances x-mx x-new) spread y))

(s/def ::x-mx
  (s/and ::mx/matrix-finite
    (fn [mx]
      (pos? (mx/columns mx)))))

(s/def ::y ::vector/vector-finite)
(s/def ::x-new ::vector/vector-finite)
(s/def ::spread ::m/finite+)

(s/fdef solve
  :args (s/with-gen
          (s/and (s/cat :args (s/keys :req [::x-mx ::y ::x-new ::spread]))
            (fn [{:keys [args]}]
              (let [{::keys [x-mx y x-new]} args]
                (and (= (mx/rows x-mx) (count y))
                  (= (mx/columns x-mx) (count x-new))))))
          #(gen/bind (s/gen ::x-mx)
             (fn [x-mx]
               (gen/fmap (fn [[y x-new spread]]
                           [{::x-mx   x-mx
                             ::y      y
                             ::x-new  x-new
                             ::spread spread}])
                 (gen/tuple
                   (gen/vector (s/gen ::m/finite) (mx/rows x-mx))
                   (gen/vector (s/gen ::m/finite) (mx/columns x-mx))
                   (s/gen ::spread))))))
  :ret ::m/number)
