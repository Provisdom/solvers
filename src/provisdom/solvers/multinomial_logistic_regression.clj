(ns provisdom.solvers.multinomial-logistic-regression
  (:require
    [clojure.spec.alpha :as s]
    [clojure.spec.gen.alpha :as gen]
    [clojure.spec.test.alpha :as st]
    [orchestra.spec.test :as ost]
    [provisdom.utility-belt.anomalies :as anomalies]
    [provisdom.math.core :as m]
    [provisdom.math.vector :as vector]
    [provisdom.math.matrix :as mx]
    [provisdom.solvers.logistic-regression :as log-regress]))

(s/def ::max-iter
  (s/with-gen (s/nilable ::m/int+)
              #(gen/one-of [(s/gen (s/int-in 100 1000))
                            (gen/return nil)])))

(s/def ::abs-accu
  (s/with-gen ::m/finite+
              #(s/gen (s/double-in :min m/tiny-dbl
                                   :max 1.0
                                   :NaN? false))))

(defn iterative-reweighted-least-squares
  ""
  ([x-mx y-mx] (iterative-reweighted-least-squares x-mx y-mx {}))
  ([x-mx y-mx {::keys [max-iter abs-accu] :or {abs-accu 1e-15}}]
   (let [yt (mx/transpose y-mx)
         sols (vec
                (pmap (fn [y]
                        (log-regress/iterative-reweighted-least-squares
                          x-mx
                          y
                          {::log-regress/abs-accu abs-accu
                           ::log-regress/max-iter max-iter}))
                      yt))]
     sols)))

(s/fdef iterative-reweighted-least-squares
        :args (s/and (s/cat :x-mx ::mx/matrix-finite
                            :y-mx ::mx/matrix-prob
                            :opts (s/? (s/keys :opt [::max-iter ::abs-accu])))
                     (fn [{:keys [x-mx y-mx]}]
                       (let [y-t (mx/transpose y-mx)]
                         (and
                           (= (mx/rows x-mx) (mx/rows y-mx))
                           (pos? (mx/columns x-mx))
                           (>= (mx/rows y-mx) 2)
                           (every? (fn [y]
                                     (and (not (every? m/one? y))
                                          (not (every? zero? y))))
                                   y-t)))))
        :ret (s/coll-of (s/or :sol ::vector/vector-finite
                              :anomaly ::anomalies/anomaly)
                        :kind vector?
                        :into []))