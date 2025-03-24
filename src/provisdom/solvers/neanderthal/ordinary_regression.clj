(ns provisdom.solvers.neanderthal.ordinary-regression
  (:require
    [clojure.spec.alpha :as s]
    [provisdom.math.core :as m]
    [provisdom.math.matrix :as mx]
    [provisdom.math.vector :as vector]
    [provisdom.neanderthal-matrix.core :as neanderthal-mx]
    [provisdom.utility-belt.anomalies :as anomalies]))

(s/def ::x-mx (mx/matrix-finite-spec {:min-columns 1 :min-rows 2}))
(s/def ::y (vector/vector-finite-spec {:min-count 2}))

(s/def ::condition-number-max ::neanderthal-mx/condition-number)

(s/def ::mean-squared-error ::m/non-)
(s/def ::standard-squared-error ::m/non-)
(s/def ::weights ::vector/vector-num)

(defn simple-ordinary-regression
  ""
  ([x-mx y]
   (simple-ordinary-regression x-mx y {}))
  ([x-mx y
    {::keys [condition-number-max]
     :or    {condition-number-max 1e15}}]
   (let [y-nmx (neanderthal-mx/matrix->neanderthal-matrix (mx/column-matrix y))
         x-nmx (neanderthal-mx/matrix->neanderthal-matrix x-mx)
         sol (neanderthal-mx/lls-with-error x-nmx y-nmx)]
     (if (anomalies/anomaly? sol)
       sol
       (let [{::neanderthal-mx/keys [mean-squared-errors
                                     standard-squared-errors
                                     solution
                                     condition-number]} sol
             mse (neanderthal-mx/neanderthal-matrix->matrix mean-squared-errors)
             sse (neanderthal-mx/neanderthal-matrix->matrix
                   standard-squared-errors)
             weights (first (neanderthal-mx/neanderthal-matrix->matrix
                              (neanderthal-mx/transpose solution)))
             mean-squared-error (max 0.0 (first (mx/diagonal mse)))
             standard-squared-error (max 0.0 (first (mx/diagonal sse)))]
         (cond
           (> condition-number condition-number-max)
           {::anomalies/message               "poorly conditioned matrix"
            ::anomalies/category              ::anomalies/exception
            ::neanderthal-mx/condition-number condition-number
            ::condition-number-max            condition-number-max
            ::anomalies/fn                    (var simple-ordinary-regression)}

           (not-any? m/nan? (concat weights [mean-squared-error
                                             standard-squared-error]))
           {::weights                weights
            ::mean-squared-error     mean-squared-error
            ::standard-squared-error standard-squared-error}

           :else
           {::anomalies/message  "Weights and errors can't be NaN."
            ::anomalies/category ::anomalies/exception
            ::anomalies/fn       (var simple-ordinary-regression)}))))))

(s/fdef simple-ordinary-regression
  :args (s/and (s/cat :x-mx ::x-mx
                      :y ::y
                      :opts (s/? (s/keys :opt [::condition-number-max])))
               (fn [{:keys [x-mx y]}]
                 (= (mx/rows x-mx) (count y))))
  :ret (s/or :anomaly ::anomalies/anomaly
             :solution (s/keys :req [::mean-squared-error
                                     ::standard-squared-error
                                     ::weights])))
