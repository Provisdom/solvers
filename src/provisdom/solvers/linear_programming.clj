(ns provisdom.solvers.linear-programming
  (:require
    [clojure.spec.alpha :as s]
    [provisdom.math.vector :as vector]
    [provisdom.solvers.internal-apache-solvers :as apache-solvers]
    [provisdom.utility-belt.anomalies :as anomalies]))

(defn linear-programming
  "Two-phase Simplex Method. Returns map of {::value, ::vector-point}.
  `objective-coefficients` should be a vector. `linear-constraints` should be
  a collection of triples of [linear-coefficients relation constraint-value],
  where linear-coefficients is a vector; relation is either `:eq`, `:leq`, or
  `:geq`; and constraint-value is a number. When '::non-negative-vars?` is set
  to true, all variables will be non-negative. `::goal` can be `:min` or
  `:max`."
  ([objective-coefficients linear-constraints]
   (linear-programming objective-coefficients linear-constraints {}))
  ([objective-coefficients linear-constraints
    {::keys [goal objective-constant non-negative-vars?]
     :or    {goal               :min
             objective-constant 0.0
             non-negative-vars? false}}]
   (let [lin-sol (apache-solvers/linear-programming
                   objective-coefficients linear-constraints
                   {::apache-solvers/goal               goal
                    ::apache-solvers/objective-constant objective-constant
                    ::apache-solvers/non-negative-vars? non-negative-vars?})]
     (if (anomalies/anomaly? lin-sol)
       lin-sol
       {::value (::apache-solvers/value lin-sol)
        ::vector-point (::apache-solvers/vector-point lin-sol)}))))

(s/def ::linear-constraints ::apache-solvers/linear-constraints)
(s/def ::goal ::apache-solvers/goal)
(s/def ::objective-constant ::apache-solvers/objective-constant)
(s/def ::non-negative-vars? ::apache-solvers/non-negative-vars?)
(s/def ::value ::apache-solvers/value)
(s/def ::vector-point ::apache-solvers/vector-point)

(s/fdef linear-programming
        :args (s/and (s/cat :objective-coefficients ::vector/vector-finite
                            :linear-constraints ::linear-constraints
                            :opts (s/? (s/keys :opts [::goal
                                                      ::objective-constant
                                                      ::non-negative-vars?])))
                     (fn [{:keys [objective-coefficients linear-constraints]}]
                       (= (count objective-coefficients)
                          (count (ffirst linear-constraints)))))
        :ret (s/or :solution (s/keys :req [::value ::vector-point])
                   :anomaly ::anomalies/anomaly))
