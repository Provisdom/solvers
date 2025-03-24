(ns provisdom.solvers.optimize-univariate
  (:require
    [clojure.spec.alpha :as s]
    [provisdom.math.core :as m]
    [provisdom.math.intervals :as intervals]
    [provisdom.solvers.internal-apache-solvers :as apache-solvers]
    [provisdom.utility-belt.anomalies :as anomalies]))

(defn optimize-univariate
  "Brent Optimizer.  Search over a `::strict-finite-interval`.
  `::guess` must be in `strict-finite-interval`."
  ([args] (optimize-univariate args {}))
  ([{::keys [univariate-f strict-finite-interval guess]}
    {::keys [max-iter goal rel-accu abs-accu]
     :or    {goal     :min
             rel-accu 1e-14
             abs-accu 1e-6}}]
   (let [max-iter (or max-iter 1000)
         sol (apache-solvers/optimize-univariate
               univariate-f
               strict-finite-interval
               guess
               {::apache-solvers/max-iter max-iter
                ::apache-solvers/goal     goal
                ::apache-solvers/rel-accu rel-accu
                ::apache-solvers/abs-accu abs-accu})]
     (if (anomalies/anomaly? sol)
       sol
       {::value (::apache-solvers/value sol)
        ::point (::apache-solvers/point sol)}))))

(s/def ::univariate-f ::apache-solvers/univariate-f)
(s/def ::strict-finite-interval ::apache-solvers/strict-finite-interval)
(s/def ::guess ::apache-solvers/initial-guess)
(s/def ::max-iter ::apache-solvers/max-iter)
(s/def ::goal ::apache-solvers/goal)
(s/def ::rel-accu ::apache-solvers/rel-accu)
(s/def ::abs-accu ::apache-solvers/abs-accu)
(s/def ::value ::m/number)
(s/def ::point ::m/number)

(s/fdef optimize-univariate
  :args (s/and (s/cat :args (s/keys :req [::guess
                                          ::strict-finite-interval
                                          ::univariate-f])
                 :opts (s/? (s/keys :opt [::abs-accu
                                          ::goal
                                          ::max-iter
                                          ::rel-accu])))
          (fn [{:keys [args]}]
            (let [{::keys [strict-finite-interval guess]} args]
              (intervals/in-bounds?
                (intervals/bounds strict-finite-interval false false)
                guess))))
  :ret (s/or :solution (s/keys :req [::value ::point])
         :anomaly ::anomalies/anomaly))
