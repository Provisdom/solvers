(ns provisdom.solvers.optimize-univariate
  (:require
    [clojure.spec.alpha :as s]
    [clojure.spec.gen.alpha :as gen]
    [clojure.spec.test.alpha :as st]
    [orchestra.spec.test :as ost]
    [provisdom.utility-belt.anomalies :as anomalies]
    [provisdom.math.intervals :as intervals]
    [provisdom.math.core :as m]
    [provisdom.solvers.internal-apache-solvers :as apache-solvers]))

(defn optimize-univariate
  "Brent Optimizer.  Search over a `::finite-interval`.
  `::guess` must be in `finite-interval` and less than the maximum."
  ([args] (optimize-univariate args {}))
  ([{::keys [univariate-f finite-interval guess]}
    {::keys [max-iter goal rel-accu abs-accu]
     :or    {goal     :min
             rel-accu 1e-14
             abs-accu 1e-6}}]
   (let [max-iter (or max-iter 1000)
         sol (apache-solvers/optimize-univariate
               univariate-f
               finite-interval
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
(s/def ::finite-interval ::apache-solvers/finite-interval)
(s/def ::guess ::apache-solvers/initial-guess)
(s/def ::max-iter ::apache-solvers/max-iter)
(s/def ::goal ::apache-solvers/goal)
(s/def ::rel-accu ::apache-solvers/rel-accu)
(s/def ::abs-accu ::apache-solvers/abs-accu)
(s/def ::value ::m/number)
(s/def ::point ::m/number)

(s/fdef optimize-univariate
        :args (s/and (s/cat :args (s/keys :req [::univariate-f ::finite-interval ::guess])
                            :opts (s/? (s/keys :opt [::max-iter ::goal ::rel-accu ::abs-accu])))
                     (fn [{:keys [args]}]
                       (let [{::keys [finite-interval guess]} args]
                         (and (< (first finite-interval)
                                 (second finite-interval))
                              (intervals/in-bounds?
                                (intervals/bounds finite-interval true false)
                                guess)))))
        :ret (s/or :solution (s/keys :req [::value ::point])
                   :anomaly ::anomalies/anomaly))