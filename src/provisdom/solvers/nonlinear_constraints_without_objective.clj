(ns provisdom.solvers.nonlinear-constraints-without-objective
  (:require
    [clojure.spec.alpha :as s]
    [clojure.spec.gen.alpha :as gen]
    [clojure.spec.test.alpha :as st]
    [orchestra.spec.test :as ost]
    [provisdom.utility-belt.anomalies :as anomalies]
    [provisdom.utility-belt.async :as async]
    [provisdom.math.core :as m]
    [provisdom.math.matrix :as mx]
    [provisdom.solvers.internal-apache-solvers :as apache-solvers]
    [provisdom.solvers.nonlinear-programming :as nlp]))

(s/def ::parallel? boolean?)

(defn- nls-selector-fn
  [met-accu]
  (fn [solutions]
    (when (first solutions)
      (let [[err pt errs]
            (reduce
              (fn [[smallest-error-sum
                    smallest-vector-point
                    smallest-weighted-constraint-errors] new-sol]
                (let [new-weighted-constraint-errors (::apache-solvers/weighted-constraint-errors new-sol)
                      new-vector-point (::apache-solvers/vector-point new-sol)
                      new-error-sum (if-not (anomalies/anomaly? new-sol)
                                      (reduce + new-weighted-constraint-errors)
                                      m/max-dbl)]
                  (if (and new-weighted-constraint-errors
                           (< new-error-sum smallest-error-sum))
                    [new-error-sum new-vector-point new-weighted-constraint-errors]
                    [smallest-error-sum smallest-vector-point smallest-weighted-constraint-errors])))
              [m/max-dbl nil nil]
              solutions)]
        (when errs
          {::weighted-constraint-errors errs
           ::vector-point               pt
           ::met-constraints?           (every? #(m/roughly? 0.0 % met-accu) errs)})))))

(defn nonlinear-least-squares
  "The default is to run `z` of the solvers, or choose `:levenberg-Marquardt`
  or `:gauss-newton`. If there is no solution, an anomaly will be returned when
  a single solver is run, otherwise nil.

  `::constraints-fn` takes and returns a vector. Each constraint function
  should return Inf+ or Inf- when out of range.

  `::constraint-jacobian-fn` is the jacobian matrix function that takes a vector
  and returns a matrix. `::constraint-jacobian-fn` should return nil when out of
  range. User may provide a numerical jacobian but performance may suffer.

  `::max-iter` - the maximum number of times to iterate in the algorithm.
  `::max-evaluations` - the maximum number of times to evaluate the model.

  `::constraint-weights` - the importance weights of each constraint. Each
   constraint error is multiplied by the square root of the weight.

   Returns map of `::vector-point`, `::weighted-constraint-errors`, and
   `::met-constraints?`."
  ([args] (nonlinear-least-squares args {}))
  ([{::keys [constraints-fn constraint-jacobian-fn vars-guess]}
    {::keys [max-iter rel-accu abs-accu max-evaluations check-by-objective?
             nls-solver-type constraint-weights met-accu parallel?]
     :or    {rel-accu            1e-14
             abs-accu            1e-6
             max-evaluations     1000
             check-by-objective? false
             nls-solver-type     :all
             met-accu            1e-6
             parallel?           false}}]
   (let [max-iter (or max-iter 1000)
         solvers (if-not (= :all nls-solver-type)
                   nls-solver-type
                   (list :levenberg-marquardt :gauss-newton))
         solver-fn (fn [solver-type]
                     #(apache-solvers/nonlinear-least-squares
                        {::apache-solvers/constraints-fn         constraints-fn
                         ::apache-solvers/constraint-jacobian-fn constraint-jacobian-fn
                         ::apache-solvers/vars-guess             vars-guess}
                        {::apache-solvers/max-iter            max-iter
                         ::apache-solvers/rel-accu            rel-accu
                         ::apache-solvers/abs-accu            abs-accu
                         ::apache-solvers/max-evaluations     max-evaluations
                         ::apache-solvers/check-by-objective? check-by-objective?
                         ::apache-solvers/use-gauss-newton?   (= solver-type :gauss-newton)
                         ::apache-solvers/constraint-weights  constraint-weights}))]
     (if (keyword? solvers)
       (let [solution ((solver-fn solvers))]
         (if (anomalies/anomaly? solution)
           solution
           {::weighted-constraint-errors (::apache-solvers/weighted-constraint-errors solution)
            ::vector-point               (::apache-solvers/vector-point solution)
            ::met-constraints?           (every? #(m/roughly? 0.0 % met-accu)
                                                 (::apache-solvers/weighted-constraint-errors solution))}))
       (async/thread-select
         (nls-selector-fn met-accu)
         (map solver-fn solvers)
         parallel?)))))

(s/def ::constraints-fn ::apache-solvers/constraints-fn)
(s/def ::constraint-jacobian-fn ::apache-solvers/constraint-jacobian-fn)
(s/def ::vars-guess ::apache-solvers/vars-guess)
(s/def ::max-iter ::apache-solvers/max-iter)
(s/def ::rel-accu ::apache-solvers/rel-accu)
(s/def ::abs-accu ::apache-solvers/abs-accu)
(s/def ::max-evaluations ::apache-solvers/max-evaluations)
(s/def ::check-by-objective? ::apache-solvers/check-by-objective?)

(s/def ::one-nls-solver-type #{:levenberg-marquardt :gauss-newton})

(s/def ::nls-solver-type
  (s/or :one ::one-nls-solver-type
        :all #{:all}
        :seq (s/coll-of ::one-nls-solver-type
                        :distinct true)))

(s/def ::constraint-weights ::apache-solvers/constraint-weights)

(s/def ::met-accu
  (s/with-gen ::m/finite+
              #(s/gen ::m/open-prob)))

(s/def ::weighted-constraint-errors ::apache-solvers/weighted-constraint-errors)
(s/def ::vector-point ::apache-solvers/vector-point)
(s/def ::met-constraints? boolean?)

(s/def ::constraints-with-jacobian-and-guess
  (s/with-gen
    (s/and (s/keys :req [::constraints-fn ::constraint-jacobian-fn ::vars-guess])
           (fn [{::keys [constraints-fn constraint-jacobian-fn vars-guess]}]
             (let [jac (constraint-jacobian-fn vars-guess)]
               (or (not jac)
                   (and (= (count (constraints-fn vars-guess)) (mx/rows jac))
                        (= (count vars-guess) (mx/columns jac)))))))
    #(gen/one-of
       (map gen/return
            (list {::constraints-fn         (fn [v]
                                              (if (= 2 (count v))
                                                (let [[a1 a2] v]
                                                  [(+ a1 (m/sq a2))
                                                   (- (m/cube a1) a2)])
                                                [m/inf+ m/inf+]))
                   ::constraint-jacobian-fn (fn [v]
                                              (when (= 2 (count v))
                                                (let [[a1 a2] v]
                                                  [[1.0 (* 3 (m/sq a1))]
                                                   [(* 2 a2) -1.0]])))
                   ::vars-guess             [2.0 -2.0]}
                  {::constraints-fn         (fn [v]
                                              (if (= 2 (count v))
                                                (let [[a1 a2] v]
                                                  [(+ a1 (m/exp a2))
                                                   (if (m/non+? a1)
                                                     m/inf-
                                                     (- (m/log a1) a2))])
                                                [m/inf+ m/inf+]))
                   ::constraint-jacobian-fn (fn [v]
                                              (when (= 2 (count v))
                                                (let [[a1 a2] v]
                                                  [[1.0 (if (not (zero? a1))
                                                          (- (/ a1))
                                                          m/inf-)]
                                                   [(m/exp a2) -1.0]])))
                   ::vars-guess             [1.9 -1.9]})))))

(s/fdef nonlinear-least-squares
  :args (s/cat :constraints-with-jacobian-and-guess ::constraints-with-jacobian-and-guess
               :opts (s/? (s/keys :opt [::max-iter ::rel-accu ::abs-accu
                                        ::max-evaluations ::check-by-objective?
                                        ::nls-solver-type ::constraint-weights
                                        ::met-accu])))
  :ret (s/nilable
         (s/or :solution (s/keys :req [::weighted-constraint-errors
                                       ::vector-point
                                       ::met-constraints?])
               :anomaly ::anomalies/anomaly)))

(defn- convert-constraints-fn-to-geq-fn
  [constraints-fn]
  (fn [da]
    (let [c-v (constraints-fn da)
          c-v-reverse (map - c-v)]
      (vec (concat c-v c-v-reverse)))))

(defn- convert-constraints-fn-to-nlp-objective
  "From last constraint."
  [constraints-fn]
  (fn [da]
    (peek (constraints-fn da))))

(defn- remove-constraint-of-constraints-fn
  [constraints-fn]
  (fn [da]
    (pop (constraints-fn da))))

(defn- remove-constraint-of-constraint-jacobian-fn
  [constraint-jacobian-fn]
  (fn [da]
    (let [jac (constraint-jacobian-fn da)]
      (when jac
        (mx/remove-row jac (dec (mx/rows jac)))))))

(defn- noc-selector-fn
  [constraint-count met-accu]
  (fn [returned-maps]
    (reduce
      (fn [best returned-map]
        (cond
          ;;test whether NLS and met all constraints
          (::met-constraints? returned-map)
          {::vector-point              (::vector-point returned-map)
           ::number-of-constraints-met constraint-count}
          ;;test whether NLP instead
          (::nlp/value returned-map)
          {::vector-point              (::nlp/vector-point returned-map)
           ::number-of-constraints-met (if (< (m/abs (::nlp/value returned-map)) met-accu)
                                         constraint-count
                                         (dec constraint-count))}
          ;;otherwise failed and returned anomaly
          :else best))
      [::vector-point nil
       ::number-of-constraints-met -1]
      returned-maps)))

(defn nonlinear-ordered-constraints
  "Each successive constraint is regarded as infinitely less important than the
  previous one.

  `::constraints-fn` takes and returns a vector. Each constraint function
  should return Inf+ or Inf- when out of range.

  `::constraint-jacobian-fn` is the jacobian matrix function that takes a vector
  and returns a matrix. `::constraint-jacobian-fn` should return nil when out of
  range. User may provide a numerical jacobian but performance may suffer.

  Returns map of `::vector-point`, `::number-of-constraints-not-met`, and
  `::squared-error-of-next-constraint`.

  Nonlinear programming and nonlinear least squares can be run in parallel.
  Nonlinear programming is used by making the square of the last constraint the
  objective. Nonlinear least squares is used by weighting each constraint
  significantly less than the previous one."
  ([args] (nonlinear-ordered-constraints args {}))
  ([{::keys [constraints-fn constraint-jacobian-fn vars-guess]}
    {::keys [max-iter rel-accu abs-accu check-by-objective? met-accu parallel?]
     :or    {rel-accu            m/dbl-close
             abs-accu            m/dbl-close
             check-by-objective? false
             met-accu            m/sgl-close
             parallel?           false}}]
   (let [max-iter (or max-iter 1000)
         constraint-count (count (constraints-fn vars-guess))
         constraint-weights (if (m/one? constraint-count)
                              (vector 1.0)
                              (vec (concat (repeat (dec constraint-count)
                                                   (/ m/sgl-close))
                                           [m/sgl-close])))
         nls #(nonlinear-least-squares
                {::constraints-fn         constraints-fn
                 ::constraint-jacobian-fn constraint-jacobian-fn
                 ::vars-guess             vars-guess}
                {::max-iter            max-iter
                 ::rel-accu            rel-accu
                 ::abs-accu            abs-accu
                 ::check-by-objective? check-by-objective?
                 ::constraint-weights  constraint-weights
                 ::met-accu            met-accu})]
     (if (>= constraint-count 2)
       (let [objective (convert-constraints-fn-to-nlp-objective constraints-fn)
             sub-constraints-fn (remove-constraint-of-constraints-fn constraints-fn)
             sub-constraint-jacobian-fn (remove-constraint-of-constraint-jacobian-fn
                                          constraint-jacobian-fn)
             geq-fn (convert-constraints-fn-to-geq-fn sub-constraints-fn)
             nlp #(nlp/constrained-nonlinear-programming
                    {::nlp/objective  objective
                     ::nlp/vars-guess vars-guess
                     ::nlp/geq-fn     geq-fn}
                    {::nlp/max-iter max-iter
                     ::nlp/abs-accu abs-accu})
             next #(nonlinear-ordered-constraints
                     {::constraints-fn         sub-constraints-fn
                      ::constraint-jacobian-fn sub-constraint-jacobian-fn
                      ::vars-guess             vars-guess}
                     {::max-iter            max-iter
                      ::rel-accu            rel-accu
                      ::abs-accu            abs-accu
                      ::check-by-objective? check-by-objective?
                      ::met-accu            met-accu})
             sol (when (> constraint-count 1)
                   (async/thread-select
                     (noc-selector-fn constraint-count met-accu)
                     [nls nlp]
                     parallel?))]
         (if (::vector-point sol)
           sol
           (next)))
       (let [sol (nls)]
         (if (anomalies/anomaly? sol)
           {::anomalies/category ::anomalies/no-solve
            ::anomalies/message  (format "Could not meet first constraint within %f." met-accu)
            ::anomalies/fn       (var nonlinear-ordered-constraints)}
           {::vector-point              (::vector-point sol)
            ::number-of-constraints-met 1}))))))

(s/def ::number-of-constraints-met ::m/long-non-)

(s/fdef nonlinear-ordered-constraints
  :args (s/cat :constraints-with-jacobian-and-guess ::constraints-with-jacobian-and-guess
               :opts (s/? (s/keys :opt [::max-iter ::rel-accu ::abs-accu
                                        ::check-by-objective? ::met-accu])))
  :ret (s/nilable
         (s/or :solution (s/keys :req [::vector-point ::number-of-constraints-met])
               :anomaly ::anomalies/anomaly)))