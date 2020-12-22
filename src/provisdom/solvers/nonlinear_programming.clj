(ns provisdom.solvers.nonlinear-programming
  (:require
    [clojure.spec.alpha :as s]
    [clojure.spec.gen.alpha :as gen]
    [clojure.spec.test.alpha :as st]
    [orchestra.spec.test :as ost]
    [provisdom.utility-belt.anomalies :as anomalies]
    [provisdom.utility-belt.async :as async]
    [provisdom.math.core :as m]
    [provisdom.math.vector :as vector]
    [provisdom.math.tensor :as tensor]
    [provisdom.math.derivatives :as derivatives]
    [provisdom.math.intervals :as intervals]
    [provisdom.solvers.internal-apache-solvers :as apache-solvers]
    [provisdom.solvers.internal-wrappers :as wrap]))

(s/def ::parallel? boolean?)

(s/def ::objective ::apache-solvers/objective)
(s/def ::vars-guess ::apache-solvers/vars-guess)
(s/def ::max-iter ::apache-solvers/max-iter)
(s/def ::goal ::apache-solvers/goal)
(s/def ::rel-accu ::apache-solvers/rel-accu)
(s/def ::abs-accu ::apache-solvers/abs-accu)
(s/def ::check-by-objective? ::apache-solvers/check-by-objective?)
(s/def ::value ::m/number)
(s/def ::vector-point ::vector/vector)
(s/def ::cobyla-initial-change ::wrap/initial-change)

(defn convert-objective
  [objective goal]
  (if (= goal :min)
    objective
    (fn [da]
      (- (objective da)))))

;;;CONSTRAINED
(defn constrained-nonlinear-programming
  "Solves using Cobyla.

  Returns a map of ::vector-point and ::value.

  `::objective` should take a double-array and return a number.
  `::vars-guess` should be a vector that is not all zeros.
  `::geq-fn` is the functional constraint, which should accept a double-array of
    variable values and return a vector of constraint values that must be
    greater than or equal to zero to meet the constraints.

  `::max-iter` - (100000) the maximum number of times to iterate in the
    algorithm. This is a lower bound on the size of the trust region. Final
    accuracy in the optimization is not precisely guaranteed.

  `::goal` - `:min` (default) or `:max`.
  `::abs-accu` - (1e-8) absolute accuracy for stopping.
  `::cobyla-initial-change` is a reasonable first change in the variable values
    (default 0.5).
  `::debug-print-level` (0, 1, 2, or 3) specifies the level of output to the
    console (default 0)."
  ([args] (constrained-nonlinear-programming args {}))
  ([{::keys [objective vars-guess geq-fn]}
    {::keys [max-iter goal abs-accu cobyla-initial-change debug-print-level]
     :or    {goal                  :min
             abs-accu              1e-8
             cobyla-initial-change 1.0
             debug-print-level     0}}]
   (let [max-iter (or max-iter 100000)
         cobyla-sol (wrap/nlp-cobyla
                      {::wrap/objective    (convert-objective objective goal)
                       ::wrap/geq-fn       geq-fn
                       ::wrap/cobyla-guess vars-guess}
                      {::wrap/max-iter          max-iter
                       ::wrap/initial-change    cobyla-initial-change
                       ::wrap/abs-accu          abs-accu
                       ::wrap/debug-print-level debug-print-level})]
     (if (anomalies/anomaly? cobyla-sol)
       cobyla-sol
       {::vector-point cobyla-sol
        ::value        (objective (double-array cobyla-sol))}))))

(s/def ::geq-fn ::wrap/geq-fn)
(s/def ::debug-print-level ::wrap/debug-print-level)

(s/def ::objective-constraints-and-guess
  (s/with-gen
    (s/keys :req [::objective ::vars-guess ::geq-fn])
    #(gen/one-of
       (map
         gen/return
         (list {::objective  (fn [da]
                               (if (= (count da) 2)
                                 (let [[a b] da]
                                   (+ a b))
                                 m/inf+))
                ::vars-guess [0.0 0.0]
                ::geq-fn     (fn [da]
                               (if (= (count da) 2)
                                 (let [[a b] da]
                                   [(dec a) (inc b)])
                                 [m/inf- m/inf-]))})))))

(s/fdef constrained-nonlinear-programming
  :args (s/cat :objective-constraints-and-guess ::objective-constraints-and-guess
               :opts (s/? (s/keys :opt [::max-iter ::goal ::abs-accu
                                        ::cobyla-initial-change
                                        ::debug-print-level])))
  :ret (s/or :solution (s/keys :req [::vector-point ::value])
             :anomaly ::anomalies/anomaly))

;;;UNBOUNDED
(defn- unbounded-selector-fn
  [goal]
  (fn [results]
    (when (pos? (count results))
      (reduce (fn [best next]
                (if (and (::apache-solvers/value next)
                         ((if (= goal :min) < >)
                          (::apache-solvers/value next)
                          (::apache-solvers/value best)))
                  next
                  best))
              {::apache-solvers/value        (if (= goal :min)
                                               m/max-dbl
                                               m/min-dbl)
               ::apache-solvers/vector-point [m/nan m/nan]}
              results))))

(defn unbounded-nonlinear-programming
  "Returns a map of ::vector-point and ::value.

  `::unbounded-solver-type` options:
   The default, `:all`, runs all solvers (can be parallel on different threads).
   Alternatively, choose one of the following, or a collection containing one or
   more of the following:
     `:cobyla`
     `:powell`
     `:nelder-mead`
     `:multi-directional-simplex`
     `:conjugate-gradient-using-polak-ribiere`
     `:conjugate-gradient-using-fletcher-reeves`.

  `::objective` should take a double-array and return a number.
  `::vars-guess` should be a vector.

  The `::gradient` function should take a double-array and return a vector. It
  is used with the conjugate gradient methods, otherwise a numerical derivative
  is used.

  `::max-iter` - (1000) the maximum number of times to iterate in the algorithm.
  `::goal` - `:min` (default) or `:max`.
  `::rel-accu` - (1e-14) relative accuracy for stopping.
  `::abs-accu` - (1e-6) absolute accuracy for stopping.
  `::check-by-objective?` - (false) Whether convergence checker uses objective
  function or point value for calculating error.
  `::initial-step-size-for-conjugate-gradient` - default is 1.0.
  `::cobyla-initial-change` is a reasonable first change in the variable values
  (default 0.5)."
  ([args] (unbounded-nonlinear-programming args {}))
  ([{::keys [objective vars-guess]}
    {::keys [gradient max-iter goal rel-accu abs-accu check-by-objective?
             unbounded-solver-type initial-step-size-for-conjugate-gradient
             cobyla-initial-change parallel?]
     :or    {goal                                     :min
             rel-accu                                 1e-14
             abs-accu                                 1e-6
             check-by-objective?                      false
             unbounded-solver-type                    :all
             initial-step-size-for-conjugate-gradient 1.0
             cobyla-initial-change                    0.5
             parallel?                                false}}]
   (let [max-iter (or max-iter 1000)
         gradient (or gradient (fn [da]
                                 ((derivatives/gradient-fn
                                    (fn [v]
                                      (objective (double-array v))))
                                  (vec da))))
         solvers (if-not (= :all unbounded-solver-type)
                   (if (keyword? unbounded-solver-type)
                     (list unbounded-solver-type)
                     unbounded-solver-type)
                   (list :cobyla :multi-directional-simplex :nelder-mead :powell
                         :conjugate-gradient-using-polak-ribiere
                         :conjugate-gradient-using-fletcher-reeves))
         solver-fn (fn [solver-type]
                     (case solver-type
                       :cobyla #(let [cobyla-sol (wrap/nlp-cobyla
                                                   {::wrap/objective    (convert-objective objective goal)
                                                    ::wrap/cobyla-guess vars-guess
                                                    ::wrap/geq-fn       (constantly [])}
                                                   {::wrap/max-iter       max-iter
                                                    ::wrap/initial-change cobyla-initial-change
                                                    ::wrap/abs-accu       abs-accu})]
                                  (if (anomalies/anomaly? cobyla-sol)
                                    cobyla-sol
                                    {::apache-solvers/vector-point cobyla-sol
                                     ::apache-solvers/value        (objective (double-array cobyla-sol))}))

                       (:conjugate-gradient-using-polak-ribiere
                         :conjugate-gradient-using-fletcher-reeves)
                       #(apache-solvers/optimize-without-constraints-and-with-gradient
                          {::apache-solvers/objective  objective
                           ::apache-solvers/gradient   gradient
                           ::apache-solvers/vars-guess vars-guess}
                          {::apache-solvers/max-iter                    max-iter
                           ::apache-solvers/goal                        goal
                           ::apache-solvers/rel-accu                    rel-accu
                           ::apache-solvers/abs-accu                    abs-accu
                           ::apache-solvers/check-by-objective?         check-by-objective?
                           ::apache-solvers/update-using-polak-ribiere? (= solver-type
                                                                           :conjugate-gradient-using-polak-ribiere)
                           ::apache-solvers/initial-step-size           initial-step-size-for-conjugate-gradient})

                       #(apache-solvers/optimize-without-constraints-or-gradient
                          objective
                          vars-guess
                          {::apache-solvers/max-iter                                        max-iter
                           ::apache-solvers/goal                                            goal
                           ::apache-solvers/rel-accu                                        rel-accu
                           ::apache-solvers/abs-accu                                        abs-accu
                           ::apache-solvers/without-constraints-or-gradient-nlp-solver-type solver-type})))
         apache-solution (if (m/one? (count solvers))
                           ((solver-fn (first solvers)))
                           (async/thread-select
                             (unbounded-selector-fn goal)
                             (map solver-fn solvers)
                             parallel?))]
     (cond (anomalies/anomaly? apache-solution) apache-solution
           (nil? apache-solution) {::anomalies/category ::anomalies/no-solve
                                   ::anomalies/message  "No solution."
                                   ::anomalies/fn       (var unbounded-nonlinear-programming)}
           :else {::value        (::apache-solvers/value apache-solution)
                  ::vector-point (::apache-solvers/vector-point apache-solution)}))))

(s/def ::one-unbounded-solver-type
  (s/or :apache-solvers ::apache-solvers/without-constraints-or-gradient-nlp-solver-type
        :not-apache-solvers #{:conjugate-gradient-using-polak-ribiere
                              :conjugate-gradient-using-fletcher-reeves
                              :cobyla}))

(s/def ::unbounded-solver-type
  (s/or :one ::one-unbounded-solver-type
        :all #{:all}
        :seq (s/coll-of ::one-unbounded-solver-type
                        :distinct true)))

(s/def ::gradient ::apache-solvers/gradient)
(s/def ::initial-step-size-for-conjugate-gradient ::apache-solvers/initial-step-size)

(s/def ::objective-with-guess
  (s/with-gen
    (s/keys :req [::objective ::vars-guess])
    #(gen/one-of
       (map
         gen/return
         (list {::objective  (fn [da]
                               (if (= (count da) 2)
                                 (let [[a b] da]
                                   (+ (m/sq a) (m/sq b)))
                                 m/inf+))
                ::vars-guess [0.0 0.0]})))))

(s/fdef unbounded-nonlinear-programming
  :args (s/cat :objective-with-guess ::objective-with-guess
               :opts (s/? (s/keys :opt [::gradient ::max-iter ::goal ::rel-accu
                                        ::abs-accu ::check-by-objective?
                                        ::unbounded-solver-type
                                        ::initial-step-size-for-conjugate-gradient
                                        ::cobyla-initial-change])))
  :ret (s/or :solution (s/keys :req [::vector-point ::value])
             :anomaly ::anomalies/anomaly))

;;;BOUNDED
(s/def ::bobyqa-interpolation-points ::apache-solvers/bobyqa-interpolation-points)
(s/def ::var-intervals (s/coll-of ::intervals/finite-interval :kind vector? :into []))

(s/def ::objective-with-guess-and-intervals
  (s/with-gen
    (s/and (s/keys :req [::objective ::vars-guess ::var-intervals])
           (fn [{::keys [vars-guess var-intervals]}]
             (and (= (count vars-guess) (count var-intervals))
                  (every? true? (map (fn [vi vg]
                                       (intervals/in-interval? vi vg))
                                     var-intervals
                                     vars-guess)))))
    #(gen/one-of
       (map
         gen/return
         (list {::objective     (fn [da]
                                  (if (= (count da) 2)
                                    (let [[a b] da]
                                      (+ a b))
                                    m/inf+))
                ::vars-guess    [0.0 0.0]
                ::var-intervals [[-1.0 1.0] [-1.0 1.0]]})))))

(s/def ::met-bounds-accu
  (s/with-gen ::m/finite-non-
              #(s/gen ::m/prob)))

(defn- convert-bounds-to-geq-fn
  [var-lower-bounds var-upper-bounds]
  (fn [da]
    (vec (concat (map (fn [d lb]
                        (- d lb))
                      da
                      var-lower-bounds)
                 (map (fn [d ub]
                        (- ub d))
                      da
                      var-upper-bounds)))))

(defn- bounded-selector-fn
  [goal met-bounds-accu var-lower-bounds var-upper-bounds]
  (fn [results]
    (when (pos? (count results))
      (reduce (fn [best next]
                (if (and (::apache-solvers/value next)
                         ((if (= goal :min) < >)
                          (::apache-solvers/value next)
                          (::apache-solvers/value best))
                         (every? (fn [[variable lower-bound upper-bound]]
                                   (intervals/in-interval-roughly?
                                     [lower-bound upper-bound]
                                     variable
                                     met-bounds-accu))
                                 (map vector
                                      (::apache-solvers/vector-point next)
                                      var-lower-bounds
                                      var-upper-bounds)))
                  next
                  best))
              {::apache-solvers/value        (if (= goal :min) m/max-dbl m/min-dbl)
               ::apache-solvers/vector-point [m/nan m/nan]}
              results))))

(defn bounded-nonlinear-programming-without-evolutionary
  "To keep it pure, this function does not include an evolutionary solver, like
  [[bounded-nonlinear-programming-including-evolutionary!]].

  Returns a map of ::vector-point and ::value.

  `::bounded-without-evolutionary-solver-type` options:
   The default is `:cobyla`. Alternatively, can choose `:bobyqa`. Or, choose
   `:all` to run all solvers in parallel on different threads.

  `::objective` should take a double-array and return a number.
  `::vars-guess` should be a vector.

  `::max-iter` - (1000) the maximum number of times to iterate in the algorithm.
  `::goal` - `:min` (default) or `:max`.
  `::abs-accu` - (1e-8) absolute accuracy for stopping.
  `::met-bounds-accu` - (1e-6) accuracy for treating bounds as met during
    threading.
  `::cobyla-initial-change` is a reasonable first change in the variable values
  (default 0.5).

  Powell's BOBYQA algorithm (Bound Optimization BY Quadratic Approximation).
  `:bobyqa` is particularly well suited for high dimensional problems where
  derivatives are not available. In most cases it outperforms the
  PowellOptimizer significantly. `:bobyqa` could also be considered as a
  replacement of any derivative-based optimizer when the derivatives are
  approximated by finite differences. Recommend setting the
  `::bobyqa-interpolation-points` in [n+2, (n+1)(n+2)/2], where 'n' is the
  number of dimensions."
  ([args] (bounded-nonlinear-programming-without-evolutionary args {}))
  ([{::keys [objective vars-guess var-intervals]}
    {::keys [max-iter goal abs-accu bounded-without-evolutionary-solver-type
             met-bounds-accu bobyqa-interpolation-points cobyla-initial-change
             parallel?]
     :or    {goal                                     :min
             abs-accu                                 1e-8
             bounded-without-evolutionary-solver-type :cobyla
             met-bounds-accu                          1e-6
             cobyla-initial-change                    0.5
             parallel?                                false}}]
   (let [max-iter (or max-iter 1000)
         solvers (if-not (= :all bounded-without-evolutionary-solver-type)
                   (if (keyword? bounded-without-evolutionary-solver-type)
                     (list bounded-without-evolutionary-solver-type)
                     bounded-without-evolutionary-solver-type)
                   (list :bobyqa :cobyla))
         var-lower-bounds (mapv first var-intervals)
         var-upper-bounds (mapv second var-intervals)
         solver-fn (fn [solver-type]
                     (case solver-type
                       :cobyla #(let [cobyla-sol (wrap/nlp-cobyla
                                                   {::wrap/objective    (convert-objective objective goal)
                                                    ::wrap/geq-fn       (convert-bounds-to-geq-fn
                                                                          var-lower-bounds var-upper-bounds)
                                                    ::wrap/cobyla-guess vars-guess}
                                                   {::wrap/max-iter       max-iter
                                                    ::wrap/initial-change cobyla-initial-change
                                                    ::wrap/abs-accu       abs-accu})]
                                  (if (anomalies/anomaly? cobyla-sol)
                                    cobyla-sol
                                    {::apache-solvers/vector-point cobyla-sol
                                     ::apache-solvers/value        (objective (double-array cobyla-sol))}))
                       :bobyqa #(apache-solvers/optimize-bobyqa
                                  objective var-lower-bounds var-upper-bounds vars-guess
                                  {::apache-solvers/bobyqa-interpolation-points bobyqa-interpolation-points
                                   ::apache-solvers/max-iter                    max-iter
                                   ::apache-solvers/goal                        goal})))
         apache-solution (if (m/one? (count solvers))
                           ((solver-fn (first solvers)))
                           (async/thread-select
                             (bounded-selector-fn goal met-bounds-accu var-lower-bounds var-upper-bounds)
                             (map solver-fn solvers)
                             parallel?))]
     (cond (anomalies/anomaly? apache-solution) apache-solution
           (nil? apache-solution) {::anomalies/category ::anomalies/no-solve
                                   ::anomalies/message  "No solution."
                                   ::anomalies/fn       (var bounded-nonlinear-programming-without-evolutionary)}
           :else {::value        (::apache-solvers/value apache-solution)
                  ::vector-point (::apache-solvers/vector-point apache-solution)}))))

(s/def ::one-bounded-without-evolutionary-solver-type #{:bobyqa :cobyla})

(s/def ::bounded-without-evolutionary-solver-type
  (s/or :one ::one-bounded-without-evolutionary-solver-type
        :all #{:all}
        :seq (s/coll-of ::one-bounded-without-evolutionary-solver-type
                        :distinct true)))

(s/fdef bounded-nonlinear-programming-without-evolutionary
  :args (s/cat :objective-with-guess-and-intervals ::objective-with-guess-and-intervals
               :opts (s/? (s/keys :opt [::max-iter ::goal ::abs-accu
                                        ::bounded-without-evolutionary-solver-type
                                        ::met-bounds-accu ::bobyqa-interpolation-points
                                        ::cobyla-initial-change])))
  :ret (s/or :solution (s/keys :req [::vector-point ::value])
             :anomaly ::anomalies/anomaly))

(defn bounded-nonlinear-programming-including-evolutionary$
  "This function always includes an evolutionary solver, unlike
  [[bounded-without-evolutionary-nonlinear-programming]].

  Returns a map of ::vector-point and ::value.

  `::bounded-including-evolutionary-solver-type` options:
  The default, `:all`, runs all solvers simultaneously on different threads.
  Alternatively, choose one of the following, or a collection containing one or
  more of the following:
  `:cobyla`
  `:bobyqa`
  `:cma-es-only`.

  `::objective` should take a double-array and return a number.
  `::vars-guess` should be a vector.

  `::max-iter` - (1000) the maximum number of times to iterate in the algorithm.
  `::goal` - `:min` (default) or `:max`.
  `::abs-accu` - (1e-8) absolute accuracy for stopping.
  `::met-bounds-accu` - (1e-6) accuracy for treating bounds as met during
    threading.
  `::cobyla-initial-change` is a reasonable first change in the variable values
    (default 0.5).

  An implementation of the active Covariance Matrix Adaptation Evolution
  Strategy (CMA-ES) for non-linear, non-convex, non-smooth, global function
  minimization. See https://www.lri.fr/~hansen/cmaesintro.html for more on
  CMA-ES. The CMA-Evolution Strategy (CMA-ES) is a reliable stochastic
  optimization method which should be applied if derivative-based methods,
  e.g. quasi-Newton BFGS or conjugate gradient, fail due to a rugged search
  landscape (e.g. noise, local optima, outlier, etc.) of the `objective`
  function. Like a quasi-Newton method, the CMA-ES learns and applies a
  variable metric on the underlying search space. Unlike a quasi-Newton method,
  the CMA-ES neither estimates nor uses gradients, making it considerably more
  reliable in terms of finding a good, or even close to optimal, solution. In
  general, on smooth objective functions the CMA-ES is roughly ten times slower
  than BFGS (counting objective function evaluations, no gradient provided). For
  up to N=10 variables, the derivative-free simplex direct search method (Nelder
  and Mead) can be faster, but it is far less reliable than CMA-ES. The CMA-ES
  is particularly well suited for non-separable and/or badly conditioned
  problems. To observe the advantage of CMA compared to a conventional evolution
  strategy, it will usually take about 30 N function evaluations. On difficult
  problems, the complete optimization (a single run) is expected to take roughly
  between 30 N and 300 N^2 function evaluations.

  `::sigmas` - Input sigma values. They define the initial coordinate-wise
  standard deviations for sampling new search points around the initial
  `vars-guess`. It is suggested to set them to the estimated distance from the
  initial to the desired optimum. Small values induce the search to be more
  local (and very small values are more likely to find a local optimum close to
  the `vars-guess`). Too small values might however lead to early termination.

  `::population-size` - Number of offspring and is the primary strategy
  parameter. In the absence of better clues, a good default could be an integer
  close to 4 + 3 ln(n), where 'n' is the number of optimized parameters.
  Increasing the `::population-size` improves global search properties at the
  expense of speed (which in general decreases at most linearly with increasing
  `::population-size`).

  `::stop-fitness` - Stop if `objective` function value is smaller than
  `::stop-fitness`.

  `::active-cma?` - Chooses the covariance matrix update method.

  `::diagonal-only` - Number of initial iterations, where the covariance matrix
  remains diagonal.

  `::check-feasible-count` - Determines how often new random objective variables
  are generated in case they are out of bounds.

  Powell's BOBYQA algorithm (Bound Optimization BY Quadratic Approximation).
  Faster but less robust than CMA-ES. `:bobyqa` is particularly well suited for
  high dimensional problems where derivatives are not available. In most cases
  it outperforms the PowellOptimizer significantly. `:bobyqa` could also be
  considered as a replacement of any derivative-based optimizer when the
  derivatives are approximated by finite differences. Recommend setting the
  `::bobyqa-interpolation-points` in [n+2, (n+1)(n+2)/2], where 'n' is the
  number of dimensions."
  ([args] (bounded-nonlinear-programming-including-evolutionary$ args {}))
  ([{::keys [objective vars-guess var-intervals]}
    {::keys [max-iter goal rel-accu abs-accu check-by-objective?
             bounded-including-evolutionary-solver-type met-bounds-accu sigmas
             population-size stop-fitness active-cma? diagonal-only
             check-feasible-count bobyqa-interpolation-points cobyla-initial-change
             parallel?]
     :or    {goal                                       :min
             rel-accu                                   1e-14
             abs-accu                                   1e-8
             check-by-objective?                        false
             bounded-including-evolutionary-solver-type :all
             met-bounds-accu                            1e-6
             stop-fitness                               0.0
             active-cma?                                true
             diagonal-only                              0
             check-feasible-count                       0
             cobyla-initial-change                      0.5
             parallel?                                  false}}]
   (let [max-iter (or max-iter 1000)
         solvers (if-not (= :all bounded-including-evolutionary-solver-type)
                   (if (keyword? bounded-including-evolutionary-solver-type)
                     (list bounded-including-evolutionary-solver-type)
                     bounded-including-evolutionary-solver-type)
                   (list :bobyqa :cobyla))
         solvers (distinct (conj solvers :cma-es-only))
         var-lower-bounds (mapv first var-intervals)
         var-upper-bounds (mapv second var-intervals)
         solver-fn (fn [solver-type]
                     (case solver-type
                       :cobyla #(let [cobyla-sol (wrap/nlp-cobyla
                                                   {::wrap/objective    (convert-objective objective goal)
                                                    ::wrap/geq-fn       (convert-bounds-to-geq-fn
                                                                          var-lower-bounds var-upper-bounds)
                                                    ::wrap/cobyla-guess vars-guess}
                                                   {::wrap/max-iter       max-iter
                                                    ::wrap/initial-change cobyla-initial-change
                                                    ::wrap/abs-accu       abs-accu})]
                                  (if (anomalies/anomaly? cobyla-sol)
                                    cobyla-sol
                                    {::apache-solvers/vector-point cobyla-sol
                                     ::apache-solvers/value        (objective (double-array cobyla-sol))}))
                       :bobyqa #(apache-solvers/optimize-bobyqa
                                  objective var-lower-bounds var-upper-bounds vars-guess
                                  {::apache-solvers/bobyqa-interpolation-points bobyqa-interpolation-points
                                   ::apache-solvers/max-iter                    max-iter
                                   ::apache-solvers/goal                        goal})
                       :cma-es-only #(apache-solvers/optimize-cma-evolution$
                                       objective var-lower-bounds var-upper-bounds vars-guess
                                       {::apache-solvers/max-iter             max-iter
                                        ::apache-solvers/goal                 goal
                                        ::apache-solvers/rel-accu             rel-accu
                                        ::apache-solvers/abs-accu             abs-accu
                                        ::apache-solvers/check-by-objective?  check-by-objective?
                                        ::apache-solvers/sigmas               sigmas
                                        ::apache-solvers/population-size      population-size
                                        ::apache-solvers/stop-fitness         stop-fitness
                                        ::apache-solvers/active-cma?          active-cma?
                                        ::apache-solvers/diagonal-only        diagonal-only
                                        ::apache-solvers/check-feasible-count check-feasible-count})))
         apache-solution (if (m/one? (count solvers))
                           ((solver-fn (first solvers)))
                           (async/thread-select
                             (bounded-selector-fn goal met-bounds-accu var-lower-bounds var-upper-bounds)
                             (map solver-fn solvers)
                             parallel?))]
     (cond (anomalies/anomaly? apache-solution) apache-solution
           (nil? apache-solution) {::anomalies/category ::anomalies/no-solve
                                   ::anomalies/message  "No solution."
                                   ::anomalies/fn       (var bounded-nonlinear-programming-including-evolutionary$)}
           :else {::value        (::apache-solvers/value apache-solution)
                  ::vector-point (::apache-solvers/vector-point apache-solution)}))))

(s/def ::sigmas ::apache-solvers/sigmas)
(s/def ::population-size ::apache-solvers/population-size)
(s/def ::stop-fitness ::apache-solvers/stop-fitness)
(s/def ::active-cma? ::apache-solvers/active-cma?)
(s/def ::diagonal-only ::apache-solvers/diagonal-only)
(s/def ::check-feasible-count ::apache-solvers/check-feasible-count)
(s/def ::one-bounded-including-evolutionary-solver-type #{:bobyqa :cobyla :cma-es-only})

(s/def ::bounded-including-evolutionary-solver-type
  (s/or :one ::one-bounded-including-evolutionary-solver-type
        :all #{:all}
        :seq (s/coll-of ::one-bounded-including-evolutionary-solver-type
                        :distinct true)))

(s/fdef bounded-nonlinear-programming-including-evolutionary$
  :args (s/cat :objective-with-guess-and-intervals ::objective-with-guess-and-intervals
               :opts (s/? (s/keys :opt [::max-iter ::goal ::rel-accu
                                        ::abs-accu ::check-by-objective?
                                        ::bounded-including-evolutionary-solver-type
                                        ::met-bounds-accu ::sigmas ::population-size
                                        ::stop-fitness ::active-cma? ::diagonal-only
                                        ::check-feasible-count
                                        ::bobyqa-interpolation-points
                                        ::cobyla-initial-change])))
  :ret (s/or :solution (s/keys :req [::vector-point ::value])
             :anomaly ::anomalies/anomaly))