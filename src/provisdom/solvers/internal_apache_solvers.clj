(ns provisdom.solvers.internal-apache-solvers
  (:require
    [clojure.spec.alpha :as s]
    [clojure.spec.gen.alpha :as gen]
    [clojure.spec.test.alpha :as st]
    [orchestra.spec.test :as ost]
    [provisdom.utility-belt.anomalies :as anomalies]
    [provisdom.math.core :as m]
    [provisdom.math.arrays :as arrays]
    [provisdom.math.vector :as vector]
    [provisdom.math.matrix :as mx]
    [provisdom.math.tensor :as tensor]
    [provisdom.math.random :as random]
    [provisdom.math.intervals :as intervals]
    [provisdom.apache-math.apache-vector :as apache-v]
    [provisdom.apache-math.apache-matrix :as apache-mx]
    [provisdom.apache-math.alternative-random :as alt-random])
  (:import
    [java.util ArrayList]
    [org.apache.commons.math3.exception
     TooManyEvaluationsException TooManyIterationsException]
    [org.apache.commons.math3.analysis
     UnivariateFunction MultivariateFunction MultivariateVectorFunction
     MultivariateMatrixFunction]
    [org.apache.commons.math3.analysis.differentiation
     UnivariateDifferentiableFunction FiniteDifferencesDifferentiator]
    [org.apache.commons.math3.analysis.solvers
     BaseUnivariateSolver BisectionSolver BracketingNthOrderBrentSolver
     BrentSolver IllinoisSolver MullerSolver MullerSolver2 PegasusSolver
     RegulaFalsiSolver RiddersSolver SecantSolver NewtonRaphsonSolver]
    [org.apache.commons.math3.optim
     OptimizationData InitialGuess SimpleBounds PointValuePair
     SimpleValueChecker SimplePointChecker SimpleVectorValueChecker MaxEval]
    [org.apache.commons.math3.optim.univariate
     BrentOptimizer SearchInterval UnivariateObjectiveFunction
     UnivariatePointValuePair]
    [org.apache.commons.math3.optim.linear
     SimplexSolver LinearObjectiveFunction LinearConstraintSet LinearConstraint
     Relationship NonNegativeConstraint]
    [org.apache.commons.math3.optim.nonlinear.scalar
     GoalType ObjectiveFunction ObjectiveFunctionGradient]
    [org.apache.commons.math3.optim.nonlinear.scalar.noderiv
     SimplexOptimizer MultiDirectionalSimplex NelderMeadSimplex PowellOptimizer
     BOBYQAOptimizer CMAESOptimizer CMAESOptimizer$Sigma
     CMAESOptimizer$PopulationSize]
    [org.apache.commons.math3.optim.nonlinear.scalar.gradient
     NonLinearConjugateGradientOptimizer
     NonLinearConjugateGradientOptimizer$BracketingStep
     NonLinearConjugateGradientOptimizer$Formula]
    [org.apache.commons.math3.fitting.leastsquares
     LevenbergMarquardtOptimizer GaussNewtonOptimizer LeastSquaresFactory
     LeastSquaresProblem$Evaluation]
    [org.apache.commons.math3.analysis.interpolation
     BicubicInterpolator BicubicInterpolatingFunction
     PiecewiseBicubicSplineInterpolator
     PiecewiseBicubicSplineInterpolatingFunction TricubicInterpolator
     TricubicInterpolatingFunction DividedDifferenceInterpolator
     HermiteInterpolator LinearInterpolator LoessInterpolator
     MicrosphereInterpolator NevilleInterpolator SplineInterpolator
     AkimaSplineInterpolator UnivariatePeriodicInterpolator]
    [org.apache.commons.math3.linear
     RealMatrix RealVector ConjugateGradient SymmLQ
     PreconditionedIterativeLinearSolver RealLinearOperator]))

;;;TODO:
;;;CREATE FUNCTIONS FROM THESE APACHE PACKAGES
;;;org.apache.commons.math3.genetics (GeneticAlgorithm)
;;;org.apache.commons.math3.ml.clustering (Clustering algorithms)
;;;org.apache.commons.math3.ode + subcategories (ODE)

;;max-dim-length for generators
(def mdl 6)

(s/def ::max-iter
  (s/with-gen (s/nilable ::m/int+)
              #(gen/one-of [(s/gen (s/int-in 100 1000))
                            (gen/return nil)])))

(s/def ::max-evaluations
  (s/with-gen ::m/int+
              #(s/gen (s/int-in 100 1000))))

(s/def ::goal #{:min :max})
(s/def ::value ::m/number)
(s/def ::point ::m/number)
(s/def ::initial-guess ::m/finite)
(s/def ::vector-point ::vector/vector)
(s/def ::weighted-constraint-errors ::vector/vector)
(s/def ::finite-interval ::intervals/finite-interval)
(s/def ::strict-finite-interval (intervals/strict-interval-spec ::m/finite))

(s/def ::weighted-constraint-errors-and-vector-point
  (s/keys :req [::weighted-constraint-errors ::vector-point]))

(s/def ::value-and-point (s/keys :req [::value ::point]))
(s/def ::value-and-vector-point (s/keys :req [::vector-point ::value]))
(s/def ::var-lower-bounds ::vector/vector-finite)
(s/def ::var-upper-bounds ::vector/vector-finite)
(s/def ::check-by-objective? boolean?)

(s/def ::vars-guess
  (s/and ::vector/vector-finite
         (fn [v]
           (pos? (count v)))))

;;note that accu's must be >= 9e-16 (or at least 3e-16),
;; otherwise Apache thinks it's zero
(def accu-precision 9e-16)

(s/def ::rel-accu
  (s/with-gen (s/and ::m/finite+
                     (fn [x]
                       (>= x accu-precision)))
              #(s/gen (s/double-in :min accu-precision
                                   :max 1.0
                                   :NaN? false))))

(s/def ::abs-accu
  (s/with-gen (s/and ::m/finite+
                     (fn [x]
                       (>= x accu-precision)))
              #(s/gen (s/double-in :min accu-precision
                                   :max 1.0
                                   :NaN? false))))

(s/def ::number->number
  (s/fspec :args (s/cat :a ::m/number)
           :ret ::m/number))

(s/def ::array->number
  (s/fspec :args (s/cat :da ::arrays/double-array)
           :ret ::m/number))

(s/def ::array->vector
  (s/fspec :args (s/cat :da ::arrays/double-array)
           :ret ::vector/vector))

(s/def ::array->nilable-matrix
  (s/fspec :args (s/cat :da ::arrays/double-array)
           :ret (s/nilable ::mx/matrix)))

(s/def ::univariate-f ::number->number)
(s/def ::objective ::array->number)
(s/def ::gradient ::array->vector)

;;;APACHE FUNCTIONS
(defn- ^UnivariateFunction univariate-function
  "`f` is a function that takes and returns a number."
  [f]
  (reify UnivariateFunction (value [_ x] (f x))))

(defn- ^MultivariateFunction multivariate-function
  "`f` is a function that takes a vector and returns a number."
  [f]
  (reify MultivariateFunction (value [_ v] (f v))))

(defn- ^MultivariateVectorFunction multivariate-vector-function
  "`f` is a function that takes and returns a vector."
  [f]
  (reify MultivariateVectorFunction (value [_ v] (double-array (f v)))))

(defn- ^MultivariateMatrixFunction multivariate-matrix-function
  "`f` is a function that takes a vector and returns a matrix."
  [f]
  (reify MultivariateMatrixFunction
    (value [_ v] (arrays/array2D :double (f v)))))

(defn- ^UnivariateDifferentiableFunction univariate-differentiable-function
  "`deriv-fn` is a function that takes and returns a number."
  ([deriv-fn ^long points ^double step-size]
   (.differentiate
     (FiniteDifferencesDifferentiator. points step-size)
     (univariate-function deriv-fn)))
  ([deriv-fn points step-size var-low-bound var-high-bound]
   (.differentiate
     (FiniteDifferencesDifferentiator.
       points step-size var-low-bound var-high-bound)
     (univariate-function deriv-fn))))

;;;OPTIMIZATION HELPERS
(defn- ->value-and-point-map
  [^UnivariatePointValuePair pv]
  {::value (.getValue pv)
   ::point (.getPoint pv)})

(defn- ->value-and-vector-point-map
  [^PointValuePair pv]
  {::value        (.getValue pv)
   ::vector-point (vec (.getPoint pv))})

(defn- ->weighted-constraint-errors-and-vector-point-map
  [^LeastSquaresProblem$Evaluation e]
  {::vector-point               (apache-v/apache-vector->vector
                                  (.getPoint e))
   ::weighted-constraint-errors (apache-v/apache-vector->vector
                                  (.getResiduals e))})

(defn- goal-fn
  [goal]
  (if (= goal :min)
    (GoalType/MINIMIZE)
    (GoalType/MAXIMIZE)))

(defn- checker-fn
  [check-by-objective? rel-accu abs-accu]
  (if check-by-objective?
    (SimpleValueChecker. rel-accu abs-accu)
    (SimplePointChecker. rel-accu abs-accu)))

(defn- vector-checker-fn
  [check-by-objective? rel-accu abs-accu]
  (if check-by-objective?
    (SimpleVectorValueChecker. rel-accu abs-accu)
    (SimplePointChecker. rel-accu abs-accu)))

;;;ROOT SOLVERS
(s/def ::univariate-f-with-interval-and-guess
  (s/with-gen
    (s/and (s/keys :req [::univariate-f ::finite-interval ::initial-guess])
           (fn [{::keys [finite-interval initial-guess]}]
             (intervals/in-interval? finite-interval initial-guess)))
    #(gen/one-of
       (map
         gen/return
         (list {::univariate-f    identity
                ::finite-interval [-5.0 5.0]
                ::initial-guess   3.0}
               {::univariate-f    (fn [v]
                                    (- (m/cube v) (* 3 v)))
                ::finite-interval [-50.0 50.0]
                ::initial-guess   3.0}
               {::univariate-f    (fn [v]
                                    (- (m/exp v) (* 5 v)))
                ::finite-interval [-50.0 50.0]
                ::initial-guess   3.0})))))

(s/def ::root-solver-type
  #{:bisection :bracketing-nth-order-brent :brent :illinois :muller :muller2
    :pegasus :regula-falsi :ridders :secant :newton-raphson})

(defn root-solver
  "`root-solver-type` options:
  :bisection, :bracketing-nth-order-brent, :brent, :illinois, :muller, :muller2,
  :newton-raphson, :pegasus, :regula-falsi, :ridders, :secant.
  `univariate-f` should return NaN when out of range."
  ([args] (root-solver args {}))
  ([{::keys [univariate-f finite-interval initial-guess]}
    {::keys [max-iter root-solver-type rel-accu abs-accu]
     :or    {root-solver-type :brent
             rel-accu         1e-14
             abs-accu         1e-6}}]
   (let [max-iter (or max-iter 1000)
         [lower upper] (mapv double finite-interval)]
     (if (= root-solver-type :newton-raphson)
       (try (.solve (NewtonRaphsonSolver. abs-accu)
                    max-iter
                    (univariate-differentiable-function univariate-f 2 0.25)
                    lower
                    upper)
            (catch Exception e
              {::anomalies/message  (str (.getMessage e))
               ::anomalies/fn       (var root-solver)
               ::anomalies/category ::anomalies/third-party}))
       (let [^BaseUnivariateSolver s
             (case root-solver-type
               :bisection (BisectionSolver. rel-accu abs-accu)
               :bracketing-nth-order-brent (BracketingNthOrderBrentSolver.
                                             rel-accu abs-accu 5)
               :brent (BrentSolver. rel-accu abs-accu)
               :illinois (IllinoisSolver. rel-accu abs-accu)
               :muller (MullerSolver. rel-accu abs-accu)
               :muller2 (MullerSolver2. rel-accu abs-accu)
               :pegasus (PegasusSolver. rel-accu abs-accu)
               :regula-falsi (RegulaFalsiSolver. rel-accu abs-accu)
               :ridders (RiddersSolver. rel-accu abs-accu)
               :secant (SecantSolver. rel-accu abs-accu)
               nil)
             uni-fn (univariate-function univariate-f)]
         (when s (try (.solve s max-iter uni-fn lower upper initial-guess)
                      (catch Exception e
                        {::anomalies/message  (str (.getMessage e))
                         ::anomalies/fn       (var root-solver)
                         ::anomalies/category ::anomalies/third-party}))))))))

(s/fdef root-solver
  :args (s/cat :univariate-f-with-interval-and-guess
               ::univariate-f-with-interval-and-guess
               :opts (s/? (s/keys :opt [::max-iter
                                        ::root-solver-type
                                        ::rel-accu
                                        ::abs-accu])))
  :ret (s/nilable (s/or :finite ::m/finite
                        :anomaly ::anomalies/anomaly)))

;;;UNIVARIATE OPTIMIZE
(defn optimize-univariate
  "Brent Optimizer.  Search over a `strict-finite-interval`.
  `initial-guess` must be in `strict-finite-interval`."
  ([univariate-f strict-finite-interval initial-guess]
   (optimize-univariate univariate-f strict-finite-interval initial-guess {}))
  ([univariate-f [low high] initial-guess
    {::keys [max-iter goal rel-accu abs-accu]
     :or    {goal     :min
             rel-accu 1e-14
             abs-accu 1e-6}}]
   (let [max-iter (or max-iter 1000)
         data [(SearchInterval. low high initial-guess)
               (MaxEval. max-iter)
               (UnivariateObjectiveFunction. (univariate-function univariate-f))
               (goal-fn goal)]
         s (BrentOptimizer. rel-accu abs-accu)]
     (try
       (let [point-value-pair (.optimize s (into-array OptimizationData data))]
         (->value-and-point-map point-value-pair))
       (catch TooManyEvaluationsException _
         {::anomalies/message  (format "Max iterations (%d) exceeded" max-iter)
          ::anomalies/fn       (var optimize-univariate)
          ::anomalies/category ::anomalies/third-party})
       (catch Exception e
         {::anomalies/message  (str (.getMessage e))
          ::anomalies/fn       (var optimize-univariate)
          ::anomalies/category ::anomalies/third-party})))))

(s/fdef optimize-univariate
  :args (s/and (s/cat :univariate-f ::univariate-f
                      :strict-finite-interval ::strict-finite-interval
                      :initial-guess ::initial-guess
                      :opts (s/? (s/keys :opts [::max-iter
                                                ::goal
                                                ::rel-accu
                                                ::abs-accu])))
               (fn [{:keys [strict-finite-interval initial-guess]}]
                 (intervals/in-bounds?
                   (intervals/bounds strict-finite-interval false false)
                   initial-guess)))
  :ret (s/or :solution ::value-and-point
             :anomaly ::anomalies/anomaly))

;;;LINEAR PROGRAMMING
(defn- linear-constraint->apache-linear-constraint
  [[linear-coefficients relation constraint-value]]
  (let [r (case relation
            :eq (Relationship/EQ)
            :leq (Relationship/LEQ)
            :geq (Relationship/GEQ))]
    (LinearConstraint. (double-array linear-coefficients)
                       ^Relationship r
                       (double constraint-value))))

(defn- add-linear-constraint
  [array-list constraint]
  (doto ^ArrayList array-list
    (.add (linear-constraint->apache-linear-constraint constraint))))

(defn linear-programming
  "`objective-coefficients` should be a vector. `linear-constraints` should be
  a collection of triples of [linear-coefficients relation constraint-value],
  where linear-coefficients is a vector; relation is either `:eq`, `:leq`, or
  `:geq`; and constraint-value is a number. When '::non-negative-vars?` is set
  to true, all variables will be non-negative."
  ([objective-coefficients linear-constraints]
   (linear-programming objective-coefficients linear-constraints {}))
  ([objective-coefficients linear-constraints
    {::keys [goal objective-constant non-negative-vars?]
     :or    {goal               :min
             objective-constant 0.0
             non-negative-vars? false}}]
   (let [array-list (ArrayList.)
         _ (doseq [constraint linear-constraints]
             (add-linear-constraint array-list constraint))
         data [(LinearObjectiveFunction.
                 (double-array objective-coefficients)
                 (double objective-constant))
               (LinearConstraintSet. array-list)
               (goal-fn goal)
               (NonNegativeConstraint. non-negative-vars?)]]
     (try
       (let [pv (.optimize (SimplexSolver.)
                           (into-array OptimizationData data))]
         (->value-and-vector-point-map pv))
       (catch Exception e
         {::anomalies/message  (str (.getMessage e))
          ::anomalies/fn       (var linear-programming)
          ::anomalies/category ::anomalies/third-party})))))

(s/def ::linear-coefficients ::vector/vector-finite)
(s/def ::relation #{:eq :leq :geq})
(s/def ::constraint-value ::m/finite)

(s/def ::linear-constraint
  (s/tuple ::linear-coefficients ::relation ::constraint-value))

(s/def ::linear-constraints
  (s/with-gen
    (s/and (s/coll-of ::linear-constraint)
           (fn [lcs]
             (every? #(= (count (ffirst lcs)) (count (first %)))
                     (rest lcs))))
    #(gen/bind (gen/large-integer* {:min 1 :max mdl})
               (fn [i]
                 (gen/vector
                   (gen/tuple (gen/vector (s/gen ::m/finite) i)
                              (s/gen ::relation)
                              (s/gen ::constraint-value))
                   1
                   mdl)))))

(s/def ::objective-constant ::m/finite)
(s/def ::non-negative-vars? boolean?)

(s/fdef linear-programming
  :args (s/and (s/cat :objective-coefficients ::vector/vector-finite
                      :linear-constraints ::linear-constraints
                      :opts (s/? (s/keys :opts [::goal
                                                ::objective-constant
                                                ::non-negative-vars?])))
               (fn [{:keys [objective-coefficients linear-constraints]}]
                 (= (count objective-coefficients)
                    (count (ffirst linear-constraints)))))
  :ret (s/or :solution ::value-and-vector-point
             :anomaly ::anomalies/anomaly))

;;;ITERATIVE LINEAR LEAST SQUARES
(defn iterative-linear-least-squares
  "Normally, use an LLS algorithm from a matrix library instead of this one.
  `symmetric-m` × y = `v`. Returns the vector, 'y'.

  A default stopping criterion is implemented. The iterations stop when
  || r || ≤ `rel-accu` || `v` ||, where `v` is the right-hand side vector,
  where 'r' is the current estimate of the residual, and `rel-accu` a
  user-specified tolerance. It should be noted that 'r' is the so-called updated
  residual, which might differ from the true residual due to rounding-off errors
  (see e.g. Strakos and Tichy, 2002).

  By default, uses the SYMMLQ algorithm. Alternatively, set
  `::use-conjugate-gradient?` to true. The implementation of the SYMMLQ
  iterative linear solver was proposed by Paige and Saunders (1975). This
  implementation is largely based on the FORTRAN code by
  Pr. Michael A. Saunders. SYMMLQ is designed to solve the system of linear
  equations A × y = b where 'A' is an n × n self-adjoint linear operator
  (defined as a RealLinearOperator), and 'b' is a given vector. The operator 'A'
  is not required to be positive definite. If 'A' is known to be definite, the
  method of conjugate gradients might be preferred, since it will require about
  the same number of iterations as SYMMLQ but slightly less work per iteration.
  SYMMLQ is designed to solve the system (A - shift × I) × y = b, where 'shift'
  is a specified scalar value. If 'shift' and 'b' are suitably chosen, the
  computed vector 'y' may approximate an (unnormalized) eigenvector of 'A', as
  in the methods of inverse iteration and/or Rayleigh-quotient iteration. Again,
  the linear operator (A - shift × I) need not be positive definite (but must be
  self-adjoint). The work per iteration is very slightly less if 'shift' = 0."
  ([symmetric-m v] (iterative-linear-least-squares symmetric-m v {}))
  ([symmetric-m v {::keys [max-iter rel-accu use-conjugate-gradient?
                           vector-initial-guess check-for-positive-definiteness?]
                   :or    {rel-accu                         m/dbl-close
                           use-conjugate-gradient?          false
                           check-for-positive-definiteness? true}}]
   (if (m/one? (count v))
     [(m/div (first v) (ffirst symmetric-m))]
     (let [max-iter (or max-iter 10000)
           ^RealMatrix a (apache-mx/apache-matrix symmetric-m)
           ^RealVector b (apache-v/apache-vector v)
           ^RealVector g (when vector-initial-guess
                           (apache-v/apache-vector vector-initial-guess))
           ^PreconditionedIterativeLinearSolver s
           (if use-conjugate-gradient?
             (ConjugateGradient. (long max-iter)
                                 (double rel-accu)
                                 (boolean check-for-positive-definiteness?))
             (SymmLQ. (long max-iter)
                      (double rel-accu)
                      (boolean check-for-positive-definiteness?)))]
       (try (vec (.toArray (if g
                             ^RealVector (.solve s
                                                 ^RealLinearOperator a
                                                 b
                                                 g)
                             ^RealVector (.solve s a b))))
            (catch Exception e
              {::anomalies/message  (str (.getMessage e))
               ::anomalies/fn       (var iterative-linear-least-squares)
               ::anomalies/category ::anomalies/third-party}))))))

(s/def ::use-conjugate-gradient? boolean?)
(s/def ::vector-initial-guess (s/nilable ::vector/vector-finite))
(s/def ::check-for-positive-definiteness? boolean?)

(s/fdef iterative-linear-least-squares
  :args (s/and (s/cat :symmetric-m ::mx/symmetric-matrix
                      :v ::vector/vector-finite
                      :opts (s/? (s/keys :opt
                                         [::max-iter
                                          ::rel-accu
                                          ::use-conjugate-gradient?
                                          ::vector-initial-guess
                                          ::check-for-positive-definiteness?])))
               (fn [{:keys [symmetric-m v]}]
                 (and (= (mx/rows symmetric-m) (count v))
                      (pos? (count v)))))
  :ret (s/or :solution ::vector/vector
             :anomaly ::anomalies/anomaly))

;;;NONLINEAR LEAST SQUARES and SYSTEMS
(defn nonlinear-least-squares
  "`::constraints-fn` takes and returns a vector. Each constraint function
  should return Inf+ or Inf- when out of range.

  `::constraint-jacobian-fn` is the jacobian matrix function that takes a
  vector and returns a matrix. `::constraint-jacobian-fn` should return nil when
  out of range.

  `::max-iter` - the maximum number of times to iterate in the algorithm.
  `::max-evaluations` - the maximum number of times to evaluate the model.

  Solver is Levenberg-Marquardt by default. Alternatively, set
  `::use-gauss-newton?` to true.

  `::constraint-weights` - the importance weights of each constraint. Each
  constraint error is multiplied by the square root of the weight.

  Returns map of `::vector-point` and `::weighted-constraint-errors`."
  ([args] (nonlinear-least-squares args {}))
  ([{::keys [constraints-fn constraint-jacobian-fn vars-guess]}
    {::keys [max-iter rel-accu abs-accu max-evaluations
             check-by-objective? use-gauss-newton? constraint-weights]
     :or    {rel-accu            1e-14
             abs-accu            1e-6
             max-evaluations     1000
             check-by-objective? false
             use-gauss-newton?   false}}]
   (let [max-iter (or max-iter 1000)
         constraint-count (count (constraints-fn vars-guess))
         c (multivariate-vector-function constraints-fn)
         j (multivariate-matrix-function constraint-jacobian-fn)
         s (if use-gauss-newton?
             (GaussNewtonOptimizer.)
             (LevenbergMarquardtOptimizer.))
         checker (LeastSquaresFactory/evaluationChecker
                   (vector-checker-fn check-by-objective? rel-accu abs-accu))
         observed (apache-v/apache-vector (vec (repeat constraint-count 0.0)))
         start (apache-v/apache-vector vars-guess)
         weights (if (and (some? constraint-weights)
                          (= (count constraint-weights) constraint-count))
                   (apache-mx/apache-matrix
                     (mx/diagonal-matrix constraint-weights))
                   (apache-mx/apache-matrix
                     (mx/diagonal-matrix (repeat constraint-count 1.0))))]
     (try
       (when s
         (let [multivariate-jacobian-fn (LeastSquaresFactory/model c j)
               problem (LeastSquaresFactory/create
                         multivariate-jacobian-fn observed start
                         weights checker max-evaluations max-iter)
               e (.optimize s problem)]
           (->weighted-constraint-errors-and-vector-point-map e)))
       (catch TooManyEvaluationsException _
         {::anomalies/message  (format "Max evaluations (%d) exceeded."
                                       max-evaluations)
          ::anomalies/fn       (var nonlinear-least-squares)
          ::anomalies/category ::anomalies/third-party})
       (catch TooManyIterationsException _
         {::anomalies/message  (format "Max iterations (%d) exceeded."
                                       max-iter)
          ::anomalies/fn       (var nonlinear-least-squares)
          ::anomalies/category ::anomalies/third-party})
       (catch Exception e
         {::anomalies/message  (str (.getMessage e))
          ::anomalies/fn       (var nonlinear-least-squares)
          ::anomalies/category ::anomalies/third-party})))))

(s/def ::constraint-jacobian-fn ::array->nilable-matrix)
(s/def ::constraints-fn ::array->vector)

(s/def ::constraints-with-jacobian-and-guess
  (s/with-gen
    (s/and (s/keys :req [::constraints-fn
                         ::constraint-jacobian-fn
                         ::vars-guess])
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

(s/def ::use-gauss-newton? boolean?)
(s/def ::constraint-weights (s/nilable ::vector/vector-finite+))

(s/fdef nonlinear-least-squares
  :args (s/cat :constraints-with-jacobian-and-guess
               ::constraints-with-jacobian-and-guess
               :opts (s/? (s/keys :opt [::max-iter
                                        ::rel-accu
                                        ::abs-accu
                                        ::max-evaluations
                                        ::check-by-objective?
                                        ::use-gauss-newton?
                                        ::constraint-weights])))
  :ret (s/or :solution ::weighted-constraint-errors-and-vector-point
             :anomaly ::anomalies/anomaly))

;;;NLP
;;OPTIMIZE WITHOUT CONSTRAINTS AND WITH GRADIENT
(defn optimize-without-constraints-and-with-gradient
  "The Conjugate Gradient algorithm is best when you can calculate the
  `gradient` of the `objective`. By default, updating is done using
  Fletcher-Reeves. Alternatively, set `update-using-polak-ribiere?` to true."
  ([args] (optimize-without-constraints-and-with-gradient args {}))
  ([{::keys [objective gradient vars-guess]}
    {::keys [max-iter goal rel-accu abs-accu check-by-objective?
             update-using-polak-ribiere? initial-step-size]
     :or    {goal                        :min
             rel-accu                    1e-14
             abs-accu                    1e-6
             check-by-objective?         false
             update-using-polak-ribiere? false
             initial-step-size           1.0}}]
   (let [max-iter (or max-iter 1000)
         updater (if update-using-polak-ribiere?
                   NonLinearConjugateGradientOptimizer$Formula/POLAK_RIBIERE
                   NonLinearConjugateGradientOptimizer$Formula/FLETCHER_REEVES)
         checker (checker-fn check-by-objective? rel-accu abs-accu)
         data [(NonLinearConjugateGradientOptimizer$BracketingStep.
                 initial-step-size)
               (MaxEval. max-iter)
               (ObjectiveFunction. (multivariate-function objective))
               (ObjectiveFunctionGradient.
                 (multivariate-vector-function gradient))
               (goal-fn goal)
               (InitialGuess. (double-array vars-guess))]]
     (try
       (let [pv (.optimize
                  (NonLinearConjugateGradientOptimizer. updater checker)
                  (into-array OptimizationData data))]
         (->value-and-vector-point-map pv))
       (catch TooManyEvaluationsException _
         {::anomalies/message  (format "Max iterations (%d) exceeded." max-iter)
          ::anomalies/fn       (var
                                 optimize-without-constraints-and-with-gradient)
          ::anomalies/category ::anomalies/third-party})
       (catch Exception e
         {::anomalies/message  (str (.getMessage e))
          ::anomalies/fn       (var
                                 optimize-without-constraints-and-with-gradient)
          ::anomalies/category ::anomalies/third-party})))))

(s/def ::update-using-polak-ribiere? boolean?)

(s/def ::initial-step-size
  (s/with-gen ::m/finite+
              #(s/gen (s/double-in :infinite? false
                                   :NaN? false
                                   :min 1e-2
                                   :max 100.0))))

(s/def ::objective-with-gradient-and-guess
  (s/with-gen
    (s/keys :req [::objective ::gradient ::vars-guess])
    #(gen/one-of
       (map gen/return
            (list {::objective  (fn [da]
                                  (if (= 2 (count da))
                                    (let [[a1 a2] da]
                                      (+ a1 (m/sq a2)))
                                    m/nan))
                   ::gradient   (fn [da]
                                  (if (= 2 (count da))
                                    (let [[a1 a2] da]
                                      [1.0 (* 2 a2)])
                                    [m/nan m/nan]))
                   ::vars-guess [2.0 -2.0]}
                  {::objective  (fn [da]
                                  (if (= 2 (count da))
                                    (let [[a1 a2] da]
                                      (+ a1 (m/exp a2)))
                                    m/nan))
                   ::gradient   (fn [da]
                                  (if (= 2 (count da))
                                    (let [[a1 a2] da]
                                      [1.0 (m/exp a2)])
                                    [m/nan m/nan]))
                   ::vars-guess [1.9 -1.9]})))))

(s/fdef optimize-without-constraints-and-with-gradient
  :args (s/cat :objective-with-gradient-and-guess
               ::objective-with-gradient-and-guess
               :opts (s/? (s/keys :opt [::max-iter
                                        ::goal
                                        ::rel-accu
                                        ::abs-accu
                                        ::check-by-objective?
                                        ::update-using-polak-ribiere?
                                        ::initial-step-size])))
  :ret (s/or :solution ::value-and-vector-point
             :anomaly ::anomalies/anomaly))

;;OPTIMIZE NO-DERIVATIVE NO-CONSTRAINTS
(defn optimize-without-constraints-or-gradient
  "`::without-constraints-or-gradient-nlp-solver-type` options: :powell,
  :nelder-mead, :multi-directional-simplex. These algorithms may work well when
  you can't or don't want to calculate the gradient of the `objective`
  function."
  ([objective vars-guess]
   (optimize-without-constraints-or-gradient objective vars-guess {}))
  ([objective vars-guess
    {::keys [max-iter goal rel-accu abs-accu
             without-constraints-or-gradient-nlp-solver-type]
     :or    {goal                                            :min
             without-constraints-or-gradient-nlp-solver-type :multi-directional-simplex
             rel-accu                                        1e-14
             abs-accu                                        1e-6}}]
   (let [max-iter (or max-iter 1000)
         dim-count (count vars-guess)
         initial (double-array vars-guess)
         data [(MaxEval. max-iter)
               (ObjectiveFunction. (multivariate-function objective))
               (goal-fn goal)
               (InitialGuess. initial)]
         data (case without-constraints-or-gradient-nlp-solver-type
                :powell data
                :nelder-mead (conj data
                                   (doto (NelderMeadSimplex. dim-count)
                                     (.build initial)))
                :multi-directional-simplex (conj
                                             data
                                             (doto (MultiDirectionalSimplex.
                                                     dim-count)
                                               (.build initial))))
         s (if (= without-constraints-or-gradient-nlp-solver-type :powell)
             (PowellOptimizer. rel-accu abs-accu)
             (SimplexOptimizer. rel-accu abs-accu))]
     (try
       (let [pv (.optimize s (into-array OptimizationData data))]
         (->value-and-vector-point-map pv))
       (catch TooManyEvaluationsException _
         {::anomalies/message  (format "Max iterations (%d) exceeded." max-iter)
          ::anomalies/fn       (var optimize-without-constraints-or-gradient)
          ::anomalies/category ::anomalies/third-party})
       (catch Exception e
         {::anomalies/message  (str (.getMessage e))
          ::anomalies/fn       (var optimize-without-constraints-or-gradient)
          ::anomalies/category ::anomalies/third-party})))))

(s/def ::without-constraints-or-gradient-nlp-solver-type
  #{:powell :nelder-mead :multi-directional-simplex})

(s/fdef optimize-without-constraints-or-gradient
  :args (s/cat :objective ::objective
               :vars-guess ::vars-guess
               :opts (s/? (s/keys :opt [::max-iter
                                        ::goal
                                        ::rel-accu
                                        ::abs-accu
                                        ::without-constraints-or-gradient-nlp-solver-type])))
  :ret (s/or :solution ::value-and-vector-point
             :anomaly ::anomalies/anomaly))

;;OPTIMIZE NO-DERIVATIVE WITH SIMPLE BOUNDS ONLY
(defn optimize-cma-evolution$
  "An implementation of the active Covariance Matrix Adaptation Evolution
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
  than BFGS (counting objective function evaluations, no gradients provided).
  For up to N=10 variables, the derivative-free simplex direct search method
  (Nelder and Mead) can be faster, but it is far less reliable than CMA-ES. The
  CMA-ES is particularly well suited for non-separable and/or badly conditioned
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

  `::check-by-objective?` - Whether convergence checker uses objective function
  or point value for calculating error."
  ([objective var-lower-bounds var-upper-bounds vars-guess]
   (optimize-cma-evolution$
     objective var-lower-bounds var-upper-bounds vars-guess {}))
  ([objective var-lower-bounds var-upper-bounds vars-guess
    {::keys [max-iter goal rel-accu abs-accu check-by-objective? sigmas
             population-size stop-fitness active-cma? diagonal-only
             check-feasible-count]
     :or    {goal                 :min
             rel-accu             1e-14
             abs-accu             1e-6
             check-by-objective?  false
             stop-fitness         0.0
             active-cma?          true
             diagonal-only        0
             check-feasible-count 0}}]
   (let [dim-count (count vars-guess)
         p (or population-size
               (+ 4 (m/floor' (* 3 (m/log dim-count)))))
         iter (or max-iter
                  (m/floor' (* 1000
                               (m/sq (+ 5 dim-count))
                               (m/pow p -0.5))))
         sigma (if sigmas
                 (double-array sigmas)
                 (double-array dim-count 0.3))
         data [(MaxEval. iter)
               (ObjectiveFunction. (multivariate-function objective))
               (goal-fn goal)
               (InitialGuess. (double-array vars-guess))
               (CMAESOptimizer$Sigma. sigma)
               (CMAESOptimizer$PopulationSize. p)
               (SimpleBounds. (double-array var-lower-bounds)
                              (double-array var-upper-bounds))]
         checker (checker-fn check-by-objective? rel-accu abs-accu)
         ;for debugging on `s` below
         ; (.getStatisticsDHistory, .getStatisticsFitnessHistory,
         ;  .getStatisticsMeanHistory, .getStatisticsSigmaHistory)
         generate-stats? false
         s (CMAESOptimizer. iter stop-fitness active-cma? diagonal-only
                            check-feasible-count (alt-random/mersenne-rng
                                                   (random/rnd-long!))
                            generate-stats? checker)]
     (try
       ;:value has to be calculated below because it's
       ; not consistent with :point values otherwise
       (let [pv (.optimize s (into-array OptimizationData data))
             pt (::vector-point (->value-and-vector-point-map pv))]
         {::vector-point pt
          ::value        (objective (double-array pt))})
       (catch TooManyEvaluationsException _
         {::anomalies/message  (format "Max iterations (%d) exceeded." max-iter)
          ::anomalies/fn       (var optimize-cma-evolution$)
          ::anomalies/category ::anomalies/third-party})
       (catch Exception e
         {::anomalies/message  (str (.getMessage e))
          ::anomalies/fn       (var optimize-cma-evolution$)
          ::anomalies/category ::anomalies/third-party})))))

(s/def ::sigmas
  (s/with-gen (s/nilable ::vector/vector-finite+)
              #(gen/one-of
                 [(gen/vector
                    (s/gen
                      (s/double-in :infinite? false
                                   :NaN? false
                                   :min 1e-2
                                   :max 100.0))
                    1
                    mdl)
                  (gen/return nil)])))

(s/def ::population-size
  (s/with-gen (s/nilable ::m/int+)
              #(gen/one-of
                 [(s/gen (s/int-in 4 10))
                  (gen/return nil)])))

(s/def ::stop-fitness
  (s/with-gen ::m/non-
              #(s/gen (s/double-in :infinite? false
                                   :NaN? false
                                   :min 1e-2
                                   :max 100.0))))

(s/def ::active-cma? boolean?)

(s/def ::diagonal-only
  (s/with-gen ::m/int-non-
              #(s/gen (s/int-in 0 1000))))

(s/def ::check-feasible-count
  (s/with-gen ::m/int-non-
              #(s/gen (s/int-in 0 1000))))

(s/fdef optimize-cma-evolution$
  :args (s/cat :objective ::objective
               :var-lower-bounds ::var-lower-bounds
               :var-upper-bounds ::var-upper-bounds
               :vars-guess ::vars-guess
               :opts (s/? (s/keys :opt [::max-iter
                                        ::goal
                                        ::rel-accu
                                        ::abs-accu
                                        ::check-by-objective?
                                        ::sigmas
                                        ::population-size
                                        ::stop-fitness
                                        ::active-cma?
                                        ::diagonal-only
                                        ::check-feasible-count])))
  :ret (s/or :solution ::value-and-vector-point
             :anomaly ::anomalies/anomaly))

(defn optimize-bobyqa
  "Powell's BOBYQA algorithm (Bound Optimization BY Quadratic Approximation).
  Faster but less robust than [[optimize-cma-evolution]]. BOBYQA is particularly
  well suited for high dimensional problems where derivatives are not available.
  In most cases it outperforms the PowellOptimizer significantly. BOBYQA could
  also be considered as a replacement of any derivative-based optimizer when the
  derivatives are approximated by finite differences. Recommend setting the
  `::bobyqa-interpolation-points` in [n+2, (n+1)(n+2)/2], where 'n' is the
  number of dimensions."
  ([objective var-lower-bounds var-upper-bounds vars-guess]
   (optimize-bobyqa objective var-lower-bounds var-upper-bounds vars-guess {}))
  ([objective var-lower-bounds var-upper-bounds vars-guess
    {::keys [max-iter goal bobyqa-interpolation-points]
     :or    {goal :min}}]
   (let [max-iter (or max-iter 1000)
         data [(MaxEval. max-iter)
               (ObjectiveFunction. (multivariate-function objective))
               (goal-fn goal)
               (InitialGuess. (double-array vars-guess))
               (SimpleBounds. (double-array var-lower-bounds)
                              (double-array var-upper-bounds))]
         p (or bobyqa-interpolation-points
               (m/ceil' (* 1.5 (inc (count vars-guess)))))
         s (BOBYQAOptimizer. p)]
     (try
       (let [pv (.optimize s (into-array OptimizationData data))]
         (->value-and-vector-point-map pv))
       (catch TooManyEvaluationsException _
         {::anomalies/message  (format "Max iterations (%d) exceeded." max-iter)
          ::anomalies/fn       (var optimize-bobyqa)
          ::anomalies/category ::anomalies/third-party})
       (catch Exception e
         {::anomalies/message  (str (.getMessage e))
          ::anomalies/fn       (var optimize-bobyqa)
          ::anomalies/category ::anomalies/third-party})))))

(s/def ::bobyqa-interpolation-points
  (s/with-gen
    (s/nilable ::m/int+)
    #(gen/one-of
       [(s/gen (s/int-in 3 8))
        (gen/return nil)])))

(s/fdef optimize-bobyqa
  :args (s/cat :objective ::objective
               :var-lower-bounds ::var-lower-bounds
               :var-upper-bounds ::var-upper-bounds
               :vars-guess ::vars-guess
               :opts (s/? (s/keys :opt [::max-iter
                                        ::goal
                                        ::bobyqa-interpolation-points])))
  :ret (s/or :solution ::value-and-vector-point
             :anomaly ::anomalies/anomaly))

;;;INTERPOLATION
(s/def ::strictly-ascending-vector-finite
  (s/with-gen
    (s/and (s/coll-of ::m/finite
                      :kind clojure.core/vector?
                      :into []
                      :min-count 2)
           (fn [v]
             (let [dv (map double v)]
               (= dv (sort (distinct dv))))))
    #(gen/bind
       (gen/vector-distinct
         (s/gen ::m/finite)
         {:min-elements 2
          :max-elements 5})
       (fn [v]
         (gen/return (vec (sort v)))))))

(s/def ::x-vals ::strictly-ascending-vector-finite)
(s/def ::y-vals ::strictly-ascending-vector-finite)
(s/def ::z-vals ::strictly-ascending-vector-finite)

(s/def ::f-vals
  (s/with-gen
    (s/coll-of ::m/finite
               :kind clojure.core/vector?
               :into []
               :min-count 2)
    #(gen/vector (s/gen ::m/finite) 2 5)))

(defn interpolation-1D
  "Implements an algorithm for interpolation of real univariate functions.
  `::x-vals` - All the x-coordinates of the interpolation points, in strictly
  ascending order.
  `::f-vals` - The values of the interpolation points on all the grid knots:
  `f-vals`[i] = f(`x-vals`[i]).

  If a positive `::period` length is entered, then the data is assumed to be
  periodic.

  `::interpolation-1D-type` can be:
    `:polynomial` - (default) Polynomial Divided Difference Algorithm. For
      reference, see Introduction to Numerical Analysis, ISBN 038795452X,
      chapter 2.
    `:neville` - Neville's Algorithm. For reference, see
      Introduction to Numerical Analysis, ISBN 038795452X, chapter 2.
    `:loess` - The Local Regression Algorithm (also Loess, Lowess). Can not
      extrapolate.
    `:cubic` - Computes a natural (also known as 'free', 'unclamped')
      cubic spline interpolation for the data set. Can not extrapolate.
    `:linear` - Linear function. Can not extrapolate.
    `:akima` - Also a cubic spline interpolation. Can not extrapolate. The Akima
      algorithm requires that the number of points >= 5.

  For the `:loess` interpolation only:
   `::bandwidth-for-loess` when computing the loess fit at a particular point,
     this fraction of source points closest to the current point is taken into
     account for computing a least-squares regression. A sensible value is
     usually 0.25 to 0.5, the default value is 0.3. For reference, see
     William S. Cleveland - Robust Locally Weighted Regression and Smoothing
     Scatterplots.

  Returns a function that accepts a point 'x' and returns the interpolated
  value."
  ([args] (interpolation-1D args {}))
  ([{::keys [x-vals f-vals]}
    {::keys [interpolation-1D-type period bandwidth-for-loess]
     :or    {interpolation-1D-type :polynomial
             period                0.0
             bandwidth-for-loess   0.3}}]
   (let [s (case interpolation-1D-type
             :polynomial (DividedDifferenceInterpolator.)
             :neville (NevilleInterpolator.)
             :loess (LoessInterpolator. bandwidth-for-loess 2)
             :cubic (SplineInterpolator.)
             :linear (LinearInterpolator.)
             :akima (AkimaSplineInterpolator.))]
     (try (let [^UnivariateFunction f (.interpolate
                                        (if (pos? period)
                                          (UnivariatePeriodicInterpolator.
                                            s period)
                                          s)
                                        (double-array x-vals)
                                        (double-array f-vals))]
            (fn [x]
              (try (.value f (double x))
                   (catch Exception e
                     {::anomalies/message  (str "Returned Function: "
                                                (.getMessage e))
                      ::anomalies/fn       (var interpolation-1D)
                      ::anomalies/category ::anomalies/third-party}))))
          (catch Exception e {::anomalies/message  (str (.getMessage e))
                              ::anomalies/fn       (var interpolation-1D)
                              ::anomalies/category ::anomalies/third-party})))))

(s/def ::interpolation-1D-type
  #{:polynomial :neville :loess :cubic :linear :akima})

(s/def ::period
  (s/with-gen ::m/finite-non-
              #(s/gen (s/double-in :infinite? false
                                   :NaN? false
                                   :min 1e-3
                                   :max 1e3))))

(s/def ::bandwidth-for-loess
  (s/with-gen ::m/open-prob
              #(s/gen (s/double-in :infinite? false
                                   :NaN? false
                                   :min 0.25
                                   :max 0.5))))

(s/def ::x-vals-with-f-vals
  (s/with-gen
    (s/and (s/keys :req [::x-vals ::f-vals])
           (fn [{::keys [x-vals f-vals]}]
             (= (count f-vals) (count x-vals))))
    #(gen/bind
       (s/gen ::x-vals)
       (fn [x]
         (gen/bind
           (gen/vector (s/gen ::m/num)
                       (count x))
           (fn [f-v]
             (gen/return {::x-vals x
                          ::f-vals f-v})))))))

(s/fdef interpolation-1D
  :args (s/cat :x-vals-with-f-vals ::x-vals-with-f-vals
               :opts (s/? (s/keys :opt [::interpolation-1D-type
                                        ::period
                                        ::bandwidth-for-loess])))
  :ret (s/or :solution (s/fspec :args (s/cat :x ::m/num)
                                :ret (s/or :inner-solution ::m/number
                                           :inner-anomaly ::anomalies/anomaly))
             :anomaly ::anomalies/anomaly))

(defn interpolation-1D-using-derivatives
  "`data` should be a collection of seqs, where each seq contains
  [x value & first derivative, the second derivative, and so ...]
  Uses Hermite Interpolator. Returns a function that accepts a point 'x' and
  returns the interpolated value."
  [data]
  (let [hi (HermiteInterpolator.)
        _ (doseq [d data]
            (.addSamplePoint hi
                             (double (first d))
                             (arrays/array2D :double (partition 1 (rest d)))))]
    (fn [x]
      (try
        (first (.value hi (double x)))
        (catch Exception e
          {::anomalies/message  (str "Returned Function: " (.getMessage e))
           ::anomalies/fn       (var interpolation-1D-using-derivatives)
           ::anomalies/category ::anomalies/third-party})))))

(s/def ::point-datum
  (s/with-gen
    (s/coll-of ::m/num :kind vector? :into [] :min-count 3)
    #(gen/vector (s/gen ::m/num) 3 5)))

(s/def ::point-data
  (s/with-gen
    (s/and (s/coll-of ::point-datum)
           (fn [p-data]
             (apply distinct? (map (comp double first) p-data))))
    #(gen/vector (s/gen ::point-datum) 3 5)))

(s/fdef interpolation-1D-using-derivatives
  :args (s/cat :data ::point-data)
  :ret (s/fspec :args (s/cat :x ::m/num)
                :ret (s/or :inner-solution ::m/number
                           :inner-anomaly ::anomalies/anomaly)))

(defn interpolation-2D
  "Generates a bicubic interpolation function. Prior to generating the
  interpolating function, the input is smoothed by default using polynomial
  fitting. Alternatively, set `::smooth?` to false. The non-smoothed version
  requires at least 5 data points.
  `::x-vals` - All the x-coordinates of the interpolation points, in strictly
  ascending order.
  `::y-vals` - All the y-coordinates of the interpolation points, in strictly
  ascending order.
  `::f-matrix` - The values of the interpolation points on all the grid knots:
   `f-matrix`[i][j] = f(`x-vals`[i], `y-vals`[j]).
  Returns a function that accepts the points 'x' and 'y' and returns the
  interpolated value."
  ([args] (interpolation-2D args {}))
  ([{::keys [x-vals y-vals f-matrix]}
    {::keys [smooth?]
     :or    {smooth? true}}]
   (if smooth?
     (try (let [^BicubicInterpolatingFunction f
                (.interpolate
                  (BicubicInterpolator.)
                  (double-array x-vals)
                  (double-array y-vals)
                  (arrays/array2D :double f-matrix))]
            (fn [x y]
              (if (.isValidPoint f x y)
                (.value f x y)
                {::anomalies/message  (format (str "Returned Function: point"
                                                   " %.5f, %.5f is not valid.")
                                              (double x)
                                              (double y))
                 ::anomalies/fn       (var interpolation-2D)
                 ::anomalies/category ::anomalies/third-party})))
          (catch Exception e
            {::anomalies/message  (str (.getMessage e))
             ::anomalies/fn       (var interpolation-2D)
             ::anomalies/category ::anomalies/third-party}))
     (try (let [^PiecewiseBicubicSplineInterpolatingFunction f
                (.interpolate
                  (PiecewiseBicubicSplineInterpolator.)
                  (double-array x-vals)
                  (double-array y-vals)
                  (arrays/array2D :double f-matrix))]
            (fn [x y]
              (if (.isValidPoint f x y)
                (.value f x y)
                {::anomalies/message  (format (str "Returned Function: point"
                                                   " %.5f, %.5f is not valid.")
                                              (double x)
                                              (double y))
                 ::anomalies/fn       (var interpolation-2D)
                 ::anomalies/category ::anomalies/third-party})))
          (catch Exception e
            {::anomalies/message  (str (.getMessage e))
             ::anomalies/fn       (var interpolation-2D)
             ::anomalies/category ::anomalies/third-party})))))

(s/def ::smooth? boolean?)
(s/def ::f-matrix ::mx/matrix-num)

(s/def ::vals-with-matrix
  (s/with-gen
    (s/and (s/keys :req [::x-vals ::y-vals ::f-matrix])
           (fn [{::keys [x-vals y-vals f-matrix]}]
             (and (= (mx/rows f-matrix) (count x-vals))
                  (= (mx/columns f-matrix) (count y-vals)))))
    #(gen/bind
       (gen/tuple (s/gen ::x-vals)
                  (s/gen ::y-vals))
       (fn [[x y]]
         (gen/bind (gen/vector
                     (gen/vector (s/gen ::m/num)
                                 (count y))
                     (count x))
                   (fn [f-m]
                     (gen/return {::x-vals   x
                                  ::y-vals   y
                                  ::f-matrix f-m})))))))

(s/fdef interpolation-2D
  :args (s/and (s/cat :vals-with-matrix ::vals-with-matrix
                      :opts (s/? (s/keys :opt [::smooth?])))
               (fn [{:keys [vals-with-matrix opts]}]
                 (let [{::keys [f-matrix]} vals-with-matrix
                       {::keys [smooth?]} opts]
                   (or smooth? (>= (tensor/ecount f-matrix) 5)))))
  :ret (s/or :solution (s/fspec :args (s/cat :x ::m/num
                                             :y ::m/num)
                                :ret (s/or :inner-solution ::m/number
                                           :inner-anomaly ::anomalies/anomaly))
             :anomaly ::anomalies/anomaly))

(defn interpolation-3D
  "`::x-vals` - All the x-coordinates of the interpolation points, in strictly
  ascending order.
  `::y-vals` - All the y-coordinates of the interpolation points, in strictly
  ascending order.
  `::z-vals` - All the z-coordinates of the interpolation points, in strictly
  ascending order.
  `::f-tensor` - the values of the interpolation points on all the grid knots:
  `f-tensor`[i][j][k] = f(`x-vals`[i], `y-vals`[j], `z-vals`[k]).
  Returns a function that accepts the points 'x', 'y', and 'z' and returns the
  interpolated value."
  [{::keys [x-vals y-vals z-vals f-tensor]}]
  (try (let [^TricubicInterpolatingFunction f
             (.interpolate
               (TricubicInterpolator.)
               (double-array x-vals)
               (double-array y-vals)
               (double-array z-vals)
               (arrays/array3D :double f-tensor))]
         (fn [x y z]
           (try (.value f x y z)
                (catch Exception e
                  {::anomalies/message  (str "Returned Function: "
                                             (.getMessage e))
                   ::anomalies/fn       (var interpolation-3D)
                   ::anomalies/category ::anomalies/third-party}))))
       (catch Exception e
         {::anomalies/message  (str (.getMessage e))
          ::anomalies/fn       (var interpolation-3D)
          ::anomalies/category ::anomalies/third-party})))

(s/def ::f-tensor
  (s/and ::tensor/tensor3D
         (fn [t]
           (every? m/num? (flatten t)))))

(s/def ::vals-with-tensor
  (s/with-gen
    (s/and (s/keys :req [::x-vals ::y-vals ::z-vals ::f-tensor])
           (fn [{::keys [x-vals y-vals z-vals f-tensor]}]
             (let [[xc yc zc] (tensor/shape f-tensor)]
               (and (= xc (count x-vals))
                    (= yc (count y-vals))
                    (= zc (count z-vals))))))
    #(gen/bind
       (gen/tuple (s/gen ::x-vals)
                  (s/gen ::y-vals)
                  (s/gen ::z-vals))
       (fn [[x y z]]
         (gen/bind (gen/vector
                     (gen/vector
                       (gen/vector (s/gen ::m/num)
                                   (count z))
                       (count y))
                     (count x))
                   (fn [f-t]
                     (gen/return {::x-vals   x
                                  ::y-vals   y
                                  ::z-vals   z
                                  ::f-tensor f-t})))))))

(s/fdef interpolation-3D
  :args (s/cat :vals-with-tensor ::vals-with-tensor)
  :ret (s/or :solution (s/fspec :args (s/cat :x ::m/num
                                             :y ::m/num
                                             :z ::m/num)
                                :ret (s/or :inner-solution ::m/number
                                           :inner-anomaly ::anomalies/anomaly))
             :anomaly ::anomalies/anomaly))

(defn interpolation-ND-microsphere$
  "Interpolator that implements the algorithm described in William Dudziak's MS
  thesis. Results are randomized. `::i-matrix` - the arguments for the
  interpolation points. `i-matrix`[i][0] is the first component of interpolation
  point 'i', `i-matrix`[i][1] is the second component, and so on until
  `i-matrix`[i][d-1], the last component of that interpolation point (where 'd'
  is thus the dimension of the space). `::f-vals` - the values for the
  interpolation points. Returns a function that takes a vector and returns a
  number."
  [{::keys [i-matrix f-vals]}]
  (let [^MultivariateFunction f (.interpolate
                                  (MicrosphereInterpolator.)
                                  (arrays/array2D :double i-matrix)
                                  (double-array f-vals))]
    (fn [v]
      (try (.value f (double-array v))
           (catch Exception e
             {::anomalies/message  (str "Returned Function: " (.getMessage e))
              ::anomalies/fn       (var interpolation-ND-microsphere$)
              ::anomalies/category ::anomalies/third-party})))))

(s/def ::i-matrix ::mx/matrix-num)

(s/def ::i-matrix-with-f-vals
  (s/with-gen
    (s/and (s/keys :req [::i-matrix ::f-vals])
           (fn [{::keys [i-matrix f-vals]}]
             (= (count f-vals) (count i-matrix))))
    #(gen/bind
       (gen/tuple (s/gen ::f-vals) (s/gen (s/int-in 1 5)))
       (fn [[f-v c]]
         (gen/bind
           (gen/vector
             (gen/vector (s/gen ::m/num)
                         c)
             (count f-v))
           (fn [i-mx]
             (gen/return {::i-matrix i-mx
                          ::f-vals   f-v})))))))

(s/fdef interpolation-ND-microsphere$
  :args (s/cat :i-matrix-with-f-vals ::i-matrix-with-f-vals)
  :ret (s/fspec :args (s/cat :v ::vector/vector-num)
                :ret (s/or :inner-solution ::m/number
                           :inner-anomaly ::anomalies/anomaly)))
