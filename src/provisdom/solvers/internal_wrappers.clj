(ns provisdom.solvers.internal-wrappers
  (:require
    [clojure.spec.alpha :as s]
    [clojure.spec.gen.alpha :as gen]
    [provisdom.apache-math.apache-matrix :as apache-mx]
    [provisdom.math.arrays :as arrays]
    [provisdom.math.core :as m]
    [provisdom.math.matrix :as mx]
    [provisdom.math.tensor :as tensor]
    [provisdom.math.vector :as vector]
    [provisdom.utility-belt.anomalies :as anomalies])
  (:import [de.xypron.jcobyla Cobyla Calcfc CobylaExitStatus]
           [com.joptimizer.optimizers OptimizationRequest
                                      JOptimizer
                                      OptimizationResponse]
           [com.joptimizer.functions PDQuadraticMultivariateRealFunction
                                     ConvexMultivariateRealFunction
                                     LinearMultivariateRealFunction]))

;;;JOPTIMIZER
;;http://www.joptimizer.com/
;;JOPTIMIZER can make a 2D version that has quad-constraints too,
;;cone constraints, linear programming, other stuff...
(s/def ::linear-objective (s/nilable ::vector/vector-finite))
(s/def ::quadratic-guess (s/nilable ::vector/vector-finite))

(s/def ::equal
  (s/and (s/nilable (s/tuple ::mx/matrix-finite ::vector/vector-finite))
         (fn [[lhs rhs]]
           (or (not lhs)
               (= (mx/rows lhs) (count rhs))))))

(s/def ::less-than
  (s/and (s/nilable (s/tuple ::mx/matrix-finite ::vector/vector-finite))
         (fn [[lhs rhs]]
           (or (not lhs)
               (= (mx/rows lhs) (count rhs))))))

(s/def ::apache-quadratic-objective ::apache-mx/positive-definite-apache-matrix-finite)

(defn quadratic-programming-joptimizer
  "Minimizes (1/2) x^T P x + q^T x over x, where 'P' is the
  `apache-quadratic-objective` matrix and 'q' is the `::linear-objective`
  vector. `apache-quadratic-objective` should be a positive (symmetrical) Apache
  matrix. `::equal` and `::less-than` should be tuples of matrices and vectors,
  where the matrix is the left-hand side and the vector is the right-hand-side.
  `::quadratic-guess` must meet constraints, and strictly meet the `less-than`
  constraints."
  ([apache-quadratic-objective]
   (quadratic-programming-joptimizer apache-quadratic-objective {}))
  ([apache-quadratic-objective
    {::keys [linear-objective equal less-than quadratic-guess]}]
   (if (apache-mx/empty-apache-matrix? apache-quadratic-objective)
     []
     (let [[eq-mx eq-v] equal
           [lt-mx lt-v] less-than
           obj (PDQuadraticMultivariateRealFunction.
                 ^"[[D" (arrays/array2D :double (apache-mx/apache-matrix->matrix apache-quadratic-objective))
                 (when linear-objective (double-array linear-objective))
                 0.0)
           equal-count (if equal
                         (count eq-v)
                         0)
           less-than-count (if less-than
                             (count lt-v)
                             0)
           req (OptimizationRequest.)
           opt (JOptimizer.)
           less-than-arr (into-array
                           ConvexMultivariateRealFunction
                           (map (fn [i]
                                  (LinearMultivariateRealFunction.
                                    (double-array (mx/get-row lt-mx i))
                                    (double (- (get lt-v i m/nan)))))
                                (range less-than-count)))]
       (.setF0 req obj)
       (when quadratic-guess (.setInitialPoint req (double-array quadratic-guess)))
       (when (pos? less-than-count) (.setFi req less-than-arr))
       (when (pos? equal-count) (.setA req ^"[[D" (arrays/array2D :double eq-mx)))
       (when (pos? equal-count) (.setB req (double-array eq-v)))
       (.setToleranceFeas req 1E-12)
       (.setTolerance req 1E-12)
       (.setOptimizationRequest opt req)
       (try (let [sol (.optimize opt)]
              (if (or (not sol) (<= sol 1))
                (vec (.getSolution
                       ^OptimizationResponse (.getOptimizationResponse opt)))
                {::anomalies/message  "Quadratic Programming Fail"
                 ::anomalies/fn       (var quadratic-programming-joptimizer)
                 ::anomalies/category ::anomalies/third-party}))
            (catch Exception e {::anomalies/message  (str (.getMessage e))
                                ::anomalies/fn       (var quadratic-programming-joptimizer)
                                ::anomalies/category ::anomalies/third-party}))))))

(s/fdef quadratic-programming-joptimizer
        :args (s/and (s/cat :apache-quadratic-objective ::apache-quadratic-objective
                            :opts (s/? (s/keys :opt [::linear-objective
                                                     ::equal
                                                     ::less-than
                                                     ::quadratic-guess])))
                     (fn [{:keys [apache-quadratic-objective opts]}]
                       (let [{::keys [linear-objective equal less-than quadratic-guess]} opts
                             vars (apache-mx/rows apache-quadratic-objective)
                             [eq-mx eq-v] equal
                             [lt-mx lt-v] less-than]
                         (and (or (not linear-objective)
                                  (= vars (count linear-objective)))
                              (or (not quadratic-guess)
                                  (= vars (count quadratic-guess)))
                              (or (not equal)
                                  (and (= vars (mx/columns eq-mx))
                                       (>= vars (mx/rows eq-mx))))
                              (or (not less-than)
                                  (= vars (mx/columns lt-mx))))))
                     (fn [{:keys [opts]}]
                       (let [{::keys [equal less-than quadratic-guess]} opts
                             [eq-mx eq-v] equal
                             [lt-mx lt-v] less-than]
                         (and (or (not equal)
                                  (not quadratic-guess)
                                  (every?
                                    zero?
                                    (tensor/subtract
                                      (mx/get-column
                                        (mx/mx* eq-mx (mx/column-matrix quadratic-guess))
                                        0)
                                      eq-v)))
                              (or (not less-than)
                                  (not quadratic-guess)
                                  (every?
                                    m/non+?
                                    (tensor/subtract
                                      (mx/get-column
                                        (mx/mx* lt-mx (mx/column-matrix quadratic-guess))
                                        0)
                                      lt-v)))))))
        :ret (s/or :anomaly ::anomalies/anomaly
                   :solution ::vector/vector))

;;;COBYLA ;https://github.com/cureos/jcobyla
(s/def ::double-array->number
  (s/fspec :args (s/cat :a ::arrays/double-array)
           :ret ::m/number))

(s/def ::double-array->vector
  (s/fspec :args (s/cat :a ::arrays/double-array)
           :ret ::vector/vector))

(s/def ::max-iter
  (s/with-gen (s/nilable ::m/int+)
              #(gen/one-of [(s/gen (s/int-in 100 1000))
                            (gen/return nil)])))

(s/def ::initial-change
  (s/with-gen
    ::m/finite+
    #(s/gen (s/double-in :infinite? false :NaN? false :min 1e-3 :max 1e3))))

(s/def ::abs-accu
  (s/with-gen
    ::m/finite+
    #(s/gen (s/double-in :infinite? false :NaN? false :min 9e-16 :max 1.0))))

(s/def ::debug-print-level
  (s/with-gen #{0 1 2 3}
              #(gen/return 0)))

(defn- cobyla-calc
  [f cons-fn]
  (reify Calcfc
    (compute [_ _ m x con]
      (do (when cons-fn
            (doseq [i (range m)]
              (aset-double con i (get (cons-fn x) i m/nan))))
          (f x)))))

(defn nlp-cobyla
  "`::objective` should accept a double-array of variable values and return the
  objective value. `::geq-fn` should accept a double-array of variable values
  and return a vector of constraint values that must be greater than or equal to
  zero to meet constraint. `::cobyla-guess` should be a double-array that is not
  all zeros. `::max-iter` is the maximum allowed number of function evaluations
  (default 100000). `::initial-change` is a reasonable first change in the
  variable values (default 0.5). `::abs-accu` is the precision (default 1e-8).
  This is a lower bound on the size of the trust region. Final accuracy in the
  optimization is not precisely guaranteed. `::debug-print-level`
  (0, 1, 2, or 3) specifies the level of output to the console (default 0)."
  ([args] (nlp-cobyla args {}))
  ([{::keys [objective geq-fn cobyla-guess]}
    {::keys [max-iter initial-change abs-accu debug-print-level]
     :or    {initial-change    0.5
             abs-accu          1e-8
             debug-print-level 0}}]
   (let [max-iter (or max-iter 1000)
         cobyla-guess (double-array cobyla-guess)
         ^CobylaExitStatus s (Cobyla/findMinimum
                               (cobyla-calc objective geq-fn)
                               (count cobyla-guess)
                               (count (geq-fn cobyla-guess))
                               cobyla-guess
                               (double initial-change)
                               (double abs-accu)
                               (long debug-print-level)
                               (long max-iter))]
     (condp = s
       CobylaExitStatus/NORMAL (vec cobyla-guess)
       CobylaExitStatus/MAX_ITERATIONS_REACHED {::anomalies/message  "Max Iterations Reached"
                                                ::anomalies/fn       (var nlp-cobyla)
                                                ::anomalies/category ::anomalies/third-party}
       CobylaExitStatus/DIVERGING_ROUNDING_ERRORS {::anomalies/message  "Diverging Rounding Errors"
                                                   ::anomalies/fn       (var nlp-cobyla)
                                                   ::anomalies/category ::anomalies/third-party}
       {::anomalies/message  "Cobyla Error"
        ::anomalies/fn       (var nlp-cobyla)
        ::anomalies/category ::anomalies/third-party}))))

(s/def ::objective ::double-array->number)
(s/def ::geq-fn ::double-array->vector)

;should this be no zeros? or just not all zeros? -- doesn't seem to be needed
(s/def ::cobyla-guess
  (s/and ::vector/vector-finite
         (fn [v]
           (pos? (count v)))
         (constantly true) #_(not (every? zero? %))))

(s/def ::objective-constraints-and-guess
  (s/with-gen
    (s/keys :req [::objective ::geq-fn ::cobyla-guess])
    #(gen/one-of
       (map
         gen/return
         (list {::objective    (fn [da]
                                 (if (= (count da) 2)
                                   (let [[a b] da]
                                     (+ a b))
                                   m/inf+))
                ::geq-fn       (fn [da]
                                 (if (= (count da) 2)
                                   (let [[a b] da]
                                     [(dec a) (inc b)])
                                   [m/inf- m/inf-]))
                ::cobyla-guess [0.0 0.0]})))))

(s/fdef nlp-cobyla
        :args (s/cat :objective-constraints-and-guess ::objective-constraints-and-guess
                     :opts (s/? (s/keys :opt [::max-iter ::initial-change
                                              ::abs-accu ::debug-print-level])))
        :ret (s/or :solution ::vector/vector
                   :anomaly ::anomalies/anomaly))


