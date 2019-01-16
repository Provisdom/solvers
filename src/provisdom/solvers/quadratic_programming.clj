(ns provisdom.solvers.quadratic-programming
  (:require
    [clojure.spec.alpha :as s]
    [clojure.spec.gen.alpha :as gen]
    [clojure.spec.test.alpha :as st]
    [orchestra.spec.test :as ost]
    [provisdom.utility-belt.anomalies :as anomalies]
    [provisdom.math.core :as m]
    [provisdom.math.vector :as vector]
    [provisdom.math.matrix :as mx]
    [provisdom.math.tensor :as tensor]
    [provisdom.math.apache-matrix :as apache-mx]
    [provisdom.solvers.internal-wrappers :as wrap]
    [provisdom.solvers.nonlinear-programming :as nlp]))

(defn- qp-to-objective
  [apache-quadratic-objective linear-objective constant-objective]
  (fn [da]
    (if (zero? (count da))
      m/nan
      (let [v (vec da)
            obj-mx (apache-mx/apache-matrix->matrix apache-quadratic-objective)
            new-mx (mx/mx* (mx/row-matrix v) obj-mx)
            final-mx (mx/mx* new-mx (mx/column-matrix v))]
        (+ constant-objective
           (* 0.5 (ffirst final-mx))
           (if linear-objective
             (vector/dot-product linear-objective v)
             0.0))))))

(comment "unused"
         (defn- qp-to-nlp-constraints
           [equal less-than]
           (fn [da]
             (let [v (vec da)]
               (concat
                 (map (fn [[lhs-mx rhs-v]]
                        (tensor/subtract rhs-v (mx/mx* lhs-mx v)))
                      (concat less-than equal))
                 (map (fn [[lhs-mx rhs-v]]
                        (tensor/subtract (mx/mx* lhs-mx v) rhs-v))
                      equal))))))

(defn quadratic-programming
  "Minimizes (1/2) x^T P x + q^T x over x, where 'P' is the
  `apache-quadratic-objective` matrix and 'q' is the `::linear-objective`
  vector. `apache-quadratic-objective` should be a positive (symmetrical) Apache
  matrix. `::equal` and `::less-than` should be tuples of matrices and vectors,
  where the matrix is the left-hand side and the vector is the right-hand-side.
  `::quadratic-guess` must meet constraints and must strictly meet the `less-than`
  constraints. Returns a map with `::vector-point` and `::value`."
  ([apache-quadratic-objective]
   (quadratic-programming apache-quadratic-objective {}))
  ([apache-quadratic-objective
    {::keys [linear-objective constant-objective equal less-than guess]
     :or    {constant-objective 0.0}}]
   (let [solution (wrap/quadratic-programming-joptimizer
                    apache-quadratic-objective
                    {::wrap/linear-objective linear-objective
                     ::wrap/equal            equal
                     ::wrap/less-than        less-than
                     ::wrap/quadratic-guess  guess})
         objective (qp-to-objective apache-quadratic-objective linear-objective constant-objective)]
     (if (anomalies/anomaly? solution)
       solution
       {::vector-point solution
        ::value        (objective solution)}))))

(s/def ::apache-quadratic-objective ::wrap/apache-quadratic-objective)
(s/def ::linear-objective ::wrap/linear-objective)
(s/def ::constant-objective ::m/finite)
(s/def ::equal ::wrap/equal)
(s/def ::less-than ::wrap/less-than)
(s/def ::guess ::wrap/quadratic-guess)
(s/def ::vector-point ::nlp/vector-point)
(s/def ::value ::nlp/value)

(s/fdef quadratic-programming
        :args (s/and (s/cat :apache-quadratic-objective ::apache-quadratic-objective
                            :opts (s/? (s/keys :opt [::linear-objective
                                                     ::equal
                                                     ::less-than
                                                     ::guess])))
                     (fn [{:keys [apache-quadratic-objective opts]}]
                       (let [{::keys [linear-objective equal less-than guess]} opts
                             vars (apache-mx/rows apache-quadratic-objective)
                             [eq-mx eq-v] equal
                             [lt-mx lt-v] less-than]
                         (and (or (not linear-objective)
                                  (= vars (count linear-objective)))
                              (or (not guess)
                                  (= vars (count guess)))
                              (or (not equal)
                                  (and (= vars (mx/columns eq-mx))
                                       (>= vars (mx/rows eq-mx))))
                              (or (not less-than)
                                  (= vars (mx/columns lt-mx))))))
                     (fn [{:keys [opts]}]
                       (let [{::keys [equal less-than guess]} opts
                             [eq-mx eq-v] equal
                             [lt-mx lt-v] less-than]
                         (and (or (not equal)
                                  (not guess)
                                  (every?
                                    zero?
                                    (tensor/subtract
                                      (mx/get-column
                                        (mx/mx* eq-mx (mx/column-matrix guess))
                                        0)
                                      eq-v)))
                              (or (not less-than)
                                  (not guess)
                                  (every?
                                    m/non+?
                                    (tensor/subtract
                                      (mx/get-column
                                        (mx/mx* lt-mx (mx/column-matrix guess))
                                        0)
                                      lt-v)))))))
        :ret (s/or :solution (s/keys :req [::vector-point ::value])
                   :anomaly ::anomalies/anomaly))