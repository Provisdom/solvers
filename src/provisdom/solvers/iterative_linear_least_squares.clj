(ns provisdom.solvers.iterative-linear-least-squares
  (:require
    [clojure.spec.alpha :as s]
    [provisdom.math.core :as m]
    [provisdom.math.matrix :as mx]
    [provisdom.math.vector :as vector]
    [provisdom.solvers.internal-apache-solvers :as apache-solvers]
    [provisdom.utility-belt.anomalies :as anomalies]))

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
                   :or    {max-iter                         10000
                           rel-accu                         m/dbl-close
                           use-conjugate-gradient?          false
                           check-for-positive-definiteness? true}}]
   (apache-solvers/iterative-linear-least-squares
     symmetric-m v
     {::apache-solvers/max-iter                         max-iter
      ::apache-solvers/rel-accu                         rel-accu
      ::apache-solvers/use-conjugate-gradient?          use-conjugate-gradient?
      ::apache-solvers/vector-initial-guess             vector-initial-guess
      ::apache-solvers/check-for-positive-definiteness? check-for-positive-definiteness?})))

(s/def ::max-iter ::apache-solvers/max-iter)
(s/def ::rel-accu ::apache-solvers/rel-accu)
(s/def ::use-conjugate-gradient? ::apache-solvers/use-conjugate-gradient?)
(s/def ::vector-initial-guess ::apache-solvers/vector-initial-guess)
(s/def ::check-for-positive-definiteness? ::apache-solvers/check-for-positive-definiteness?)

(s/fdef iterative-linear-least-squares
        :args (s/and (s/cat :symmetric-m ::mx/symmetric-matrix
                            :v ::vector/vector-finite
                            :opts (s/? (s/keys :opt [::max-iter
                                                     ::rel-accu
                                                     ::use-conjugate-gradient?
                                                     ::vector-initial-guess
                                                     ::check-for-positive-definiteness?])))
                     (fn [{:keys [symmetric-m v]}]
                       (and (= (mx/rows symmetric-m) (count v))
                            (pos? (count v)))))
        :ret (s/or :solution ::vector/vector
                   :anomaly ::anomalies/anomaly))
