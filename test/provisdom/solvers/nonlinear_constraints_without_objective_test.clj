(ns provisdom.solvers.nonlinear-constraints-without-objective-test
  (:require [clojure.test :refer :all]
            [provisdom.test.core :refer :all]
            [provisdom.math.core :as m]
            [provisdom.solvers.nonlinear-constraints-without-objective :as nlc]
            [clojure.spec.test.alpha :as st]
            [orchestra.spec.test :as ost]))

;17 seconds

(set! *warn-on-reflection* true)

(ost/instrument)

(defn constraints-fn12
  [da]
  (if (= (count da) 2)
    (let [[a b] da]
      [(- a b 0.5)
       (- (m/sq a) (m/sq b) 0.5)])
    [m/inf- m/inf-]))

(defn jac12
  [da]
  (when (= (count da) 2)
    (let [[a b] da]
      [[1.0 -1.0]
       [(* 2.0 a) (* -2.0 b)]])))

(defn constraints-fn123
  [da]
  (if (= (count da) 2)
    (let [[a b] da]
      [(- a b 0.5)
       (- (m/sq a) (m/sq b) 0.5)
       (- (m/cube a) (m/cube b) 0.5)])
    [m/inf- m/inf-]))

(defn jac123
  [da]
  (when (= (count da) 2)
    (let [[a b] da]
      [[1.0 -1.0]
       [(* 2.0 a) (* -2.0 b)]
       [(* 3.0 (m/sq a)) (* -3.0 (m/sq b))]])))

(def vars-guess [1.5 0.5])

(deftest nonlinear-least-squares-test
  (is (spec-check nlc/nonlinear-least-squares
                  {:coll-check-limit 10
                   :coll-error-limit 10
                   :fspec-iterations 10
                   :recursion-limit  1
                   :test-check       {:num-tests 400}}))
  (is= {::nlc/met-constraints?           true
        ::nlc/vector-point               [0.75 0.25]
        ::nlc/weighted-constraint-errors [0.0 0.0]}
       (nlc/nonlinear-least-squares
         {::nlc/constraints-fn         constraints-fn12
          ::nlc/constraint-jacobian-fn jac12
          ::nlc/vars-guess             vars-guess}
         {::nlc/nls-solver-type :levenberg-marquardt}))
  (is= {::nlc/met-constraints?           true
        ::nlc/vector-point               [0.7499999999999999 0.2499999999999999]
        ::nlc/weighted-constraint-errors [0.0 5.551115123125783E-17]}
       (nlc/nonlinear-least-squares
         {::nlc/constraints-fn         constraints-fn12
          ::nlc/constraint-jacobian-fn jac12
          ::nlc/vars-guess             vars-guess}
         {::nlc/nls-solver-type :gauss-newton}))
  (is= {::nlc/met-constraints?           true
        ::nlc/vector-point               [0.75 0.25]
        ::nlc/weighted-constraint-errors [0.0 0.0]}
       (nlc/nonlinear-least-squares
         {::nlc/constraints-fn         constraints-fn12
          ::nlc/constraint-jacobian-fn jac12
          ::nlc/vars-guess             vars-guess}))
  (is= {::nlc/met-constraints?           false
        ::nlc/vector-point               [0.7982990600415211 0.3153940258519974]
        ::nlc/weighted-constraint-errors [0.017094965810476315 -0.037807997720045616 0.022632179393191842]}
       (nlc/nonlinear-least-squares
         {::nlc/constraints-fn         constraints-fn123
          ::nlc/constraint-jacobian-fn jac123
          ::nlc/vars-guess             vars-guess}))
  (is= {::nlc/met-constraints?           true
        ::nlc/vector-point               [0.75 0.24999999999999992]
        ::nlc/weighted-constraint-errors [-1.1102230246251565E-10 0.0]}
       (nlc/nonlinear-least-squares
         {::nlc/constraints-fn         constraints-fn12
          ::nlc/constraint-jacobian-fn jac12
          ::nlc/vars-guess             vars-guess}
         {::nlc/constraint-weights [1e12 1e12]}))
  (is= {::nlc/met-constraints?           false
        ::nlc/vector-point               [0.7704812980835266 0.27205964469589844]
        ::nlc/weighted-constraint-errors [0.0035292903173199186 -0.019624980424417737 0.028061431796740102]}
       (nlc/nonlinear-least-squares
         {::nlc/constraints-fn         constraints-fn123
          ::nlc/constraint-jacobian-fn jac123
          ::nlc/vars-guess             vars-guess}
         {::nlc/constraint-weights [5.0 1.0 0.2]}))
  (is= {::nlc/met-constraints?           false
        ::nlc/vector-point               [0.7500001406248287 0.25000014062488135]
        ::nlc/weighted-constraint-errors [5.268008251846368E-11 -1.4062480235832453E-7 9.374978906273723E-5]}
       (nlc/nonlinear-least-squares
         {::nlc/constraints-fn         constraints-fn123
          ::nlc/constraint-jacobian-fn jac123
          ::nlc/vars-guess             vars-guess}
         {::nlc/constraint-weights [1e6 1.0 1e-6]})))

(deftest nonlinear-ordered-constraints-test
  (is (spec-check nlc/nonlinear-ordered-constraints
                  {:coll-check-limit 10
                   :coll-error-limit 10
                   :fspec-iterations 10
                   :recursion-limit  1
                   :test-check       {:num-tests 50}}))
  (is= {::nlc/vector-point              [0.75 0.25000000000000006]
        ::nlc/number-of-constraints-met 2}
       (nlc/nonlinear-ordered-constraints
         {::nlc/constraints-fn         constraints-fn123
          ::nlc/constraint-jacobian-fn jac123
          ::nlc/vars-guess             vars-guess})))

#_(ost/unstrument)