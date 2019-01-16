(ns provisdom.solvers.quadratic-programming-test
  (:require [clojure.test :refer :all]
            [provisdom.test.core :refer :all]
            [provisdom.solvers.quadratic-programming :as qp]
            [provisdom.math.apache-matrix :as apache-mx]
            [clojure.spec.test.alpha :as st]
            [orchestra.spec.test :as ost]))

;10 seconds

(set! *warn-on-reflection* true)

(ost/instrument)

(def eq1 [[[1.0 1.0]] [3.0]])
(def lt2 [[[1.0 1.0]] [3.0]])
(def obj1 (apache-mx/apache-matrix [[1.0 0.5] [0.5 2.0]]))

(deftest quadratic-programming-test
  (is (spec-check qp/quadratic-programming
                  {:coll-check-limit 10
                   :coll-error-limit 10
                   :fspec-iterations 10
                   :recursion-limit  1
                   :test-check       {:num-tests 150}}))
  (is= {::qp/value        0.0
        ::qp/vector-point [0.0 0.0]}
       (qp/quadratic-programming obj1))
  (is= {::qp/value        1.4375000000000018
        ::qp/vector-point [3.250000000000001 -0.2500000000000002]}
       (qp/quadratic-programming obj1
                                 {::qp/linear-objective [-1.0 1.0]
                                  ::qp/equal            eq1
                                  ::qp/guess            [2.0 1.0]}))
  (is= {::qp/value        1.4375000000000013
        ::qp/vector-point [3.250000000000001 -0.25000000000000044]}
       (qp/quadratic-programming obj1
                                 {::qp/linear-objective [-1.0 1.0]
                                  ::qp/equal            eq1}))
  (is= {::qp/value        -1.1428571428571428
        ::qp/vector-point [1.4285714285713031 -0.857142857142899]}
       (qp/quadratic-programming obj1
                                 {::qp/linear-objective [-1.0 1.0]
                                  ::qp/less-than        lt2}))
  (is= {::qp/value        -0.5257142857142857
        ::qp/vector-point [1.0857142857142095 -0.1714285714285968]}
       (qp/quadratic-programming
         obj1
         {::qp/linear-objective [-1.0 -0.2]
          ::qp/less-than        lt2
          ::qp/guess            [-2.4210526315789473 0.4210526315789473]}))
  (is= {::qp/value        -0.5257142857142858
        ::qp/vector-point [1.085714285714155 -0.17142857142861498]}
       (qp/quadratic-programming obj1
                                 {::qp/linear-objective [-1.0 -0.2]
                                  ::qp/less-than        lt2})))

#_(ost/unstrument)