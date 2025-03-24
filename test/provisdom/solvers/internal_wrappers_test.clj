(ns provisdom.solvers.internal-wrappers-test
  (:require [clojure.test :refer :all]
            [provisdom.test.core :refer :all]
            [provisdom.solvers.internal-wrappers :as wrap]
            [provisdom.apache-math.apache-matrix :as apache-mx]
            [clojure.spec.test.alpha :as st]
            [orchestra.spec.test :as ost]
            [provisdom.math.core :as m]))

;;;20 seconds

(set! *warn-on-reflection* true)

(ost/instrument)

(deftest quadratic-programming-joptimizer-test
  (is (spec-check wrap/quadratic-programming-joptimizer
                  {:coll-check-limit 10
                   :coll-error-limit 10
                   :fspec-iterations 10
                   :recursion-limit  1
                   :test-check       {:num-tests 200}}))
  (is= [0.0 0.0]
       (wrap/quadratic-programming-joptimizer
         (apache-mx/->apache-matrix [[1.0 0.5] [0.5 2.0]])))
  (is= [1.0857142857142856 -0.17142857142857143]
       (wrap/quadratic-programming-joptimizer
         (apache-mx/->apache-matrix [[1.0 0.5] [0.5 2.0]])
         {::wrap/linear-objective [-1.0 -0.2]}))
  (is= [1.4285714285714286 -0.8571428571428571]
       (wrap/quadratic-programming-joptimizer
         (apache-mx/->apache-matrix [[1.0 0.5] [0.5 2.0]])
         {::wrap/linear-objective [-1.0 1.0]}))
  (is= [1.4285714285713031 -0.857142857142899]
       (wrap/quadratic-programming-joptimizer
         (apache-mx/->apache-matrix [[1.0 0.5] [0.5 2.0]])
         {::wrap/linear-objective [-1.0 1.0]
          ::wrap/less-than        [[[1.0 1.0]] [3.0]]}))
  (is= [3.250000000000001 -0.2500000000000002]
       (wrap/quadratic-programming-joptimizer
         (apache-mx/->apache-matrix [[1.0 0.5] [0.5 2.0]])
         {::wrap/linear-objective [-1.0 1.0]
          ::wrap/equal            [[[1.0 1.0]] [3.0]]
          ::wrap/quadratic-guess  [2.0 1.0]})))

(defn obj-fn1
  [da]
  (if (= (count da) 2)
    (let [[a b] da]
      (+ a b))
    m/inf+))

(defn geq-fn1
  [da]
  (if (= (count da) 2)
    (let [[a b] da]
      [(dec a) (inc b)])
    [m/inf- m/inf-]))

(deftest nlp-cobyla-test
  (is (spec-check wrap/nlp-cobyla))
  (is= [1.0 -1.0]
       (wrap/nlp-cobyla {::wrap/objective    obj-fn1
                         ::wrap/geq-fn       geq-fn1
                         ::wrap/cobyla-guess [0.0 0.0]})))

#_(ost/unstrument)