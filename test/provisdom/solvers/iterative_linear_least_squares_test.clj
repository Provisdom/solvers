(ns provisdom.solvers.iterative-linear-least-squares-test
  (:require [clojure.test :refer :all]
            [provisdom.test.core :refer :all]
            [provisdom.solvers.iterative-linear-least-squares :as iterative-lls]
            [clojure.spec.test.alpha :as st]
            [orchestra.spec.test :as ost]))

;8 seconds

(set! *warn-on-reflection* true)

(ost/instrument)

(deftest iterative-linear-least-squares-test
  (is (spec-check iterative-lls/iterative-linear-least-squares
                  {:coll-check-limit 10
                   :coll-error-limit 10
                   :fspec-iterations 10
                   :recursion-limit  1
                   :test-check       {:num-tests 200}}))
  (is= [1.3333333333333333 0.33333333333333326]
       (iterative-lls/iterative-linear-least-squares
         [[1 2] [2 1]]
         [2 3]))
  (is= [5.999999999999998 1.9999999999999984]
       (iterative-lls/iterative-linear-least-squares
         [[1.0 0.5] [0.5 3.0]]
         [7.0 9.0]))
  (is= [6.26666666666666 1.4666666666666666]
       (iterative-lls/iterative-linear-least-squares
         [[1.0 0.5] [0.5 4.0]]
         [7.0 9.0])))

#_(ost/unstrument)