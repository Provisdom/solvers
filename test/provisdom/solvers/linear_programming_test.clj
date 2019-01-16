(ns provisdom.solvers.linear-programming-test
  (:require [clojure.test :refer :all]
            [provisdom.test.core :refer :all]
            [provisdom.solvers.linear-programming :as lp]
            [clojure.spec.test.alpha :as st]
            [orchestra.spec.test :as ost]))

;3 seconds

(set! *warn-on-reflection* true)

(ost/instrument)

(def c1 [[1.0 0.0 0.0] :geq 2.0])
(def c2 [[1.0 2.0 3.0] :geq 5.0])
(def c3 [[9.0 4.0 14.0] :eq 49.0])

(deftest linear-programming-test
  (is (spec-check lp/linear-programming))
  (is= {::lp/value        2.7500000000000018
        ::lp/vector-point [2.0 -3.187499999999999 3.125]}
       (lp/linear-programming [3.0 2.0 1.0] [c1 c2 c3]))
  (is= {::lp/value        4.000000000000002
        ::lp/vector-point [2.0 -3.187499999999999 3.125]}
       (lp/linear-programming [3.0 2.0 1.0] [c1 c2 c3]
                              {::lp/objective-constant 1.25}))
  (is= {::lp/value        8.214285714285714
        ::lp/vector-point [2.0 0.0 2.214285714285714]}
       (lp/linear-programming [3.0 2.0 1.0] [c1 c2 c3]
                              {::lp/non-negative-vars? true}))
  (is= {::lp/value        5.000000000000002
        ::lp/vector-point [2.0 -3.187499999999999 3.125]}
       (lp/linear-programming [1.0 2.0 3.0]
                              [[[1.0 0.0 0.0] :leq 2.0]
                               [[1.0 2.0 3.0] :leq 5.0]
                               c3]
                              {::lp/goal :max}))
  (is= {::lp/value        2.7500000000000018
        ::lp/vector-point [2.0 -3.187499999999999 3.125]}
       (lp/linear-programming [3.0 2.0 1.0] [c1 c1 c2 c3])))

#_(ost/unstrument)