(ns provisdom.solvers.generalized-regression-neural-network-test
  (:require
    [clojure.spec.test.alpha :as st]
    [clojure.test :refer :all]
    [provisdom.solvers.generalized-regression-neural-network :as gen-reg-nn]
    [provisdom.test.core :refer :all]))

;;? seconds

(set! *warn-on-reflection* true)

(deftest solve-for-spread-test
  (with-instrument `gen-reg-nn/solve-for-spread
    (is (spec-check gen-reg-nn/solve-for-spread)))
  (with-instrument (st/instrumentable-syms)
    (is= 1.7988160558737314
      (gen-reg-nn/solve-for-spread
        {::gen-reg-nn/x-training-mx   [[1.0 2.0 3.0]
                                       [3.0 2.0 5.0]
                                       [8.0 3.0 6.0]]
         ::gen-reg-nn/y-training      [5.0 2.0 6.0]
         ::gen-reg-nn/x-test-mx       [[4.0 2.0 6.0]
                                       [2.0 4.0 6.0]]
         ::gen-reg-nn/y-test          [1.0 5.0]
         ::gen-reg-nn/spread-guess    1.0
         ::gen-reg-nn/spread-interval [1e-14 10.0]}))))

(deftest solve-test
  (with-instrument `gen-reg-nn/solve
    (is (spec-check gen-reg-nn/solve)))
  (with-instrument (st/instrumentable-syms)
  (is= 3.676761977142881
    (gen-reg-nn/solve
      {::gen-reg-nn/x-mx   [[1.0 2.0 3.0]
                            [3.0 2.0 5.0]
                            [8.0 3.0 6.0]
                            [4.0 2.0 6.0]
                            [2.0 4.0 6.0]]
       ::gen-reg-nn/y      [5.0 2.0 6.0 1.0 5.0]
       ::gen-reg-nn/x-new  [2.0 3.0 4.0]
       ::gen-reg-nn/spread 1.0}))
  (is= 3.5383753122108947
    (gen-reg-nn/solve
      {::gen-reg-nn/x-mx   [[1.0 2.0 3.0]
                            [3.0 2.0 5.0]
                            [8.0 3.0 6.0]
                            [4.0 2.0 6.0]
                            [2.0 4.0 6.0]]
       ::gen-reg-nn/y      [5.0 2.0 6.0 1.0 5.0]
       ::gen-reg-nn/x-new  [2.0 3.0 4.0]
       ::gen-reg-nn/spread 1.7988160558737314}))))
