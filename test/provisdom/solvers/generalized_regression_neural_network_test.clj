(ns provisdom.solvers.generalized-regression-neural-network-test
  (:require
    [clojure.test :refer :all]
    [provisdom.test.core :refer :all]
    [provisdom.solvers.generalized-regression-neural-network :as gen-reg-nn]
    [orchestra.spec.test :as ost]))

;;? seconds

(set! *warn-on-reflection* true)

(ost/instrument)

(deftest solve-for-spread-test
  (is (spec-check gen-reg-nn/solve-for-spread
                  {:coll-check-limit 10
                   :coll-error-limit 10
                   :fspec-iterations 10
                   :recursion-limit  1
                   :test-check       {:num-tests 300}}))
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
          ::gen-reg-nn/spread-interval [1e-14 10.0]})))

(deftest solve-test
  (is (spec-check gen-reg-nn/solve))
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
          ::gen-reg-nn/spread 1.7988160558737314})))

#_(ost/unstrument)