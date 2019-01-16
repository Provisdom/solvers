(ns provisdom.solvers.multinomial-logistic-regression-test
  (:require
    [clojure.test :refer :all]
    [provisdom.test.core :refer :all]
    [provisdom.solvers.multinomial-logistic-regression :as multi-log-regress]
    [orchestra.spec.test :as ost]))

;;6 seconds

(set! *warn-on-reflection* true)

(ost/instrument)

(deftest iterative-reweighted-least-squares-test
  (is (spec-check multi-log-regress/iterative-reweighted-least-squares
                  {:coll-check-limit 10
                   :coll-error-limit 10
                   :fspec-iterations 10
                   :recursion-limit  1
                   :test-check       {:num-tests 8}}))
  (is= [[0.26338258410495197 -0.02601474092508002]
        [0.42977931467539693 -0.05858562072786558]]
       (multi-log-regress/iterative-reweighted-least-squares
         [[1 2] [3 5] [7 23]] [[0.5 0.56] [0.7 0.8] [0.8 0.88]])))

#_(ost/unstrument)