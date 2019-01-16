(ns provisdom.solvers.logistic-regression-test
  (:require
    [clojure.test :refer :all]
    [provisdom.test.core :refer :all]
    [provisdom.solvers.logistic-regression :as log-regress]
    [orchestra.spec.test :as ost]))

;;8 seconds

(set! *warn-on-reflection* true)

(ost/instrument)

(deftest iterative-reweighted-least-squares-test
  (is (spec-check log-regress/iterative-reweighted-least-squares
                  {:coll-check-limit 10
                   :coll-error-limit 10
                   :fspec-iterations 10
                   :recursion-limit  1
                   :test-check       {:num-tests 25}}))
  (is= [0.26338258410495197 -0.02601474092508002]
       (log-regress/iterative-reweighted-least-squares
         [[1 2] [3 5] [7 23]] [0.5 0.7 0.8]))
  (is= [0.42977931467539693 -0.05858562072786558]
       (log-regress/iterative-reweighted-least-squares
         [[1 2] [3 5] [7 23]] [0.56 0.8 0.88]))
  (is= [0.7607052740510958 -0.15751338511807894]
       (log-regress/iterative-reweighted-least-squares
         [[0 30] [1 2] [3 5] [7 23] [10 10]] [0.0 0.56 0.8 0.88 1.0]))
  (is= [0.760705275242646 -0.15751338542141294]
       (log-regress/iterative-reweighted-least-squares
         [[0 30] [1 2] [3 5] [7 23] [10 10]]
         [0.0 0.56 0.8 0.88 1.0]
         {::log-regress/ridge-lambda 1e-9})))

(deftest log-likelihood-of-data-given-parameters-test
  ;;gen/return -- could use gen/bind to create gen tests -- simple function
  (is (spec-check log-regress/log-likelihood-of-data-given-parameters
                  {:coll-check-limit 10
                   :coll-error-limit 10
                   :fspec-iterations 10
                   :recursion-limit  1
                   :test-check       {:num-tests 1}}))
  (is= -1.8153910326882052
       (log-regress/log-likelihood-of-data-given-parameters
         {::log-regress/x-mx [[1 2] [3 5] [7 23]]
          ::log-regress/y [0.5 0.7 0.8]
          ::log-regress/parameters [0.26338258410495197 -0.02601474092508002]})))

#_(ost/unstrument)