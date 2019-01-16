(ns provisdom.solvers.internal-metaheuristics-test
  (:require
    [clojure.test :refer :all]
    [provisdom.test.core :refer :all]
    [provisdom.solvers.internal-metaheuristics :as meta]
    [orchestra.spec.test :as ost]))

;;0 seconds

(set! *warn-on-reflection* true)

(ost/instrument)

(def f identity)
(def guess 2.0)

#_(deftest swarm-test
  (is (spec-check meta/swarm))
  (is= 0.0 (meta/swarm f guess)))

#_(ost/unstrument)