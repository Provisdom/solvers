(ns provisdom.solvers.integer-root-solvers-test
  (:require
    [clojure.test :refer :all]
    [provisdom.test.core :refer :all]
    [provisdom.solvers.integer-root-solvers :as int-root]
    [orchestra.spec.test :as ost]))

;;10 seconds

(set! *warn-on-reflection* true)

(ost/instrument)

(deftest integer-root-solver-test
  (is (spec-check int-root/integer-root-solver))
  (is= 1 (int-root/integer-root-solver #(+ % -0.5) [-2000 200]))
  (is= 0 (int-root/integer-root-solver identity [-2000 200]))
  (is= 0 (int-root/integer-root-solver #(+ % 0.5) [-2000 200])))

#_(ost/unstrument)