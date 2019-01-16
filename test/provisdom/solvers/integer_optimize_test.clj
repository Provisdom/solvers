(ns provisdom.solvers.integer-optimize-test
  (:require
    [clojure.test :refer :all]
    [provisdom.test.core :refer :all]
    [provisdom.solvers.integer-optimize :as int-opt]
    [orchestra.spec.test :as ost]))

;;5 seconds

(set! *warn-on-reflection* true)

(ost/instrument)

(deftest integer-optimize-test
  (is (spec-check int-opt/integer-optimize))
  (is= 0 (int-opt/integer-optimize #(+ 25 (* % % -0.5)) 0 [-2000 200]))
  (is= 200 (int-opt/integer-optimize identity 5 [-2000 200]))
  (is= 200 (int-opt/integer-optimize identity -5 [-2000 200]))
  (is= -2000 (int-opt/integer-optimize #(- %) 5 [-2000 200]))
  (is= -2000 (int-opt/integer-optimize #(- %) -5 [-2000 200]))
  (is= 15 (int-opt/integer-optimize #(+ (* 30 %) (* % % -1)) 5 [-2000 200]))
  (is= 15 (int-opt/integer-optimize #(+ (* 30 %) (* % % -1)) -5 [-2000 200]))
  (is= -15 (int-opt/integer-optimize #(+ (* -30 %) (* % % -1)) -5 [-2000 200]))
  (is= -15 (int-opt/integer-optimize #(+ (* -30 %) (* % % -1)) 5 [-2000 200])))

#_(ost/unstrument)