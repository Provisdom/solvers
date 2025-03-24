(ns provisdom.solvers.neanderthal.curve-fitting.smoothing-splines-test
  (:require
    [clojure.spec.test.alpha :as st]
    [clojure.test :refer :all]
    [provisdom.solvers.neanderthal.curve-fitting.smoothing-splines :as splines]
    [provisdom.test.core :refer :all]))

;9 seconds

(set! *warn-on-reflection* true)

(deftest smoothing-cubic-splines-dof-test
  (with-instrument `splines/smoothing-cubic-splines-dof
    (is (spec-check splines/smoothing-cubic-splines-dof)))
  (with-instrument (st/instrumentable-syms)
    (is= 2.34173669467787
      (splines/smoothing-cubic-splines-dof
        {::splines/x-vals              [0.0 1.0 2.0 3.0]
         ::splines/variances           [1.0 1.0 1.0 1.0]
         ::splines/smoothing-parameter 0.5}))))

(deftest smoothing-cubic-splines-test
  (with-instrument `splines/smoothing-cubic-splines
    (is (spec-check splines/smoothing-cubic-splines)))
  (with-instrument (st/instrumentable-syms)
    (let [ret (splines/smoothing-cubic-splines
                {::splines/x-vals              [0.0 1.0 2.0 3.0]
                 ::splines/f-vals              [1.0 4.0 25.0 81.0]
                 ::splines/smoothing-parameter 0.5})
          s-f (::splines/smoothing-cubic-splines-fn ret)]
      (is= 2.34173669467787 (::splines/degrees-of-freedom ret))
      (is= -24.587301587301592 (s-f -1.0))
      (is= -7.54341736694678 (s-f 0))
      (is= 10.92436974789916 (s-f 1))
      (is= 36.78151260504202 (s-f 2.0))
      (is= 106.5873015873016 (s-f 4.0)))))
