(ns provisdom.solvers.neanderthal.ordinary-regression-test
  (:require
    [clojure.test :refer :all]
    [provisdom.test.core :refer :all]
    [provisdom.solvers.neanderthal.ordinary-regression :as ord-regress]
    [provisdom.math.core :as m]
    [clojure.spec.test.alpha :as st]
    [orchestra.spec.test :as ost]))

;;4 seconds

(set! *warn-on-reflection* true)

(ost/instrument)

(def x-mx
  [[7.0 2.8] [8.0 9.3] [3.0 9.9]
   [3.3 2.2] [3.8 5.9] [2.2 8.9]
   [4.6 2.5] [2.4 3.5] [8.4 2.3]])

(def y [1 2 3 4 5 6 7 8 9])

(deftest simple-ordinary-regression-test
  (is (spec-check ord-regress/simple-ordinary-regression))
  (is (data-approx=
        {::ord-regress/weights                [0.624874327896168
                                               0.22810879461175296]
         ::ord-regress/mean-squared-error     12.177883256018546
         ::ord-regress/standard-squared-error 15.65727847202385}
        (ord-regress/simple-ordinary-regression x-mx y))))

#_(ost/unstrument)