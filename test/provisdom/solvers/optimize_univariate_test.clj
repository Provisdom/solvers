(ns provisdom.solvers.optimize-univariate-test
  (:require [clojure.test :refer :all]
            [provisdom.test.core :refer :all]
            [provisdom.math.core :as m]
            [provisdom.solvers.optimize-univariate :as opt-uni]
            [clojure.spec.test.alpha :as st]
            [orchestra.spec.test :as ost]))

;18 seconds

(set! *warn-on-reflection* true)

(ost/instrument)

(deftest optimize-univariate-test
  (is (spec-check opt-uni/optimize-univariate))
  (is= {::opt-uni/point 0.6299606261931993
        ::opt-uni/value -3.472470393710553}
       (opt-uni/optimize-univariate
         {::opt-uni/univariate-f    (fn [x]
                                      (- (m/pow x 4.0) x 3.0))
          ::opt-uni/finite-interval [-10.0 10.0]
          ::opt-uni/guess           0.0}))
  (is= {::opt-uni/point 1.9554987882436676
        ::opt-uni/value 3.07905097874895}
       (opt-uni/optimize-univariate
         {::opt-uni/univariate-f    (fn [x]
                                      (- (* 5.5 x) (m/pow x 2.3) 3.0))
          ::opt-uni/finite-interval [-10.0 10.0]
          ::opt-uni/guess           0.0}
         {::opt-uni/goal :max})))

#_(ost/unstrument)