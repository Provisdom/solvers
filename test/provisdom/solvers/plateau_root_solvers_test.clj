(ns provisdom.solvers.plateau-root-solvers-test
  (:require [clojure.test :refer :all]
            [orchestra.spec.test :as ost]
            [provisdom.test.core :refer :all]
            [provisdom.math.core :as m]
            [provisdom.math.intervals :as intervals]
            [provisdom.solvers.plateau-root-solvers :as plateau]))

;3 seconds

(set! *warn-on-reflection* true)

(ost/instrument)

(deftest plateau-root-solver-test
  (is (spec-check plateau/plateau-root-solver))
  (is= #::plateau{:lower-plateau 5
                  :lower-value   5.0
                  :upper-plateau 6
                  :upper-value   6.25}
    (plateau/plateau-root-solver
      {::plateau/interval  [0.0 10.0]
       ::plateau/plateau-f (fn [x]
                             [(long (intervals/bound-by-interval
                                      [m/min-long m/max-long]
                                      x))
                              (> x 5.1)])}))
  (is= #::plateau{:lower-plateau 4
                  :lower-value   4.375
                  :upper-plateau 5
                  :upper-value   5.0}
    (plateau/plateau-root-solver
      {::plateau/interval  [0.0 10.0]
       ::plateau/plateau-f (fn [x]
                             [(long (intervals/bound-by-interval
                                      [m/min-long m/max-long]
                                      x))
                              (> x 4.9)])})))

(deftest multi-dimensional-plateau-root-solver-test
  (is (spec-check plateau/multi-dimensional-plateau-root-solver))
  (is= {::plateau/plateau-solutions [#::plateau{:lower-plateau 5
                                                :lower-value   5.0
                                                :upper-plateau 6
                                                :upper-value   6.25}]}
    (plateau/multi-dimensional-plateau-root-solver
      {::plateau/intervals       [[0.0 10.0]]
       ::plateau/multi-plateau-f (fn [v]
                                   (if (empty? v)
                                     []
                                     (let [x (first v)
                                           x (intervals/bound-by-interval
                                               [m/min-long m/max-long]
                                               x)]
                                       [[(long x) (> x 5.1)]])))}))
  (is= {::plateau/plateau-solutions [#::plateau{:lower-plateau 4
                                                :lower-value   4.375
                                                :upper-plateau 5
                                                :upper-value   5.0}]}
    (plateau/multi-dimensional-plateau-root-solver
      {::plateau/intervals       [[0.0 10.0]]
       ::plateau/multi-plateau-f (fn [v]
                                   (if (empty? v)
                                     []
                                     (let [x (first v)
                                           x (intervals/bound-by-interval
                                               [m/min-long m/max-long]
                                               x)]
                                       [[(long x) (> x 4.9)]])))}))
  (is= {::plateau/plateau-solutions [#::plateau{:lower-plateau 7
                                                :lower-value   7.5
                                                :upper-plateau 8
                                                :upper-value   8.75}
                                     #::plateau{:lower-plateau 5
                                                :lower-value   5.0
                                                :upper-plateau 6
                                                :upper-value   6.25}]}
    (plateau/multi-dimensional-plateau-root-solver
      {::plateau/intervals       [[0.0 10.0] [0.0 10.0]]
       ::plateau/multi-plateau-f (fn [v]
                                   (if (= 2 (count v))
                                     (let [[x y] v
                                           x (intervals/bound-by-interval
                                               [m/min-long m/max-long]
                                               x)
                                           y (intervals/bound-by-interval
                                               [m/min-long m/max-long]
                                               y)]
                                       [[(long x) (> (- x (* 0.5 y)) 4.9)]
                                        [(long y) (> y 5.1)]])
                                     []))})))

#_(ost/unstrument)