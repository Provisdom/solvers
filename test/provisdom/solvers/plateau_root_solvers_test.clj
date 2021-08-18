(ns provisdom.solvers.plateau-root-solvers-test
  (:require [clojure.test :refer :all]
            [orchestra.spec.test :as ost]
            [provisdom.test.core :refer :all]
            [provisdom.math.core :as m]
            [provisdom.math.intervals :as intervals]
            [provisdom.solvers.plateau-root-solvers :as plateau]))

;96 seconds

(set! *warn-on-reflection* true)

(ost/instrument)

;;;UNIVARIATE
(deftest plateau-update-test
  (is (spec-check plateau/plateau-update))
  (is= {::plateau/lower {::plateau/plateau 2
                         ::plateau/value   2.0}}
    (plateau/plateau-update
      {::plateau/plateau-f        (fn [x]
                                    [(long (intervals/bound-by-interval
                                             [m/min-long m/max-long]
                                             x))
                                     (>= x 6.0)])
       ::plateau/plateau-solution {::plateau/lower {::plateau/plateau 0
                                                    ::plateau/value   0.0}}
       ::plateau/value            2.0})))

(deftest plateau-root-solver-test
  (is (spec-check plateau/plateau-root-solver))
  (is= {::plateau/lower {::plateau/plateau 5
                         ::plateau/value   5.0}
        ::plateau/upper {::plateau/plateau 6
                         ::plateau/value   6.25}}
    (plateau/plateau-root-solver
      {::plateau/interval  [0.0 10.0]
       ::plateau/plateau-f (fn [x]
                             [(long (intervals/bound-by-interval
                                      [m/min-long m/max-long]
                                      x))
                              (>= x 6.0)])}))
  (is= {::plateau/lower {::plateau/plateau 4
                         ::plateau/value   4.375}
        ::plateau/upper {::plateau/plateau 5
                         ::plateau/value   5.0}}
    (plateau/plateau-root-solver
      {::plateau/interval  [0.0 10.0]
       ::plateau/plateau-f (fn [x]
                             [(long (intervals/bound-by-interval
                                      [m/min-long m/max-long]
                                      x))
                              (>= x 5.0)])}))
  (is= {::plateau/upper {::plateau/plateau 0
                         ::plateau/value   0.0}}
    (plateau/plateau-root-solver
      {::plateau/interval  [0.0 10.0]
       ::plateau/plateau-f (fn [x]
                             [(long (intervals/bound-by-interval
                                      [m/min-long m/max-long]
                                      x))
                              (>= x -1.0)])}))
  (is= {::plateau/lower {::plateau/plateau 10
                         ::plateau/value   10.0}}
    (plateau/plateau-root-solver
      {::plateau/interval  [0.0 10.0]
       ::plateau/plateau-f (fn [x]
                             [(long (intervals/bound-by-interval
                                      [m/min-long m/max-long]
                                      x))
                              (>= x 11.0)])})))

(deftest tighten-plateau-solution-test
  (is (spec-check plateau/tighten-plateau-solution))
  (is= {::plateau/lower {::plateau/plateau 5
                         ::plateau/value   5.999999642372131}
        ::plateau/upper {::plateau/plateau 6
                         ::plateau/value   6.000000238418579}}
    (plateau/tighten-plateau-solution
      {::plateau/lower-bound      5.0
       ::plateau/plateau-f        (fn [x]
                                    [(long (intervals/bound-by-interval
                                             [m/min-long m/max-long]
                                             x))
                                     (>= x 6.0)])
       ::plateau/plateau-solution {::plateau/lower {::plateau/plateau 5
                                                    ::plateau/value   5.0}
                                   ::plateau/upper {::plateau/plateau 6
                                                    ::plateau/value   6.25}}}))
  (is= {::plateau/lower {::plateau/plateau 4
                         ::plateau/value   4.999999403953552}
        ::plateau/upper {::plateau/plateau 5
                         ::plateau/value   5.0}}
    (plateau/tighten-plateau-solution
      {::plateau/lower-bound      4.375
       ::plateau/plateau-f        (fn [x]
                                    [(long (intervals/bound-by-interval
                                             [m/min-long m/max-long]
                                             x))
                                     (>= x 5.0)])
       ::plateau/plateau-solution {::plateau/lower {::plateau/plateau 4
                                                    ::plateau/value   4.375}
                                   ::plateau/upper {::plateau/plateau 5
                                                    ::plateau/value   5.0}}})))

;;;MULTIVARIATE
(deftest single-pass-test
  (is (spec-check plateau/single-pass))
  (is= {::plateau/last-values  [1.0]
        ::plateau/lower-values [1.0]}
    (plateau/single-pass
      #::plateau{:abs-accu        1e-1
                 :guesses         [0.7]
                 :last-values     [0.6]
                 :lower-values    [0.4]
                 :max-iter        3
                 :multi-plateau-f (fn [v]
                                    (if (empty? v)
                                      []
                                      (let [x (first v)
                                            x (intervals/bound-by-interval
                                                [m/min-long m/max-long]
                                                x)]
                                        [[(long x) (>= x 6.0)]])))
                 :upper-bounds    [1.0]})))

(deftest multi-dimensional-plateau-root-solver-test
  (is (spec-check plateau/multi-dimensional-plateau-root-solver))
  (is= {::plateau/plateau-solutions
        [{::plateau/lower {::plateau/plateau 5
                           ::plateau/value   5.0}
          ::plateau/upper {::plateau/plateau 6
                           ::plateau/value   6.25}}]
        ::plateau/single-pass-solutions
        [{::plateau/last-values  [6.25]
          ::plateau/lower-values [5.0]}
         {::plateau/last-values  [6.25]
          ::plateau/lower-values [5.0]}]}
    (plateau/multi-dimensional-plateau-root-solver
      {::plateau/intervals       [[0.0 10.0]]
       ::plateau/multi-plateau-f (fn [v]
                                   (if (empty? v)
                                     []
                                     (let [x (first v)
                                           x (intervals/bound-by-interval
                                               [m/min-long m/max-long]
                                               x)]
                                       [[(long x) (>= x 6.0)]])))}))
  (is= {::plateau/plateau-solutions
        [{::plateau/lower {::plateau/plateau 4
                           ::plateau/value   4.375}
          ::plateau/upper {::plateau/plateau 5
                           ::plateau/value   5.0}}]
        ::plateau/single-pass-solutions
        [{::plateau/last-values  [5.0]
          ::plateau/lower-values [4.375]}
         {::plateau/last-values  [5.0]
          ::plateau/lower-values [4.375]}]}
    (plateau/multi-dimensional-plateau-root-solver
      {::plateau/intervals       [[0.0 10.0]]
       ::plateau/multi-plateau-f (fn [v]
                                   (if (empty? v)
                                     []
                                     (let [x (first v)
                                           x (intervals/bound-by-interval
                                               [m/min-long m/max-long]
                                               x)]
                                       [[(long x) (>= x 5.0)]])))}))
  (is= {::plateau/plateau-solutions
        [{::plateau/lower {::plateau/plateau 4
                           ::plateau/value   7.5}
          ::plateau/upper {::plateau/plateau 5
                           ::plateau/value   8.75}}
         {::plateau/lower {::plateau/plateau 5
                           ::plateau/value   5.0}
          ::plateau/upper {::plateau/plateau 6
                           ::plateau/value   6.25}}]
        ::plateau/single-pass-solutions
        [{::plateau/last-values  [5.0 6.25]
          ::plateau/lower-values [4.375 5.0]}
         {::plateau/last-values  [8.75 6.25]
          ::plateau/lower-values [7.5 5.0]}
         {::plateau/last-values  [8.75 6.25]
          ::plateau/lower-values [7.5 5.0]}]}
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
                                       [[(long (- x (* 0.5 y)))
                                         (>= (- x (* 0.5 y)) 5.0)]
                                        [(long y) (>= y 6.0)]])
                                     []))}
      {::plateau/partition-size 1})))

(deftest tighten-multi-plateau-solution-test
  (is (spec-check plateau/tighten-multi-plateau-solution))
  (is= {::plateau/plateau-solutions
        [{::plateau/lower {::plateau/plateau 5
                           ::plateau/value   5.999999642372131}
          ::plateau/upper {::plateau/plateau 6
                           ::plateau/value   6.000000238418579}}]
        ::plateau/single-pass-solutions
        [{::plateau/last-values  [6.25]
          ::plateau/lower-values [5.0]}
         {::plateau/last-values  [6.25]
          ::plateau/lower-values [5.0]}]}
    (plateau/tighten-multi-plateau-solution
      {::plateau/multi-plateau-f
       (fn [v]
         (if (empty? v)
           []
           (let [x (first v)
                 x (intervals/bound-by-interval
                     [m/min-long m/max-long]
                     x)]
             [[(long x) (>= x 6.0)]])))
       ::plateau/multi-plateau-solution
       {::plateau/plateau-solutions
        [{::plateau/lower {::plateau/plateau 5
                           ::plateau/value   5.0}
          ::plateau/upper {::plateau/plateau 6
                           ::plateau/value   6.25}}]
        ::plateau/single-pass-solutions
        [{::plateau/last-values  [6.25]
          ::plateau/lower-values [5.0]}
         {::plateau/last-values  [6.25]
          ::plateau/lower-values [5.0]}]}}))
  (is= {::plateau/plateau-solutions
        [{::plateau/lower {::plateau/plateau 4
                           ::plateau/value   8.124999403953552}
          ::plateau/upper {::plateau/plateau 5
                           ::plateau/value   8.125}}
         {::plateau/lower {::plateau/plateau 5
                           ::plateau/value   5.999999642372131}
          ::plateau/upper {::plateau/plateau 6
                           ::plateau/value   6.000000238418579}}]
        ::plateau/single-pass-solutions
        [{::plateau/last-values  [5.0 6.25]
          ::plateau/lower-values [4.375 5.0]}
         {::plateau/last-values  [8.75 6.25]
          ::plateau/lower-values [7.5 5.0]}
         {::plateau/last-values  [8.75 6.25]
          ::plateau/lower-values [7.5 5.0]}]}
    (plateau/tighten-multi-plateau-solution
      {::plateau/multi-plateau-f
       (fn [v]
         (if (= 2 (count v))
           (let [[x y] v
                 x (intervals/bound-by-interval
                     [m/min-long m/max-long]
                     x)
                 y (intervals/bound-by-interval
                     [m/min-long m/max-long]
                     y)]
             [[(long (- x (* 0.5 y)))
               (>= (- x (* 0.5 y)) 5.0)]
              [(long y) (>= y 6.0)]])
           []))
       ::plateau/multi-plateau-solution
       {::plateau/plateau-solutions
        [{::plateau/lower {::plateau/plateau 4
                           ::plateau/value   7.5}
          ::plateau/upper {::plateau/plateau 5
                           ::plateau/value   8.75}}
         {::plateau/lower {::plateau/plateau 5
                           ::plateau/value   5.0}
          ::plateau/upper {::plateau/plateau 6
                           ::plateau/value   6.25}}]
        ::plateau/single-pass-solutions
        [{::plateau/last-values  [5.0 6.25]
          ::plateau/lower-values [4.375 5.0]}
         {::plateau/last-values  [8.75 6.25]
          ::plateau/lower-values [7.5 5.0]}
         {::plateau/last-values  [8.75 6.25]
          ::plateau/lower-values [7.5 5.0]}]}})))

#_(ost/unstrument)
