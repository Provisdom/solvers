(ns provisdom.solvers.internal-apache-solvers-test
  (:require
    [clojure.test :refer :all]
    [provisdom.test.core :refer :all]
    [provisdom.solvers.internal-apache-solvers :as apache-solvers]
    [provisdom.utility-belt.anomalies :as anomalies]
    [provisdom.math.core :as m]
    [provisdom.math.matrix :as mx]
    [provisdom.math.tensor :as tensor]
    [provisdom.math.random :as random]
    [clojure.spec.test.alpha :as st]
    [orchestra.spec.test :as ost]))

;;130 seconds

(set! *warn-on-reflection* true)

(ost/instrument)

;;ROOT SOLVER TESTS
(def fgb0 {::apache-solvers/univariate-f    identity
           ::apache-solvers/initial-guess   -5.0
           ::apache-solvers/finite-interval [-100.0 100.0]})

(def fgb1 {::apache-solvers/univariate-f    (fn [x]
                                              (- (m/sq x) (* 20.0 x)))
           ::apache-solvers/initial-guess   -5.0
           ::apache-solvers/finite-interval [-10.0 10.0]})

(def fgb2 {::apache-solvers/univariate-f    (fn [x]
                                              (- (m/exp x) (* 2000.0 x)))
           ::apache-solvers/initial-guess   -3.0
           ::apache-solvers/finite-interval [-5.0 5.0]})

(def fgb3 {::apache-solvers/univariate-f    (comp inc m/sq)
           ::apache-solvers/initial-guess   -5.0
           ::apache-solvers/finite-interval [-10.0 10.0]})

(deftest root-solver-test
  (is (spec-check apache-solvers/root-solver))
  ;;test 0
  (is= 0.0 (apache-solvers/root-solver fgb0))
  ;;test 1 -- 0.0
  (is= -4.992965308883406E-7
       (apache-solvers/root-solver
         fgb1 {::apache-solvers/root-solver-type :brent}))
  (is= -2.980232238769531E-7
       (apache-solvers/root-solver
         fgb1 {::apache-solvers/root-solver-type :bisection}))
  (is= 0.0
       (apache-solvers/root-solver
         fgb1 {::apache-solvers/root-solver-type :bracketing-nth-order-brent}))
  (is= 0.0
       (apache-solvers/root-solver
         fgb1 {::apache-solvers/root-solver-type :illinois}))
  (is= 0.0
       (apache-solvers/root-solver
         fgb1 {::apache-solvers/root-solver-type :muller}))
  (is= 2.0E-6
       (apache-solvers/root-solver
         fgb1 {::apache-solvers/root-solver-type :muller2}))
  (is= -2.7947626310537565E-16
       (apache-solvers/root-solver
         fgb1 {::apache-solvers/root-solver-type :pegasus}))
  (is= 4.441642203583786E-17
       (apache-solvers/root-solver
         fgb1 {::apache-solvers/root-solver-type :regula-falsi}))
  (is= 0.0
       (apache-solvers/root-solver
         fgb1 {::apache-solvers/root-solver-type :ridders}))
  (is= 1.1992433948198355E-15
       (apache-solvers/root-solver
         fgb1 {::apache-solvers/root-solver-type :secant}))
  (is= 7.812805199624811E-4
       (apache-solvers/root-solver
         fgb1 {::apache-solvers/root-solver-type :newton-raphson}))
  ;;test 2 -- 5.002501876668296E-4
  (is= 5.002501514765274E-4
       (apache-solvers/root-solver
         fgb2 {::apache-solvers/root-solver-type :brent}))
  (is= 5.003809928894043E-4
       (apache-solvers/root-solver
         fgb2 {::apache-solvers/root-solver-type :bisection}))
  (is= 5.002502003152998E-4
       (apache-solvers/root-solver
         fgb2 {::apache-solvers/root-solver-type :bracketing-nth-order-brent}))
  (is= 5.002501846562727E-4
       (apache-solvers/root-solver
         fgb2 {::apache-solvers/root-solver-type :illinois}))
  (is= 5.002501876668296E-4
       (apache-solvers/root-solver
         fgb2 {::apache-solvers/root-solver-type :muller}))
  (is= 5.002501876668296E-4
       (apache-solvers/root-solver
         fgb2 {::apache-solvers/root-solver-type :muller2}))
  (is= 5.002501876668296E-4
       (apache-solvers/root-solver
         fgb2 {::apache-solvers/root-solver-type :pegasus}))
  (is= 5.002501876668299E-4
       (apache-solvers/root-solver
         fgb2 {::apache-solvers/root-solver-type :regula-falsi}))
  (is= 5.002501876668296E-4
       (apache-solvers/root-solver
         fgb2 {::apache-solvers/root-solver-type :ridders}))
  (is= 5.002501876675032E-4
       (apache-solvers/root-solver
         fgb2 {::apache-solvers/root-solver-type :secant}))
  (is= 5.04165457653621E-4
       (apache-solvers/root-solver
         fgb2 {::apache-solvers/root-solver-type :newton-raphson}))
  ;;test 3
  (is= {::anomalies/message  (str "function values at endpoints do not have different signs, "
                                  "endpoints: [-10, 10], values: [101, 101]")
        ::anomalies/fn       #'provisdom.solvers.internal-apache-solvers/root-solver
        ::anomalies/category ::anomalies/third-party}
       (apache-solvers/root-solver fgb3)))

;;OPTIMIZE UNIVARIATE TESTS
(defn uni-fn1
  [x]
  (m/sq (inc x)))

(deftest optimize-univariate-test
  (is (spec-check apache-solvers/optimize-univariate))
  (is= {::apache-solvers/value 0.0
        ::apache-solvers/point -1.0}
       (apache-solvers/optimize-univariate uni-fn1 [-5 5] 0.0)))

;;LINEAR PROGRAMMING TESTS
(deftest linear-programming-test
  (is (spec-check apache-solvers/linear-programming))
  (is= {::apache-solvers/value        -14.0
        ::apache-solvers/vector-point [4.0 -1.0 -3.0]}
       (apache-solvers/linear-programming [1 3 5]
                                          [[[1 2 0] :eq 2]
                                           [[0 1 0] :geq -1]
                                           [[0 0 1] :geq -3]])))

;;ITERATIVE LINEAR LEAST SQUARES
(deftest iterative-linear-least-squares-test
  (is (spec-check apache-solvers/iterative-linear-least-squares
                  {:coll-check-limit 10
                   :coll-error-limit 10
                   :fspec-iterations 10
                   :recursion-limit  1
                   :test-check       {:num-tests 200}}))
  (is= [1.3333333333333333 0.33333333333333326]
       (apache-solvers/iterative-linear-least-squares
         [[1 2] [2 1]]
         [2 3]))
  (is= [5.999999999999998 1.9999999999999984]
       (apache-solvers/iterative-linear-least-squares
         [[1.0 0.5] [0.5 3.0]]
         [7.0 9.0]))
  (is= [6.26666666666666 1.4666666666666666]
       (apache-solvers/iterative-linear-least-squares
         [[1.0 0.5] [0.5 4.0]]
         [7.0 9.0])))

;;NONLINEAR LEAST SQUARES TESTS
(defn- constraints-fn
  [da]
  (if (= (count da) 2)
    (let [[a b] da]
      [(- a b 0.5)])
    [m/nan]))

(defn- jacobian-fn
  [da]
  (when (= (count da) 2)
    [[1.0 -1.0]]))

(defn- constraints-fn2
  [da]
  (if (= (count da) 2)
    (let [[a b] da]
      [(- a b 0.5)
       (- (m/sq a) (m/sq b) 0.5)])
    [m/nan m/nan]))

(defn- jacobian-fn2
  [da]
  (when (= (count da) 2)
    (let [[a b] da]
      [[1.0 -1.0]
       [(* 2 a) (* 2 b)]])))

(defn- constraints-fn3
  [da]
  (if (= (count da) 2)
    (let [[a b] da]
      [(- a b 0.5)
       (- (m/sq a) (m/sq b) 0.5)
       a])
    [m/nan m/nan m/nan]))

(defn- jacobian-fn3
  [da]
  (when (= (count da) 2)
    (let [[a b] da]
      [[1.0 -1.0]
       [(* 2 a) (* 2 b)]
       [1.0 0.0]])))

(def vars-guess [1.0 1.0])

(deftest nonlinear-least-squares-test
  (is (spec-check apache-solvers/nonlinear-least-squares))
  (is= {::apache-solvers/weighted-constraint-errors [0.0]   ;under-constrained
        ::apache-solvers/vector-point               [1.5 1.0]}
       (apache-solvers/nonlinear-least-squares
         {::apache-solvers/constraint-jacobian-fn jacobian-fn
          ::apache-solvers/constraints-fn         constraints-fn
          ::apache-solvers/vars-guess             vars-guess}))
  (is= {::apache-solvers/weighted-constraint-errors [0.0 -6.801038283654748E-7]
        ::apache-solvers/vector-point               [0.7500006801038284 0.25000068010382837]}
       (apache-solvers/nonlinear-least-squares
         {::apache-solvers/constraint-jacobian-fn jacobian-fn2
          ::apache-solvers/constraints-fn         constraints-fn2
          ::apache-solvers/vars-guess             vars-guess}))
  (is= {::apache-solvers/weighted-constraint-errors [0.0 -3.419459554209554E-7]
        ::apache-solvers/vector-point               [0.7500006243050441 0.2500006243050441]}
       (apache-solvers/nonlinear-least-squares
         {::apache-solvers/constraint-jacobian-fn jacobian-fn2
          ::apache-solvers/constraints-fn         constraints-fn2
          ::apache-solvers/vars-guess             vars-guess}
         {::apache-solvers/constraint-weights [0.7 0.3]}))
  (is= {::apache-solvers/weighted-constraint-errors [0.5 0.5 -1.0]
        ::apache-solvers/vector-point               [1.0 1.0]} ;over-constrained
       (apache-solvers/nonlinear-least-squares
         {::apache-solvers/constraint-jacobian-fn jacobian-fn3
          ::apache-solvers/constraints-fn         constraints-fn3
          ::apache-solvers/vars-guess             vars-guess}))
  (is= {::apache-solvers/weighted-constraint-errors [6.437267259157873E-7 0.1581124310822946 -0.22360782441803675]
        ::apache-solvers/vector-point               [0.5000022956995676 2.839747948836489E-6]} ;over-constrained
       (apache-solvers/nonlinear-least-squares
         {::apache-solvers/constraint-jacobian-fn jacobian-fn3
          ::apache-solvers/constraints-fn         constraints-fn3
          ::apache-solvers/vars-guess             vars-guess}
         {::apache-solvers/constraint-weights [1.4 0.4 0.2]})))

;;NONLINEAR PROGRAMMING TESTS
(defn obj-fn1
  [da]
  (if (= (count da) 2)
    (let [[a b] da]
      (+ (m/sq (inc a)) (m/sq (dec b))))
    m/inf+))

(defn grad-fn1
  [da]
  (if (= (count da) 2)
    (let [[a b] da]
      [(* 2 (inc a)) (* 2 (dec b))])
    [m/inf+ m/inf+]))

(deftest optimize-without-constraints-and-with-gradient-test
  (is (spec-check apache-solvers/optimize-without-constraints-and-with-gradient))
  (is= {::apache-solvers/value        0.0
        ::apache-solvers/vector-point [-1.0 1.0]}
       (apache-solvers/optimize-without-constraints-and-with-gradient
         {::apache-solvers/objective  obj-fn1
          ::apache-solvers/gradient   grad-fn1
          ::apache-solvers/vars-guess [0.5 0.5]}))
  (is= {::apache-solvers/value        0.0
        ::apache-solvers/vector-point [-1.0 1.0]}
       (apache-solvers/optimize-without-constraints-and-with-gradient
         {::apache-solvers/objective  #(- (obj-fn1 %))
          ::apache-solvers/gradient   grad-fn1
          ::apache-solvers/vars-guess [0.5 0.5]}
         {::apache-solvers/goal :max})))

(deftest optimize-without-constraints-or-gradient-test
  (is (spec-check apache-solvers/optimize-without-constraints-or-gradient
                  {:coll-check-limit 10
                   :coll-error-limit 10
                   :fspec-iterations 10
                   :recursion-limit  1
                   :test-check       {:num-tests 300}}))
  (is= {::apache-solvers/value        0.0
        ::apache-solvers/vector-point [-1.0 1.0]}
       (apache-solvers/optimize-without-constraints-or-gradient
         obj-fn1 [0.5 0.5]))
  (is= {::apache-solvers/value        0.0
        ::apache-solvers/vector-point [-1.0 1.0]}
       (apache-solvers/optimize-without-constraints-or-gradient
         #(- (obj-fn1 %)) [0.5 0.5] {::apache-solvers/goal :max})))

(deftest optimize-cma-evolution!-test
  (is (spec-check apache-solvers/optimize-cma-evolution!))
  (is= {::apache-solvers/value        0.3125
        ::apache-solvers/vector-point [-0.5 0.75]}
       (random/bind-seed
         3
         (apache-solvers/optimize-cma-evolution!
           obj-fn1 [-0.5 -0.5] [0.75 0.75] [0.5 0.5])))
  (is= {::apache-solvers/value        -0.3125
        ::apache-solvers/vector-point [-0.5 0.75]}
       (random/bind-seed
         3
         (apache-solvers/optimize-cma-evolution!
           #(- (obj-fn1 %)) [-0.5 -0.5] [0.75 0.75] [0.5 0.5]
           {::apache-solvers/goal :max}))))

(deftest optimize-bobyqa-test
  (is (spec-check apache-solvers/optimize-bobyqa))
  (is= {::apache-solvers/value        0.3125
        ::apache-solvers/vector-point [-0.5 0.75]}
       (apache-solvers/optimize-bobyqa
         obj-fn1 [-0.5 0.5] [0.75 0.75] [0.5 0.5]))
  (is= {::apache-solvers/value        -0.3125
        ::apache-solvers/vector-point [-0.5 0.75]}
       (apache-solvers/optimize-bobyqa
         #(- (obj-fn1 %)) [-0.5 0.5] [0.75 0.75] [0.5 0.5]
         {::apache-solvers/goal :max})))

;;INTERPOLATION TESTS
(def x-vals [1 2 3 4])
(def f-vals [8 2 4 3])

(deftest interpolation-1D-test
  (is (spec-check apache-solvers/interpolation-1D))
  (let [f-1D (apache-solvers/interpolation-1D
               {::apache-solvers/x-vals x-vals
                ::apache-solvers/f-vals f-vals})]
    (is= 4.5625 (f-1D 3.5))
    (is= 17.4375 (f-1D 0.5)))
  (is= {::anomalies/message  "number of points (4)"
        ::anomalies/fn       (var apache-solvers/interpolation-1D)
        ::anomalies/category ::anomalies/third-party}
       (apache-solvers/interpolation-1D
         {::apache-solvers/x-vals x-vals
          ::apache-solvers/f-vals f-vals}
         {::apache-solvers/interpolation-1D-type :akima}))
  (is= {::anomalies/message  "4 is smaller than the minimum (5)"
        ::anomalies/fn       (var apache-solvers/interpolation-1D)
        ::anomalies/category ::anomalies/third-party}
       (apache-solvers/interpolation-1D
         {::apache-solvers/x-vals x-vals
          ::apache-solvers/f-vals f-vals}
         {::apache-solvers/period                0.5
          ::apache-solvers/interpolation-1D-type :polynomial}))
  (let [f-1D (apache-solvers/interpolation-1D
               {::apache-solvers/x-vals x-vals
                ::apache-solvers/f-vals f-vals}
               {::apache-solvers/interpolation-1D-type :linear})]
    (is= 3.5 (f-1D 3.5))
    (is= {::anomalies/message  "Returned Function: 0.5 out of [1, 4] range"
          ::anomalies/fn       (var apache-solvers/interpolation-1D)
          ::anomalies/category ::anomalies/third-party}
         (f-1D 0.5))))

(deftest interpolation-1D-using-derivatives-test
  (is (spec-check apache-solvers/interpolation-1D-using-derivatives))
  (let [f-1D-d (apache-solvers/interpolation-1D-using-derivatives
                 [[1.0 3.0 2.0]
                  [2.0 4.0 1.5]])]
    (is= 9.0 (f-1D-d 3.0))))

(def xs [2.0 4.0 5.0 6.0 7.0])
(def ys [18.0 20.0 21.0 23.0 26.0 45.0])

(defn inner-fxy
  [x y]
  (if (m/pos? y)
    (+ (m/log y)
       x
       (* (m/exp (- x))
          (+ y (m/cube y))))
    m/nan))

(def fxy
  (mx/compute-matrix (count xs)
                     (count ys)
                     (fn [r c]
                       (if (and (< r (count xs))
                                (< c (count ys)))
                         (inner-fxy (nth xs r) (nth ys c))
                         m/nan))))

(deftest interpolation-2D-test
  (is (spec-check apache-solvers/interpolation-2D))
  (let [f2d-t (apache-solvers/interpolation-2D
                {::apache-solvers/x-vals   xs
                 ::apache-solvers/y-vals   ys
                 ::apache-solvers/f-matrix fxy}
                {::apache-solvers/smooth? true})]
    (is= {::anomalies/message  "Returned Function: point 3.00000, 20.00000 is not valid."
          ::anomalies/fn       #'provisdom.solvers.internal-apache-solvers/interpolation-2D
          ::anomalies/category ::anomalies/third-party}
         (f2d-t 3 20))
    (is= 264.1818628669024 (f2d-t 4.1 26.0))
    (is {::anomalies/message  "Returned Function: point 2.10000, 26.00000 is not valid."
         ::anomalies/fn       #'provisdom.solvers.internal-apache-solvers/interpolation-2D
         ::anomalies/category ::anomalies/third-party}
        (f2d-t 2.1 20.0)))
  (let [f2d-f (apache-solvers/interpolation-2D
                {::apache-solvers/x-vals   xs
                 ::apache-solvers/y-vals   ys
                 ::apache-solvers/f-matrix fxy}
                {::apache-solvers/smooth? false})]
    (is= 496.6707017028581 (f2d-f 3 20))
    (is= 285.9697557586808 (f2d-f 4.1 26.0))))

(def xv [0 2 4 6 8 10])
(def yv [-25 -15 -5 5 15 25])
(def zv [-10 -8 -6 -4 -2 0])
(def f-xyz
  (tensor/compute-tensor [6 6 6]
                         (fn [indices]
                           (if (= (count indices) 3)
                             (let [[x y z] indices]
                               (if (and (<= 0 x 6)
                                        (<= 0 y 6)
                                        (<= 0 z 6))
                                 (+ (get xv x m/nan)
                                    (get yv y m/nan)
                                    (get zv z m/nan))
                                 m/nan))
                             m/nan))))

(deftest interpolation-3D-test
  (is (spec-check apache-solvers/interpolation-3D))
  (let [f3d (apache-solvers/interpolation-3D
              {::apache-solvers/x-vals   xv
               ::apache-solvers/y-vals   yv
               ::apache-solvers/z-vals   zv
               ::apache-solvers/f-tensor f-xyz})]
    (is= 0.5234375 (f3d 1.5 2.5 -3.5))
    (is= -35.0 (f3d 0 -25 -10))
    (is= 35.0 (f3d 10 25 0))
    (is= -33.0 (f3d 2 -25 -10))
    (is= -25.0 (f3d 0 -15 -10))
    (is= -33.0 (f3d 0 -25 -8))
    (is= -25.0 (f3d 10 -25 -10))
    (is= 15.0 (f3d 0 25 -10))
    (is= -25.0 (f3d 0 -25 0))))

(def i-nd [[0 0 0 0] [4 1 2 3] [9 9 9 9]])
(def f-nd [0 5 9])

(deftest interpolation-ND-microsphere$-test
  (is (spec-check apache-solvers/interpolation-ND-microsphere$))
  (let [fnd (apache-solvers/interpolation-ND-microsphere$
              {::apache-solvers/i-matrix i-nd
               ::apache-solvers/f-vals   f-nd})]
    (comment "Uncontrollably randomized"
             (is= 0.9906283925553352 (fnd [1 1 1 1]))
             (is= 4.890693642100985 (fnd [10 -1 1 1])))))

#_(ost/unstrument)