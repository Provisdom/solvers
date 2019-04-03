(ns provisdom.solvers.curve-fitting-test
  (:require [clojure.test :refer :all]
            [provisdom.test.core :refer :all]
            [provisdom.solvers.curve-fitting :as curve-fitting]
            [provisdom.math.core :as m]
            [provisdom.math.series :as series]
            [clojure.spec.test.alpha :as st]
            [orchestra.spec.test :as ost]))

;6 seconds

(set! *warn-on-reflection* true)

(ost/instrument)

(def points-data [[2.0 26.0] [4.0 173.0] [5.0 281.0] [6.0 322.0] [7.0 37.0]])

(def xs
  (mapv first points-data))

(def g-xs
  (mapv second points-data))

(deftest linear-least-squares-line-fitting-test
  (is (spec-check curve-fitting/linear-least-squares-line-fitting
                  {:coll-check-limit 10
                   :coll-error-limit 10
                   :fspec-iterations 10
                   :recursion-limit  1
                   :test-check       {:num-tests 200}}))
  (let [lls-map (curve-fitting/linear-least-squares-line-fitting
                  {::curve-fitting/x-vals   xs
                   ::curve-fitting/f-vals   g-xs
                   ::curve-fitting/basis-fn (fn [x]
                                              [1.0 x (m/cube x)])})]
    (is (data-approx= [-374.1841730714473 198.65888604231853 -2.7220798937280524]
                      (::curve-fitting/line-fitting-weights lls-map)))
    (is (data-approx= 265.55365784229997
                      ((::curve-fitting/line-fitting-fn lls-map) 5.5))))
  (let [lls-map-x (curve-fitting/linear-least-squares-line-fitting
                    {::curve-fitting/x-vals   (conj xs 6.0) ;can have duplicate x-vals
                     ::curve-fitting/f-vals   (conj g-xs 3.23432)
                     ::curve-fitting/basis-fn (fn [x]
                                                [1.0 x (m/cube x)])})]
    (is (data-approx= [-287.27029929684994 162.806101134645 -2.3851099982080766]
                      (::curve-fitting/line-fitting-weights lls-map-x)))
    (is (data-approx= 211.3405809918287
                      ((::curve-fitting/line-fitting-fn lls-map-x) 5.5))))
  (let [lls-map2 (curve-fitting/linear-least-squares-line-fitting
                   {::curve-fitting/x-vals   xs
                    ::curve-fitting/f-vals   g-xs
                    ::curve-fitting/basis-fn (fn [x]
                                               [1.0 (m/sq x) (m/cube x)])})]
    (is (data-approx= [-147.53975088364837 49.710219695486124 -6.462473288353114]
                      (::curve-fitting/line-fitting-weights lls-map2)))
    (is (data-approx= 281.0004015550576
                      ((::curve-fitting/line-fitting-fn lls-map2) 5.5))))
  (let [ch1-4 (curve-fitting/linear-least-squares-line-fitting
                {::curve-fitting/x-vals   xs
                 ::curve-fitting/f-vals   g-xs
                 ::curve-fitting/basis-fn (series/polynomial-fn
                                            4 {::series/chebyshev-kind 1})})]
    (is (data-approx= [-1377.7020833331776 1360.6833333331895 -267.49999999997624
                       24.316666666665014 -0.7979166666666255]
                      (::curve-fitting/line-fitting-weights ch1-4)))
    (is (data-approx= 322.47187499999745
                      ((::curve-fitting/line-fitting-fn ch1-4) 5.5))))
  (let [ch1-3 (curve-fitting/linear-least-squares-line-fitting
                {::curve-fitting/x-vals   xs
                 ::curve-fitting/f-vals   g-xs
                 ::curve-fitting/basis-fn (series/polynomial-fn
                                            3 {::series/chebyshev-kind 1})})]
    (is (data-approx=
          [1125.4528301886708 -953.7825696316196 136.6630727762796 -5.672955974842742]
          (::curve-fitting/line-fitting-weights ch1-3)))
    (is (data-approx= 329.3530997304597
                      ((::curve-fitting/line-fitting-fn ch1-3) 5.5))))
  (let [ch1-2 (curve-fitting/linear-least-squares-line-fitting
                {::curve-fitting/x-vals   xs
                 ::curve-fitting/f-vals   g-xs
                 ::curve-fitting/basis-fn (series/polynomial-fn
                                            2 {::series/chebyshev-kind 1})})]
    (is (data-approx= [-547.6597938144306 337.4874815905734 -17.7349042709867]
                      (::curve-fitting/line-fitting-weights ch1-2)))
    (is (data-approx= 253.29455081001447
                      ((::curve-fitting/line-fitting-fn ch1-2) 5.5))))
  (let [ch2-3 (curve-fitting/linear-least-squares-line-fitting
                {::curve-fitting/x-vals   xs
                 ::curve-fitting/f-vals   g-xs
                 ::curve-fitting/basis-fn (series/polynomial-fn
                                            3 {::series/chebyshev-kind 2})})]
    (is (data-approx=
          [1057.1212938004992 -474.05480682837816 68.3315363881387 -2.8364779874213335]
          (::curve-fitting/line-fitting-weights ch2-3)))
    (is (data-approx= 329.35309973045923
                      ((::curve-fitting/line-fitting-fn ch2-3) 5.5))))
  (let [p4 (curve-fitting/linear-least-squares-line-fitting
             {::curve-fitting/x-vals   xs
              ::curve-fitting/f-vals   g-xs
              ::curve-fitting/basis-fn (series/polynomial-fn 4)})]
    (is (data-approx= [-1110.999999999948 1287.733333333297 -528.6166666666575
                       97.26666666666567 -6.383333333333295]
                      (::curve-fitting/line-fitting-weights p4)))
    (is (data-approx= 322.4718750000002
                      ((::curve-fitting/line-fitting-fn p4) 5.5))))
  (let [p1-4 (curve-fitting/linear-least-squares-line-fitting
               {::curve-fitting/x-vals   xs
                ::curve-fitting/f-vals   g-xs
                ::curve-fitting/basis-fn (series/polynomial-fn
                                           4 {::series/start-degree 1})})]
    (is (data-approx=
          [116.73530443936708 -108.88714394337724 34.8383049657826 -3.0789904661636016]
          (::curve-fitting/line-fitting-weights p1-4)))
    (is (data-approx= 326.9623453676063
                      ((::curve-fitting/line-fitting-fn p1-4) 5.5)))))

(deftest b-spline-line-fitting-test
  (is (spec-check curve-fitting/b-spline-line-fitting))
  (let [b2 (curve-fitting/b-spline-line-fitting
             {::curve-fitting/distinct-x-vals xs
              ::curve-fitting/f-vals          g-xs
              ::curve-fitting/degree          2})]
    (is= 330.11053109396266 (b2 5.5))
    (is= 290.3318202451203 (b2 6.0)))
  (let [b3 (curve-fitting/b-spline-line-fitting
             {::curve-fitting/distinct-x-vals xs
              ::curve-fitting/f-vals          g-xs
              ::curve-fitting/degree          3})]
    (is= 329.35309973045815 (b3 5.5))
    (is= 306.51482479784363 (b3 6.0))))

(def points-2D
  [[[2.0 18.0] 796.0]
   [[2.0 20.0] 1090.0]
   [[2.0 21.0] 1261.0]
   [[2.0 23.0] 1654.0]
   [[2.0 26.0] 2387.0]
   [[2.0 45.0] 12344.0]
   [[4.0 18.0] 114.0]
   [[4.0 20.0] 153.0]
   [[4.0 21.0] 177.0]
   [[4.0 23.0] 230.0]
   [[4.0 26.0] 329.0]
   [[4.0 45.0] 1677.0]
   [[5.0 18.0] 47.0]
   [[5.0 20.0] 62.0]
   [[5.0 21.0] 70.0]
   [[5.0 23.0] 90.0]
   [[5.0 26.0] 126.0]
   [[5.0 45.0] 623.0]
   [[6.0 18.0] 23.0]
   [[6.0 20.0] 28.0]
   [[6.0 21.0] 32.0]
   [[6.0 23.0] 39.0]
   [[6.0 26.0] 52.0]
   [[6.0 45.0] 235.0]
   [[7.0 18.0] 15.0]
   [[7.0 20.0] 17.0]
   [[7.0 21.0] 18.0]
   [[7.0 23.0] 21.0]
   [[7.0 26.0] 26.0]
   [[7.0 45.0] 93.0]])

(def xys
  (mapv first points-2D))

(def g-xys
  (mapv second points-2D))

(deftest linear-least-squares-curve-fitting-test
  (is (spec-check curve-fitting/linear-least-squares-curve-fitting
                  {:coll-check-limit 10
                   :coll-error-limit 10
                   :fspec-iterations 10
                   :recursion-limit  1
                   :test-check       {:num-tests 200}}))
  (let [lls-z (curve-fitting/linear-least-squares-curve-fitting
                {::curve-fitting/x-matrix       [[2.0 30.0 5.0]
                                                 [3.0 40.0 6.0]
                                                 [7.0 20.0 4.0]
                                                 [6.0 15.0 8.0]
                                                 [4.0 36.0 4.0]]
                 ::curve-fitting/f-vals         [4.0 3.0 2.0 12.0 15.0]
                 ::curve-fitting/curve-basis-fn (fn [v]
                                                  (if (= 3 (count v))
                                                    (let [[x y z] v]
                                                      [1 x y z (m/sq x)])
                                                    []))})]
    (is (data-approx= 14.303501506024148
                      ((::curve-fitting/curve-fitting-fn lls-z) [5.5 28.7 4.3])))
    (is (data-approx=
          [2.858433734940034 25.911897590361473 -0.8659638554216929
           -2.4051204819277294 -3.169427710843378]
          (::curve-fitting/curve-fitting-weights lls-z))))
  (let [lls-curve (curve-fitting/linear-least-squares-curve-fitting
                    {::curve-fitting/x-matrix       xys
                     ::curve-fitting/f-vals         g-xys
                     ::curve-fitting/curve-basis-fn (fn [v]
                                                      (if (= 2 (count v))
                                                        (let [[x y] v]
                                                          [1 x y (m/cube x)])
                                                        []))})]
    (is (data-approx= 225.39971828421585
                      ((::curve-fitting/curve-fitting-fn lls-curve) [5.5 28.7])))
    (is (data-approx=
          [3576.2055568677993 -1644.3254459626135 107.97588652482267 15.591742575196873]
          (::curve-fitting/curve-fitting-weights lls-curve))))
  (let [poly-2 (curve-fitting/linear-least-squares-curve-fitting
                 {::curve-fitting/x-matrix       xys
                  ::curve-fitting/f-vals         g-xys
                  ::curve-fitting/curve-basis-fn (fn [v]
                                                   (if (= 2 (count v))
                                                     (let [[x y] v
                                                           val-v ((series/polynomial-2D-fn-by-degree 2) x y)]
                                                       val-v)
                                                     []))})]
    (is (data-approx= -261.66971859835485
                      ((::curve-fitting/curve-fitting-fn poly-2) [5.5 28.7])))
    (is (data-approx=
          [-3189.545588323287 334.41408500495845 -368.8008791074165
           2.9019124617507708 -86.23982036747992 217.56283750613647]
          (::curve-fitting/curve-fitting-weights poly-2))))
  (let [poly-3 (curve-fitting/linear-least-squares-curve-fitting
                 {::curve-fitting/x-matrix       xys
                  ::curve-fitting/f-vals         g-xys
                  ::curve-fitting/curve-basis-fn (fn [v]
                                                   (if (= 2 (count v))
                                                     (let [[x y] v
                                                           val-v ((series/polynomial-2D-fn-by-degree 3) x y)]
                                                       val-v)
                                                     []))})]
    (is (data-approx= -109.31465603752986
                      ((::curve-fitting/curve-fitting-fn poly-3) [5.5 28.7])))
    (is (data-approx=
          [-2964.644064914418 197.5044338683106 1334.2066850436406
           11.103399823769012 -202.70447843509646 174.80915849973607
           0.03301654985858392 -2.317900060917706 29.851202599925646
           -52.79533542976914]
          (::curve-fitting/curve-fitting-weights poly-3))))
  (let [poly-4 (curve-fitting/linear-least-squares-curve-fitting
                 {::curve-fitting/x-matrix       xys
                  ::curve-fitting/f-vals         g-xys
                  ::curve-fitting/curve-basis-fn (fn [v]
                                                   (if (= 2 (count v))
                                                     (let [[x y] v
                                                           val-v ((series/polynomial-2D-fn-by-degree 4) x y)]
                                                       val-v)
                                                     []))})]
    (is (data-approx= 176.26088371885817
                      ((::curve-fitting/curve-fitting-fn poly-4) [5.5 28.7])))
    (is (data-approx=
          [-2050.3968757541825 106.29722391047446 1981.8669830053382
           14.064819444200818 -211.08730725698555 -259.6457040821559
           0.14328304225982386 -7.152479556335808 76.65972517042242
           -69.41038436784879 1.3854417314597466E-4 -0.026159257745385543
           0.8018493606563802 -7.247196255726908 10.718055555555372]
          (::curve-fitting/curve-fitting-weights poly-4))))
  (let [poly-20 (curve-fitting/linear-least-squares-curve-fitting
                  {::curve-fitting/x-matrix       xys
                   ::curve-fitting/f-vals         g-xys
                   ::curve-fitting/curve-basis-fn (fn [v]
                                                    (if (= 2 (count v))
                                                      (let [[x y] v
                                                            val-v ((series/polynomial-2D-fn-by-basis-count
                                                                     20
                                                                     {::series/chebyshev-kind 1}) x y)]
                                                        val-v)
                                                      []))})]
    (is (data-approx= 123.65870528567757
                      ((::curve-fitting/curve-fitting-fn poly-20) [5.5 28.7])))
    (is (data-approx=
          [-4570.710184016991 -348.5196618375316 7126.080272273072 29.4624197555226
           -239.66381568606565 -1377.9959457490452 -0.37825629375083264
           -5.214477757191738 56.624414816925864 114.13616877081664
           0.0045612061972588325 -0.03113641892213491 0.6647388129750418
           -5.569874056570459 -3.345550489699451 -1.7560473159920388E-5
           2.0814908582451206E-5 0.001121204887739387 -0.024356456660173074
           0.18373754643701698]
          (::curve-fitting/curve-fitting-weights poly-20)))))

(deftest smoothing-cubic-spline-dof-test
  (is (spec-check curve-fitting/smoothing-cubic-spline-dof
                  {:fspec-iterations 0
                   :test-check       {:num-tests 1000}}))
  (is= 2.34173669467787
       (curve-fitting/smoothing-cubic-spline-dof
         {::curve-fitting/x-vals              [0.0 1.0 2.0 3.0]
          ::curve-fitting/variances           [1.0 1.0 1.0 1.0]
          ::curve-fitting/smoothing-parameter 0.5})))

(deftest smoothing-cubic-spline-test
  (is (spec-check curve-fitting/smoothing-cubic-spline
                  {:fspec-iterations 0
                   :test-check       {:num-tests 1000}}))
  (let [ret (curve-fitting/smoothing-cubic-spline
              {::curve-fitting/x-vals              [0.0 1.0 2.0 3.0]
               ::curve-fitting/f-vals              [1.0 4.0 25.0 81.0]
               ::curve-fitting/smoothing-parameter 0.5})
        s-f (::curve-fitting/smoothing-cubic-spline-fn ret)]
    (is= 2.34173669467787 (::curve-fitting/dof ret))
    (is= -24.587301587301592 (s-f -1.0))
    (is= -7.54341736694678 (s-f 0))
    (is= 10.92436974789916 (s-f 1))
    (is= 36.78151260504202 (s-f 2.0))
    (is= 106.5873015873016 (s-f 4.0))))

#_(ost/unstrument)