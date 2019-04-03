(ns provisdom.solvers.interpolation-test
  (:require [clojure.test :refer :all]
            [provisdom.test.core :refer :all]
            [provisdom.solvers.interpolation :as interpolation]
            [provisdom.utility-belt.anomalies :as anomalies]
            [provisdom.math.core :as m]
            [provisdom.math.matrix :as mx]
            [provisdom.math.tensor :as tensor]
            [provisdom.solvers.internal-apache-solvers :as apache-solvers]
            [clojure.spec.test.alpha :as st]
            [orchestra.spec.test :as ost]))

;43 seconds

(set! *warn-on-reflection* true)

(ost/instrument)

(defn- g
  [x]
  (+ x
     (* 2.0 (m/sq x))
     (* 3.0 (m/cube x))
     (- (m/exp x))))

(defn g'
  [x]
  (+ 1.0
     (* 4.0 x)
     (* 9.0 (m/sq x))
     (- (m/exp x))))

(defn g''
  [x]
  (+ 4.0
     (* 18.0 x)
     (- (m/exp x))))

(defn g'''
  [x]
  (+ 18.0
     (- (m/exp x))))

(defn g''''
  [x]
  (- (m/exp x)))

(defn data-fn
  [x]
  (vector x (g x) (g' x) (g'' x) (g''' x) (g'''' x)))

(def xs [2.0 4.0 5.0 6.0 7.0])
(def xs2 [2.0 4.0 5.0 6.0 7.0 8.0 9.0])
(def ys [18.0 20.0 21.0 23.0 26.0 45.0])
(def zs [-3.0 -2.0 0.0 5.0])
(def f-vals (mapv g xs))
(def f-vals2 (mapv g xs2))
(def data (mapv data-fn xs))

(def interpolation-1D-map
  {::interpolation/x-vals xs
   ::interpolation/f-vals f-vals})

(def interpolation-1D-map2
  {::interpolation/x-vals xs2
   ::interpolation/f-vals f-vals2})

(deftest g-comparison-testing
  (is= -0.9815306597126334 (g -0.5))
  (is= 49.692506039296525 (g 2.5))
  (is= 228.3578686994782 (g 4.5))
  (is= 320.4330677357796 (g 5.5))
  (is= 322.5712065072649 (g 6.0)))

(deftest interpolation-1D-quadratic-test
  (is (spec-check interpolation/interpolation-1D-quadratic))
  (let [zero-start-slope (interpolation/interpolation-1D-quadratic
                           {::interpolation/x-vals xs
                            ::interpolation/f-vals f-vals
                            ::interpolation/slope  0.0
                            ::interpolation/end-slope? false})]
    (is= {::anomalies/category ::anomalies/no-solve
          ::anomalies/message  "Out of range: -0.5"
          ::anomalies/fn       (var interpolation/interpolation-1D-quadratic)}
         (zero-start-slope -0.5))
    (is= 322.5712065072649 (zero-start-slope 6.0))
    (is= 309.22770124872096 (zero-start-slope 5.5)))
  (is= 322.5712065072649
       ((interpolation/interpolation-1D-quadratic
          {::interpolation/x-vals xs
           ::interpolation/f-vals f-vals
           ::interpolation/slope  2.0
           ::interpolation/end-slope? false})
         6.0))
  (is= 455.4272975726662
       ((interpolation/interpolation-1D-quadratic
          {::interpolation/x-vals xs
           ::interpolation/f-vals f-vals
           ::interpolation/slope  2.0
           ::interpolation/end-slope? true})
         5.5)))

(deftest interpolation-1D-cubic-clamped-test
  (is (spec-check interpolation/interpolation-1D-cubic-clamped))
  (is= 331.9882749613879
       ((interpolation/interpolation-1D-cubic-clamped
          {::interpolation/x-vals             xs
           ::interpolation/f-vals             f-vals
           ::interpolation/start-acceleration 0.0
           ::interpolation/end-acceleration   0.0})
         5.5))
  (is= 322.5712065072649
       ((interpolation/interpolation-1D-cubic-clamped
          {::interpolation/x-vals             xs
           ::interpolation/f-vals             f-vals
           ::interpolation/start-acceleration 0.0
           ::interpolation/end-acceleration   0.0})
         6.0))
  (is= 57.93523103667523
       ((interpolation/interpolation-1D-cubic-clamped
          {::interpolation/x-vals             xs
           ::interpolation/f-vals             f-vals
           ::interpolation/start-acceleration 2.0
           ::interpolation/end-acceleration   10.0})
         2.5))
  (is= 224.51782190079774
       ((interpolation/interpolation-1D-cubic-clamped
          {::interpolation/x-vals             xs
           ::interpolation/f-vals             f-vals
           ::interpolation/start-acceleration 1.0
           ::interpolation/end-acceleration   20.0})
         4.5)))

;;only :cubic-hermite, :polynomial, :neville, and :cubic-closed extrapolate
(deftest interpolation-1D-test
  (is (spec-check interpolation/interpolation-1D))
  (let [cubic (interpolation/interpolation-1D
                interpolation-1D-map)]
    (is= {::anomalies/category ::anomalies/third-party
          ::anomalies/message  "Returned Function: -0.5 out of [2, 7] range"
          ::anomalies/fn       (var apache-solvers/interpolation-1D)}
         (cubic -0.5))
    (is= 322.5712065072649 (cubic 6.0))
    (is= 331.98827496138784 (cubic 5.5)))
  (let [cubic-closed (interpolation/interpolation-1D
                       interpolation-1D-map
                       {::interpolation/interpolation-1D-type :cubic-closed})]
    (is= 4963.942526378081 (cubic-closed -0.5))
    (is= 322.5712065072649 (cubic-closed 6.0))
    (is= 341.24773908406263 (cubic-closed 5.5)))
  (let [cubic-hermite (interpolation/interpolation-1D
                        interpolation-1D-map
                        {::interpolation/interpolation-1D-type :cubic-hermite})]
    (is= -279.18465785267523 (cubic-hermite -0.5))
    (is= 322.5712065072649 (cubic-hermite 6.0))
    (is= 326.6658584439873 (cubic-hermite 5.5)))
  (let [akima (interpolation/interpolation-1D
                interpolation-1D-map
                {::interpolation/interpolation-1D-type :akima})]
    (is= {::anomalies/category ::anomalies/third-party
          ::anomalies/message  "Returned Function: -0.5 out of [2, 7] range"
          ::anomalies/fn       (var apache-solvers/interpolation-1D)}
         (akima -0.5))
    (is= 322.5712065072649 (akima 6.0))
    (is= 330.05633400646775 (akima 5.5)))
  (let [linear (interpolation/interpolation-1D
                 interpolation-1D-map
                 {::interpolation/interpolation-1D-type :linear})]
    (is= {::anomalies/category ::anomalies/third-party
          ::anomalies/message  "Returned Function: -0.5 out of [2, 7] range"
          ::anomalies/fn       (var apache-solvers/interpolation-1D)}
         (linear -0.5))
    (is= 322.5712065072649 (linear 6.0))
    (is= 302.0790237023441 (linear 5.5)))
  (let [neville (interpolation/interpolation-1D
                  interpolation-1D-map
                  {::interpolation/interpolation-1D-type :neville})]
    (is= -1890.457907360669 (neville -0.5))
    (is= 322.5712065072649 (neville 6.0))
    (is= 323.080992182385 (neville 5.5)))
  (let [polynomial (interpolation/interpolation-1D
                     interpolation-1D-map
                     {::interpolation/interpolation-1D-type :polynomial})]
    (is= -1890.457907360669 (polynomial -0.5))
    (is= 322.5712065072649 (polynomial 6.0))
    (is= 323.080992182385 (polynomial 5.5)))
  ;;loess requires more points but can be more accurate
  (is= {::anomalies/category ::anomalies/third-party
        ::anomalies/message  "bandwidth (1)"
        ::anomalies/fn       (var apache-solvers/interpolation-1D)}
       (interpolation/interpolation-1D
         interpolation-1D-map
         {::interpolation/interpolation-1D-type :loess}))
  (let [loess (interpolation/interpolation-1D
                interpolation-1D-map2
                {::interpolation/interpolation-1D-type :loess})]
    (is= {::anomalies/category ::anomalies/third-party
          ::anomalies/message  "Returned Function: -0.5 out of [2, 9] range"
          ::anomalies/fn       (var apache-solvers/interpolation-1D)}
         (loess -0.5))
    (is= 322.5712065072649 (loess 6.0))
    (is= 327.5922578043044 (loess 5.5)))
  (is= 331.98827496138784
       ((interpolation/interpolation-1D
          interpolation-1D-map
          {::interpolation/interpolation-1D-type :loess
           ::interpolation/bandwidth-for-loess   0.5})
         5.5))
  ;;can't use a period that results in repeated x-vals
  (is= {::anomalies/category ::anomalies/third-party
        ::anomalies/message  "points 0 and 1 are not strictly increasing (-0.5 >= -0.5)"
        ::anomalies/fn       (var apache-solvers/interpolation-1D)}
       (interpolation/interpolation-1D
         interpolation-1D-map
         {::interpolation/interpolation-1D-type :cubic
          ::interpolation/period                0.5}))
  ;;can't use a period that results in repeated x-vals
  (is= {::anomalies/category ::anomalies/third-party
        ::anomalies/message  "points 0 and 1 are not strictly increasing (-3 >= -3)"
        ::anomalies/fn       (var apache-solvers/interpolation-1D)}
       (interpolation/interpolation-1D
         interpolation-1D-map
         {::interpolation/period                3.0}))
  (is= 105.43512756101777
       ((interpolation/interpolation-1D
          interpolation-1D-map
          {::interpolation/period                0.44})
         -0.5))
  ;;two points close together (e.g., 5.0 and 5.0001) creates huge slope
  (is= -355576.843194702
       ((interpolation/interpolation-1D
          interpolation-1D-map
          {::interpolation/period                3.0001})
         -0.5))
  (is= -130.88605099330766
       ((interpolation/interpolation-1D
          interpolation-1D-map
          {::interpolation/period                8.5})
         -0.5))
  (is= 281.5868408974234
       ((interpolation/interpolation-1D
          interpolation-1D-map
          {::interpolation/period                5.5})
         -0.5))
  (is= 322.5712065072649
       ((interpolation/interpolation-1D
          interpolation-1D-map
          {::interpolation/period                5.5})
         6.0))
  (is= 341.15899851629416
       ((interpolation/interpolation-1D
          interpolation-1D-map
          {::interpolation/period                5.5})
         5.5)))

(deftest interpolation-1D-using-derivatives-test
  (is (spec-check interpolation/interpolation-1D-using-derivatives))
  (let [interpolation' (interpolation/interpolation-1D-using-derivatives data)]
    (is= -4.737638660160948 (interpolation' -0.5))
    (is= 49.69250604003082 (interpolation' 2.5))
    (is= 320.4330677357795 (interpolation' 5.5))
    (is= 322.57120650726495 (interpolation' 6.0))))

;;2D
(defn f-xy
  [x y]
  (+ (m/log y)
     x
     (* (m/exp (- x))
        (+ y (m/cube y)))))

(def f-matrix
  (mx/compute-matrix
    (count xs)
    (count ys)
    (fn [r c]
      (let [x (get xs r m/nan)
            y (get ys c m/nan)]
        (f-xy x y)))))

(def interpolation-2D-map
  {::interpolation/x-vals   xs
   ::interpolation/y-vals   ys
   ::interpolation/f-matrix f-matrix})

(deftest interpolation-2D-test
  (is (spec-check interpolation/interpolation-2D))
  (is= 113.66774636055577 (f-xy 5.0 25.0))                  ;actual point
  (is= 126.85943961592386 (f-xy 5.0 26.0))
  (let [polynomial (interpolation/interpolation-2D interpolation-2D-map)]
    (is= 113.66775071229222 (polynomial 5.0 25.0))
    (is= 126.85943961592386 (polynomial 5.0 26.0)))
  (let [bicubic-natural (interpolation/interpolation-2D
                          interpolation-2D-map
                          {::interpolation/interpolation-2D-type :bicubic})]
    (is= 113.3328166226598 (bicubic-natural 5.0 25.0))
    (is= 126.85943961592389 (bicubic-natural 5.0 26.0)))
  (let [bicubic-apache (interpolation/interpolation-2D
                         interpolation-2D-map
                         {::interpolation/interpolation-2D-type     :bicubic
                          ::interpolation/bicubic-spline-boundaries :apache})]
    (is= 113.34088323003708 (bicubic-apache 5.0 25.0))
    (is= 126.85943961592386 (bicubic-apache 5.0 26.0)))
  (let [bicubic-apache-smooth (interpolation/interpolation-2D
                                interpolation-2D-map
                                {::interpolation/interpolation-2D-type     :bicubic
                                 ::interpolation/bicubic-spline-boundaries :apache-with-polynomial-input-smoothing})]
    (is= 109.1102730561051 (bicubic-apache-smooth 5.0 25.0))
    (is= 126.85943961592386 (bicubic-apache-smooth 5.0 26.0)))
  (let [bicubic-closed (interpolation/interpolation-2D
                         interpolation-2D-map
                         {::interpolation/interpolation-2D-type     :bicubic
                          ::interpolation/bicubic-spline-boundaries :closed})]
    (is= 113.13417979331703 (bicubic-closed 5.0 25.0))
    (is= 126.85943961592383 (bicubic-closed 5.0 26.0)))
  (let [bilinear (interpolation/interpolation-2D
                   interpolation-2D-map
                   {::interpolation/interpolation-2D-type :bilinear})]
    (is= 114.66331578887625 (bilinear 5.0 25.0))
    (is= 126.85943961592386 (bilinear 5.0 26.0)))
  (let [bicubic-hermite (interpolation/interpolation-2D
                          interpolation-2D-map
                          {::interpolation/interpolation-2D-type :bicubic-hermite})]
    (is= 111.30800237859938 (bicubic-hermite 5.0 25.0))
    (is= 126.85943961592386 (bicubic-hermite 5.0 26.0))))

;;3D
(defn xyz-fn
  [x y z]
  (+ x
     z
     (m/log y)
     (* (m/exp (- x))
        (+ y (m/cube y)))))

(def f-xyz
  (tensor/compute-tensor
    [(count xs) (count ys) (count zs)]
    (fn [[i j k]]
      (let [x (get xs i m/nan)
            y (get ys j m/nan)
            z (get zs k m/nan)]
        (xyz-fn x y z)))))

(deftest interpolation-3D-test
  (is (spec-check interpolation/interpolation-3D))
  (is= 114.66774636055577 (xyz-fn 5.0 25.0 1.0))            ;actual points
  (is= 185.38415118074502 (xyz-fn 4.5 25.2 -0.4))
  (is= 68.58614648323473 (xyz-fn 5.0 21.0 -2.0))
  (let [f-3D (interpolation/interpolation-3D {::interpolation/x-vals   xs
                                              ::interpolation/y-vals   ys
                                              ::interpolation/z-vals   zs
                                              ::interpolation/f-tensor f-xyz})]
    (is= 111.12965547320519 (f-3D 5.0 25.0 1.0))
    (is= 132.712320895993 (f-3D 4.5 25.2 -0.4))
    (is= 68.58614648323407 (f-3D 5.0 21.0 -2.0))))

;;N-Dim
(deftest interpolation-ND-microsphere$-test
  (is (spec-check interpolation/interpolation-ND-microsphere$))
  (let [f-nd (interpolation/interpolation-ND-microsphere$
               {::interpolation/i-matrix [[0 0 0 0] [4 1 2 3] [9 9 9 9]]
                ::interpolation/f-vals   [0 5 9]})]
    ;(is= 1.1876105128225518 (f-nd [1 1 1 1]))
    ;(is= 2.0692505039789597 (f-nd [1 1 1 2]))
    ;(is= 1.849973654285914 (f-nd [1 1 2 1]))
    ;(is= 1.6021974126075473 (f-nd [1 2 1 1]))
    ;(is= 2.333487054025956 (f-nd [2 1 1 1]))
    (is= 0.0 (f-nd [0 0 0 0]))
    (is= 5.0 (f-nd [4 1 2 3]))
    (is= 9.0 (f-nd [9 9 9 9]))))

;;;SLOPE INTERPOLATION
(deftest g'-comparison-testing
  (is= 0.6434693402873666 (g' -0.5))
  (is= 55.067506039296525 (g' 2.5))
  (is= 111.23286869947819 (g' 4.5))
  (is= 50.558067735779616 (g' 5.5))
  (is= -54.42879349273511 (g' 6.0)))

(deftest slope-interpolation-1D-linear-test
  (is (spec-check interpolation/slope-interpolation-1D-linear))
  (let [linear-slope (interpolation/slope-interpolation-1D-linear
                       {::interpolation/x-vals xs
                        ::interpolation/f-vals f-vals})]
    (is= {::anomalies/category ::anomalies/no-solve
          ::anomalies/message  "Out of range: -0.5"
          ::anomalies/fn       (var interpolation/slope-interpolation-1D-linear)}
         (linear-slope -0.5))
    (is= -285.2043649357234 (linear-slope 7.0))
    (is= -285.2043649357234 (linear-slope 6.9))
    (is= -285.2043649357234 (linear-slope 6.0))
    (is= 40.98436560984152 (linear-slope 5.5))
    (is= 40.98436560984152 (linear-slope 5.0))
    (is= 73.39545303289322 (linear-slope 2.1))
    (is= 73.39545303289322 (linear-slope 2.0))))

(deftest slope-interpolation-1D-quadratic-test
  (is (spec-check interpolation/slope-interpolation-1D-quadratic))
  (let [zero-start-slope (interpolation/slope-interpolation-1D-quadratic
                           {::interpolation/x-vals xs
                            ::interpolation/f-vals f-vals
                            ::interpolation/slope  0.0
                            ::interpolation/end-slope? false})]
    (is= {::anomalies/category ::anomalies/no-solve
          ::anomalies/message  "Out of range: -0.5"
          ::anomalies/fn       (var interpolation/slope-interpolation-1D-quadratic)}
         (zero-start-slope -0.5))
    (is= -582.798385295781 (zero-start-slope 7.0))
    (is= -523.2795812237697 (zero-start-slope 6.9))
    (is= 12.38965542433425 (zero-start-slope 6.0))
    (is= 40.98436560984152 (zero-start-slope 5.5))
    (is= 69.57907579534879 (zero-start-slope 5.0))
    (is= 7.339545303289328 (zero-start-slope 2.1))
    (is= 0.0 (zero-start-slope 2.0)))
  (let [set-start-slope (interpolation/slope-interpolation-1D-quadratic
                          {::interpolation/x-vals xs
                           ::interpolation/f-vals f-vals
                           ::interpolation/slope  2.0
                           ::interpolation/end-slope? false})]
    (is= -580.798385295781 (set-start-slope 7.0))
    (is= -521.6795812237697 (set-start-slope 6.9))
    (is= 10.38965542433425 (set-start-slope 6.0))
    (is= 40.98436560984152 (set-start-slope 5.5))
    (is= 71.57907579534879 (set-start-slope 5.0))
    (is= 9.139545303289328 (set-start-slope 2.1))
    (is= 2.0 (set-start-slope 2.0)))
  (let [set-end-slope (interpolation/slope-interpolation-1D-quadratic
                        {::interpolation/x-vals xs
                         ::interpolation/f-vals f-vals
                         ::interpolation/slope  0.0
                         ::interpolation/end-slope? true})]
    (is= 0.0 (set-end-slope 7.0))
    (is= -57.040872987144475 (set-end-slope 6.9))
    (is= -570.4087298714468 (set-end-slope 6.0))
    (is= 40.984365609841575 (set-end-slope 5.5))
    (is= 652.3774610911298 (set-end-slope 5.0))
    (is= 531.8580920694922 (set-end-slope 2.1))
    (is= 582.798385295781 (set-end-slope 2.0))))

(deftest slope-interpolation-1D-cubic-clamped-test
  (is (spec-check interpolation/slope-interpolation-1D-cubic-clamped))
  (let [zero-acceleration (interpolation/slope-interpolation-1D-cubic-clamped
                            {::interpolation/x-vals             xs
                             ::interpolation/f-vals             f-vals
                             ::interpolation/start-acceleration 0.0
                             ::interpolation/end-acceleration   0.0})]
    (is= {::anomalies/category ::anomalies/no-solve
          ::anomalies/message  "Out of range: -0.5"
          ::anomalies/fn       (var interpolation/slope-interpolation-1D-cubic-clamped)}
         (zero-acceleration -0.5))
    (is= -367.34794066509505 (zero-acceleration 7.0))
    (is= -364.8836333932139 (zero-acceleration 6.9))
    (is= -120.91721347698004 (zero-acceleration 6.0))
    (is= 62.116652635164854 (zero-acceleration 5.5))
    (is= 118.35679659536976 (zero-acceleration 5.0))
    (is= 62.67514110479001 (zero-acceleration 2.1))
    (is= 62.594131190975624 (zero-acceleration 2.0)))
  (let [set-closer-accelerations (interpolation/slope-interpolation-1D-cubic-clamped
                                   {::interpolation/x-vals             xs
                                    ::interpolation/f-vals             f-vals
                                    ::interpolation/start-acceleration 2.0
                                    ::interpolation/end-acceleration   10.0})]
    (is= -364.4680957038547 (set-closer-accelerations 7.0))
    (is= -362.9401837808108 (set-closer-accelerations 6.9))
    (is= -121.67690339946068 (set-closer-accelerations 6.0))
    (is= 62.26684643361447 (set-closer-accelerations 5.5))
    (is= 118.51571132405196 (set-closer-accelerations 5.0))
    (is= 61.80609071719313 (set-closer-accelerations 2.1))
    (is= 61.532115687099676 (set-closer-accelerations 2.0)))
  (let [set-further-accelerations (interpolation/slope-interpolation-1D-cubic-clamped
                                    {::interpolation/x-vals             xs
                                     ::interpolation/f-vals             f-vals
                                     ::interpolation/start-acceleration 1.0
                                     ::interpolation/end-acceleration   20.0})]
    (is= -361.5766228356377 (set-further-accelerations 7.0))
    (is= -360.9854550986403 (set-further-accelerations 6.9))
    (is= -122.45984913589476 (set-further-accelerations 6.0))
    (is= 62.402505348343155 (set-further-accelerations 5.5))
    (is= 118.75602140157135 (set-further-accelerations 5.0))
    (is= 62.29831939936368 (set-further-accelerations 2.1))
    (is= 62.12126297392139 (set-further-accelerations 2.0))))

#_(ost/unstrument)