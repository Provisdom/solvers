(ns provisdom.solvers.interpolation
  (:require
    [clojure.spec.alpha :as s]
    [clojure.spec.gen.alpha :as gen]
    [clojure.spec.test.alpha :as st]
    [orchestra.spec.test :as ost]
    [provisdom.utility-belt.anomalies :as anomalies]
    [provisdom.math.core :as m]
    [provisdom.math.matrix :as mx]
    [provisdom.math.apache-matrix :as apache-mx]
    [provisdom.math.vector :as vector]
    [provisdom.math.tensor :as tensor]
    [provisdom.solvers.internal-apache-solvers :as apache-solvers]
    [incanter.interpolation :as incanter]
    [clojure.core.matrix :as ccm]))

(s/def ::x-vals ::apache-solvers/x-vals)
(s/def ::y-vals ::apache-solvers/y-vals)
(s/def ::z-vals ::apache-solvers/z-vals)
(s/def ::f-vals ::apache-solvers/f-vals)

(s/def ::x-vals-with-f-vals
  (s/with-gen
    (s/and (s/keys :req [::x-vals ::f-vals])
           (fn [{::keys [x-vals f-vals]}]
             (= (count f-vals) (count x-vals))))
    #(gen/bind
       (s/gen ::x-vals)
       (fn [x]
         (gen/bind
           (gen/vector (s/gen ::m/num)
                       (count x))
           (fn [f-v]
             (gen/return {::x-vals x
                          ::f-vals f-v})))))))

(s/def ::slope ::m/finite)
(s/def ::end-slope? boolean?)
(s/def ::start-acceleration ::m/finite)
(s/def ::end-acceleration ::m/finite)

(s/def ::x-vals-with-f-vals-and-slope
  (s/merge ::x-vals-with-f-vals
           (s/keys :req [::slope ::end-slope?])))

(s/def ::x-vals-with-f-vals-and-acceleration
  (s/merge ::x-vals-with-f-vals
           (s/keys :req [::start-acceleration ::end-acceleration])))

;(s/def ::interpolation-1D-type
;#{:polynomial :neville :loess :cubic :linear :akima})

(defn interpolation-1D-quadratic
  "`::x-vals` must be strictly increasing. `::slope` is the slope at the first
  'x-val' by default. Setting `::end-slope?` to true sets the slope at the last
  'x-val' instead."
  [{::keys [x-vals f-vals slope end-slope?]}]
  (if end-slope?
    (fn [x]
      ((interpolation-1D-quadratic {::x-vals     (mapv - (reverse x-vals))
                                    ::f-vals     (vec (reverse f-vals))
                                    ::slope      (- slope)
                                    ::end-slope? false})
        (- x)))
    (fn [x]
      (cond (or (< x (first x-vals))
                (> x (peek x-vals))) {::anomalies/category ::anomalies/no-solve
                                      ::anomalies/message  (str "Out of range: " x)
                                      ::anomalies/fn       (var interpolation-1D-quadratic)}
            (== x (peek x-vals)) (peek f-vals)
            (== x (first x-vals)) (first f-vals)
            :else (let [index (ffirst
                                (drop-while (fn [[_ xn]]
                                              (<= xn x))
                                            (map vector (range) x-vals)))
                        index (or index (dec (count x-vals)))
                        x1 (get x-vals index m/nan)
                        y1 (get f-vals index m/nan)]
                    (if (== x x1)
                      y1
                      (let [x0 (get x-vals (dec index) m/nan)
                            y0 (get f-vals (dec index) m/nan)
                            z (vec (reductions
                                     (fn [tot i]
                                       (- (* 2.0 (/ (- (get f-vals (inc i) m/nan)
                                                       (get f-vals i m/nan))
                                                    (- (get x-vals (inc i) m/nan)
                                                       (get x-vals i m/nan))))
                                          tot))
                                     (double slope)
                                     (range index)))
                            [z0 z1] (subvec z (dec index) (inc index))]
                        (+ y0 (* z0 (- x x0))
                           (* 0.5
                              (- z1 z0)
                              (m/sq (- x x0))
                              (/ (- x1 x0)))))))))))

(s/fdef interpolation-1D-quadratic
        :args (s/cat :x-vals-with-f-vals-and-slope ::x-vals-with-f-vals-and-slope)
        :ret (s/fspec :args (s/cat :x ::m/num)
                      :ret (s/or :inner-solution ::m/number
                                 :inner-anomaly ::anomalies/anomaly)))

;;e.g., this function is a generalized form of what can be found here: 
;;http://online.redwoods.edu/instruct/darnold/laproj/Fall98/SkyMeg/Proj.PDF
(defn interpolation-1D-cubic-clamped
  "`::start-acceleration` and `::end-acceleration` should be the second
  derivatives at the first and last data points. `::x-vals` must be strictly
  increasing."
  [{::keys [x-vals f-vals start-acceleration end-acceleration]}]
  (let [data-count (count x-vals)
        rhs (vector/compute-vector
              data-count
              (fn [index]
                (cond (zero? index) start-acceleration
                      (>= index (dec data-count)) end-acceleration
                      :else (let [x0 (get x-vals (dec index) m/nan)
                                  x1 (get x-vals index m/nan)
                                  x2 (get x-vals (inc index) m/nan)
                                  y0 (get f-vals (dec index) m/nan)
                                  y1 (get f-vals index m/nan)
                                  y2 (get f-vals (inc index) m/nan)]
                              (* 6.0 (- (/ (- y2 (double y1)) (- x2 (double x1)))
                                        (/ (- y1 (double y0)) (- x1 (double x0)))))))))
        coefficient-mx (mx/compute-matrix
                         data-count
                         data-count
                         (fn [r c]
                           (if (or (zero? r) (= (dec data-count) r))
                             (if (= r c) 1.0 0.0)
                             (cond (= r c) (* 2.0 (- (get x-vals (inc r) m/nan)
                                                     (get x-vals (dec r) m/nan)))
                                   (= r (inc c)) (- (get x-vals r m/nan)
                                                    (get x-vals (dec r) m/nan))
                                   (= r (dec c)) (- (get x-vals (inc r) m/nan)
                                                    (get x-vals r m/nan))
                                   :else 0.0))))
        z-apache-mx (::apache-mx/LLS-solution
                      (apache-mx/qr-decomposition-with-linear-least-squares
                        (apache-mx/apache-matrix coefficient-mx)
                        (apache-mx/apache-matrix (mx/column-matrix rhs))))]
    (fn [x]
      (if (anomalies/anomaly? z-apache-mx)
        z-apache-mx
        (cond (or (< x (first x-vals))
                  (> x (peek x-vals))) {::anomalies/category ::anomalies/no-solve
                                        ::anomalies/message  (str "Out of range: " x)
                                        ::anomalies/fn       (var interpolation-1D-cubic-clamped)}
              (== x (peek x-vals)) (peek f-vals)
              (== x (first x-vals)) (first f-vals)
              :else (let [index (ffirst
                                  (drop-while (fn [[_ xn]]
                                                (<= xn x))
                                              (map vector (range) x-vals)))
                          x1 (get x-vals (dec index) m/nan)
                          y1 (get f-vals (dec index) m/nan)
                          x2 (get x-vals index m/nan)
                          y2 (get f-vals index m/nan)
                          z1 (apache-mx/get-entry z-apache-mx (dec index) 0)
                          z2 (apache-mx/get-entry z-apache-mx index 0)
                          h (- x2 x1)
                          s (/ 6.0)
                          ih (/ h)]
                      (+ (* (+ (* z2 (m/cube (- x x1)))
                               (* z1 (m/cube (- x2 x))))
                            s
                            ih)
                         (* (- x x1) (- (* y2 ih) (* h z2 s)))
                         (* (- x2 x) (- (* y1 ih) (* h z1 s))))))))))

(s/fdef interpolation-1D-cubic-clamped
        :args (s/cat :x-vals-with-f-vals-and-acceleration ::x-vals-with-f-vals-and-acceleration)
        :ret (s/fspec :args (s/cat :x ::m/num)
                      :ret (s/or :inner-solution ::m/number
                                 :inner-anomaly ::anomalies/anomaly)))

(defn interpolation-1D
  "Implements an algorithm for interpolation of real univariate functions.
  `::x-vals` - All the x-coordinates of the interpolation points, in strictly
    ascending order.
  `::f-vals` - The values of the interpolation points on all the grid knots:
    `f-vals`[i] = f(`x-vals`[i]).

  `::interpolation-1D-type` can be:
    `:cubic` - (default) Computes a natural (also known as 'free', 'unclamped')
      cubic spline interpolation for the data set. Can not extrapolate.
    `:cubic-closed` - S'(start) = S'(end), S''(start) = S''(end). This is
      similar to having a full-length period, except there can be a trend.
    `:cubic-hermite` - Similar to `:cubic-closed` but tries to use derivatives.
    `:akima` - Also a cubic spline interpolation. Can not extrapolate. The Akima
      algorithm requires that the number of points >= 5.
    `:loess` - The Local Regression Algorithm (also Loess, Lowess). Can not
      extrapolate.
    `:linear` - Linear function. Can not extrapolate.
    `:polynomial` - Polynomial Divided Difference Algorithm. For reference, see
      Introduction to Numerical Analysis, ISBN 038795452X, chapter 2.
    `:neville` - Neville's Algorithm. Similar results to `:polynomial`. For
      reference, see Introduction to Numerical Analysis, ISBN 038795452X,
      chapter 2.

  Only `:cubic-hermite`, `:cubic-closed`, `:polynomial`, and `:neville` can
  extrapolate.

  If a positive `::period` length is entered, then the data is assumed to be
  periodic. Doesn't apply to `:cubic-hermite` or `:cubic-closed`.

  For the `:loess` interpolation only:
    `::bandwidth-for-loess` when computing the loess fit at a particular point,
      this fraction of source points closest to the current point is taken into
      account for computing a least-squares regression. A sensible value is
      usually 0.25 to 0.5, the default value is 0.3. For reference, see
      William S. Cleveland - Robust Locally Weighted Regression and Smoothing
      Scatterplots.

  Returns a function that accepts a point 'x' and returns the interpolated
  value."
  ([args] (interpolation-1D args {}))
  ([{::keys [x-vals f-vals]}
    {::keys [interpolation-1D-type period bandwidth-for-loess]
     :or    {interpolation-1D-type :cubic
             bandwidth-for-loess   0.3
             period                0.0}}]
   (if (or (= :cubic-hermite interpolation-1D-type)
           (= :cubic-closed interpolation-1D-type))
     (try
       (let [ii-fn (incanter/interpolate
                     (partition 2 (interleave x-vals f-vals))
                     (if (= interpolation-1D-type :cubic-closed)
                       :cubic
                       :cubic-hermite)
                     :boundaries :closed)]
         (fn [x]
           (or (ii-fn x) {::anomalies/message  "Inner Interpolation Function returned nil."
                          ::anomalies/fn       (var interpolation-1D)
                          ::anomalies/category ::anomalies/third-party})))
       (catch Exception e {::anomalies/message  (or (.getMessage e) "")
                           ::anomalies/fn       (var interpolation-1D)
                           ::anomalies/category ::anomalies/third-party}))
     (apache-solvers/interpolation-1D
       {::apache-solvers/x-vals x-vals
        ::apache-solvers/f-vals f-vals}
       {::apache-solvers/interpolation-1D-type interpolation-1D-type
        ::apache-solvers/period                period
        ::apache-solvers/bandwidth-for-loess   bandwidth-for-loess}))))

(s/def ::interpolation-1D-type
  (s/or :apache ::apache-solvers/interpolation-1D-type
        :incanter #{:cubic-closed :cubic-hermite}))

(s/def ::period ::apache-solvers/period)
(s/def ::bandwidth-for-loess ::apache-solvers/bandwidth-for-loess)

(s/fdef interpolation-1D
        :args (s/cat :x-vals-with-f-vals ::x-vals-with-f-vals
                     :opts (s/? (s/keys :opt [::interpolation-1D-type
                                              ::period
                                              ::bandwidth-for-loess])))
        :ret (s/or :solution (s/fspec :args (s/cat :x ::m/num)
                                      :ret (s/or :inner-solution ::m/number
                                                 :inner-anomaly ::anomalies/anomaly))
                   :anomaly ::anomalies/anomaly))

(defn interpolation-1D-using-derivatives
  "`data` should be a collection of seqs, where each seq contains
  [x value & first derivative, the second derivative, and so ...]
  Uses Hermite Interpolator. Returns a function that accepts a point x and
  returns the interpolated value."
  [data]
  (apache-solvers/interpolation-1D-using-derivatives data))

(s/def ::point-data ::apache-solvers/point-data)

(s/fdef interpolation-1D-using-derivatives
        :args (s/cat :data ::point-data)
        :ret (s/fspec :args (s/cat :x ::m/num)
                      :ret (s/or :inner-solution ::m/number
                                 :inner-anomaly ::anomalies/anomaly)))

(defn interpolation-2D
  "Generates an interpolation function over a rectangular grid.
  `::x-vals` - All the x-coordinates of the interpolation points, in strictly
    ascending order.
  `::y-vals` - All the y-coordinates of the interpolation points, in strictly
    ascending order.
  `::f-matrix` - The values of the interpolation points on all the grid knots:
   `f-matrix`[i][j] = f(`x-vals`[i], `y-vals`[j]).

  `::interpolation-2D-type`:
    `:polynomial` (default), `:bicubic`, `:bicubic-hermite`, `:bilinear`.

  `::bicubic-spline-boundaries` define the boundary conditions for the :bicubic
    spline. Possible values are: `:natural` (default), `:closed`, `:apache`,
    and `:apache-with-polynomial-input-smoothing`. The `:apache` requires at
    least 5 data points. Apache didn't specify the algorithm used.

  Returns a function that accepts the points 'x' and 'y' and returns the
  interpolated value."
  ([args] (interpolation-2D args {}))
  ([{::keys [x-vals y-vals f-matrix]}
    {::keys [interpolation-2D-type bicubic-spline-boundaries]
     :or    {interpolation-2D-type     :polynomial
             bicubic-spline-boundaries :natural}}]
   (if (and (= interpolation-2D-type :bicubic)
            (or (= bicubic-spline-boundaries :apache)
                (= bicubic-spline-boundaries :apache-with-polynomial-input-smoothing)))
     (apache-solvers/interpolation-2D
       {::apache-solvers/x-vals   x-vals
        ::apache-solvers/y-vals   y-vals
        ::apache-solvers/f-matrix f-matrix}
       {::apache-solvers/smooth? (= bicubic-spline-boundaries
                                    :apache-with-polynomial-input-smoothing)})
     (try (fn [x y]
            (try
              (ccm/mget ((incanter/interpolate-grid
                           (mx/transpose f-matrix)
                           interpolation-2D-type
                           :xs (map double x-vals)
                           :ys (map double y-vals)
                           :boundaries bicubic-spline-boundaries)
                          x
                          y))
              (catch Exception e
                {::anomalies/message  (str "Returned Function: " (.getMessage e))
                 ::anomalies/fn       (var interpolation-2D)
                 ::anomalies/category ::anomalies/third-party})))
          (catch Exception e
            {::anomalies/message  (.getMessage e)
             ::anomalies/fn       (var interpolation-2D)
             ::anomalies/category ::anomalies/third-party})))))

(s/def ::f-matrix ::mx/matrix-num)

(s/def ::interpolation-2D-type
  #{:polynomial :bicubic :bicubic-hermite :bilinear})

(s/def ::bicubic-spline-boundaries
  #{:natural :closed :apache :apache-with-polynomial-input-smoothing})

(s/def ::vals-with-matrix
  (s/with-gen
    (s/and (s/keys :req [::x-vals ::y-vals ::f-matrix])
           (fn [{::keys [x-vals y-vals f-matrix]}]
             (and (= (mx/rows f-matrix) (count x-vals))
                  (= (mx/columns f-matrix) (count y-vals)))))
    #(gen/bind
       (gen/tuple (s/gen ::x-vals)
                  (s/gen ::y-vals))
       (fn [[x y]]
         (gen/bind (gen/vector
                     (gen/vector (s/gen ::m/num)
                                 (count y))
                     (count x))
                   (fn [f-m]
                     (gen/return {::x-vals   x
                                  ::y-vals   y
                                  ::f-matrix f-m})))))))

(s/fdef interpolation-2D
        :args (s/and (s/cat :vals-with-matrix ::vals-with-matrix
                            :opts (s/? (s/keys :opt [::interpolation-2D-type
                                                     ::bicubic-spline-boundaries])))
                     (fn [{:keys [vals-with-matrix opts]}]
                       (let [{::keys [f-matrix]} vals-with-matrix
                             {::keys [interpolation-2D-type bicubic-spline-boundaries]} opts]
                         (or (not (and (= interpolation-2D-type :bicubic)
                                       (= bicubic-spline-boundaries :apache)))
                             (>= (tensor/ecount f-matrix) 5)))))
        :ret (s/or :solution (s/fspec :args (s/cat :x ::m/num
                                                   :y ::m/num)
                                      :ret (s/or :inner-solution ::m/number
                                                 :inner-anomaly ::anomalies/anomaly))
                   :anomaly ::anomalies/anomaly))

(defn interpolation-3D
  "`::x-vals` - All the x-coordinates of the interpolation points, in strictly
    ascending order.
  `::y-vals` - All the y-coordinates of the interpolation points, in strictly
    ascending order.
  `::z-vals` - All the z-coordinates of the interpolation points, in strictly
    ascending order.
  `::f-tensor` - the values of the interpolation points on all the grid knots:
    `f-tensor`[i][j][k] = f(`x-vals`[i], `y-vals`[j], `z-vals`[k]).
  Returns a function that accepts the points 'x', 'y', and 'z' and returns the
  interpolated value."
  [{::keys [x-vals y-vals z-vals f-tensor]}]
  (apache-solvers/interpolation-3D
    {::apache-solvers/x-vals   x-vals
     ::apache-solvers/y-vals   y-vals
     ::apache-solvers/z-vals   z-vals
     ::apache-solvers/f-tensor f-tensor}))

(s/def ::f-tensor ::apache-solvers/f-tensor)

(s/def ::vals-with-tensor
  (s/with-gen
    (s/and (s/keys :req [::x-vals ::y-vals ::z-vals ::f-tensor])
           (fn [{::keys [x-vals y-vals z-vals f-tensor]}]
             (let [[xc yc zc] (tensor/shape f-tensor)]
               (and (= xc (count x-vals))
                    (= yc (count y-vals))
                    (= zc (count z-vals))))))
    #(gen/bind
       (gen/tuple (s/gen ::x-vals)
                  (s/gen ::y-vals)
                  (s/gen ::z-vals))
       (fn [[x y z]]
         (gen/bind (gen/vector
                     (gen/vector
                       (gen/vector (s/gen ::m/num)
                                   (count z))
                       (count y))
                     (count x))
                   (fn [f-t]
                     (gen/return {::x-vals   x
                                  ::y-vals   y
                                  ::z-vals   z
                                  ::f-tensor f-t})))))))

(s/fdef interpolation-3D
        :args (s/cat :vals-with-tensor ::vals-with-tensor)
        :ret (s/or :solution (s/fspec :args (s/cat :x ::m/num
                                                   :y ::m/num
                                                   :z ::m/num)
                                      :ret (s/or :inner-solution ::m/number
                                                 :inner-anomaly ::anomalies/anomaly))
                   :anomaly ::anomalies/anomaly))

(defn interpolation-ND-microsphere$
  "Interpolator that implements the algorithm described in William Dudziak's MS
  thesis. Results are randomized.
  `::i-matrix` - the arguments for the interpolation points.
  `i-matrix`[i][0] is the first component of interpolation point 'i',
  `i-matrix`[i][1] is the second component, and so on until `i-matrix`[i][d-1],
  the last component of that interpolation point (where 'd' is thus the
  dimension of the space).
  `::f-vals` - the values for the interpolation points.
  Returns a function that takes a vector and returns a number."
  [{::keys [i-matrix f-vals]}]
  (apache-solvers/interpolation-ND-microsphere$
    {::apache-solvers/i-matrix i-matrix
     ::apache-solvers/f-vals   f-vals}))

(s/def ::i-matrix ::apache-solvers/i-matrix)

(s/def ::i-matrix-with-f-vals
  (s/with-gen
    (s/and (s/keys :req [::i-matrix ::f-vals])
           (fn [{::keys [i-matrix f-vals]}]
             (= (count f-vals) (count i-matrix))))
    #(gen/bind
       (gen/tuple (s/gen ::f-vals) (s/gen (s/int-in 1 5)))
       (fn [[f-v c]]
         (gen/bind
           (gen/vector
             (gen/vector (s/gen ::m/num)
                         c)
             (count f-v))
           (fn [i-mx]
             (gen/return {::i-matrix i-mx
                          ::f-vals   f-v})))))))

(s/fdef interpolation-ND-microsphere$
        :args (s/cat :i-matrix-with-f-vals ::i-matrix-with-f-vals)
        :ret (s/fspec :args (s/cat :v ::vector/vector-num)
                      :ret (s/or :inner-solution ::m/number
                                 :inner-anomaly ::anomalies/anomaly)))

;;;SLOPE INTERPOLATION
(defn slope-interpolation-1D-linear
  "`::x-vals` must be strictly increasing. Slopes exactly at `x-vals` take on
  the slope on the greater side of the 'x-val'."
  [{::keys [x-vals f-vals]}]
  (let [data-count (count x-vals)]
    (fn [x]
      (cond (or (< x (first x-vals))
                (> x (peek x-vals))) {::anomalies/category ::anomalies/no-solve
                                      ::anomalies/message  (str "Out of range: " x)
                                      ::anomalies/fn       (var slope-interpolation-1D-linear)}
            (== x (peek x-vals)) (/ (- (double (peek f-vals)) (nth f-vals (- data-count 2)))
                                    (- x (double (get x-vals (- data-count 2) m/nan))))
            (== x (first x-vals)) (/ (- (double (second f-vals)) (first f-vals))
                                     (- (double (second x-vals)) (first x-vals)))
            :else (let [index (ffirst
                                (drop-while (fn [[_ xn]]
                                              (<= xn x))
                                            (map vector (range) x-vals)))
                        x0 (double (get x-vals (dec index) m/nan))
                        x1 (double (get x-vals index m/nan))
                        y0 (double (get f-vals (dec index) m/nan))
                        y1 (double (get f-vals index m/nan))
                        slope (/ (- y1 y0) (- x1 x0))]
                    (if (== x x1)
                      (let [x2 (get x-vals (inc index) m/nan)
                            y2 (get f-vals (inc index) m/nan)]
                        (* 0.5 (+ slope (/ (- y2 y1) (- x2 x1)))))
                      slope))))))

(s/fdef slope-interpolation-1D-linear
        :args (s/cat :x-vals-with-f-vals ::x-vals-with-f-vals)
        :ret (s/fspec :args (s/cat :x ::m/num)
                      :ret (s/or :inner-solution ::m/number
                                 :inner-anomaly ::anomalies/anomaly)))

(defn slope-interpolation-1D-quadratic
  "`::x-vals` must be strictly increasing. `::slope` is the slope at the first
  'x-val' by default. Setting `::end-slope?` to true sets the slope at the last
  'x-val' instead."
  [{::keys [x-vals f-vals slope end-slope?]}]
  (if end-slope?
    (fn [x]
      (let [sol-fn (slope-interpolation-1D-quadratic
                     {::x-vals     (mapv - (reverse x-vals))
                      ::f-vals     (vec (reverse f-vals))
                      ::slope      (- slope)
                      ::end-slope? false})
            sol (sol-fn (- (double x)))]
        (if (anomalies/anomaly? sol)
          sol
          (- sol))))
    (fn [x]
      (cond (or (< x (first x-vals))
                (> x (peek x-vals))) {::anomalies/category ::anomalies/no-solve
                                      ::anomalies/message  (str "Out of range: " x)
                                      ::anomalies/fn       (var slope-interpolation-1D-quadratic)}
            (== x (first x-vals)) slope
            :else (let [index (ffirst
                                (drop-while (fn [[_ xn]]
                                              (<= xn x))
                                            (map vector (range) x-vals)))
                        index (or index (dec (count x-vals)))
                        x0 (double (get x-vals (dec index) m/nan))
                        x1 (double (get x-vals index m/nan))
                        midpoint? (and (== x x1) (not (== x (peek x-vals))))
                        slope-fn #(+ %3 (* (- (double %4) %3)
                                           (- (double x) %1)
                                           (/ (- (double %2) %1))))
                        z (vec (reductions
                                 (fn [tot e]
                                   (- (* 2.0
                                         (/ (- (double (get f-vals (inc e) m/nan))
                                               (get f-vals e m/nan))
                                            (- (double (get x-vals (inc e) m/nan))
                                               (get x-vals e m/nan))))
                                      tot))
                                 (double slope)
                                 (if midpoint?
                                   (range (inc index))
                                   (range index))))]
                    (if midpoint?
                      (let [[z0 z1 z2] (subvec z (dec index) (+ index 2))
                            x2 (double (get x-vals (inc index) m/nan))]
                        (* 0.5 (+ (slope-fn x0 x1 z0 z1)
                                  (slope-fn x1 x2 z1 z2))))
                      (let [[z0 z1] (subvec z (dec index) (inc index))]
                        (slope-fn x0 x1 z0 z1))))))))

(s/fdef slope-interpolation-1D-quadratic
        :args (s/cat :x-vals-with-f-vals-and-slope ::x-vals-with-f-vals-and-slope)
        :ret (s/fspec :args (s/cat :x ::m/num)
                      :ret (s/or :inner-solution ::m/number
                                 :inner-anomaly ::anomalies/anomaly)))

(defn slope-interpolation-1D-cubic-clamped
  "`::start-acceleration` and `::end-acceleration` should be the second
  derivatives at the first and last data points. `::x-vals` must be strictly
  increasing. Returns a function that returns interpolated slope instead of
  interpolated value."
  [{::keys [x-vals f-vals start-acceleration end-acceleration]}]
  (let [data-count (count x-vals)
        rhs (vector/compute-vector
              data-count
              (fn [index]
                (cond (zero? index) start-acceleration
                      (>= index (dec data-count)) end-acceleration
                      :else (let [x0 (get x-vals (dec index) m/nan)
                                  x1 (get x-vals index m/nan)
                                  x2 (get x-vals (inc index) m/nan)
                                  y0 (get f-vals (dec index) m/nan)
                                  y1 (get f-vals index m/nan)
                                  y2 (get f-vals (inc index) m/nan)]
                              (* 6.0
                                 (- (/ (- (double y2) y1) (- (double x2) x1))
                                    (/ (- (double y1) y0) (- (double x1) x0))))))))
        coefficient-mx (mx/compute-matrix
                         data-count
                         data-count
                         (fn [r c]
                           (if (or (zero? r) (= (dec data-count) r))
                             (if (= r c) 1.0 0.0)
                             (cond (= r c) (* 2.0
                                              (- (double (get x-vals (inc r) m/nan))
                                                 (get x-vals (dec r) m/nan)))
                                   (= r (inc c)) (- (double (get x-vals r m/nan))
                                                    (get x-vals (dec r) m/nan))
                                   (= r (dec c)) (- (double (get x-vals (inc r) m/nan))
                                                    (get x-vals r m/nan))
                                   :else 0.0))))
        z-apache-mx (::apache-mx/LLS-solution
                      (apache-mx/qr-decomposition-with-linear-least-squares
                        (apache-mx/apache-matrix coefficient-mx)
                        (apache-mx/apache-matrix (mx/column-matrix rhs))))]
    (fn [x]
      (if (anomalies/anomaly? z-apache-mx)
        z-apache-mx
        (cond (or (< x (first x-vals))
                  (> x (peek x-vals))) {::anomalies/category ::anomalies/no-solve
                                        ::anomalies/message  (str "Out of range: " x)
                                        ::anomalies/fn       (var slope-interpolation-1D-cubic-clamped)}
              :else (let [index (ffirst
                                  (drop-while (fn [[_ xn]]
                                                (<= xn x))
                                              (map vector (range) x-vals)))
                          index (or index (dec (count x-vals)))
                          x1 (get x-vals (dec index) m/nan)
                          y1 (get f-vals (dec index) m/nan)
                          x2 (double (get x-vals index m/nan))
                          y2 (double (get f-vals index m/nan))
                          z1 (apache-mx/get-entry z-apache-mx (dec index) 0)
                          z2 (double (apache-mx/get-entry z-apache-mx index 0))
                          h (- x2 x1)
                          s (/ 6.0)
                          ih (/ h)]
                      (+ (/ (- y2 y1) h)
                         (- (* h
                               (- z2 z1)
                               s))
                         (* 0.5
                            ih
                            (- (* z2 (m/sq (- x x1)))
                               (* z1 (m/sq (- x x2))))))))))))

(s/fdef slope-interpolation-1D-cubic-clamped
        :args (s/cat :x-vals-with-f-vals-and-acceleration ::x-vals-with-f-vals-and-acceleration)
        :ret (s/fspec :args (s/cat :x ::m/num)
                      :ret (s/or :inner-solution ::m/number
                                 :inner-anomaly ::anomalies/anomaly)))
