(ns provisdom.solvers.curve-fitting
  (:require
    [clojure.spec.alpha :as s]
    [clojure.spec.gen.alpha :as gen]
    [clojure.spec.test.alpha :as st]
    [orchestra.spec.test :as ost]
    [provisdom.utility-belt.anomalies :as anomalies]
    [provisdom.math.core :as m]
    [provisdom.math.vector :as vector]
    [provisdom.math.matrix :as mx]
    [incanter.interpolation :as incanter]
    [uncomplicate.neanderthal.native :as native]
    [uncomplicate.neanderthal.core :as neanderthal]
    [uncomplicate.neanderthal.linalg :as linear-algebra]
    [provisdom.math.neanderthal-matrix :as neanderthal-mx]))

;;;LINE FITTING
(s/def ::f-vals
  (s/coll-of ::m/finite
             :kind clojure.core/vector?
             :into []
             :min-count 2))

(defn linear-least-squares-line-fitting
  "Returns a map of a `::line-fitting-weights` vector and a function
  `::line-fitting-fn` that takes and returns a finite. `basis-fn` takes a
  finite and returns a vector. To use B-splines as the basis, use
  [[b-spline-line-fitting]]."
  [{::keys [x-vals f-vals basis-fn]}]
  (let [bf-mx (mapv basis-fn x-vals)]
    (if-not (mx/matrix? bf-mx)
      {::anomalies/message  "basis-fn returns are inconsistent"
       ::anomalies/category ::anomalies/exception
       ::anomalies/fn       (var linear-least-squares-line-fitting)}
      (let [lhs (neanderthal-mx/matrix->neanderthal-matrix bf-mx)
            rhs (neanderthal-mx/matrix->neanderthal-matrix
                  (mx/column-matrix f-vals))
            solution (neanderthal-mx/lls lhs rhs)
            lls (cond (anomalies/anomaly? solution) solution
                      (neanderthal-mx/empty-neanderthal-matrix? solution) []
                      :else (mx/get-column
                              (neanderthal-mx/neanderthal-matrix->matrix
                                solution)
                              0))]
        (if (anomalies/anomaly? lls)
          lls
          {::line-fitting-fn      (fn [x]
                                    (let [bf-v (basis-fn x)]
                                      (if (= (count lls) (count bf-v))
                                        (vector/dot-product lls bf-v)
                                        m/nan)))
           ::line-fitting-weights lls})))))

(s/def ::line-fitting-fn
  (s/fspec :args (s/cat :finite ::m/finite)
           :ret ::m/number))

(s/def ::line-fitting-weights ::vector/vector)

(s/def ::basis-fn
  (s/fspec :args (s/cat :finite ::m/finite)
           :ret ::vector/vector-finite))

(s/def ::x-vals
  (s/coll-of ::m/finite
             :kind clojure.core/vector?
             :into []
             :min-count 2))

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

(s/def ::vals-with-basis
  (s/merge ::x-vals-with-f-vals (s/keys :req [::basis-fn])))

(s/fdef linear-least-squares-line-fitting
        :args (s/cat :vals-with-basis ::vals-with-basis)
        :ret (s/or :sol (s/keys :req [::line-fitting-fn ::line-fitting-weights])
                   :anomaly ::anomalies/anomaly))

(defn b-spline-line-fitting
  "Normally returns a function that takes a finite and returns a finite or nil.
  Uses linear least squares."
  [{::keys [distinct-x-vals f-vals degree]}]
  (try (incanter/interpolate (partition 2 (interleave distinct-x-vals f-vals))
                             :linear-least-squares
                             :basis :b-spline
                             :degree degree)
       (catch Exception e
         {::anomalies/fn       (var b-spline-line-fitting)
          ::anomalies/message  (.getMessage e)
          ::anomalies/category ::anomalies/third-party})))

(s/def ::distinct-x-vals
  (s/coll-of ::m/finite
             :kind clojure.core/vector?
             :distinct true
             :into []
             :min-count 2))

(s/def ::distinct-x-vals-with-f-vals
  (s/with-gen
    (s/and (s/keys :req [::distinct-x-vals ::f-vals])
           (fn [{::keys [distinct-x-vals f-vals]}]
             (= (count f-vals) (count distinct-x-vals))))
    #(gen/bind
       (s/gen ::distinct-x-vals)
       (fn [x]
         (gen/bind
           (gen/vector (s/gen ::m/num)
                       (count x))
           (fn [f-v]
             (gen/return {::distinct-x-vals x
                          ::f-vals          f-v})))))))

(s/def ::degree ::m/non-)

(s/def ::vals-with-degree
  (s/merge ::distinct-x-vals-with-f-vals (s/keys :req [::degree])))

(s/fdef b-spline-line-fitting
        :args (s/cat :vals-with-degree ::vals-with-degree)
        :ret (s/or :sol (s/fspec :args (s/cat :finite ::m/finite)
                                 :ret (s/nilable ::m/number))
                   :anomaly ::anomalies/anomaly))

;;;CURVE FITTING
(defn linear-least-squares-curve-fitting
  "Returns a map of a `::curve-fitting-weights` vector and a function
  `::curve-fitting-fn` that takes a vector and returns a finite.
  `curve-basis-fn` takes and returns a vector."
  [{::keys [x-matrix f-vals curve-basis-fn]}]
  (let [bf-mx (mapv curve-basis-fn x-matrix)]
    (if-not (mx/matrix? bf-mx)
      {::anomalies/message  "basis-fn returns are inconsistent"
       ::anomalies/category ::anomalies/exception
       ::anomalies/fn       (var linear-least-squares-curve-fitting)}
      (let [lhs (neanderthal-mx/matrix->neanderthal-matrix bf-mx)
            rhs (neanderthal-mx/matrix->neanderthal-matrix
                  (mx/column-matrix f-vals))
            solution (neanderthal-mx/lls lhs rhs)
            lls (cond (anomalies/anomaly? solution) solution
                      (neanderthal-mx/empty-neanderthal-matrix? solution) []
                      :else (mx/get-column
                              (neanderthal-mx/neanderthal-matrix->matrix solution)
                              0))]
        (if (anomalies/anomaly? lls)
          lls
          {::curve-fitting-fn      (fn [v]
                                     (let [bf-v (curve-basis-fn v)]
                                       (if (= (count lls) (count bf-v))
                                         (vector/dot-product lls bf-v)
                                         m/nan)))
           ::curve-fitting-weights lls})))))

(s/def ::curve-fitting-fn
  (s/fspec :args (s/cat :v ::vector/vector-finite)
           :ret ::m/number))

(s/def ::curve-fitting-weights ::vector/vector)

(s/def ::curve-basis-fn
  (s/fspec :args (s/cat :v ::vector/vector-finite)
           :ret ::vector/vector-finite))

(s/def ::x-matrix
  (s/and ::mx/matrix-finite
         (fn [m]
           (>= (mx/rows m) 2))))

(s/def ::x-matrix-with-f-vals
  (s/with-gen
    (s/and (s/keys :req [::x-matrix ::f-vals])
           (fn [{::keys [x-matrix f-vals]}]
             (= (count f-vals) (mx/rows x-matrix))))
    #(gen/bind
       (s/gen ::x-matrix)
       (fn [x-mx]
         (gen/bind
           (gen/vector (s/gen ::m/num)
                       (mx/rows x-mx))
           (fn [f-v]
             (gen/return {::x-matrix x-mx
                          ::f-vals   f-v})))))))

(s/def ::matrix-and-vals-with-basis
  (s/merge ::x-matrix-with-f-vals (s/keys :req [::curve-basis-fn])))

(s/fdef linear-least-squares-curve-fitting
        :args (s/cat :matrix-and-vals-with-basis ::matrix-and-vals-with-basis)
        :ret (s/or :sol (s/keys :req [::curve-fitting-fn ::curve-fitting-weights])
                   :anomaly ::anomalies/anomaly))

;;;SMOOTHING CUBIC SPLINES
(s/def ::smoothing-parameter ::m/prob)

(s/def ::coefficients
  (s/coll-of
    (s/double-in :infinite? false
                 :NaN? false)))

(s/def ::variances ::vector/vector-finite+)

(defprotocol Function1
  "Function of one argument"
  (evaluate [this x]
    "Evaluate function at a point."))

(defprotocol Derivative1
  "Derivative of function of one argument"
  (derivative [this]
    "Return instance of function representing derivative"))

(defprotocol DegreesOfFreedom
  (dof [this]))

(defrecord Polynomial [^doubles coefficients n]
  Function1
  (evaluate [_ x]
    (apply + (map *
                  coefficients
                  (take n (iterate #(* x %) 1.0)))))

  Derivative1
  (derivative [_]
    (Polynomial. (->> coefficients
                      (map-indexed *)
                      (drop 1)
                      vec)
                 (dec n))))

(defn polynomial
  ""
  [coefficients]
  (->Polynomial coefficients (count coefficients)))

(s/fdef polynomial
        :args (s/cat :coefficients ::coefficients)
        :ret #(= Polynomial (type %)))

(defn- quincunx
  [^doubles u ^doubles v ^doubles w ^doubles q]
  (let [J (alength u)]
    (aset u 0 0.0)
    (aset v 1 (/ (aget v 1) (aget u 1)))
    (aset w 1 (/ (aget w 1) (aget u 1)))

    (doseq [j (range 2 J)]
      (let [j-1 (- j 1)
            j-2 (- j 2)]
        (aset u j (- (aget u j)
                     (* (aget u j-2) (aget w j-2) (aget w j-2))
                     (* (aget u j-1) (aget v j-1) (aget v j-1))))
        (aset v j (/ (- (aget v j)
                        (* (aget u j-1) (aget v j-1) (aget w j-1)))
                     (aget u j)))
        (aset w j (/ (aget w j) (aget u j)))))

    (aset q 1 (- (aget q 1) (* (aget v 0) (aget q 0))))
    (doseq [j (range 2 J)]
      (let [j-1 (- j 1)
            j-2 (- j 2)]
        (aset q j (- (aget q j)
                     (* (aget v j-1) (aget q j-1))
                     (* (aget w j-2) (aget q j-2))))))
    (doseq [j (range 1 J)]
      (aset q j (/ (aget q j) (aget u j))))

    (aset q (dec J) 0.0)
    (doseq [j (range (- J 3) 0 -1)]
      (aset q j (- (aget q j)
                   (* (aget v j) (aget q (+ j 1)))
                   (* (aget w j) (aget q (+ j 2))))))

    q))
#_(import '(umontreal.ssj.functionfit SmoothingCubicSpline))

(defn- resolver
  [x-array y-array v-array smoothing-parameter]
  (let [N (alength x-array)
        spline-vector (make-array Polynomial (inc N))
        n (dec N)
        h (double-array N)
        r (double-array N)
        u (double-array N)
        v (double-array N)
        w (double-array N)
        q (double-array (inc N))
        sigma (double-array (map #(if (m/non+? %)
                                    1e100
                                    (/ %))
                                 v-array))
        mu (if (m/non+? smoothing-parameter)
             1e100
             (* 2.0
                (/ (m/one- smoothing-parameter)
                   (* 3.0 smoothing-parameter))))]
    (aset h 0 (- (aget x-array 1)
                 (aget x-array 0)))
    (aset r 0 (/ 3.0
                 (aget h 0)))
    (doseq [i (range 1 n)]
      (aset h i (- (aget x-array (inc i))
                   (aget x-array i)))
      (aset r i (/ 3.0 (aget h i)))
      (aset q i (- (* 3.0 (/ (- (aget y-array (inc i))
                                (aget y-array i))
                             (aget h i)))
                   (* 3.0
                      (/ (- (aget y-array i)
                            (aget y-array (dec i)))
                         (aget h (dec i)))))))

    (doseq [i (range 1 n)]
      (aset u i (+ (* (aget r (dec i))
                      (aget r (dec i))
                      (aget sigma (dec i)))
                   (* (+ (aget r (dec i))
                         (aget r i))
                      (+ (aget r (dec i))
                         (aget r i))
                      (aget sigma i))
                   (* (aget r i)
                      (aget r i)
                      (aget sigma (inc i)))))
      (aset u i (+ (* mu
                      (aget u i))
                   (* 2.0
                      (- (aget x-array (inc i))
                         (aget x-array (dec i))))))
      (aset v i (- (- (* (+ (aget r (dec i))
                            (aget r i))
                         (aget r i)
                         (aget sigma i)))
                   (* (aget r i)
                      (+ (aget r i)
                         (aget r (inc i)))
                      (aget sigma (inc i)))))
      (aset v i (+ (* mu
                      (aget v i))
                   (aget h i)))
      (aset w i (* mu
                   (aget r i)
                   (aget r (inc i))
                   (aget sigma (inc i)))))

    (let [q (quincunx u v w q)
          params (double-array 4)
          dd (- (aget y-array 1)
                (* mu
                   (+ (* (- (+ (aget r 0)
                               (aget r 1)))
                         (aget q 1))
                      (* (aget r 1)
                         (aget q 2)))
                   (aget sigma 1)))]
      (aset params 0 (- (aget y-array 0)
                        (* mu
                           (aget r 0)
                           (aget q 1)
                           (aget sigma 0))))
      (aset params 1 (- (/ (- dd
                              (aget params 0))
                           (aget h 0))
                        (/ (* (aget q 1)
                              (aget h 0))
                           3.0)))
      (aset spline-vector 0 (polynomial (vec params)))

      (aset params 3 (/ (aget q 1)
                        (* 3.0
                           (aget h 0))))
      (aset params 2 0.0)
      (aset spline-vector 1 (polynomial (vec params)))

      (doseq [j (range 1 n)]
        (aset params 3 (/ (- (aget q (inc j))
                             (aget q j))
                          (* 3.0
                             (aget h j))))
        (aset params 2 (aget q j))
        (aset params 1 (+ (* (+ (aget q j)
                                (aget q (dec j)))
                             (aget h (dec j)))
                          (get-in (aget spline-vector j) [:coefficients 1])))
        (aset params 0 (+ (* (aget r (dec j))
                             (aget q (dec j)))
                          (* (- (- (aget r (dec j)))
                                (aget r j))
                             (aget q j))
                          (* (aget r j)
                             (aget q (inc j)))))
        (aset params 0 (- (aget y-array j)
                          (* mu
                             (aget params 0)
                             (aget sigma j))))
        (aset spline-vector (inc j) (polynomial (vec params))))

      (aset params 3 0.0)
      (aset params 2 0.0)
      (aset params 1 (-> spline-vector
                         (aget n)
                         derivative
                         (evaluate (- (aget x-array (dec N))
                                      (aget x-array (- N 2))))))
      (aset params 0 (-> spline-vector
                         (aget n)
                         (evaluate (- (aget x-array (dec N))
                                      (aget x-array (- N 2))))))
      (aset spline-vector (inc n) (polynomial (vec params))))

    spline-vector))

(defn- get-fit-polynomial-index
  [x x-vals]
  (or (->> x-vals
           (map-indexed vector)
           (filter (fn [[_ x']]
                     (< x x')))
           ffirst)
      (count x-vals)))

(defn- smoothing-cubic-spline-eigenvalues*
  [x-vals variances]
  (let [x-count (count x-vals)
        n (dec x-count)
        x-diffs (vec (cons (first x-vals)
                           (map -
                                (rest x-vals)
                                x-vals)))
        x2-skip-diffs (vec (cons 0.0
                                 (mapv #(* 2.0 (+ %1 %2))
                                       x-diffs
                                       (rest x-diffs))))
        r (mapv (partial / 3.0) x-diffs)
        x-rev-skip-diffs (vec (cons 0.0
                                    (map #(- (+ %1 %2))
                                         x-diffs
                                         (rest x-diffs))))
        R (native/dsb (dec n)
                      1
                      (interleave (rest x2-skip-diffs) (rest x-diffs)))
        S (linear-algebra/tri
            (native/dtr (:lu (linear-algebra/ptrf! (neanderthal/copy R)))))
        Q' (native/dgb (dec n)
                       x-count
                       0
                       2
                       (interleave r (rest x-rev-skip-diffs) (rest r)))
        S-1*Q' (neanderthal/mm (neanderthal/view-ge S) (native/dge Q'))
        K (neanderthal/mm (native/dgd x-count
                                      (map (partial * (/ 2.0 3.0))
                                           variances))
                          (neanderthal/trans S-1*Q')
                          S-1*Q')
        d (neanderthal-mx/singular-values K)]
    d))

(def smoothing-cubic-spline-eigenvalues
  (memoize smoothing-cubic-spline-eigenvalues*))

(defn smoothing-cubic-spline-dof
  ""
  [{::keys [x-vals variances smoothing-parameter]}]
  (let [d (smoothing-cubic-spline-eigenvalues x-vals variances)]
    (apply +
           (map #(/ (inc (* (m/one- smoothing-parameter)
                            (/ smoothing-parameter)
                            %)))
                d))))

(s/fdef smoothing-cubic-spline-dof
        :args (s/cat :args (s/keys :req [::x-vals
                                         ::variances
                                         ::smoothing-parameter]))
        :ret ::m/finite-non-)

(defrecord SmoothingCubicSpline
  [polynomials x-vals variances smoothing-parameter]
  Function1
  (evaluate [_ x]
    (let [i (get-fit-polynomial-index x x-vals)]
      (if (zero? i)
        (evaluate (polynomials i) (- x (first x-vals)))
        (evaluate (polynomials i) (- x (x-vals (dec i)))))))

  DegreesOfFreedom
  (dof [_]
    (smoothing-cubic-spline-dof
      {::x-vals              x-vals
       ::variances           variances
       ::smoothing-parameter smoothing-parameter})))

(defn smoothing-cubic-spline
  ""
  ([args] (smoothing-cubic-spline args {}))
  ([{::keys [x-vals f-vals smoothing-parameter]}
    {::keys [variances]}]
   (let [variances (or variances (repeat (count x-vals) 1.0))
         polynomials (vec (resolver (double-array x-vals)
                                    (double-array f-vals)
                                    (double-array variances)
                                    smoothing-parameter))]
     (->SmoothingCubicSpline polynomials x-vals variances smoothing-parameter))))

(s/fdef smoothing-cubic-spline
        :args (s/and (s/cat :args (s/merge ::x-vals-with-f-vals
                                           (s/keys :req [::smoothing-parameter]))
                            :opts (s/? (s/keys :opt [::variances])))
                     (fn [{:keys [args opts]}]
                       (let [{::keys [x-vals]} args
                             variances (::variances opts)]
                         (or (nil? variances)
                             (= (count x-vals) (count variances))))))
        :ret #(= SmoothingCubicSpline (type %)))

;;do the x-vals need to be strictly ascending?
;;In DOF, why does the x(0) value affect the DOF? -- currently divide by zero if x(0)=0.0
;;add documentation
;;maybe use spec instead of defprotocols




