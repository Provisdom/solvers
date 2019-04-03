(ns provisdom.solvers.curve-fitting.smoothing-spline
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

;;;SMOOTHING CUBIC SPLINES
;;; Mostly ported from https://github.com/umontreal-simul/ssj/blob/master/src/main/java/umontreal/ssj/functionfit/SmoothingCubicSpline.java
(s/def ::smoothing-parameter ::m/open-prob)

(s/def ::coefficients (s/coll-of ::m/double-finite))

(s/def ::variances ::vector/vector-finite+)

(s/def ::x-vals
  (s/with-gen
    (s/and (s/coll-of ::m/finite
                      :kind clojure.core/vector?
                      :into []
                      :min-count 4)
           (fn [x] (->> x
                        (partition 2)
                        (every? (fn [[x1 x2]] (> x2 x1))))))
    #(gen/bind
       (s/gen (s/coll-of (s/double-in :infinite? false :NaN? false) :min-count 4))
       (fn [x]
         (gen/return (->> x
                          set
                          sort
                          vec))))))

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
  (s/assert ::coefficients coefficients)
  (->Polynomial coefficients (count coefficients)))

#_(s/fdef polynomial
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

;;; See http://www.physics.muni.cz/~jancely/NM/Texty/Numerika/CubicSmoothingSpline.pdf
;;; pages 16-21. All bindings are named to match the paper, please don't change them.
(defn- smoothing-cubic-spline-eigenvalues*
  [xs vs]
  (let [N (count xs)
        n (dec N)
        h (->> xs
               (map - (->> xs (drop 1)))
               (#(cons (first %) %))
               vec)
        p (->> h
               (drop 1)
               (map #(* 2.0 (+ %1 %2)) h)
               (cons 0.0)
               vec)
        r (map (partial / 3.0) h)
        f (->> r
               (drop 1)
               (map #(- (+ %1 %2)) r)
               (cons 0.0)
               vec)
        R (native/dsb (- n 1) 1
                      (interleave (drop 1 p) (drop 1 h)))
        S (linear-algebra/tri
            (native/dtr (:lu (linear-algebra/ptrf! (neanderthal/copy R)))))
        Q' (native/dgb (dec n)
                       N
                       0
                       2
                       (interleave r (drop 1 f) (drop 1 r)))
        S-1*Q' (neanderthal/mm (neanderthal/view-ge S) (native/dge Q'))
        K (neanderthal/mm (native/dgd N
                                      (map (partial * (/ 2.0 3.0)) vs))
                          (neanderthal/trans S-1*Q')
                          S-1*Q')
        d (neanderthal-mx/singular-values K)]
    d))

(def smoothing-cubic-spline-eigenvalues
  (memoize smoothing-cubic-spline-eigenvalues*))

(defn smoothing-cubic-spline-dof
  ""
  [{::keys [x-vals variances smoothing-parameter]}]
  (anomalies/anomalous-let [d (smoothing-cubic-spline-eigenvalues x-vals variances)]
    (apply +
           (map #(/ (inc (* (m/one- smoothing-parameter)
                            (/ smoothing-parameter)
                            %)))

                d))))

(s/def ::dof (s/or :val ::m/finite-non-
                   :nan ::m/nan
                   :anomaly ::anomalies/anomaly))

(s/fdef smoothing-cubic-spline-dof
  :args (s/cat :args (s/keys :req [::x-vals
                                   ::variances
                                   ::smoothing-parameter]))
  :ret ::dof)

#_(defrecord SmoothingCubicSpline
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
  "Returns a map of a function `::smoothing-cubic-spline-fn` that takes and
   returns a finite and a finite+ `::dof`.
  `curve-basis-fn` takes and returns a vector."
  ([args] (smoothing-cubic-spline args {}))
  ([{::keys [x-vals f-vals smoothing-parameter]}
    {::keys [variances]}]
   (let [variances (or variances (vec (repeat (count x-vals) 1.0)))
         polynomials (vec (resolver (double-array x-vals)
                                    (double-array f-vals)
                                    (double-array variances)
                                    smoothing-parameter))]
     {::smoothing-cubic-spline-fn
      (fn [x]
        (let [i (get-fit-polynomial-index x x-vals)]
          (if (zero? i)
            (evaluate (polynomials i) (- x (first x-vals)))
            (evaluate (polynomials i) (- x (x-vals (dec i)))))))
      ::dof
      ((memoize
         (fn []
           (smoothing-cubic-spline-dof
             {::x-vals              x-vals
              ::variances           variances
              ::smoothing-parameter smoothing-parameter}))))})))

(s/def ::smoothing-cubic-spline-fn
  (s/fspec :args (s/cat :x ::m/finite)
           :ret ::m/finite))

(s/fdef smoothing-cubic-spline
  :args (s/and (s/cat :args (s/merge ::x-vals-with-f-vals
                                     (s/keys :req [::smoothing-parameter]))
                      :opts (s/? (s/keys :opt [::variances])))
               (fn [{:keys [args opts]}]
                 (let [{::keys [x-vals]} args
                       variances (::variances opts)]
                   (or (nil? variances)
                       (= (count x-vals) (count variances))))))
  :ret (s/keys :req [::smoothing-cubic-spline-fn ::dof])
  #_#(= SmoothingCubicSpline (type %)))




