(ns provisdom.solvers.curve-fitting.curve-fitting
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
    [provisdom.neanderthal-matrix.neanderthal-matrix :as neanderthal-mx]))

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

(s/def ::x-vals (s/coll-of ::m/finite
                           :kind clojure.core/vector?
                           :into []
                           :min-count 4))

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

