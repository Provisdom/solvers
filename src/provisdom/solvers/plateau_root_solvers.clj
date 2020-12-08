(ns provisdom.solvers.plateau-root-solvers
  (:require
    [clojure.spec.alpha :as s]
    [clojure.spec.gen.alpha :as gen]
    [provisdom.math.core :as m]
    [provisdom.math.intervals :as intervals]
    [provisdom.math.vector :as vector]
    [provisdom.utility-belt.anomalies :as anomalies]))

(s/def ::plateau ::m/long)
(s/def ::lower-plateau ::plateau)
(s/def ::upper-plateau ::plateau)
(s/def ::value ::m/finite)
(s/def ::lower-value ::value)
(s/def ::upper-value ::value)
(s/def ::values ::vector/vector-finite)
(s/def ::guess ::value)
(s/def ::lower-guess ::value)
(s/def ::upper-guess ::value)
(s/def ::guesses ::values)
(s/def ::stopped-early? boolean?)
(s/def ::lowered-upper-plateau? boolean?)
(s/def ::reached-max-passes? boolean?)
(s/def ::max-iter (s/spec m/int+? :gen #(s/gen (s/int-in 1 100))))
(s/def ::max-passes (s/spec m/int+? :gen #(s/gen (s/int-in 1 3))))
(s/def ::abs-accu ::m/finite+)
(s/def ::interval ::intervals/finite-interval)
(s/def ::lower-bound ::value)
(s/def ::lower-bounds ::values)

(s/def ::intervals
  (s/coll-of ::interval
    :gen-max 5))

(s/def ::plateau-output (s/tuple ::plateau boolean?))

(s/def ::plateau-f
  (s/fspec :args (s/cat :value ::value)
    :ret ::plateau-output))

(s/def ::plateau-and-value
  (s/keys :req [::plateau ::value]))

(s/def ::upper ::plateau-and-value)
(s/def ::lower ::plateau-and-value)

(defn- plateau-solution-contains-lower-or-upper?
  [{::keys [lower upper]}]
  (or lower upper))

(s/def ::plateau-solution
  (s/and (s/keys :opt [::lower ::stopped-early? ::upper])
    plateau-solution-contains-lower-or-upper?))

(s/def ::plateau-solutions
  (s/coll-of ::plateau-solution
    :gen-max 5))

(s/def ::multi-plateau-f
  (s/with-gen (s/fspec :args (s/cat :values ::values)
                :ret (s/coll-of ::plateau-output))
    #(gen/return
       (fn [v]
         (if (empty? v)
           []
           (let [x (first v)]
             (vec (repeat (count v) [(long (intervals/bound-by-interval
                                             [m/min-long m/max-long]
                                             x))
                                     (> x 5.1)]))))))))

(s/def ::multi-plateau-solution
  (s/keys :req [::plateau-solutions]
    :opt [::reached-max-passes?]))

(defn- plateau-update
  [plateau-f plateau-solution-or-anom x-v]
  (if (anomalies/anomaly? plateau-solution-or-anom)
    plateau-solution-or-anom
    (let [ps plateau-solution-or-anom
          {::keys [lower upper]} ps
          x-v (when x-v (double x-v))
          {l-v ::value
           l   ::plateau} lower
          {u-v ::value
           u   ::plateau} upper]
      (if (and x-v
            (or (not lower) (> x-v l-v))
            (or (not upper) (< x-v u-v)))
        (let [[x bool] (plateau-f x-v)
              #_#_blah (println "A: " bool x x-v ps)]
          (if bool
            (if (and lower (<= x l))
              {::anomalies/message  "bad function"
               ::anomalies/category ::anomalies/error
               ::anomalies/fn       (var plateau-update)}
              (assoc ps ::upper {::plateau x
                                 ::value   x-v}))
            (if (and upper (>= x u))
              {::anomalies/message  "bad function"
               ::anomalies/category ::anomalies/error
               ::anomalies/fn       (var plateau-update)}
              (assoc ps ::lower {::plateau x
                                 ::value   x-v}))))
        ps))))

(defn plateau-root-solver
  "Plateau function takes a finite and returns tuple of [long, boolean], where
  the long represents the plateau. The plateau function must be monotonically
  increasing in the `::interval`. The function must return true when the plateau
  is at or above the root, and false otherwise. Solver stops when it finds the
  minimum plateau that returns true, or when the difference between inputs is
  less than or equal to `::abs-accu`. Returns the last evaluated plateau range
  and input range."
  ([args] (plateau-root-solver args {}))
  ([{::keys [interval plateau-f]}
    {::keys [abs-accu lower-guess guess max-iter upper-guess]
     :or    {abs-accu 1e-6
             max-iter 1000}}]
   (let [f (partial plateau-update plateau-f)
         [l-v u-v] interval
         current {}
         current (f current (when guess
                              (intervals/bound-by-interval interval guess)))
         current (f current (when upper-guess
                              (intervals/bound-by-interval interval upper-guess)))
         current (f current (when lower-guess
                              (intervals/bound-by-interval interval lower-guess)))
         current (f current u-v)
         current (f current l-v)]
     (if (anomalies/anomaly? current)
       current
       (let [{::keys [lower upper]} current]
         (if (or (not lower) (not upper))
           current
           (loop [current current
                  i 0]
             (let [{::keys [lower upper]} current
                   {l-v ::value
                    l   ::plateau} lower
                   {u-v ::value
                    u   ::plateau} upper
                   done? (m/one? (- (double u) l))]
               (if (or done?
                     (>= i max-iter)
                     (<= (- u-v l-v) abs-accu))
                 (cond-> current
                   (not done?) (assoc ::stopped-early? true))
                 (let [new (f current (* 0.5 (+ l-v u-v)))]
                   (if (anomalies/anomaly? new)
                     new
                     (recur new (inc i)))))))))))))

(s/fdef plateau-root-solver
  :args (s/cat :args (s/keys :req [::interval ::plateau-f])
          :opts (s/? (s/keys :opt [::abs-accu
                                   ::lower-guess
                                   ::guess
                                   ::max-iter
                                   ::upper-guess])))
  :ret (s/or :plateau-solution ::plateau-solution
         :anomaly ::anomalies/anomaly))

(defn tighten-plateau-solution
  "Tightens the `::plateau-solution` within `::abs-accu`. If a new solution is
  found in which the upper plateau can be decreased, the optional
  `::lowered-upper-plateau?` will be true."
  ([args] (tighten-plateau-solution args {}))
  ([{::keys [plateau-f plateau-solution]}
    {::keys [abs-accu max-iter]
     :or    {abs-accu 1e-6
             max-iter 1000}}]
   (let [f (partial plateau-update plateau-f)
         {::keys [lower upper]} plateau-solution
         orig-u (::plateau upper)]
     (if (or (not lower) (not upper))
       plateau-solution
       (loop [current plateau-solution
              i 0]
         (let [{::keys [lower upper]} current
               {l-v ::value} lower
               {u-v ::value
                u   ::plateau} upper
               done? (<= (- u-v l-v) abs-accu)]
           (if (or done? (>= i max-iter))
             (cond-> current
               (not done?) (assoc ::stopped-early? true)
               (not= orig-u u) (assoc ::lowered-upper-plateau? true))
             (let [new (f current (* 0.5 (+ l-v u-v)))]
               (if (anomalies/anomaly? new)
                 new
                 (recur new (inc i)))))))))))

(s/fdef tighten-plateau-solution
  :args (s/cat :args (s/keys :req [::plateau-f ::plateau-solution])
          :opts (s/? (s/keys :opt [::abs-accu ::max-iter])))
  :ret (s/or :plateau-solution ::plateau-solution
         :anomaly ::anomalies/anomaly))

(defn single-pass
  [multi-plateau-f
   lower-values
   last-values
   upper-bounds
   abs-accu
   max-iter
   guesses]
  (reduce (fn [[low-values values] dim]
            (let [p-f #(nth (multi-plateau-f (assoc values dim %)) dim)
                  last-value (nth values dim)
                  guess (nth guesses dim)
                  sol (plateau-root-solver
                        {::interval  [last-value (nth upper-bounds dim)]
                         ::plateau-f p-f}
                        {::abs-accu abs-accu
                         ::guess    guess
                         ::max-iter max-iter})]
              (if (anomalies/anomaly? sol)
                (reduced sol)
                (let [l-v (or (::value (::lower sol)) (nth lower-values dim))
                      v (or (::value (::upper sol)) (nth upper-bounds dim))]
                  [(assoc low-values dim l-v) (assoc values dim v)]))))
    [lower-values last-values]
    (range (count last-values))))

(defn multi-dimensional-plateau-root-solver
  "Multi-dimensional plateau function takes a vector-finite and returns vector
  of tuples of [long, boolean], where the long represents the plateau for that
  dimension. The plateau function must be monotonically increasing in the
  `::intervals`, even across dimensions. For an input vector, a change to a
  dimension should only change the output of other dimensions if the output in
  the original dimension changes plateaus. Each tuple must return true when the
  plateau is at or above the root for that dimension, and false otherwise.
  Solver stops when, for every dimension, it finds the minimum plateau that
  returns true, or when the difference between inputs is less than or equal to
  `::abs-accu`. Returns the last evaluated plateau range and input range for
  each dimension."
  ([args] (multi-dimensional-plateau-root-solver args {}))
  ([{::keys [intervals multi-plateau-f]}
    {::keys [abs-accu guesses max-iter max-passes]
     :or    {abs-accu   1e-6
             max-iter   1000
             max-passes 10}}]
   (let [lower-bounds (mapv (comp double first) intervals)
         upper-bounds (mapv (comp double second) intervals)
         guesses (if (and guesses
                       (= (count guesses) (count lower-bounds)))
                   guesses
                   lower-bounds)
         sol (single-pass
               multi-plateau-f
               lower-bounds
               lower-bounds
               upper-bounds
               abs-accu
               max-iter
               guesses)
         var-f (var multi-dimensional-plateau-root-solver)]
     (loop [sol sol
            previous-last-values nil
            i 1]
       (let [[lower-values last-values] sol
             done? (= last-values previous-last-values)]
         (if (or done? (>= i max-passes))
           (let [pss (mapv
                       (fn [dim]
                         (let [p-f #(nth (multi-plateau-f
                                           (assoc last-values dim %))
                                      dim)
                               l-v (nth lower-values dim)
                               u-v (nth last-values dim)
                               [l l-bool] (when l-v (p-f l-v))
                               [u u-bool] (when u-v (p-f u-v))]
                           (if (or l-bool (false? u-bool))
                             {::anomalies/message  "bad function"
                              ::anomalies/category ::anomalies/error
                              ::anomalies/fn       var-f}
                             (cond-> {}
                               l (assoc ::lower {::plateau l
                                                 ::value   l-v})
                               u (assoc ::upper {::plateau u
                                                 ::value   u-v})))))
                       (range (count last-values)))]
             (or (first (filter anomalies/anomaly? pss))
               (cond-> {::plateau-solutions pss}
                 (not done?) (assoc ::reached-max-passes? true))))
           (let [sol (single-pass
                       multi-plateau-f
                       lower-values
                       last-values
                       upper-bounds
                       abs-accu
                       max-iter
                       last-values)]
             (if (anomalies/anomaly? sol)
               sol
               (recur sol last-values (inc i))))))))))

(s/fdef multi-dimensional-plateau-root-solver
  :args (s/cat :args (s/keys :req [::intervals ::multi-plateau-f])
          :opts (s/? (s/keys :opt [::abs-accu
                                   ::guesses
                                   ::max-iter
                                   ::max-passes])))
  :ret (s/or :multi-plateau-solution ::multi-plateau-solution
         :anomaly ::anomalies/anomaly))

(defn tighten-multi-plateau-solution
  "Tightens the `::multi-plateau-solution` within `::abs-accu`. If tightening
  lowers the upper plateau of any dimension, an anomaly will be returned because
  changing the upper plateau could drastically change the solution in the other
  dimensions. To avoid the possibility of this, ensure the `::abs-accu`
  is not smaller than that used for the original solution."
  ([args] (tighten-multi-plateau-solution args {}))
  ([{::keys [multi-plateau-f multi-plateau-solution]}
    {::keys [abs-accu max-iter]
     :or    {abs-accu 1e-6
             max-iter 1000}}]
   (let [pss (::plateau-solutions multi-plateau-solution)
         values (mapv (fn [ps] (or (::value (::upper ps))
                                 (::value (::lower ps))))
                  pss)
         message (str "Tightening lowered the upper plateau. There is a better"
                   " solution than the original multi-plateau-solution.")
         npss (mapv
                (fn [dim]
                  (let [p-f #(nth (multi-plateau-f (assoc values dim %)) dim)
                        plateau-solution (nth pss dim)
                        new-ps (tighten-plateau-solution
                                 {::plateau-f        p-f
                                  ::plateau-solution plateau-solution}
                                 {::abs-accu abs-accu
                                  ::max-iter max-iter})]
                    (if (::lowered-upper-plateau? new-ps)
                      {::anomalies/category ::anomalies/error
                       ::anomalies/message  message
                       ::anomalies/fn       (var tighten-multi-plateau-solution)}
                      new-ps)))
                (range (count pss)))]
     (or (first (filter anomalies/anomaly? npss))
       (assoc multi-plateau-solution ::plateau-solutions npss)))))

(s/fdef tighten-multi-plateau-solution
  :args (s/cat :args (s/keys :req [::multi-plateau-f ::multi-plateau-solution])
          :opts (s/? (s/keys :opt [::abs-accu ::max-iter])))
  :ret (s/or :multi-plateau-solution ::multi-plateau-solution
         :anomaly ::anomalies/anomaly))