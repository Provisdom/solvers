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
(s/def ::lower-values ::values)
(s/def ::last-values ::values)
(s/def ::guess ::value)
(s/def ::lower-guess ::value)
(s/def ::upper-guess ::value)
(s/def ::guesses ::values)
(s/def ::lower-guesses ::guesses)
(s/def ::upper-guesses ::guesses)
(s/def ::stopped-early? boolean?)
(s/def ::lowered-upper-plateau? boolean?)
(s/def ::reached-max-passes? boolean?)
(s/def ::max-iter (s/spec m/int+? :gen #(s/gen (s/int-in 1 100))))
(s/def ::max-passes (s/spec m/int+? :gen #(s/gen (s/int-in 1 3))))
(s/def ::abs-accu ::m/finite+)
(s/def ::interval ::intervals/finite-interval)
(s/def ::lower-bound ::value)
(s/def ::lower-bounds ::values)
(s/def ::upper-bounds ::values)

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

(s/def ::plateau-solution
  (s/keys :req [(or ::lower ::upper)]
    :opt [::stopped-early?]))

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

(s/def ::single-pass-solution
  (s/keys :req [::last-values ::lower-values]))

(s/def ::single-pass-solutions
  (s/coll-of ::single-pass-solution
    :gen-max 5))

(s/def ::multi-plateau-solution
  (s/keys :req [::plateau-solutions ::single-pass-solutions]
    :opt [::reached-max-passes?]))

(s/def ::dimension ::m/int-non-)
(s/def ::partition-size ::m/int+)

;;;UNIVARIATE
(defn plateau-update
  "Updates the `::plateau-solution` by testing the `::value`."
  [{::keys [plateau-f plateau-solution value]}]
  (let [{::keys [lower upper]} plateau-solution
        x-v (double value)
        {l-v ::value
         l   ::plateau} lower
        {u-v ::value
         u   ::plateau} upper
        anom {::anomalies/message  "bad function"
              ::anomalies/category ::anomalies/error
              ::anomalies/fn       (var plateau-update)}]
    (if (and (or (not lower) (> x-v l-v))
          (or (not upper) (< x-v u-v)))
      (let [[x bool] (plateau-f x-v)]
        (if bool
          (if (and lower (<= x l))
            anom
            (assoc plateau-solution ::upper {::plateau x
                                             ::value   x-v}))
          (if (and upper (>= x u))
            anom
            (assoc plateau-solution ::lower {::plateau x
                                             ::value   x-v}))))
      plateau-solution)))

(s/fdef plateau-update
  :args (s/cat :args (s/keys :req [::plateau-f ::plateau-solution ::value]))
  :ret (s/or :plateau-solution ::plateau-solution
         :anomaly ::anomalies/anomaly))

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
   (let [f (fn [z value]
             (cond (anomalies/anomaly? z) z
                   (nil? value) z
                   (empty? z) (let [[p bool] (plateau-f value)]
                                (assoc z (if bool ::upper ::lower)
                                         {::plateau p
                                          ::value   value}))
                   :else (plateau-update
                           {::plateau-f        plateau-f
                            ::plateau-solution z
                            ::value            value})))
         [l-v u-v] interval
         current (-> {}
                   (f (when guess
                        (intervals/bound-by-interval interval guess)))
                   (f (when upper-guess
                        (intervals/bound-by-interval interval upper-guess)))
                   (f (when lower-guess
                        (intervals/bound-by-interval interval lower-guess)))
                   (f u-v)
                   (f l-v))]
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
   (let [{::keys [lower upper]} plateau-solution
         orig-u (::plateau upper)]
     (if (or (not lower) (not upper))
       plateau-solution
       (loop [current plateau-solution
              i 0]
         (let [{::keys [lower upper]} current
               {l-v ::value} lower
               {u-v ::value
                u   ::plateau} upper
               done? (<= (- u-v (double l-v)) abs-accu)]
           (if (or done? (>= i max-iter))
             (cond-> current
               (not done?) (assoc ::stopped-early? true)
               (not= orig-u u) (assoc ::lowered-upper-plateau? true))
             (let [new (plateau-update
                         {::plateau-f        plateau-f
                          ::plateau-solution current
                          ::value            (* 0.5 (+ l-v u-v))})]
               (if (anomalies/anomaly? new)
                 new
                 (recur new (inc i)))))))))))

(s/fdef tighten-plateau-solution
  :args (s/cat :args (s/keys :req [::plateau-f ::plateau-solution])
          :opts (s/? (s/keys :opt [::abs-accu ::max-iter])))
  :ret (s/or :plateau-solution ::plateau-solution
         :anomaly ::anomalies/anomaly))

;;;MULTIVARIATE
(defn single-dimension
  "Calculates the changes for just a single dimension."
  [{::keys [abs-accu
            guesses
            last-values
            lower-guesses
            max-iter
            multi-plateau-f
            upper-bounds
            upper-guesses]}
   dim]
  (let [p-f #(get (multi-plateau-f (assoc last-values dim %)) dim [0 true])
        upper-bound (get upper-bounds dim m/max-dbl)
        last-value (min (get last-values dim 0.0) upper-bound)
        guess (get guesses dim 0.0)
        lower-guess (when lower-guesses (get lower-guesses dim nil))
        upper-guess (when upper-guesses (get upper-guesses dim nil))]
    (plateau-root-solver
      {::interval  [last-value upper-bound]
       ::plateau-f p-f}
      (cond-> {::abs-accu abs-accu
               ::guess    guess
               ::max-iter max-iter}
        lower-guess (assoc ::lower-guess lower-guess)
        upper-guess (assoc ::upper-guess upper-guess)))))

(s/fdef single-dimension
  :args (s/cat :args (s/keys :req [::abs-accu
                                   ::guesses
                                   ::max-iter
                                   ::multi-plateau-f
                                   ::upper-bounds]
                       :opt [::lower-guesses ::upper-guesses])
          :dim ::dimension)
  :ret (s/or :plateau-solution ::plateau-solution
         :anomaly ::anomalies/anomaly))

(defn single-pass
  "Runs one pass to update the `::last-values` and `::lower-values`. Can set
  `::partition-size` to greater than 1 to run chunks in parallel."
  [{::keys [last-values lower-values partition-size upper-bounds]
    :or    {partition-size 1}
    :as    args}]
  (let [dims (count last-values)
        dim-blocks (partition-all partition-size (range dims))
        fmap (if (= partition-size 1) map pmap)]
    (reduce
      (fn [{values ::last-values
            :as    acc}
           dim-block]
        (if (anomalies/anomaly? acc)
          (reduced acc)
          (let [dims+sols (fmap (fn [dim]
                                  [dim
                                   (single-dimension
                                     (assoc args ::last-values values)
                                     dim)])
                            dim-block)]
            (reduce
              (fn [acc2 [dim sol]]
                (if (anomalies/anomaly? sol)
                  (reduced sol)
                  (let [l-v (or (::value (::lower sol))
                              (get lower-values dim 0.0))
                        v (or (::value (::upper sol))
                            (get upper-bounds dim 0.0))]
                    (-> acc2
                      (assoc-in [::last-values dim] v)
                      (assoc-in [::lower-values dim] l-v)))))
              acc
              dims+sols))))
      {::last-values  last-values
       ::lower-values lower-values}
      dim-blocks)))

(s/fdef single-pass
  :args (s/cat :args (s/keys :req [::abs-accu
                                   ::guesses
                                   ::last-values
                                   ::lower-values
                                   ::max-iter
                                   ::multi-plateau-f
                                   ::upper-bounds]
                       :opt [::lower-guesses ::partition-size ::upper-guesses]))
  :ret (s/or :single-pass-solution ::single-pass-solution
         :anomaly ::anomalies/anomaly))

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
    {::keys [abs-accu
             lower-guesses
             guesses
             max-iter
             max-passes
             partition-size
             upper-guesses]
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
               (cond-> {::abs-accu        abs-accu
                        ::guesses         guesses
                        ::last-values     lower-bounds
                        ::lower-values    lower-bounds
                        ::max-iter        max-iter
                        ::multi-plateau-f multi-plateau-f
                        ::upper-bounds    upper-bounds}
                 lower-guesses (assoc ::lower-guesses lower-guesses)
                 partition-size (assoc ::partition-size partition-size)
                 upper-guesses (assoc ::upper-guesses upper-guesses)))]
     (if (anomalies/anomaly? sol)
       sol
       (loop [sols [sol]
              previous-last-values nil
              i 1]
         (let [{::keys [lower-values last-values]} (peek sols)
               done? (= last-values previous-last-values)]
           (if (or done? (>= i max-passes))
             (let [pss (mapv (fn [dim]
                               (let [p-f #(get (multi-plateau-f
                                                 (assoc last-values dim %))
                                            dim
                                            [0 true])
                                     l-v (get lower-values dim)
                                     u-v (get last-values dim)
                                     [l l-bool] (when l-v (p-f l-v))
                                     [u u-bool] (when (and u-v (not= u-v l-v))
                                                  (p-f u-v))
                                     [l l-v u u-v] (if l-bool
                                                     [nil nil l l-v]
                                                     [l l-v u u-v])]
                                 (cond-> {}
                                   l (assoc ::lower {::plateau l
                                                     ::value   l-v})
                                   u (assoc ::upper {::plateau u
                                                     ::value   u-v}))))
                         (range (count last-values)))]
               (cond-> {::plateau-solutions     pss
                        ::single-pass-solutions sols}
                 (not done?) (assoc ::reached-max-passes? true)))
             (let [sol (single-pass
                         {::abs-accu        abs-accu
                          ::guesses         last-values
                          ::last-values     last-values
                          ::lower-values    lower-values
                          ::max-iter        max-iter
                          ::multi-plateau-f multi-plateau-f
                          ::upper-bounds    upper-bounds})]
               (if (anomalies/anomaly? sol)
                 sol
                 (recur (conj sols sol) last-values (inc i)))))))))))

(s/fdef multi-dimensional-plateau-root-solver
  :args (s/cat :args (s/keys :req [::intervals ::multi-plateau-f])
          :opts (s/? (s/keys :opt [::abs-accu
                                   ::guesses
                                   ::lower-guesses
                                   ::max-iter
                                   ::max-passes
                                   ::partition-size
                                   ::upper-guesses])))
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
                  (let [p-f #(get (multi-plateau-f (assoc values dim %))
                               dim
                               [0 true])
                        plateau-solution (get pss dim {::lower {::plateau 0
                                                                ::value   0.0}})
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
