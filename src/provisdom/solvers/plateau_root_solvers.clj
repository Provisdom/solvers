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

(s/def ::plateau-solution
  (s/keys :req [::upper-plateau ::upper-value]
    :opt [::lower-plateau ::lower-value ::stopped-early?]))

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
  [plateau-f args x-v]
  (if (anomalies/anomaly? args)
    args
    (let [[l-v l u-v u] args
          x-v (when x-v (double x-v))]
      (if (and x-v
            (or (> x-v l-v) (and (not l) (== x-v l-v)))
            (or (< x-v u-v) (and (not u) (== x-v u-v))))
        (let [[x bool] (plateau-f x-v)
              #_#_blah (println bool x x-v args)]
          (if bool
            (if (and l (<= x l))
              {::anomalies/message  "bad function"
               ::anomalies/category ::anomalies/error
               ::anomalies/fn       (var plateau-update)}
              [l-v l x-v x])
            (if (and u (>= x u))
              {::anomalies/message  "bad function"
               ::anomalies/category ::anomalies/error
               ::anomalies/fn       (var plateau-update)}
              [x-v x u-v u])))
        args))))

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
         current [l-v nil u-v nil]
         current (f current guess)
         current (f current upper-guess)
         current (f current lower-guess)
         current (f current u-v)
         current (f current l-v)]
     (if (anomalies/anomaly? current)
       current
       (let [[l-v l u-v u] current]
         (cond (not l) {::upper-plateau u
                        ::upper-value   u-v}
               (not u) {::upper-plateau l
                        ::upper-value   l-v}
               :else (loop [[l-v l u-v u] current
                            i 0]
                       (if (or (m/one? (- (double u) l))
                             (>= i max-iter)
                             (<= (- u-v l-v) abs-accu))
                         (cond-> {::lower-plateau l
                                  ::lower-value   l-v
                                  ::upper-plateau u
                                  ::upper-value   u-v}

                           (not (m/one? (- (double u) l)))
                           (assoc ::stopped-early? true))
                         (let [a (f [l-v l u-v u] (* 0.5 (+ l-v u-v)))]
                           (if (anomalies/anomaly? a)
                             a
                             (recur a (inc i))))))))))))

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
  "Tightens the `::plateau-solution` within `::abs-accu`."
  ([args] (tighten-plateau-solution args {}))
  ([{::keys [lower-bound plateau-f plateau-solution]}
    {::keys [abs-accu max-iter]
     :or    {abs-accu 1e-6
             max-iter 1000}}]
   (let [f (partial plateau-update plateau-f)
         ps plateau-solution
         {::keys [lower-plateau lower-value upper-plateau upper-value]} ps
         [l-v l] (if (or (not lower-value) (< lower-bound lower-value))
                   [lower-bound nil]
                   [lower-value lower-plateau])
         current [l-v l upper-value upper-plateau]
         current (f current l-v)]
     (if-not l
       {::upper-plateau upper-plateau
        ::upper-value   upper-value}
       (loop [[l-v l u-v u] current
              i 0]
         (if (or (>= i max-iter) (<= (- u-v l-v) abs-accu))
           (cond-> {::lower-plateau l
                    ::lower-value   l-v
                    ::upper-plateau u
                    ::upper-value   u-v}
             (= i max-iter) (assoc ::stopped-early? true))
           (let [a (f [l-v l u-v u] (* 0.5 (+ l-v u-v)))]
             (if (anomalies/anomaly? a)
               a
               (recur a (inc i))))))))))

(s/fdef tighten-plateau-solution
  :args (s/cat :args (s/keys :req [::lower-bound ::plateau-f ::plateau-solution])
          :opts (s/? (s/keys :opt [::abs-accu ::max-iter])))
  :ret (s/or :plateau-solution ::plateau-solution
         :anomaly ::anomalies/anomaly))

(defn single-pass
  [multi-plateau-f
   guesses
   uppers
   abs-accu
   max-iter
   start-plateau-solutions]
  (reduce (fn [plateau-solutions dim]
            (let [v (mapv ::upper-value plateau-solutions)
                  p-f #(nth (multi-plateau-f (assoc v dim %)) dim)
                  last-value (nth v dim)
                  guess (get guesses dim last-value)
                  sol (plateau-root-solver
                        {::interval  [last-value (nth uppers dim)]
                         ::plateau-f p-f}
                        {::abs-accu abs-accu
                         ::guess    guess
                         ::max-iter max-iter})]
              (if (anomalies/anomaly? sol)
                (reduced sol)
                (let [lv (or (::lower-value sol)
                           (::lower-value (nth plateau-solutions dim)))
                      lp (or (::lower-plateau sol)
                           (::lower-plateau (nth plateau-solutions dim)))
                      sol (cond-> sol
                            lv (assoc ::lower-value lv)
                            lp (assoc ::lower-plateau lp))]
                  (assoc plateau-solutions dim sol)))))
    start-plateau-solutions
    (range (count start-plateau-solutions))))

(defn multi-dimensional-plateau-root-solver
  "Multi-dimensional plateau function takes a vector-finite and returns vector
  of tuples of [long, boolean], where the long represents the plateau for that
  dimension. The plateau function must be monotonically increasing in the
  `::intervals`, even across dimensions. Each tuple must return true when the
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
   (let [lowers (mapv (comp double first) intervals)
         uppers (mapv (comp double second) intervals)
         plateau-solutions (single-pass
                             multi-plateau-f
                             guesses
                             uppers
                             abs-accu
                             max-iter
                             (mapv (fn [l] {::upper-value l}) lowers))]
     (loop [plateau-solutions plateau-solutions
            last-upper-values nil
            i 1]
       (let [upper-values (mapv ::upper-value plateau-solutions)]
         (if (or (= upper-values last-upper-values)
               (>= i max-passes))
           (cond-> {::plateau-solutions plateau-solutions}

             (not= upper-values last-upper-values)
             (assoc ::reached-max-passes? true))

           (let [new-plateau-solutions (single-pass
                                         multi-plateau-f
                                         upper-values
                                         uppers
                                         abs-accu
                                         max-iter
                                         plateau-solutions)]
             (if (anomalies/anomaly? new-plateau-solutions)
               new-plateau-solutions
               (recur new-plateau-solutions upper-values (inc i))))))))))

(s/fdef multi-dimensional-plateau-root-solver
  :args (s/cat :args (s/keys :req [::intervals ::multi-plateau-f])
          :opts (s/? (s/keys :opt [::abs-accu
                                   ::guesses
                                   ::max-iter
                                   ::max-passes])))
  :ret (s/or :multi-plateau-solution ::multi-plateau-solution
         :anomaly ::anomalies/anomaly))

(defn tighten-single-pass
  [multi-plateau-f
   lower-bounds
   abs-accu
   max-iter
   start-plateau-solutions]
  (reduce (fn [plateau-solutions dim]
            (let [v (mapv ::upper-value plateau-solutions)
                  p-f #(nth (multi-plateau-f (assoc v dim %)) dim)
                  plateau-solution (nth plateau-solutions dim)
                  sol (tighten-plateau-solution
                        {::lower-bound      (get lower-bounds dim m/min-dbl)
                         ::plateau-f        p-f
                         ::plateau-solution plateau-solution}
                        {::abs-accu abs-accu
                         ::max-iter max-iter})]
              (if (anomalies/anomaly? sol)
                (reduced sol)
                (let [lv (or (::lower-value sol)
                           (::lower-value (nth plateau-solutions dim)))
                      lp (or (::lower-plateau sol)
                           (::lower-plateau (nth plateau-solutions dim)))
                      sol (cond-> sol
                            lv (assoc ::lower-value lv)
                            lp (assoc ::lower-plateau lp))]
                  (assoc plateau-solutions dim sol)))))
    start-plateau-solutions
    (range (count start-plateau-solutions))))

(defn tighten-multi-plateau-solution
  "Tightens the `::multi-plateau-solution` within `::abs-accu`."
  ([args] (tighten-multi-plateau-solution args {}))
  ([{::keys [lower-bounds multi-plateau-f multi-plateau-solution]}
    {::keys [abs-accu max-iter max-passes]
     :or    {abs-accu   1e-6
             max-iter   1000
             max-passes 10}}]
   (loop [plateau-solutions (::plateau-solutions multi-plateau-solution)
          last-upper-values nil
          i 0]
     (let [upper-values (mapv ::upper-value plateau-solutions)]
       (if (or (= upper-values last-upper-values)
             (>= i max-passes))
         (cond-> {::plateau-solutions plateau-solutions}

           (not= upper-values last-upper-values)
           (assoc ::reached-max-passes? true))

         (let [new-plateau-solutions (tighten-single-pass
                                       multi-plateau-f
                                       lower-bounds
                                       abs-accu
                                       max-iter
                                       plateau-solutions)]
           (if (anomalies/anomaly? new-plateau-solutions)
             new-plateau-solutions
             (recur new-plateau-solutions upper-values (inc i)))))))))

(s/fdef tighten-multi-plateau-solution
  :args (s/cat :args (s/keys :req [::lower-bounds
                                   ::multi-plateau-f
                                   ::multi-plateau-solution])
          :opts (s/? (s/keys :opt [::abs-accu ::max-iter])))
  :ret (s/or :multi-plateau-solution ::multi-plateau-solution
         :anomaly ::anomalies/anomaly))