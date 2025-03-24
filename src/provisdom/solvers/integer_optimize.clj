(ns provisdom.solvers.integer-optimize
  (:require
    [clojure.spec.alpha :as s]
    [clojure.spec.gen.alpha :as gen]
    [provisdom.math.core :as m]
    [provisdom.math.intervals :as intervals]
    [provisdom.utility-belt.anomalies :as anomalies]))

(s/def ::optimize-f
  (s/with-gen
    (s/fspec :args (s/cat :n ::m/long)
      :ret ::m/finite)
    #(gen/one-of
       (map gen/return
         (list (fn [n]
                 (let [v (m/pow n -2)]
                   (if (m/finite? v)
                     v
                     m/min-dbl))))))))

(s/def ::point-and-val
  (s/and (s/tuple ::m/int ::m/finite)
    (fn [[point _]]
      (>= point 2))))

(defn- recursive-bounded
  "`best` should be between `worse` and `better` but have a higher `best-val`
  than either. `best` can be nil."
  [adj-f [worse worse-val] [better better-val] [best best-val]]
  (loop [[worse worse-val] [worse worse-val]
         [better better-val] [better better-val]
         [best best-val] [best best-val]]
    (cond (and (nil? best)
            (m/one? (m/abs (- better worse)))) better

          (nil? best)
          (let [new (max (m/floor' (m/cbrt (* (double worse) better better)))
                      (inc (min worse better)))
                new-val (adj-f new)]
            (cond (<= new-val worse-val)
                  {::anomalies/category        ::anomalies/error
                   ::anomalies/fn              (var recursive-bounded)
                   ::anomalies/solver-category ::anomalies/bad-supplied-function
                   ::anomalies/data            [[new new-val]
                                                [worse worse-val]
                                                [better better-val]]}

                  (> new-val better-val) (recur
                                           [worse worse-val]
                                           [better better-val]
                                           [new new-val])
                  :else (recur
                          [new new-val]
                          [better better-val]
                          nil)))

          (m/one? (m/abs (- better best))) (recur
                                             [worse worse-val]
                                             [best best-val]
                                             nil)
          (m/one? (m/abs (- worse best))) (recur
                                            [better better-val]
                                            [best best-val]
                                            nil)
          :else (let [new (max (m/floor' (m/cbrt (* (double better) best best)))
                            (inc (min better best)))
                      new-val (adj-f new)]
                  (cond (= new-val best-val) (recur
                                               [new new-val]
                                               [best best-val]
                                               nil)

                        (<= new-val worse-val)
                        {::anomalies/category        ::anomalies/error
                         ::anomalies/fn              (var recursive-bounded)
                         ::anomalies/solver-category ::anomalies/bad-supplied-function
                         ::anomalies/data            [[new new-val]
                                                      [worse worse-val]
                                                      [better better-val]
                                                      [best best-val]]}

                        (> new-val best-val) (recur
                                               [better better-val]
                                               [best best-val]
                                               [new new-val])
                        :else (recur
                                [worse worse-val]
                                [new new-val]
                                [best best-val]))))))

(s/fdef recursive-bounded
  :args (s/and (s/cat :adj-f ::optimize-f
                 :worse-tuple ::point-and-val
                 :better-tuple ::point-and-val
                 :best-tuple (s/nilable ::point-and-val))
          (fn [{:keys [worse-tuple better-tuple best-tuple]}]
            (let [[worse worse-val] worse-tuple
                  [better better-val] better-tuple
                  [best best-val] best-tuple]
              (and (>= better-val worse-val)
                (not= worse better)
                (or (not best)
                  (and (not= best worse)
                    (not= best better)
                    (> best-val better-val)
                    (or (and (> best worse) (< best better))
                      (and (< best worse) (> best better)))))))))
  :ret (s/or :sol ::m/int+
         :anomaly ::anomalies/anomaly))

(defn- recursive-unbounded
  "Only `smaller` tuple can be nil.  `bigger` must be >= 2 and >= `smaller` + 2,
  and `bigger-val` must be greater than `smaller-val`."
  [adj-f [smaller smaller-val] [bigger bigger-val] cap]
  (loop [[smaller smaller-val] [smaller smaller-val]
         [bigger bigger-val] [bigger bigger-val]]
    (let [new (int (min (* 2.0 bigger) cap))
          new-val (adj-f new)]
      (cond (and (= new cap)
              (> new-val bigger-val)) (recursive-bounded
                                        adj-f
                                        [bigger bigger-val]
                                        [new new-val]
                                        nil)
            (> new-val bigger-val) (recur [bigger bigger-val] [new new-val])
            (or (nil? smaller)
              (= new-val bigger-val)) (recursive-bounded
                                        adj-f
                                        [new new-val]
                                        [bigger bigger-val]
                                        nil)
            (> new-val smaller-val) (recursive-bounded
                                      adj-f
                                      [smaller smaller-val]
                                      [new new-val]
                                      [bigger bigger-val])
            :else (recursive-bounded
                    adj-f
                    [new new-val]
                    [smaller smaller-val]
                    [bigger bigger-val])))))

(s/fdef recursive-unbounded
  :args (s/and (s/cat :adj-f ::optimize-f
                 :smaller-tuple (s/nilable ::point-and-val)
                 :bigger-tuple ::point-and-val
                 :cap ::m/int+)
          (fn [{:keys [smaller-tuple bigger-tuple cap]}]
            (let [[smaller smaller-val] smaller-tuple
                  [bigger bigger-val] bigger-tuple]
              (and (> cap bigger)
                (or (not smaller)
                  (and (> bigger-val smaller-val)
                    (>= bigger (+ smaller 2))))))))
  :ret (s/or :sol ::m/int+
         :anomaly ::anomalies/anomaly))

(defn integer-optimize
  "Custom integer maximizer that exponentially focuses search around `guess`.
  For the given range of integers, `f` must either:
  (1) have no change of sign of the derivative, or
  (2) have one change of sign of the derivative, which must be at the maximum."
  [f guess [lower upper]]
  (let [zero-val (f guess)
        one-val (if (>= upper (inc guess))
                  (f (inc guess))
                  m/min-dbl)
        best-val (cond (= zero-val one-val) guess

                       (> zero-val one-val)
                       (let [m1-val (if (<= lower (dec guess))
                                      (f (dec guess))
                                      m/min-dbl)]
                         (if (>= zero-val m1-val)
                           guess
                           (let [m2-val (if (<= lower (- guess 2))
                                          (f (- guess 2))
                                          m/min-dbl)]
                             (if (>= m1-val m2-val)
                               (dec guess)
                               (if (= lower (- guess 2))
                                 (- guess 2)
                                 (let [v (recursive-unbounded
                                           (fn [g]
                                             (f (- guess g)))
                                           nil
                                           [2 m2-val]
                                           (int (- guess lower)))]
                                   (if (anomalies/anomaly? v)
                                     v
                                     (- guess v))))))))

                       :else
                       (let [two-val (if (>= upper (+ guess 2))
                                       (f (+ guess 2))
                                       m/min-dbl)]
                         (if (>= one-val two-val)
                           (inc guess)
                           (if (= upper (+ guess 2))
                             (+ guess 2)
                             (let [v (recursive-unbounded
                                       (fn [g]
                                         (f (+ g guess)))
                                       nil
                                       [2 two-val]
                                       (int (- upper guess)))]
                               (if (anomalies/anomaly? v)
                                 v
                                 (+ guess v)))))))]
    best-val))

(s/fdef integer-optimize
  :args (s/with-gen
          (s/and (s/cat :f ::optimize-f
                   :guess ::m/long
                   :interval ::intervals/long-interval)
            (fn [{:keys [guess interval]}]
              (and (intervals/in-interval? interval guess)
                (m/int? (- guess (first interval)))
                (m/int? (- (second interval) guess)))))
          #(gen/bind (gen/tuple
                       (s/gen ::m/long)
                       (s/gen ::intervals/int-interval))
             (fn [[g ii]]
               (gen/tuple (s/gen ::optimize-f)
                 (gen/return g)
                 (gen/tuple
                   (gen/return
                     (long (max m/min-long
                             (+ (double g) (first ii)))))
                   (gen/return
                     (long (min m/max-long
                             (+ (double g) (second ii))))))))))
  :ret (s/or :sol ::m/long
         :anomaly ::anomalies/anomaly))
