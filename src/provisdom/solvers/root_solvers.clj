(ns provisdom.solvers.root-solvers
  (:require
    [clojure.spec.alpha :as s]
    [clojure.spec.gen.alpha :as gen]
    [clojure.spec.test.alpha :as st]
    [orchestra.spec.test :as ost]
    [provisdom.utility-belt.async :as async]
    [provisdom.utility-belt.anomalies :as anomalies]
    [provisdom.math.core :as m]
    [provisdom.math.derivatives :as derivatives]
    [provisdom.math.intervals :as intervals]
    [provisdom.solvers.internal-apache-solvers :as apache-solvers]))

(s/def ::parallel? boolean?)

(s/def ::one-root-solver-type
  (s/or :apache-solvers ::apache-solvers/root-solver-type
    :not-apache-solvers #{:brent-dekker :modified-newton-raphson}))

(s/def ::root-solver-type
  (s/or :one ::one-root-solver-type
    :all #{:all}
    :seq (s/coll-of ::one-root-solver-type :distinct true)))

(s/def ::max-iter ::apache-solvers/max-iter)
(s/def ::rel-accu ::apache-solvers/rel-accu)
(s/def ::abs-accu ::apache-solvers/abs-accu)
(s/def ::guess ::apache-solvers/initial-guess)
(s/def ::univariate-f ::apache-solvers/univariate-f)
(s/def ::interval ::apache-solvers/finite-interval)

(s/def ::uni-f-with-guess-and-interval
  (s/with-gen
    (s/and (s/keys :req [::univariate-f ::guess ::interval])
      (fn [{::keys [guess interval]}]
        (intervals/in-interval? interval guess)))
    #(gen/one-of (map gen/return
                   (list {::univariate-f identity
                          ::guess        3.0
                          ::interval     [-5.0 5.0]}
                     {::univariate-f (fn [v]
                                       (- (m/cube v) (* 3 v)))
                      ::guess        3.0
                      ::interval     [-50.0 50.0]}
                     {::univariate-f (fn [v]
                                       (- (m/exp v) (* 5 v)))
                      ::guess        3.0
                      ::interval     [-50.0 50.0]})))))

(s/def ::univariate-f-with-interval
  (s/with-gen
    (s/keys :req [::univariate-f ::interval])
    #(gen/one-of (map gen/return
                   (list {::univariate-f identity
                          ::interval     [-5.0 5.0]}
                     {::univariate-f (fn [v]
                                       (+ (m/cube v) (* 3 v)))
                      ::interval     [-5.0 5.0]}
                     {::univariate-f (fn [v]
                                       (+ (m/exp v) (m/cube v) (* 3 v)))
                      ::interval     [-5.0 5.0]}
                     {::univariate-f (comp inc m/sq)
                      ::interval     [-5.0 5.0]})))))

(s/def ::mnr-derivative-fn
  (s/with-gen
    (s/nilable ::univariate-f)
    #(gen/return nil)))

(s/def ::univariate-f-with-mnr-derivative-fn-and-guess
  (s/with-gen
    (s/keys :req [::univariate-f ::mnr-derivative-fn ::guess])
    #(gen/one-of
       (map
         gen/return
         (list {::univariate-f      identity
                ::mnr-derivative-fn (constantly 1.0)
                ::guess             4.0}

           {::univariate-f      (fn [v] (+ (m/cube v) (* 3 v)))
            ::mnr-derivative-fn (fn [v] (* 3 (inc (m/sq v))))
            ::guess             -1.0}

           {::univariate-f      (fn [v] (+ (m/exp v) (m/cube v) (* 3 v)))
            ::mnr-derivative-fn (fn [v] (+ (m/exp v) (* 3 (inc (m/sq v)))))
            ::guess             3.0}

           {::univariate-f      (comp inc m/sq)
            ::mnr-derivative-fn (partial * 2)
            ::guess             5.0})))))

(defn root-solver-quadratic
  "Returns a tuple with the two solutions [- +] to 'x' from the quadratic
  equation, `a`*x^2 + `b`*x + `c`. Arguments: `a`, `b`, `c`: coefficients of a
  qaudratic equation."
  [a b c]
  (let [a (double a)
        b (double b)
        d (- (m/sq b) (* 4 a c))]
    (if (and (not (zero? a)) (m/finite-non-? d))
      (let [t1 (- b)
            t2 (m/sqrt d)
            t3 (* 2 a)]
        [(/ (- t1 t2) t3)
         (/ (+ t1 t2) t3)])
      {::anomalies/message  "No solution."
       ::anomalies/fn       (var root-solver-quadratic)
       ::anomalies/category ::anomalies/no-solve})))

(s/fdef root-solver-quadratic
  :args (s/cat :a ::m/num
          :b ::m/num
          :c ::m/num)
  :ret (s/or :solution (s/tuple ::m/num ::m/num)
         :anomaly ::anomalies/anomaly))

(defn- modified-newton-raphson
  "Modification dynamically alters a change factor based on previous changes.
  Can use a numerical derivative. Depending on function, guess may have to be
  close or result will become NaN."
  [{::keys [univariate-f mnr-derivative-fn guess]} abs-accu max-iter]
  (let [df (or mnr-derivative-fn (derivatives/derivative-fn univariate-f))
        dg (double guess)
        max-iter (or max-iter 1000)]
    (loop [i 0
           x1 dg
           x2 dg
           x3 dg]
      (let [fx1 (univariate-f x1)
            i (inc i)]
        (if (m/roughly? fx1 0 abs-accu)
          x1
          (if (> i max-iter)
            {::anomalies/message  (format "Too many iterations %d at %f" i x1)
             ::anomalies/fn       (var modified-newton-raphson)
             ::anomalies/category ::anomalies/no-solve}
            (let [d (df x1)]
              (if (m/roughly? d 0 m/dbl-close)
                {::anomalies/message  (format "Derivative too close to zero at %f" x1)
                 ::anomalies/fn       (var modified-newton-raphson)
                 ::anomalies/category ::anomalies/no-solve}
                (let [den (- (* 2 x2) x1 x3)
                      c (if (or (m/roughly? x2 x3 m/dbl-close)
                              (m/roughly? den 0 m/dbl-close))
                          1
                          (m/div (- x2 x3) den))
                      x (- x1 (/ (* c fx1) d))]
                  (recur i x x1 x2))))))))))

(s/fdef modified-newton-raphson
  :args (s/cat :univariate-f-with-mnr-derivative-fn-and-guess
          ::univariate-f-with-mnr-derivative-fn-and-guess
          :abs-accu ::abs-accu
          :max-iter ::max-iter)
  :ret (s/or :finite ::m/finite
         :anomaly ::anomalies/anomaly))

(defn- brent-dekker
  "References: http://en.wikipedia.org/wiki/Brent's_method"
  [{::keys [univariate-f interval]} abs-accu max-iter]
  (let [[a b] (mapv double interval)
        fa (double (univariate-f a))
        fb (double (univariate-f b))
        d m/max-dbl
        e abs-accu
        max-iter (or max-iter 1000)
        msg "root is not bracketed at [%f %f] with output of [%f %f]"]
    (if (m/non-? (* fa fb))
      {::anomalies/message  (format msg a b fa fb)
       ::anomalies/fn       (var brent-dekker)
       ::anomalies/category ::anomalies/no-solve}
      (let [z (< (m/abs fa) (m/abs fb))
            s b
            fs fb
            b (if z a b)
            fb (if z fa fb)
            a (if z s a)
            fa (if z fs fa)]
        (loop [a a
               fa fa
               b b
               fb fb
               c a
               fc fa
               d d
               m true
               i 0]
          (if (or (m/roughly? fb 0 m/dbl-close)
                (<= (m/abs (- a b)) e))
            b
            (let [s (if (and (not= fa fc) (not= fb fc))
                      (+ (m/div (* a fb fc) (* (- fa fb) (- fa fc)))
                        (m/div (* b fa fc) (* (- fb fa) (- fb fc)))
                        ;;Inverse quadratic interpolation
                        (/ (* c fa fb)
                          (- fc fa)
                          (- fc fb)))
                      (- b (m/div (* fb (- b a)) (- fb fa)))) ;Secant Rule
                  t (* (+ (* 3 a) b) 0.25)
                  m (or (not (or (and (> s t) (< s b))
                               (and (< s t) (> s b))))
                      (and m
                        (>= (m/abs (- s b)) (* 0.5 (m/abs (- b c)))))
                      (and (not m)
                        (>= (m/abs (- s b))
                          (* 0.5 (m/abs (- c d)))))
                      (and m
                        (< (m/abs (- b c)) e))
                      (and (not m)
                        (< (m/abs (- c d)) e)))
                  s (if m
                      (* 0.5 (+ a b))
                      s)
                  fs (univariate-f s)
                  d c
                  c b
                  fc fb
                  z (neg? (* fa fs))
                  z2 (if z
                       (< (m/abs fa) (m/abs fs))
                       (< (m/abs fs) (m/abs fb)))
                  s2 b
                  fs2 fb
                  b (if z
                      (if z2 a s)
                      (if z2 s b))
                  fb (if z
                       (if z2 fa fs)
                       (if z2 fs fb))
                  a (if z
                      (if z2 s a)
                      (if z2 s2 s))
                  fa (if z
                       (if z2 fs fa)
                       (if z2 fs2 fs))
                  i (inc i)
                  d (double d)
                  msg2 "too many iterations %d at %f with output of %f"]
              (if (> i max-iter)
                {::anomalies/message  (format msg2 i b fb)
                 ::anomalies/fn       (var brent-dekker)
                 ::anomalies/category ::anomalies/no-solve}
                (recur a fa b fb c fc d m i)))))))))

(s/fdef brent-dekker
  :args (s/cat :univariate-f-with-interval ::univariate-f-with-interval
          :abs-accu ::abs-accu
          :max-iter ::max-iter)
  :ret (s/or :finite ::m/finite
         :anomaly ::anomalies/anomaly))

(defn- root-selector-fn
  [univariate-f interval]
  (fn [results]
    (second (reduce
              (fn [[smallest-value smallest-point] new-point]
                (let [new-value (if (number? new-point)
                                  (m/abs (univariate-f new-point))
                                  m/max-dbl)]
                  (if (and (< new-value smallest-value)
                        (intervals/in-interval? interval new-point))
                    [new-value new-point]
                    [smallest-value smallest-point])))
              [m/max-dbl nil]
              results))))

(defn root-solver
  "The default is to run `:all` of the solvers, or choose one of the following,
  or a collection containing one or more of the following:
   `:bisection`, `:bracketing-nth-order-brent`, `:brent`, `:brent-dekker`,
   `:illinois`, `:modified-newton-raphson`, `:muller`, `:muller2`,
   `:newton-raphson`, `:pegasus`, `:regula-falsi`, `:ridders`, `:secant`.

  `::mnr-derivative-fn` - the function derivative for the modified newton 
  raphson method.

  `:modified-newton-raphson` does not use the interval to solve and will use a
  numerical derivative if `::mnr-derivative-fn` is nil. `:brent-dekker` and
  `:newton-raphson` do not use the `::guess`. `::univariate-f` and
  `::mnr-derivative-fn` should return NaN instead of throwing errors or
  returning nil.

  Except for `modified-newton-raphson`, solvers require that the intervals are
  set such that their function value have opposite signs and the ideal solution
  is in between."
  ([args] (root-solver args {}))
  ([{::keys [univariate-f guess interval]}
    {::keys [max-iter rel-accu abs-accu root-solver-type mnr-derivative-fn
             parallel?]
     :or    {rel-accu         1e-6
             abs-accu         1e-14
             root-solver-type :all
             parallel?        false}}]
   (let [max-iter (or max-iter 1000)
         solvers (if-not (= :all root-solver-type)
                   root-solver-type
                   (list :bisection :bracketing-nth-order-brent :brent
                     :brent-dekker :illinois :modified-newton-raphson :muller
                     :muller2 :newton-raphson :pegasus :regula-falsi :ridders
                     :secant))
         solver-fn (fn [solver-type]
                     (condp = solver-type
                       :brent-dekker #(brent-dekker {::univariate-f univariate-f
                                                     ::interval     interval}
                                        abs-accu
                                        max-iter)
                       :modified-newton-raphson #(modified-newton-raphson
                                                   {::univariate-f      univariate-f
                                                    ::mnr-derivative-fn mnr-derivative-fn
                                                    ::guess             guess}
                                                   abs-accu
                                                   max-iter)
                       #(apache-solvers/root-solver
                          {::apache-solvers/univariate-f    univariate-f
                           ::apache-solvers/initial-guess   guess
                           ::apache-solvers/finite-interval interval}
                          {::apache-solvers/max-iter         max-iter
                           ::apache-solvers/root-solver-type solver-type
                           ::apache-solvers/rel-accu         rel-accu
                           ::apache-solvers/abs-accu         abs-accu})))]
     (cond (keyword? solvers) ((solver-fn solvers))
           (m/one? (count solvers)) ((solver-fn (first solvers)))
           :else (let [sol (async/thread-select
                             (root-selector-fn univariate-f interval)
                             (map solver-fn solvers)
                             parallel?)]
                   (or sol {::anomalies/message  "No solution"
                            ::anomalies/fn       (var root-solver)
                            ::anomalies/category ::anomalies/no-solve}))))))

(s/fdef root-solver
  :args (s/cat :uni-f-with-guess-and-interval ::uni-f-with-guess-and-interval
          :opts (s/? (s/keys :opt [::max-iter
                                   ::rel-accu
                                   ::abs-accu
                                   ::root-solver-type
                                   ::mnr-derivative-fn])))
  :ret (s/or :finite ::m/finite
         :anomaly ::anomalies/anomaly))