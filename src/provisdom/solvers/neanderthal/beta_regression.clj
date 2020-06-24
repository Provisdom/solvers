(ns provisdom.solvers.neanderthal.beta-regression
  (:require
    [clojure.spec.alpha :as s]
    [clojure.spec.gen.alpha :as gen]
    [clojure.spec.test.alpha :as st]
    [clojure.pprint :refer [pprint]]
    [orchestra.spec.test :as ost]
    [provisdom.utility-belt.anomalies :as anomalies]
    [provisdom.math.core :as m]
    [provisdom.math.special-functions :as special-fns]
    [provisdom.math.vector :as vector]
    [provisdom.math.matrix :as mx]
    [provisdom.math.tensor :as tensor]
    [provisdom.math.derivatives :as der]
    [provisdom.apache-math.apache-matrix :as apache-mx]
    [provisdom.solvers.nonlinear-constraints-without-objective :as nlls]
    [provisdom.solvers.nonlinear-programming :as nlp]
    [uncomplicate.neanderthal.native :as nat]
    [uncomplicate.neanderthal.core :as n]
    [uncomplicate.neanderthal.vect-math :as vec]))

(defn descent [f g xs error rate epochs]
  (let [ns (range (count xs))
        r rate #_(/ rate error 2)]
    (letfn [(delta [ys i e] (into [] (map (fn [j y] (if (= i j) (+ y e) y))
                                          ns ys)))
            (step [ys]
              (mapv #(+ %1 (* r %2)) ys (g ys)))]
      (loop [ys xs
             m 0]
        (let [zs (step ys)]
          (println ys (g ys) zs)
          #_(if (= 0 (mod m 20)) (println zs))
          (if (or (> m epochs)
                  #_(< (norm (map - zs ys)) error))
            ys
            (recur zs (inc m))))))))

(defn digamma
  [work x]
  (n/copy! x work)
  (n/alter! work (fn ^double [^long _ ^double x] (special-fns/digamma x)))
  work)

(defn digamma!
  [x]
  (n/alter! x (fn ^double [^long _ ^double x] (special-fns/digamma x)))
  x)

(defn logit
  [neg-ones work x]
  (let [one-x (n/axpy x neg-ones)
        quot (vec/div! x one-x work)]
    (vec/log! quot work)))

(defn logit!
  [ones x]
  (let [one-x (n/axpy! -1.0 x ones)
        quot (vec/div! x one-x)]
    (vec/log! quot)))

(defn logistic
  [ones neg-ones work x]
  (do
    (vec/mul! x neg-ones work)
    (vec/exp! work)
    (n/axpy! ones work)
    (vec/inv! work)))

(defn logistic-derivative
  [ones x]
  (let [f (logistic ones x)
        f2 (vec/pow f 2)]
    (n/axpby! 1.0 f -1.0 f2)))

(defn logistic!
  [ones x]
  (let [ex (vec/exp! (n/scal! -1.0 x))]
    (vec/inv! (n/axpy! 1.0 ones ex))))

(defn logistic-derivative!
  [ones x]
  (let [f (logistic! ones x)
        f2 (vec/pow f 2)]
    (n/axpby! 1.0 f -1.0 f2)))

(defn x->mu
  [beta x]
  (let [xbeta (vector/dot-product x beta)]
    (special-fns/logistic xbeta)))

(defn- y*
  [y]
  (m/log (/ y (m/one- y))))

(defn mu*
  [phi mu]
  (let [mu*phi (* mu phi)]
    (- (special-fns/digamma mu*phi)
       (special-fns/digamma (* (m/one- mu) phi)))))

(defn objective
  [x-mx y [phi & beta :as b]]
  (let [beta (vec beta)
        eta (mx/get-column (mx/mx* x-mx (mx/column-matrix beta)) 0)
        mu (mapv special-fns/logistic eta)
        l_t (mapv (fn [mu y]
                    (+ (special-fns/log-gamma phi)
                       (- (special-fns/log-gamma (* mu phi)))
                       (- (special-fns/log-gamma (* (m/one- mu) phi)))
                       (* (- (* mu phi) 1.0) (m/log y))
                       (* (- (* (m/one- mu) phi) 1.0) (m/log (m/one- y)))))
                 mu y)
        f (apply + l_t)]
    #_(println l_t)
    #_(println f)
    (if (m/nan? f) m/inf+ f)))

(defn gradient-old
  [x-mx y [phi & beta :as b]]
  (let [beta (vec beta)
        eta (mx/get-column (mx/mx* x-mx (mx/column-matrix beta)) 0)
        mu (mapv special-fns/logistic eta)
        mu* (mapv (partial mu* phi) mu)
        y* (mapv y* y)
        T (mapv #(* (special-fns/logistic %) (m/one- (special-fns/logistic %))) eta)
        Ty*-mu* (mapv (fn [t mu* y*]
                        (* phi t (- y* mu*)))
                      T mu* y*)
        u_beta (mx/get-column (mx/mx* (mx/transpose x-mx) (mx/column-matrix Ty*-mu*)) 0)
        u_phi (apply + (map (fn [mu y]
                              (+ (special-fns/digamma phi)
                                 (- (* mu (special-fns/digamma (* mu phi))))
                                 (- (* (m/one- mu) (special-fns/digamma (* (m/one- mu) phi))))
                                 (* mu (m/log y))
                                 (* (m/one- mu) (m/log (m/one- y)))))
                            mu y))
        s (vec (cons u_phi u_beta))]
    #_(println s)
    #_(throw (ex-info "" {}))
    s))

(defn gradient
  [ones -ones X y* log-y log-one-y work [phi & beta :as b]]
  (let [beta (nat/dv beta)
        eta (n/mv! X beta (work 0))
        mu (logistic ones -ones (work 1) eta)
        phis (n/scal phi ones)
        mu*phi (vec/mul! phis mu (work 4))
        one-mu (n/axpy -1.0 mu ones)
        one-mu*phi (vec/mul! phis one-mu (work 3))
        digamma-phi (digamma (work 5) phis)
        digamma-mu*phi (digamma (work 6) mu*phi)
        digamma-one-mu*phi (digamma (work 7) one-mu*phi)
        mu* (n/axpy -1.0 digamma-one-mu*phi digamma-mu*phi)
        T (nat/dgd (.n ones) (logistic-derivative! ones eta))
        y*-mu* (n/axpy! y* (vec/mul! mu* -ones))
        Ty*-mu* (n/mv! T y*-mu* (work 2))
        u_beta (n/mv phi (n/trans X) Ty*-mu*)
        u_phi (n/sum
                (n/axpy 1.0 digamma-phi
                        -1.0 (vec/mul! digamma-mu*phi mu)
                        -1.0 (vec/mul! digamma-one-mu*phi one-mu)
                        1.0 (vec/mul! mu log-y)
                        1.0 (vec/mul! one-mu log-one-y)))
        #_(apply + (map (fn [mu y]
                          (+ (special-fns/digamma phi)
                             (- (* mu (special-fns/digamma (* mu phi))))
                             (- (* (m/one- mu) (special-fns/digamma (* (m/one- mu) phi))))
                             (* mu (m/log y))
                             (* (m/one- mu) (m/log (m/one- y)))))
                        mu y))
        s (vec (cons u_phi u_beta))]
    #_(println s)
    #_(throw (ex-info "" {}))
    #_(clojure.pprint/pprint s)
    s))

(defn jacobian
  [x-mx y [phi & beta :as b]]
  (let []))

(defn solve
  ([x-mx y beta0] (solve x-mx y beta0 {}))
  ([x-mx y beta0 {::keys [old? max-iter max-evaluations abs-accu] :or {abs-accu 1e-8 max-iter 1000 max-evaluations 1000}}]
   (let [n (count y)
         X (nat/dge x-mx)
         yv (nat/dv y)
         ones (nat/dv (repeat n 1.0))
         -ones (nat/dv (repeat n -1.0))
         log-y (vec/log yv)
         log-one-y (vec/log (n/axpy -1.0 yv ones))
         y* (n/axpy -1.0 log-one-y log-y)
         work1 (vec (repeatedly 8 #(nat/dv (repeat n 0.0))))
         work2 (vec (repeatedly 8 #(nat/dv (repeat n 0.0))))
         grad1 (partial gradient ones -ones X y* log-y log-one-y work1)
         grad1-old (partial gradient-old x-mx y)
         grad2 (partial gradient ones -ones X y* log-y log-one-y work2)
         jacof (der/jacobian-fn (if old? grad1-old grad1))]
     (descent (comp - (partial objective x-mx y)) #(mapv - (gradient-old x-mx y %)) beta0 1e-3 0.00001 1000)
     #_(nlp/unbounded-nonlinear-programming {::nlp/objective (partial objective x-mx y)
                                             ::nlp/vars-guess beta0}
                                            {::nlp/gradient (partial gradient x-mx y)
                                             ::nlp/goal :max
                                             ::nlp/unbounded-solver-type :nelder-mead
                                             #_#_::nlp/initial-step-size-for-conjugate-gradient 100000000.0
                                             ::nlp/abs-accu abs-accu
                                             ::nlp/max-iter max-iter})
     #_(nlls/nonlinear-least-squares {::nlls/constraints-fn         (if old? grad1-old grad2)
                                      ::nlls/constraint-jacobian-fn (fn [x] (vec (jacof (vec x))))
                                      ::nlls/vars-guess             beta0}
                                     {::nlls/nls-solver-type :newton-raphson
                                      ::nlls/max-evaluations max-evaluations
                                      ::nlls/max-iter        max-iter
                                      ::nlls/abs-accu        abs-accu}))))
