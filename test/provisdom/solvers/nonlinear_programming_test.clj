(ns provisdom.solvers.nonlinear-programming-test
  (:require [clojure.test :refer :all]
            [provisdom.test.core :refer :all]
            [provisdom.math.core :as m]
            [provisdom.math.random :as random]
            [provisdom.solvers.nonlinear-programming :as nlp]
            [clojure.spec.test.alpha :as st]
            [orchestra.spec.test :as ost]))

;26 seconds

(set! *warn-on-reflection* true)

(ost/instrument)

(defn- objective
  [da]
  (if (= (count da) 2)
    (let [[a b] da]
      (- (m/exp a)
         (m/sq a)
         a
         (- (m/sq b))))
    m/nan))

(defn- gradient
  [da]
  (if (= (count da) 2)
    (let [[a b] da]
      [(dec (- (m/exp a) (* 2 a)))
       (* 2 b)])
    [m/nan m/nan]))

(defn geq-fn1
  [da]
  (if (= (count da) 2)
    (let [[a b] da]
      [(- a b 0.5)
       (- (- a b 0.5))])
    [m/inf- m/inf-]))

(defn geq-fn2
  [da]
  (if (= (count da) 2)
    (let [[a b] da]
      [(- a b 0.5)
       (- (- a b 0.5))
       (- (m/sq a) b 0.4)
       (- a 0.1)
       (- 1.1 a)
       (- b 0.1)
       (- 1.1 b)])
    [m/inf- m/inf- m/inf- m/inf-
     m/inf- m/inf- m/inf-]))

(defn geq-fn3
  [da]
  (if (= (count da) 2)
    (let [[a b] da]
      [(- a b 0.5)
       (- (- a b 0.5))
       (- (m/sq a) b 0.4)
       (- (m/pow a 3.0) (m/pow b 3.0) 0.9)
       (- a 0.1)
       (- 1.1 a)
       (- b 0.1)
       (- 1.1 b)])
    [m/inf- m/inf- m/inf- m/inf-
     m/inf- m/inf- m/inf- m/inf-]))

(defn geq-fn4
  [da]
  (if (= (count da) 2)
    (let [[a b] da]
      [(- (m/sq a) b 0.4)])
    [m/inf- m/inf- m/inf-]))

(def vars-guess [0.6 0.6])
(def var-intervals [[0.1 1.1] [0.1 1.1]])
(def var-intervals2 [[-0.1 2.1] [-0.1 2.1]])

(deftest constrained-nonlinear-programming-test
  (is (spec-check nlp/constrained-nonlinear-programming
                  {:coll-check-limit 10
                   :coll-error-limit 10
                   :fspec-iterations 10
                   :recursion-limit  1
                   :test-check       {:num-tests 60}}))
  (is= {::nlp/vector-point [0.6931471805091619 0.19314718050916202]
        ::nlp/value        0.8637056388801094}
       (nlp/constrained-nonlinear-programming
         {::nlp/objective  objective
          ::nlp/vars-guess vars-guess
          ::nlp/geq-fn     geq-fn1}))
  (is= {::nlp/vector-point [0.8872983338968569 0.38729833389685697]
        ::nlp/value        0.9039629552680116}
       (nlp/constrained-nonlinear-programming
         {::nlp/objective  objective
          ::nlp/vars-guess vars-guess
          ::nlp/geq-fn     geq-fn2}))
  (is= {::nlp/vector-point [1.011030003726171 0.5110300037261711]
        ::nlp/value        0.9763704437553828}
       (nlp/constrained-nonlinear-programming
         {::nlp/objective  objective
          ::nlp/vars-guess vars-guess
          ::nlp/geq-fn     geq-fn3}))
  (is= {::nlp/vector-point [1.2564312102129944 1.0759313612014115E-10]
        ::nlp/value        0.6778118266163524}
       (nlp/constrained-nonlinear-programming
         {::nlp/objective  objective
          ::nlp/vars-guess vars-guess
          ::nlp/geq-fn     geq-fn4})))

(deftest unbounded-nonlinear-programming-test
  (is (spec-check nlp/unbounded-nonlinear-programming
                  {:coll-check-limit 10
                   :coll-error-limit 10
                   :fspec-iterations 10
                   :recursion-limit  1
                   :test-check       {:num-tests 200}}))
  (is= {::nlp/value        0.6778118266175843
        ::nlp/vector-point [1.256431450577305 -1.0899505598883398E-6]}
       (nlp/unbounded-nonlinear-programming {::nlp/objective  objective
                                             ::nlp/vars-guess [0.1 0.1]}))
  (is= {::nlp/value        -0.6778118266175843
        ::nlp/vector-point [1.256431450577305 -1.0899505598883398E-6]}
       (nlp/unbounded-nonlinear-programming
         {::nlp/objective  (fn [x]
                             (- (objective x)))
          ::nlp/vars-guess [0.1 0.1]}
         {::nlp/goal :max}))
  (is= {::nlp/value        0.6778118266175384
        ::nlp/vector-point [1.2564313350002578 -1.0835075855017152E-6]}
       (nlp/unbounded-nonlinear-programming {::nlp/objective  objective
                                             ::nlp/vars-guess vars-guess}))
  (is= {::nlp/value        0.6778120040393626
        ::nlp/vector-point [1.25625 3.9062500000002637E-4]}
       (nlp/unbounded-nonlinear-programming
         {::nlp/objective  objective
          ::nlp/vars-guess [0.1 0.1]}
         {::nlp/unbounded-solver-type [:multi-directional-simplex
                                       :nelder-mead
                                       :powell]}))
  (is= {::nlp/value        0.6778118266175843
        ::nlp/vector-point [1.256431450577305 -1.0899505598883398E-6]}
       (nlp/unbounded-nonlinear-programming
         {::nlp/objective  objective
          ::nlp/vars-guess [0.1 0.1]}
         {::nlp/unbounded-solver-type :cobyla}))
  (is= {::nlp/value        0.6778120040393626
        ::nlp/vector-point [1.25625 3.9062500000002637E-4]}
       (nlp/unbounded-nonlinear-programming
         {::nlp/objective  objective
          ::nlp/vars-guess [0.1 0.1]}
         {::nlp/unbounded-solver-type :multi-directional-simplex}))
  (is= {::nlp/value        0.6778120220853294
        ::nlp/vector-point [1.2561596612082213 -3.737686114618499E-4]}
       (nlp/unbounded-nonlinear-programming
         {::nlp/objective  objective
          ::nlp/vars-guess [0.1 0.1]}
         {::nlp/unbounded-solver-type :nelder-mead}))
  (is= {::nlp/value        0.6778120099062883
        ::nlp/vector-point [1.256923363600252 1.3877787807814457E-17]}
       (nlp/unbounded-nonlinear-programming
         {::nlp/objective  objective
          ::nlp/vars-guess [0.1 0.1]}
         {::nlp/unbounded-solver-type :powell}))
  (is= {::nlp/value        0.6778118267055769
        ::nlp/vector-point [1.2564419965907032 -1.09108569361859E-6]}
       (nlp/unbounded-nonlinear-programming
         {::nlp/objective  objective
          ::nlp/vars-guess [0.1 0.1]}
         {::nlp/unbounded-solver-type :conjugate-gradient-using-fletcher-reeves
          ::nlp/gradient              gradient}))
  (is= {::nlp/value        0.677811826705496
        ::nlp/vector-point [1.2564419917200609 -1.0904278115382364E-6]}
       (nlp/unbounded-nonlinear-programming
         {::nlp/objective  objective
          ::nlp/vars-guess [0.1 0.1]}
         {::nlp/unbounded-solver-type :conjugate-gradient-using-polak-ribiere
          ::nlp/gradient              gradient})))

(deftest bounded-nonlinear-programming-without-evolutionary-test
  (is (spec-check nlp/bounded-nonlinear-programming-without-evolutionary
                  {:coll-check-limit 10
                   :coll-error-limit 10
                   :fspec-iterations 10
                   :recursion-limit  1
                   :test-check       {:num-tests 350}}))
  (is= {::nlp/value        0.7041660239464331
        ::nlp/vector-point [1.1 0.1]}
       (nlp/bounded-nonlinear-programming-without-evolutionary
         {::nlp/objective     objective
          ::nlp/vars-guess    vars-guess
          ::nlp/var-intervals var-intervals}))
  (is= {::nlp/value        0.7041660239464331
        ::nlp/vector-point [1.1 0.09999999999999998]}
       (nlp/bounded-nonlinear-programming-without-evolutionary
         {::nlp/objective     objective
          ::nlp/vars-guess    vars-guess
          ::nlp/var-intervals var-intervals}
         {::nlp/bounded-without-evolutionary-solver-type :cobyla}))
  (is= {::nlp/value        0.7041660239464331
        ::nlp/vector-point [1.1 0.1]}
       (nlp/bounded-nonlinear-programming-without-evolutionary
         {::nlp/objective     objective
          ::nlp/vars-guess    vars-guess
          ::nlp/var-intervals var-intervals}
         {::nlp/bounded-without-evolutionary-solver-type :bobyqa}))
  (is= {::nlp/value        0.677811826616352
        ::nlp/vector-point [1.2564312010957783 -6.845604985308992E-9]}
       (nlp/bounded-nonlinear-programming-without-evolutionary
         {::nlp/objective     objective
          ::nlp/vars-guess    vars-guess
          ::nlp/var-intervals var-intervals2}
         {::nlp/bounded-without-evolutionary-solver-type :bobyqa})))

(deftest bounded-nonlinear-programming-including-evolutionary!-test
  (is (spec-check nlp/bounded-nonlinear-programming-including-evolutionary!
                  {:coll-check-limit 10
                   :coll-error-limit 10
                   :fspec-iterations 10
                   :recursion-limit  1
                   :test-check       {:num-tests 150}}))
  (is= {::nlp/value        0.7041660239464331
        ::nlp/vector-point [1.1 0.1]}
       (random/bind-seed
         3
         (nlp/bounded-nonlinear-programming-including-evolutionary!
           {::nlp/objective     objective
            ::nlp/vars-guess    vars-guess
            ::nlp/var-intervals var-intervals})))
  (is= {::nlp/value        0.7041660239464331
        ::nlp/vector-point [1.1 0.1]}
       (random/bind-seed
         3
         (nlp/bounded-nonlinear-programming-including-evolutionary!
           {::nlp/objective     objective
            ::nlp/vars-guess    vars-guess
            ::nlp/var-intervals var-intervals}
           {::nlp/bounded-including-evolutionary-solver-type :bobyqa})))
  (is= {::nlp/value        0.7041660239464331
        ::nlp/vector-point [1.1 0.1]}
       (random/bind-seed
         3
         (nlp/bounded-nonlinear-programming-including-evolutionary!
           {::nlp/objective     objective
            ::nlp/vars-guess    vars-guess
            ::nlp/var-intervals var-intervals}
           {::nlp/bounded-including-evolutionary-solver-type :cobyla})))
  (is= {::nlp/value        0.7041660239464331
        ::nlp/vector-point [1.1 0.1]}
       (random/bind-seed
         3
         (nlp/bounded-nonlinear-programming-including-evolutionary!
           {::nlp/objective     objective
            ::nlp/vars-guess    vars-guess
            ::nlp/var-intervals var-intervals}
           {::nlp/bounded-including-evolutionary-solver-type :cma-es-only})))
  (is= {::nlp/value        0.677811826616352
        ::nlp/vector-point [1.2564312010957783 -6.845604985308992E-9]}
       (random/bind-seed
         3
         (nlp/bounded-nonlinear-programming-including-evolutionary!
           {::nlp/objective     objective
            ::nlp/vars-guess    vars-guess
            ::nlp/var-intervals var-intervals2})))
  (is= {::nlp/value        0.6778118266163528
        ::nlp/vector-point [1.256431192726513 -2.4399976876458034E-8]}
       (random/bind-seed
         3
         (nlp/bounded-nonlinear-programming-including-evolutionary!
           {::nlp/objective     objective
            ::nlp/vars-guess    vars-guess
            ::nlp/var-intervals var-intervals2}
           {::nlp/bounded-including-evolutionary-solver-type :cma-es-only})))
  (is= {::nlp/value        -0.7041660239464331
        ::nlp/vector-point [1.1 0.1]}
       (random/bind-seed
         3
         (nlp/bounded-nonlinear-programming-including-evolutionary!
           {::nlp/objective     (fn [x]
                                  (- (objective x)))
            ::nlp/vars-guess    vars-guess
            ::nlp/var-intervals var-intervals}
           {::nlp/goal                                       :max
            ::nlp/bounded-including-evolutionary-solver-type :cma-es-only}))))

#_(ost/unstrument)