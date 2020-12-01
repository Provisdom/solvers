(ns provisdom.solvers.root-solvers-test
  (:require [clojure.test :refer :all]
            [provisdom.test.core :refer :all]
            [provisdom.math.core :as m]
            [provisdom.solvers.root-solvers :as root-solvers]
            [provisdom.utility-belt.anomalies :as anomalies]
            [provisdom.solvers.internal-apache-solvers :as apache-solvers]
            [orchestra.spec.test :as ost]))

;54 seconds

(set! *warn-on-reflection* true)

(ost/instrument)

(deftest root-solver-quadratic-test
  (is (spec-check root-solvers/root-solver-quadratic))
  (is= [-0.3333333333333333 1.0]
    (root-solvers/root-solver-quadratic 3.0 -2.0 -1.0))
  (is= [-0.8507810593582121 2.350781059358212]
    (root-solvers/root-solver-quadratic 2.0 -3.0 -4.0))
  (is= [-0.25 0.0] (root-solvers/root-solver-quadratic 4.0 1.0 0.0))
  (is= {::anomalies/message  "No solution."
        ::anomalies/fn       (var root-solvers/root-solver-quadratic)
        ::anomalies/category ::anomalies/no-solve}
    (root-solvers/root-solver-quadratic 3.0 2.0 1.0))
  (is= {::anomalies/message  "No solution."
        ::anomalies/fn       (var root-solvers/root-solver-quadratic)
        ::anomalies/category ::anomalies/no-solve}
    (root-solvers/root-solver-quadratic 4.0 2.0 1.0)))

(defn f
  [x]
  (+ (m/pow x 2.3) (- x) 0.1))

(defn df
  [x]
  (dec (* 2.3 (m/pow x 1.3))))

(def guess 1.0)

(deftest root-solver-test
  (is (spec-check root-solvers/root-solver
        {:coll-check-limit 10
         :coll-error-limit 10
         :fspec-iterations 10
         :recursion-limit  1
         :test-check       {:num-tests 300}}))
  ;;by default, use all of the solvers
  (is= 0.9148015417397132
    (root-solvers/root-solver {::root-solvers/univariate-f f
                               ::root-solvers/guess        guess
                               ::root-solvers/interval     [-20.0 20.0]}))
  (is= 0.9148015417397137
    (root-solvers/root-solver {::root-solvers/univariate-f f
                               ::root-solvers/guess        guess
                               ::root-solvers/interval     [0.1 1.8]}))
  (is= 0.9148015417397134
    (root-solvers/root-solver {::root-solvers/univariate-f f
                               ::root-solvers/guess        guess
                               ::root-solvers/interval     [0.5 1.4]}))
  ;;can use a mix
  (is= 0.9148015417397132
    (root-solvers/root-solver
      {::root-solvers/univariate-f f
       ::root-solvers/guess        guess
       ::root-solvers/interval     [-20.0 20.0]}
      {::root-solvers/root-solver-type
       [:illinois :bracketing-nth-order-brent
        :brent-dekker :muller :muller2 :bisection]}))
  ;modified Newton-Raphson (MNR) ignores 'interval' and can use analytical
  ; derivative or will use numerical derivative by default
  ; performs well and robust even w/o derivative
  (is= 0.9148015417397137
    (root-solvers/root-solver
      {::root-solvers/univariate-f f
       ::root-solvers/guess        guess
       ::root-solvers/interval     [-20.0 20.0]}
      {::root-solvers/root-solver-type :modified-newton-raphson}))
  (is= 0.9148015417397137
    (root-solvers/root-solver
      {::root-solvers/univariate-f f
       ::root-solvers/guess        guess
       ::root-solvers/interval     [0.1 1.8]}
      {::root-solvers/root-solver-type :modified-newton-raphson}))
  (is= 0.9148015417397137
    (root-solvers/root-solver
      {::root-solvers/univariate-f f
       ::root-solvers/guess        guess
       ::root-solvers/interval     [0.5 1.4]}
      {::root-solvers/root-solver-type :modified-newton-raphson}))
  (is= 0.9148015417397137                                   ;same as w/o df
    (root-solvers/root-solver
      {::root-solvers/univariate-f f
       ::root-solvers/guess        guess
       ::root-solvers/interval     [-20.0 20.0]}
      {::root-solvers/root-solver-type  :modified-newton-raphson
       ::root-solvers/mnr-derivative-fn df}))
  (is= 0.9148015417397184                                   ;changed df->df2 to make sure
    (root-solvers/root-solver
      {::root-solvers/univariate-f f
       ::root-solvers/guess        guess
       ::root-solvers/interval     [-20.0 20.0]}
      {::root-solvers/root-solver-type  :modified-newton-raphson
       ::root-solvers/mnr-derivative-fn (fn [x]
                                          (dec (* 2.3 (m/pow x 1.4))))}))
  ;;Newton-Raphson isn't super accurate
  (is= {::anomalies/message  "illegal state: maximal count (1,000) exceeded: evaluations"
        ::anomalies/fn       (var apache-solvers/root-solver)
        ::anomalies/category ::anomalies/third-party}
    (root-solvers/root-solver
      {::root-solvers/univariate-f f
       ::root-solvers/guess        guess
       ::root-solvers/interval     [-20.0 20.0]}
      {::root-solvers/root-solver-type :newton-raphson}))
  (is= 0.8925955018349163
    (root-solvers/root-solver
      {::root-solvers/univariate-f f
       ::root-solvers/guess        guess
       ::root-solvers/interval     [0.1 1.8]}
      {::root-solvers/root-solver-type :newton-raphson}))
  (is= 0.8925955018349163
    (root-solvers/root-solver
      {::root-solvers/univariate-f f
       ::root-solvers/guess        guess
       ::root-solvers/interval     [0.5 1.4]}
      {::root-solvers/root-solver-type :newton-raphson}))
  ;;illinois
  (is (anomalies/anomaly?
        (root-solvers/root-solver {::root-solvers/univariate-f f
                                   ::root-solvers/guess        guess
                                   ::root-solvers/interval     [-20.0 20.0]}
          {::root-solvers/root-solver-type :illinois})))
  (is= {::anomalies/message  (str "function values at endpoints do not have different signs, "
                               "endpoints: [0.1, 1.8], values: [0.005, 2.165]")
        ::anomalies/fn       (var apache-solvers/root-solver)
        ::anomalies/category ::anomalies/third-party}
    (root-solvers/root-solver {::root-solvers/univariate-f f
                               ::root-solvers/guess        guess
                               ::root-solvers/interval     [0.1 1.8]}
      {::root-solvers/root-solver-type :illinois}))
  (is= 0.914801638500105
    (root-solvers/root-solver {::root-solvers/univariate-f f
                               ::root-solvers/guess        guess
                               ::root-solvers/interval     [0.5 1.4]}
      {::root-solvers/root-solver-type :illinois}))
  ;;bracketing-nth-order-brent
  (is (anomalies/anomaly?
        (root-solvers/root-solver
          {::root-solvers/univariate-f f
           ::root-solvers/guess        guess
           ::root-solvers/interval     [-20.0 20.0]}
          {::root-solvers/root-solver-type :bracketing-nth-order-brent})))
  (is= {::anomalies/message  (str "function values at endpoints do not have different signs, "
                               "endpoints: [0.1, 1.8], values: [0.005, 2.165]")
        ::anomalies/fn       (var apache-solvers/root-solver)
        ::anomalies/category ::anomalies/third-party}
    (root-solvers/root-solver
      {::root-solvers/univariate-f f
       ::root-solvers/guess        guess
       ::root-solvers/interval     [0.1 1.8]}
      {::root-solvers/root-solver-type :bracketing-nth-order-brent}))
  (is= 0.9148015511654863
    (root-solvers/root-solver
      {::root-solvers/univariate-f f
       ::root-solvers/guess        guess
       ::root-solvers/interval     [0.5 1.4]}
      {::root-solvers/root-solver-type :bracketing-nth-order-brent}))
  ;;brent
  (is (anomalies/anomaly?
        (root-solvers/root-solver
          {::root-solvers/univariate-f f
           ::root-solvers/guess        guess
           ::root-solvers/interval     [-20.0 20.0]}
          {::root-solvers/root-solver-type :brent})))
  (is= {::anomalies/message  (str "function values at endpoints do not have different signs, "
                               "endpoints: [0.1, 1.8], values: [0.005, 2.165]")
        ::anomalies/fn       (var apache-solvers/root-solver)
        ::anomalies/category ::anomalies/third-party}
    (root-solvers/root-solver
      {::root-solvers/univariate-f f
       ::root-solvers/guess        guess
       ::root-solvers/interval     [0.1 1.8]}
      {::root-solvers/root-solver-type :brent}))
  (is= 0.9148015404078098
    (root-solvers/root-solver
      {::root-solvers/univariate-f f
       ::root-solvers/guess        guess
       ::root-solvers/interval     [0.5 1.4]}
      {::root-solvers/root-solver-type :brent}))
  ;;brent-dekker
  (is= 0.9148015417397132
    (root-solvers/root-solver
      {::root-solvers/univariate-f f
       ::root-solvers/guess        guess
       ::root-solvers/interval     [-20.0 20.0]}
      {::root-solvers/root-solver-type :brent-dekker}))
  (is= {::anomalies/message  (str "root is not bracketed at [0.100000 1.800000] with "
                               "output of [0.005012 2.164798]")
        ::anomalies/fn       (var root-solvers/brent-dekker)
        ::anomalies/category ::anomalies/no-solve}
    (root-solvers/root-solver
      {::root-solvers/univariate-f f
       ::root-solvers/guess        guess
       ::root-solvers/interval     [0.1 1.8]}
      {::root-solvers/root-solver-type :brent-dekker}))
  (is= 0.9148015417397135
    (root-solvers/root-solver
      {::root-solvers/univariate-f f
       ::root-solvers/guess        guess
       ::root-solvers/interval     [0.5 1.4]}
      {::root-solvers/root-solver-type :brent-dekker}))
  ;;muller
  (is (anomalies/anomaly?
        (root-solvers/root-solver {::root-solvers/univariate-f f
                                   ::root-solvers/guess        guess
                                   ::root-solvers/interval     [-20.0 20.0]}
          {::root-solvers/root-solver-type :muller})))
  (is= {::anomalies/message  (str "function values at endpoints do not have different signs, "
                               "endpoints: [0.1, 1.8], values: [0.005, 2.165]")
        ::anomalies/fn       (var apache-solvers/root-solver)
        ::anomalies/category ::anomalies/third-party}
    (root-solvers/root-solver {::root-solvers/univariate-f f
                               ::root-solvers/guess        guess
                               ::root-solvers/interval     [0.1 1.8]}
      {::root-solvers/root-solver-type :muller}))
  (is= 0.914801541740109
    (root-solvers/root-solver {::root-solvers/univariate-f f
                               ::root-solvers/guess        guess
                               ::root-solvers/interval     [0.5 1.4]}
      {::root-solvers/root-solver-type :muller}))
  ;;muller2
  (is= {::anomalies/message  "illegal state: maximal count (1,000) exceeded: evaluations"
        ::anomalies/fn       (var apache-solvers/root-solver)
        ::anomalies/category ::anomalies/third-party}
    (root-solvers/root-solver {::root-solvers/univariate-f f
                               ::root-solvers/guess        guess
                               ::root-solvers/interval     [-20.0 20.0]}
      {::root-solvers/root-solver-type :muller2}))
  (is= {::anomalies/message  (str "function values at endpoints do not have different signs, "
                               "endpoints: [0.1, 1.8], values: [0.005, 2.165]")
        ::anomalies/fn       (var apache-solvers/root-solver)
        ::anomalies/category ::anomalies/third-party}
    (root-solvers/root-solver {::root-solvers/univariate-f f
                               ::root-solvers/guess        guess
                               ::root-solvers/interval     [0.1 1.8]}
      {::root-solvers/root-solver-type :muller2}))
  (is= 0.9148015417397134
    (root-solvers/root-solver {::root-solvers/univariate-f f
                               ::root-solvers/guess        guess
                               ::root-solvers/interval     [0.5 1.4]}
      {::root-solvers/root-solver-type :muller2}))
  ;;bisection
  (is= -19.999999999999996                                  ;essentially doesn't work w/o opposite signs
    (root-solvers/root-solver {::root-solvers/univariate-f f
                               ::root-solvers/guess        guess
                               ::root-solvers/interval     [-20.0 20.0]}
      {::root-solvers/root-solver-type :bisection}))
  (is= 1.7999999999999972                                   ;essentially doesn't work w/o opposite signs
    (root-solvers/root-solver {::root-solvers/univariate-f f
                               ::root-solvers/guess        guess
                               ::root-solvers/interval     [0.1 1.8]}
      {::root-solvers/root-solver-type :bisection}))
  (is= 0.9148015417397114
    (root-solvers/root-solver {::root-solvers/univariate-f f
                               ::root-solvers/guess        guess
                               ::root-solvers/interval     [0.5 1.4]}
      {::root-solvers/root-solver-type :bisection}))
  ;;pegasus
  (is (anomalies/anomaly?
        (root-solvers/root-solver {::root-solvers/univariate-f f
                                   ::root-solvers/guess        guess
                                   ::root-solvers/interval     [-20.0 20.0]}
          {::root-solvers/root-solver-type :pegasus})))
  (is= {::anomalies/message  (str "function values at endpoints do not have different signs, "
                               "endpoints: [0.1, 1.8], values: [0.005, 2.165]")
        ::anomalies/fn       (var apache-solvers/root-solver)
        ::anomalies/category ::anomalies/third-party}
    (root-solvers/root-solver {::root-solvers/univariate-f f
                               ::root-solvers/guess        guess
                               ::root-solvers/interval     [0.1 1.8]}
      {::root-solvers/root-solver-type :pegasus}))
  (is= 0.9148015417397134
    (root-solvers/root-solver {::root-solvers/univariate-f f
                               ::root-solvers/guess        guess
                               ::root-solvers/interval     [0.5 1.4]}
      {::root-solvers/root-solver-type :pegasus}))
  ;;regula-falsi
  (is (anomalies/anomaly?
        (root-solvers/root-solver
          {::root-solvers/univariate-f f
           ::root-solvers/guess        guess
           ::root-solvers/interval     [-20.0 20.0]}
          {::root-solvers/root-solver-type :regula-falsi})))
  (is= {::anomalies/message  (str "function values at endpoints do not have different signs, "
                               "endpoints: [0.1, 1.8], values: [0.005, 2.165]")
        ::anomalies/fn       (var apache-solvers/root-solver)
        ::anomalies/category ::anomalies/third-party}
    (root-solvers/root-solver
      {::root-solvers/univariate-f f
       ::root-solvers/guess        guess
       ::root-solvers/interval     [0.1 1.8]}
      {::root-solvers/root-solver-type :regula-falsi}))
  (is= 0.914801541739713
    (root-solvers/root-solver
      {::root-solvers/univariate-f f
       ::root-solvers/guess        guess
       ::root-solvers/interval     [0.5 1.4]}
      {::root-solvers/root-solver-type :regula-falsi}))
  ;;ridders
  (is (anomalies/anomaly?
        (root-solvers/root-solver {::root-solvers/univariate-f f
                                   ::root-solvers/guess        guess
                                   ::root-solvers/interval     [-20.0 20.0]}
          {::root-solvers/root-solver-type :ridders})))
  (is= {::anomalies/message  (str "function values at endpoints do not have different signs, "
                               "endpoints: [0.1, 1.8], values: [0.005, 2.165]")
        ::anomalies/fn       (var apache-solvers/root-solver)
        ::anomalies/category ::anomalies/third-party}
    (root-solvers/root-solver {::root-solvers/univariate-f f
                               ::root-solvers/guess        guess
                               ::root-solvers/interval     [0.1 1.8]}
      {::root-solvers/root-solver-type :ridders}))
  (is= 0.9148015417397091
    (root-solvers/root-solver {::root-solvers/univariate-f f
                               ::root-solvers/guess        guess
                               ::root-solvers/interval     [0.5 1.4]}
      {::root-solvers/root-solver-type :ridders}))
  ;;secant
  (is (anomalies/anomaly?
        (root-solvers/root-solver {::root-solvers/univariate-f f
                                   ::root-solvers/guess        guess
                                   ::root-solvers/interval     [-20.0 20.0]}
          {::root-solvers/root-solver-type :secant})))
  (is= {::anomalies/message  (str "function values at endpoints do not have different signs, "
                               "endpoints: [0.1, 1.8], values: [0.005, 2.165]")
        ::anomalies/fn       (var apache-solvers/root-solver)
        ::anomalies/category ::anomalies/third-party}
    (root-solvers/root-solver {::root-solvers/univariate-f f
                               ::root-solvers/guess        guess
                               ::root-solvers/interval     [0.1 1.8]}
      {::root-solvers/root-solver-type :secant}))
  (is= 0.9148015417354283
    (root-solvers/root-solver {::root-solvers/univariate-f f
                               ::root-solvers/guess        guess
                               ::root-solvers/interval     [0.5 1.4]}
      {::root-solvers/root-solver-type :secant})))

(deftest brent-dekker-test
  (is (spec-check root-solvers/brent-dekker)))

(deftest modified-newton-raphson-test
  (is (spec-check root-solvers/modified-newton-raphson
        {:coll-check-limit 10
         :coll-error-limit 10
         :fspec-iterations 10
         :recursion-limit  1
         :test-check       {:num-tests 300}})))

#_(ost/unstrument)