(ns provisdom.solvers.stepwise-regression-test
  (:require
    [clojure.test :refer :all]
    [provisdom.test.core :refer :all]
    [provisdom.solvers.stepwise-regression :as stepwise]
    [provisdom.math.core :as m]
    [provisdom.math.neanderthal-matrix :as neanderthal-mx]
    [clojure.spec.test.alpha :as st]
    [orchestra.spec.test :as ost]))

;;94 seconds

(set! *warn-on-reflection* true)

(ost/instrument)

(def x-mx
  [[7.0 2.8] [8.0 9.3] [3.0 9.9]
   [3.3 2.2] [3.8 5.9] [2.2 8.9]
   [4.6 2.5] [2.4 3.5] [8.4 2.3]])

(def y [1 2 3 4 5 6 7 8 9])

(defn prob-of-model-fn
  [component-group]
  1.0)

(deftest ordinary-stepwise-regression-fn-test
  (is (spec-check stepwise/ordinary-stepwise-regression-fn
                  {:coll-check-limit 10
                   :coll-error-limit 10
                   :fspec-iterations 10
                   :recursion-limit  1
                   :test-check       {:num-tests 300}}))
  (is (data-approx=
        {::stepwise/weights [0.624874327896168 0.22810879461175296]
         ::stepwise/error   12.177883256018546}
        ((stepwise/ordinary-stepwise-regression-fn)
          x-mx
          y))))

(deftest logistic-stepwise-regression-fn-test
  (is (spec-check stepwise/logistic-stepwise-regression-fn
                  {:coll-check-limit 10
                   :coll-error-limit 10
                   :fspec-iterations 10
                   :recursion-limit  1
                   :test-check       {:num-tests 250}}))
  (is= {::stepwise/weights [0.02918512018138139 -0.0516875178183345]
        ::stepwise/error   0.0}
       ((stepwise/logistic-stepwise-regression-fn)
         x-mx
         [0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9])))

(deftest one-step-solve-test
  (is (spec-check stepwise/one-step-solve
                  {:coll-check-limit 10
                   :coll-error-limit 10
                   :fspec-iterations 10
                   :recursion-limit  1
                   :test-check       {:num-tests 5}}))
  (is (data-approx=
        #::stepwise
            {:score           -27.31458022888188
             :component-group [:a :b]
             :error           12.177883256018546
             :weights         [0.624874327896168 0.22810879461175296]}
        (update
          (stepwise/one-step-solve
            #::stepwise
                {:x-mx                         x-mx
                 :y                            y
                 :regression-fn                (stepwise/ordinary-stepwise-regression-fn)
                 :selection-score-fn           stepwise/least-squares-bic-selection-score-fn
                 :prob-of-model-fn             prob-of-model-fn
                 :best-component-group-info    #::stepwise
                                                   {:score           -1000.0
                                                    :component-group {:a {::stepwise/basis-fn
                                                                          (fn [x]
                                                                            [(get x 0 0.0)])}}
                                                    :error           0.1
                                                    :weights         [0.3]}
                 :alternative-component-groups [{:a {::stepwise/basis-fn
                                                     (fn [x]
                                                       [(get x 0 0.0)])}
                                                 :b {::stepwise/basis-fn
                                                     (fn [x]
                                                       [(get x 1 0.0)])}}]})
          ::stepwise/component-group
          keys))))

(deftest one-step-forward-solve-alternative-component-groups-test
  (is (spec-check stepwise/one-step-forward-solve-alternative-component-groups))
  (is= [[:a :b]]
       (mapv keys (stepwise/one-step-forward-solve-alternative-component-groups
                    {:a {::stepwise/basis-fn
                         (fn [x]
                           [(get x 0 0.0)])}}
                    {:b {::stepwise/basis-fn
                         (fn [x]
                           [(get x 1 0.0)])}}))))

(deftest one-step-backward-solve-alternative-component-groups-test
  (is (spec-check stepwise/one-step-backward-solve-alternative-component-groups))
  (is= #{'(:a) '(:b)}
       (set (map keys (stepwise/one-step-backward-solve-alternative-component-groups
                        {:a {::stepwise/basis-fn
                             (fn [x]
                               [(get x 0 0.0)])}
                         :b {::stepwise/basis-fn
                             (fn [x]
                               [(get x 1 0.0)])}})))))

(deftest solve-test
  (is (spec-check stepwise/solve
                  {:coll-check-limit 10
                   :coll-error-limit 10
                   :fspec-iterations 10
                   :recursion-limit  1
                   :test-check       {:num-tests 8}}))
  (is= #::stepwise
           {:score           -26.534808032342788
            :component-group [:a]
            :error           13.0720297305116
            :weights         [0.81874624474264]}
       (update
         (stepwise/solve
           #::stepwise
               {:x-mx               x-mx
                :y                  y
                :component-group    {:a {::stepwise/basis-fn
                                         (fn [x]
                                           [(get x 0 0.0)])}
                                     :b {::stepwise/basis-fn
                                         (fn [x]
                                           [(get x 1 0.0)])}}
                :selection-score-fn stepwise/least-squares-bic-selection-score-fn
                :prob-of-model-fn   prob-of-model-fn})
         ::stepwise/component-group
         keys)))

#_(ost/unstrument)