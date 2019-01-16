(require
  '[clojure.test :refer :all]
  '[provisdom.test.core :refer :all]
  '[provisdom.solvers.stepwise-regression :as stepwise]
  '[provisdom.math.core :as m]
  '[clojure.spec.test.alpha :as st]
  '[orchestra.spec.test :as ost])

(set! *warn-on-reflection* true)

(defonce foo (ost/instrument))


(update
  (stepwise/solve
    {::stepwise/x-mx               [[7.0 2.8] [8.0 9.3] [3.0 9.9]
                                    [3.3 2.2] [3.8 5.9] [2.2 8.9]
                                    [4.6 2.5] [2.4 3.5] [8.4 2.3]]
     ::stepwise/y                  [1 2 3 4 5 6 7 8 9]
     ::stepwise/component-group    {:a {::stepwise/basis-fn
                                        (fn [x]
                                          [(get x 0 0.0)])}
                                    :b {::stepwise/basis-fn
                                        (fn [x]
                                          [(get x 1 0.0)])}}
     ::stepwise/selection-score-fn (fn [error
                                        data-count
                                        parameter-count]
                                     (let [log-p (* -0.5
                                                    data-count
                                                    (inc (m/log
                                                           (* m/two-pi
                                                              error))))
                                           bic (- log-p
                                                  (* 0.5
                                                     parameter-count
                                                     (m/log data-count)))]
                                       (m/exp bic)))
     ::stepwise/prob-of-model-fn   (fn [component-group]
                                     1.0)})
  ::stepwise/component-group
  keys)