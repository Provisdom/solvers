(ns provisdom.solvers.stepwise-regression
  (:require
    [clojure.spec.alpha :as s]
    [clojure.spec.gen.alpha :as gen]
    [clojure.spec.test.alpha :as st]
    [orchestra.spec.test :as ost]
    [provisdom.utility-belt.anomalies :as anomalies]
    [provisdom.math.core :as m]
    [provisdom.math.vector :as vector]
    [provisdom.math.matrix :as mx]
    [provisdom.math.neanderthal-matrix :as neanderthal-mx]
    [provisdom.solvers.logistic-regression :as log-regress]))

(declare ordinary-stepwise-regression-fn)

(s/def ::independent-variables
  (s/and ::vector/vector-finite
         (fn [v]
           (pos? (count v)))))

(s/def ::independent-variables-mx
  (s/and ::mx/matrix-finite
         (fn [mx]
           (and
             (>= (mx/rows mx) 2)
             (pos? (mx/columns mx))))))

(s/def ::weights ::vector/vector-num)

(s/def ::error ::m/non-)

(s/def ::x
  (s/and ::vector/vector-finite
         (fn [v]
           (pos? (count v)))))

(s/def ::x-mx
  (s/and ::mx/matrix-finite
         (fn [mx]
           (and
             (>= (mx/rows mx) 2)
             (pos? (mx/columns mx))))))

(s/def ::y
  (s/and ::vector/vector-finite
         (fn [y]
           (>= (count y) 2))))

(s/def ::basis-fn
  (s/fspec :args (s/cat :x ::x)
           :ret ::independent-variables))

(s/def ::extra-degrees ::m/long-non-)

(s/def ::component-name
  (s/with-gen any?
              #(s/gen (s/or :keyword keyword?
                            :number ::m/number))))

(s/def ::component-info
  (s/keys :req [::basis-fn]
          :opt [::extra-degrees]))

(s/def ::component-group
  (s/with-gen
    (s/map-of ::component-name ::component-info)
    #(gen/map (s/gen ::component-name)
              (s/gen ::component-info)
              {:min-elements 0
               :max-elements 3})))

(s/def ::component-groups
  (s/with-gen (s/coll-of ::component-group)
              #(gen/vector (s/gen ::component-group) 0 3)))

(s/def ::alternative-component-groups ::component-groups)

(s/def ::regression-fn
  (s/with-gen
    (s/fspec :args (s/and (s/cat :independent-variables-mx ::independent-variables-mx
                                 :y ::y)
                          (fn [{:keys [independent-variables-mx y]}]
                            (= (mx/rows independent-variables-mx) (count y))))
             :ret (s/or :solution (s/keys :req [::weights ::error])
                        :anomaly ::anomalies/anomaly))
    #(gen/return (ordinary-stepwise-regression-fn))))

(s/def ::prior ::m/prob)
(s/def ::data-count ::m/long+)
(s/def ::parameter-count ::m/long+)

(s/def ::selection-score-fn
  (s/fspec :args (s/cat :args (s/keys :req [::prior
                                            ::error
                                            ::data-count
                                            ::parameter-count]))
           :ret ::m/num))

(s/def ::prob-of-model-fn
  (s/with-gen
    (s/fspec :args (s/cat :component-group ::component-group)
             :ret ::m/prob)
    #(gen/return (fn [cg] 0.1))))

(s/def ::score ::m/num)

(s/def ::best-component-group-info
  (s/keys :req [::score ::component-group ::error ::weights]))

(s/def ::solve-type #{:forward :backward :both})

(s/def ::max-iter
  (s/with-gen (s/nilable ::m/int+)
              #(gen/one-of [(s/gen (s/int-in 1 3))
                            (gen/return nil)])))

(defn ordinary-stepwise-regression-fn
  ""
  []
  (fn [independent-variables-mx y]
    (let [y-nmx (neanderthal-mx/matrix->neanderthal-matrix
                  (mx/column-matrix y))
          ind-vars-nmx (neanderthal-mx/matrix->neanderthal-matrix
                         independent-variables-mx)
          soln (neanderthal-mx/lls-with-error
                 ind-vars-nmx
                 y-nmx)]
      (if (anomalies/anomaly? soln)
        soln
        (let [{::neanderthal-mx/keys [mean-squared-errors solution condition]} soln
              mse (neanderthal-mx/neanderthal-matrix->matrix mean-squared-errors)
              sol (neanderthal-mx/neanderthal-matrix->matrix solution)
              weights (vec (flatten sol))
              error (max (ffirst mse) 0.0)]
          (cond
            (< condition 1e-8)
            {::anomalies/message  "poorly conditioned matrix"
             ::anomalies/category ::anomalies/exception
             ::neanderthal-mx/condition condition
             ::anomalies/fn       (var ordinary-stepwise-regression-fn)}

            (not-any? m/nan? weights)
            {::weights weights
             ::error   error}

            :else
            {::anomalies/message  "weights can't be NaN"
             ::anomalies/category ::anomalies/exception
             ::anomalies/fn       (var ordinary-stepwise-regression-fn)}))))))

(s/fdef ordinary-stepwise-regression-fn
        :args (s/cat)
        :ret ::regression-fn)

(defn logistic-stepwise-regression-fn
  ""
  ([] (logistic-stepwise-regression-fn {}))
  ([{::log-regress/keys [max-iter abs-accu]}]
   (fn [independent-variables-mx y]
     (if (and (vector/vector-prob? y)
              (not (every? m/one? y))
              (not (every? zero? y)))
       (let [sol (log-regress/iterative-reweighted-least-squares
                   independent-variables-mx
                   y
                   {::log-regress/max-iter max-iter
                    ::log-regress/abs-accu (or abs-accu 1e-9)})]
         (if (anomalies/anomaly? sol)
           sol
           ;;(log-regress/standard-squared-error sol)
           (let [error 0.0]
             {::weights sol
              ::error   error})))
       {::anomalies/message  (str "`y` must be a 'vector-prob', and not every value"
                                  " can be zero, and not every value can be one.")
        ::anomalies/category ::anomalies/exception
        ::anomalies/fn       (var logistic-stepwise-regression-fn)}))))

(s/fdef logistic-stepwise-regression-fn
        :args (s/cat :opts (s/? (s/keys :opt [::log-regress/max-iter
                                              ::log-regress/abs-accu])))
        :ret ::regression-fn)

(defn least-squares-bic-selection-score-fn
  [{::keys [prior error data-count parameter-count]}]
  (let [log-p (* -0.5
                 data-count
                 (inc (m/log
                        (* m/two-pi
                           error))))
        bic (- log-p
               (* 0.5
                  parameter-count
                  (m/log data-count)))]
    (if (zero? prior)
      m/inf-
      (+ (m/log prior) bic))))

(s/def least-squares-bic-selection-score-fn ::selection-score-fn)

(defn one-step-solve
  ""
  [{::keys [x-mx y regression-fn selection-score-fn prob-of-model-fn
            best-component-group-info alternative-component-groups]}]
  (let [data-count (count y)
        component-group-count (count alternative-component-groups)
        new-best (loop [i 0
                        best best-component-group-info]
                   (if (>= i component-group-count)
                     best
                     (let [component-group (nth alternative-component-groups i)
                           b-fn (fn [x]
                                  (vec (mapcat (fn [[_ {::keys [basis-fn _]}]]
                                                 (basis-fn x))
                                               component-group)))
                           parameter-count (+ 1
                                              (count (b-fn (first x-mx)))
                                              (reduce
                                                (fn [tot component]
                                                  (+ tot
                                                     (get component ::extra-degrees 0)))
                                                0
                                                component-group))
                           prior (prob-of-model-fn component-group)
                           best-possible-new-score (selection-score-fn
                                                     {::prior           prior
                                                      ::error           0.0
                                                      ::data-count      data-count
                                                      ::parameter-count parameter-count})
                           best-score (::score best)]
                       (if (<= best-possible-new-score best-score)
                         (recur (inc i) best)
                         (let [ind-vars-mx (mapv b-fn x-mx)
                               sol (regression-fn ind-vars-mx y)]
                           (if (anomalies/anomaly? sol)
                             (recur (inc i) best)
                             (let [{::keys [error weights]} sol
                                   score (selection-score-fn
                                           {::prior           prior
                                            ::error           error
                                            ::data-count      data-count
                                            ::parameter-count parameter-count})]
                               (if (> score best-score)
                                 (recur (inc i) {::score           score
                                                 ::component-group component-group
                                                 ::error           error
                                                 ::weights         weights})
                                 (recur (inc i) best)))))))))]
    new-best))

(s/fdef one-step-solve
        :args (s/cat :args (s/and (s/keys :req [::x-mx
                                                ::y
                                                ::regression-fn
                                                ::selection-score-fn
                                                ::prob-of-model-fn
                                                ::best-component-group-info
                                                ::alternative-component-groups])
                                  (fn [{::keys [x-mx y]}]
                                    (= (mx/rows x-mx) (count y)))))
        :ret ::best-component-group-info)

(defn one-step-forward-solve-alternative-component-groups
  ""
  [best-component-group component-group]
  (mapv (fn [[component-name component-info]]
          (assoc best-component-group component-name component-info))
        component-group))

(s/fdef one-step-forward-solve-alternative-component-groups
        :args (s/cat :best-component-group ::component-group
                     :component-group ::component-group)
        :ret ::alternative-component-groups)

(defn one-step-backward-solve-alternative-component-groups
  ""
  [best-component-group]
  (mapv (fn [[component-name _]]
          (dissoc best-component-group component-name))
        best-component-group))

(s/fdef one-step-backward-solve-alternative-component-groups
        :args (s/cat :best-component-group ::component-group)
        :ret ::alternative-component-groups)

(defn solve
  "Optional:
  `max-iter` counts both forward and backward moves (default is 10).
  `solve-type` can be `:forward`, `:backward`, or `both` (default).
  `regression-fn` uses ordinary regression by default."
  ([args] (solve args {}))
  ([{::keys [x-mx y component-group selection-score-fn prob-of-model-fn]}
    {::keys [max-iter solve-type regression-fn]}]
   (let [regression-fn (or regression-fn (ordinary-stepwise-regression-fn))
         max-iter (or max-iter 10)
         starting-component-group (if (= solve-type :backward)
                                    component-group
                                    {})]
     (loop [i 0
            tested-component-groups [{}]
            best-component-group-info {::component-group starting-component-group
                                       ::score           m/min-dbl
                                       ::error           m/max-dbl
                                       ::weights         []}
            forward? (not= solve-type :backward)
            last-try? (or (= solve-type :backward) (= solve-type :forward))]
       (if (>= i max-iter)
         best-component-group-info
         (let [alternative-component-groups (if forward?
                                              (one-step-forward-solve-alternative-component-groups
                                                (::component-group best-component-group-info)
                                                component-group)
                                              (one-step-backward-solve-alternative-component-groups
                                                (::component-group best-component-group-info)))
               alternative-component-groups (reduce (fn [tot alt-group]
                                                      (let [alt-keys (set (keys alt-group))]
                                                        (if (some (fn [test-group]
                                                                    (= (set (keys test-group))
                                                                       alt-keys))
                                                                  tested-component-groups)
                                                          tot
                                                          (conj tot alt-group))))
                                                    []
                                                    alternative-component-groups)
               tested-component-groups (concat tested-component-groups
                                               alternative-component-groups)
               new-best (one-step-solve
                          {::x-mx                         x-mx
                           ::y                            y
                           ::regression-fn                regression-fn
                           ::selection-score-fn           selection-score-fn
                           ::prob-of-model-fn             prob-of-model-fn
                           ::best-component-group-info    best-component-group-info
                           ::alternative-component-groups alternative-component-groups})
               new-forward? (condp = solve-type
                              :forward true
                              :backward false
                              (not forward?))]
           (if (and last-try? (= new-best best-component-group-info))
             new-best
             (recur (inc i)
                    tested-component-groups
                    new-best
                    new-forward?
                    (= new-best best-component-group-info)))))))))

(s/fdef solve
        :args (s/cat :args (s/and (s/keys :req [::x-mx
                                                ::y
                                                ::component-group
                                                ::selection-score-fn
                                                ::prob-of-model-fn])
                                  (fn [{::keys [x-mx y]}]
                                    (= (mx/rows x-mx) (count y))))
                     :opts (s/? (s/keys :opt [::max-iter ::solve-type ::regression-fn])))
        :ret ::best-component-group-info)