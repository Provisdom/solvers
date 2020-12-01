(ns provisdom.solvers.integer-root-solvers
  (:require
    [clojure.spec.alpha :as s]
    [clojure.spec.gen.alpha :as gen]
    [provisdom.math.core :as m]
    [provisdom.math.intervals :as intervals]))

(defn integer-root-solver
  "Bisection algorithm. Returns the minimum integer that has a functional value
  greater than or equal to zero. Function must be non-positive at the `lower`
  bound and non-negative at the `upper` bound."
  [strictly-increasing-fn [lower upper]]
  (loop [n (m/round (+ (* 0.5 lower) (* 0.5 upper)) :down)
         l lower
         u upper]
    (let [value (strictly-increasing-fn n)]
      (cond (zero? value) n
            (pos? value) (if (== n l)
                           n
                           (recur (m/round (+ (* 0.5 n) (* 0.5 l)) :down) l n))
            :else (recur (m/round (+ (* 0.5 n) (* 0.5 u)) :up) (inc n) u)))))

(s/fdef integer-root-solver
        :args (s/and (s/cat :increasing-fn (s/with-gen
                                             (s/fspec :args (s/cat :n ::m/long)
                                                      :ret ::m/num)
                                             #(gen/one-of
                                                (map gen/return
                                                     (list m/exp m/sq))))
                            :long-interval ::intervals/long-interval)
                     (fn [{:keys [increasing-fn long-interval]}]
                       (let [[lower upper] long-interval]
                         (and (m/non+? (increasing-fn lower))
                              (m/non-? (increasing-fn upper))))))
        :ret ::m/long)