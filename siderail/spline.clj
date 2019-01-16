(require '[provisdom.solvers.curve-fitting :refer :all]
         '[oz.core :as oz]
         '[clojure.pprint :refer [pprint]]
         '[provisdom.math.apache-matrix :as amx]
         '[provisdom.math.matrix :as mx]
         '[provisdom.math.neanderthal-matrix :as nmx]
         '[uncomplicate.neanderthal.core :as n]
         '[uncomplicate.neanderthal.linalg :as la]
         '[uncomplicate.neanderthal.native :as nat]
         #_'[uncomplicate.neanderthal.opencl :as ocl])

(import '(umontreal.ssj.functionfit SmoothingCubicSpline))

(defonce plot-server (oz/start-plot-server!))

(comment
  (def rho 1.0)
  (def x (vec (range 2.0 18.0 2.0)))
  (def y (mapv * x x))
  (def x' (range 0.0 20.0 1.0))
  (def s (smoothing-cubic-spline x y rho))
  (def s' (SmoothingCubicSpline. (double-array x) (double-array y) rho))
  (evaluate s 4.0)
  (map * x' x')
  (def y' (map (partial evaluate s) x'))
  (def y'' (map #(.evaluate s' %) x'))

  :end)

(comment
  (def rho 0.9991008092716556)
  (def x (vec (range 0.0 (* 2.0 Math/PI) (/ (* 2.0 Math/PI) 100.0))))
  (def y (mapv #(+ (Math/sin %) (* (rand))) x))
  (def v (repeat (count y) (/ 1.0)))
  (def s (smoothing-cubic-spline x y v rho))
  (def s' (SmoothingCubicSpline. (double-array x) (double-array y) rho))
  (def x' (range -1 7 0.1))
  (def y' (map (partial evaluate s) x'))
  (def y'' (map #(.evaluate s' %) x'))
  :end)

(defn scatter
  [data x-field y-field color-field]
  {:layer    [{:data     {:values (filter #(= :prediction (:type %)) data)}
               :encoding {:x {:field x-field}
                          :y {:field y-field}}
               :mark     :line}
              {:data     {:values (filter #(= :data (:type %)) data)}
               :encoding {:x {:field x-field}
                          :y {:field y-field}}
               :mark     :circle}]
   :encoding {:color {:field :type :type :nominal}}
   :width    500 :height 500})

(comment
  (let [data (vec (concat (map (fn [x y] {:x x :y y :type :data}) x y)
                          (map (fn [x y] {:x x :y y :type :prediction}) x' y')))]
    (oz/v! (scatter data :x :y nil)))
  :end)


(comment
  (def rho 0.9)
  (def x (vec (sort (repeatedly 100 #(rand (* 10.0 2.0 Math/PI))))))
  (def y (mapv #(* 2.0 (+ (Math/sin %) (* (rand)))) x))
  (def s (smoothing-cubic-spline x y rho))
  (def rhos (range 0.0 1.001 0.001))
  #_(def lambdas (mapv #(Math/pow 10 %) (range -20 0)))
  (def splines (mapv (partial smoothing-cubic-spline x y) rhos))
  (let [data (mapv (fn [s] {:rho (:rho s) :dof (dof s) :type :prediction}) splines)]
    #_(pprint data)
    (oz/v! (scatter data :rho :dof nil)))

  :end)

(comment
  (def rho 0.9)
  (def x (vec (sort (repeatedly 100 #(rand (* 10.0 2.0 Math/PI))))))
  (def y (mapv #(* 2.0 (+ (Math/sin %) (* (rand)))) x))
  (def s (smoothing-cubic-spline x y rho))
  (def vs (range 0.01 10.0 0.01))
  #_(def lambdas (mapv #(Math/pow 10 %) (range -20 0)))
  (def splines (mapv #(smoothing-cubic-spline x y (repeat (count x) %) rho) vs))
  (let [data (mapv (fn [s] {:log-variance (Math/log10 (first (:vs s))) :dof (dof s) :type :prediction}) splines)]
    #_(pprint data)
    (oz/v! (scatter data :log-variance :dof nil)))

  :end)