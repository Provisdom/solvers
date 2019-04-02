(require '[provisdom.solvers.beta-regression :refer :all :as br]
         '[provisdom.math.special-functions :as sfn]
         '[provisdom.solvers.curve-fitting :as cf]
         '[provisdom.math.core :as m]
         '[uncomplicate.neanderthal.native :as nat]
         '[uncomplicate.neanderthal.core :as n]
         '[uncomplicate.neanderthal.vect-math :as vec]
         '[criterium.core :as crit])

(def data
  [[38.4 6.1 220 235 6.9]
   [40.3 4.8 231 307 14.4]
   [40.0 6.1 217 212 7.4]
   [31.8 0.2 316 365 8.5]
   [40.8 3.5 210 218 8.0]
   [41.3 1.8 267 235 2.8]
   [38.1 1.2 274 285 5.0]
   [50.8 8.6 190 205 12.2]
   [32.2 5.2 236 267 10.0]
   [38.4 6.1 220 300 15.2]
   [40.3 4.8 231 367 26.8]
   [32.2 2.4 284 351 14.0]
   [31.8 0.2 316 379 14.7]
   [41.3 1.8 267 275 6.4]
   [38.1 1.2 274 365 17.6]
   [50.8 8.6 190 275 22.3]
   [32.2 5.2 236 360 24.8]
   [38.4 6.1 220 365 26.0]
   [40.3 4.8 231 395 34.9]
   [40.0 6.1 217 272 18.2]
   [32.2 2.4 284 424 23.2]
   [31.8 0.2 316 428 18.0]
   [40.8 3.5 210 273 13.1]
   [41.3 1.8 267 358 16.1]
   [38.1 1.2 274 444 32.1]
   [50.8 8.6 190 345 34.7]
   [32.2 5.2 236 402 31.7]
   [38.4 6.1 220 410 33.6]
   [40.0 6.1 217 340 30.4]
   [40.8 3.5 210 347 26.6]
   [41.3 1.8 267 416 27.8]
   [50.8 8.6 190 407 45.7]])

(def data
  [[1 50.8 8.6 190 1 0 0 0 0 0 0 0 0 205 12.2]
   [1 50.8 8.6 190 1 0 0 0 0 0 0 0 0 275 22.3]
   [1 50.8 8.6 190 1 0 0 0 0 0 0 0 0 345 34.7]
   [1 50.8 8.6 190 1 0 0 0 0 0 0 0 0 407 45.7]
   [2 40.8 3.5 210 0 1 0 0 0 0 0 0 0 218 8.0]
   [2 40.8 3.5 210 0 1 0 0 0 0 0 0 0 273 13.1]
   [2 40.8 3.5 210 0 1 0 0 0 0 0 0 0 347 26.6]
   [3 40.0 6.1 217 0 0 1 0 0 0 0 0 0 212 7.4]
   [3 40.0 6.1 217 0 0 1 0 0 0 0 0 0 272 18.2]
   [3 40.0 6.1 217 0 0 1 0 0 0 0 0 0 340 30.4]
   [4 38.4 6.1 220 0 0 0 1 0 0 0 0 0 235 6.9]
   [4 38.4 6.1 220 0 0 0 1 0 0 0 0 0 300 15.2]
   [4 38.4 6.1 220 0 0 0 1 0 0 0 0 0 365 26.0]
   [4 38.4 6.1 220 0 0 0 1 0 0 0 0 0 410 33.6]
   [5 40.3 4.8 231 0 0 0 0 1 0 0 0 0 307 14.4]
   [5 40.3 4.8 231 0 0 0 0 1 0 0 0 0 367 26.8]
   [5 40.3 4.8 231 0 0 0 0 1 0 0 0 0 395 34.9]
   [6 32.2 5.2 236 0 0 0 0 0 1 0 0 0 267 10.0]
   [6 32.2 5.2 236 0 0 0 0 0 1 0 0 0 360 24.8]
   [6 32.2 5.2 236 0 0 0 0 0 1 0 0 0 402 31.7]
   [7 41.3 1.8 267 0 0 0 0 0 0 1 0 0 235 2.8]
   [7 41.3 1.8 267 0 0 0 0 0 0 1 0 0 275 6.4]
   [7 41.3 1.8 267 0 0 0 0 0 0 1 0 0 358 16.1]
   [7 41.3 1.8 267 0 0 0 0 0 0 1 0 0 416 27.8]
   [8 38.1 1.2 274 0 0 0 0 0 0 0 1 0 285 5.0]
   [8 38.1 1.2 274 0 0 0 0 0 0 0 1 0 365 17.6]
   [8 38.1 1.2 274 0 0 0 0 0 0 0 1 0 444 32.1]
   [9 32.2 2.4 284 0 0 0 0 0 0 0 0 1 351 14.0]
   [9 32.2 2.4 284 0 0 0 0 0 0 0 0 1 424 23.2]
   [10 31.8 0.2 316 0 0 0 0 0 0 0 0 0 365 8.5]
   [10 31.8 0.2 316 0 0 0 0 0 0 0 0 0 379 14.7]
   [10 31.8 0.2 316 0 0 0 0 0 0 0 0 0 428 18.0]])

(def rs (vec (repeatedly 2 #_1000000 #(vec (repeatedly 100 rand)))))
(def rv (mapv nat/dv rs))
(def ones (nat/dv (repeat 100 1.0)))
(def work (n/copy ones))
(def m1 (nat/dv (repeat 100 -1.0)))
(def r (first rv))
(defn f
  [x y]
  (sfn/logistic (+ x (- y) (* x y))))

#_(def data (for [x (range -1 1 0.1)
                  y(range -1 1 0.1)]
              [x y (+ (f x y) (rand 0.1))]))

(def x-mx #_(mapv #(->> % butlast (cons 1.0) vec) data)
  (mapv #(->> %
              (drop 4)
              (butlast)
              (map double)
              (cons 1.0)
              vec)
        data))
(def xm (nat/dge x-mx))
(def y (mapv #(/ (last %) 100.0) data))
(def n (count y))
(def k (count (first x-mx)))

(def z (map sfn/logit y))

(def fit (cf/linear-least-squares-curve-fitting {::cf/x-matrix x-mx ::cf/f-vals z ::cf/curve-basis-fn identity}))
(def beta0 (::cf/curve-fitting-weights fit))
(def bv (nat/dv beta0))
(def z' (mapv (::cf/curve-fitting-fn fit) x-mx))
(def sigma (/ (apply + (map m/sq (map - z z'))) (- n k)))
(def mu-tah (map sfn/logistic z'))

(def phi0 (dec (/ (apply + (map #(/ (* sigma (* % (m/one- %)))) mu-tah)) n)))

(comment
  (time (solve x-mx y (vec (cons phi0 beta0)) {::br/abs-accu 1e-8 ::br/max-evaluations 10000000 ::br/max-iter 10000000}))

  (gradient x-mx y (vec (cons phi0 beta0)))

  :end)