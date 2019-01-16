(ns provisdom.solvers.asd
  (:require
    [uncomplicate.neanderthal.core :as core]
    [uncomplicate.neanderthal.native :as native]))

(def a (native/dge 3 2 [1 2 3 4 5 6]))
(def b (native/dge 2 3 [10 20 30 40 50 60]))
