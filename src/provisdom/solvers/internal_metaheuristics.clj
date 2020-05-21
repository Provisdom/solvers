(ns provisdom.solvers.internal-metaheuristics
  (:require
    [clojure.spec.alpha :as s]
    [clojure.spec.gen.alpha :as gen]
    [clojure.spec.test.alpha :as st]
    [orchestra.spec.test :as ost]
    [provisdom.math.core :as m]
    [clojure.set :only [intersection difference]])
  (:import [java.util ArrayList Collections]))

;;;Each particle could follow different rules; 
;;;dynamic constriction and neighborhood pull, small inertia and global pull
;;;Allow some initial points to be set; grid or random start or combo 
;;;like random from grid
;;; All constraints dynamically in objective, squared power of error or so; 
;;;can not allow some particles movement to break equality dynamically 
;;;shrinking regions
;;; mutation/crossover/PSO in one algo

;;;OPT4j has tons of links to java libraries

;;;Look at .NET PSO
;;;http://www.particleswarm.info/Programs.html
;;;also probably want to write our own DoNLP2 from .NET since that can be 
;;;used as a hill climber

;;;The original versions of these came from here: 
;;;https://github.com/dermatthias/metaheuristics

(def DEBUG? false)

;; replace this with distributions.random
(defn- normal [mu sigma]
  (let [r (new java.util.Random)]
    (+ mu (* sigma (.nextGaussian r)))))

;;this makes no sense at all...
(defn- scramble [l]
  (let [items (ArrayList.),
        _ (doseq [li l] (.add items li))]
    (do (Collections/shuffle items) (seq items))))

;;Invasive weed optimization
(defstruct plant :seedlist :position :pfit :tfit)

(defstruct population
  :plantlist :gbest :gworst :minseed :maxseed :sigmaInit :sigmaFinal)

(defn- init-plant [nDimensions max]
  (agent {:seedlist (list),
          :position (double-array (for [i (range nDimensions)] (rand max))),
          :pfit     0.0, :tfit 0.0}))

(defn- init-population
  [nPlants nDimensions nSeedMin nSeedMax sigmaInit sigmaFinal max]
  (let [plants (map (fn [_] (init-plant nDimensions max)) (range nPlants))
        gbest (init-plant nDimensions max)
        gworst (init-plant nDimensions max)]
    (agent (struct population plants gbest gworst (int nSeedMin) (int nSeedMax)
                   sigmaInit sigmaFinal))))

(defn- set-bestworst [population ftype]
  (let [plants (for [p (:plantlist population)] @p)
        sorted (sort-by :pfit ftype plants)]
    (assoc population :gbest (agent (first sorted)) :gworst (agent (last sorted)))))

;; eval seeds
(defn- eval-seed [seed fitness]
  (let [^doubles pos (:position seed), fit (fitness pos)]
    (assoc seed :pfit fit)))

(defn- eval-plantseeds [plant fitness]
  (let [seeds (:seedlist @plant)]
    (doseq [s seeds] (send s eval-seed fitness))
    (apply await seeds)))

;; create new seeds
(defn- create-new-seed [^doubles pos sigma]
  (let [dim (count pos),
        pos-offset (double-array (map (fn [_] (normal 0 sigma)) (range dim))),
        newpos (amap pos i ret (+ (aget pos i) (aget pos-offset i)))]
    (agent {:seedlist (list), :position newpos, :pfit 0.0, :tfit 0.0})))

(defn- generate-seeds [plant population maxIt modulation iteration]
  "Generates new seeds"
  (let [pfit (double (:pfit plant))
        ^doubles pos (:position plant)
        gbestfit (double (:pfit @(:gbest population)))
        gworstfit (double (:pfit @(:gworst population)))
        minSeed (int (:minseed population))
        maxSeed (int (:maxseed population))
        sigmaInit (double (:sigmaInit population))
        sigmaFinal (double (:sigmaFinal population))
        nSeeds (int (+ (* (/ (- pfit gworstfit) (- gbestfit gworstfit))
                          maxSeed)
                       (* (/ (- pfit gbestfit) (- gworstfit gbestfit))
                          minSeed)))
        sigma (+ (* (/ (m/pow (- maxIt iteration) modulation)
                       (m/pow maxIt modulation))
                    (- sigmaInit sigmaFinal))
                 sigmaFinal)
        newSeeds (map (fn [_] (create-new-seed pos sigma)) (range nSeeds))]
    (assoc plant :seedlist newSeeds)))

(defn- competition [population ftype limit style]
  (let [seeds (for [p (:plantlist population)
                    s (:seedlist @p)]
                @s)
        plants (for [p (:plantlist population)]
                 @p)
        ;; best up to limit
        newpop (take limit (sort-by :pfit ftype (concat seeds plants)))
        newagents (for [np newpop]
                    (agent np))]
    (assoc population :plantlist newagents
                      :gbest (agent (first newpop))
                      :gworst (agent (last newpop)))))

(defn- grow [fitness ftype population maxIt nPlantsMax modulation iteration]
  ;; generate new seeds for each plant-agent
  (if DEBUG? (println "creating new seeds..."))
  (dorun (map #(send % generate-seeds @population maxIt modulation iteration)
              (:plantlist @population)))
  (apply await (:plantlist @population))
  (if DEBUG?
    (println "no. of seeds:" (count (for [p (:plantlist @population) seeds
                                          (:seedlist @p)] seeds))))
  ;; update seed fitness
  (if DEBUG? (println "updating seed fitness..."))
  ;; could be parallelized
  (dorun (map #(eval-plantseeds % fitness) (:plantlist @population)))
  ;; competitive exclusion
  (if DEBUG? (println "competition..."))
  (send population competition ftype nPlantsMax 1)
  (await population)
  (if DEBUG? (println "best in generation" iteration ":"
                      (:pfit @(:gbest @population)) "\n---")))

(defn iwo
  "Starts the IWO algorithm. Algorithm based on 
   http://dx.doi.org/10.1016/j.ecoinf.2006.07.003
   Parameters:
   fitness - the fitness function to be used. only one parameter: 
               the position of the particle (double-array).
   ftype - defines the type of optimization problem (minimize (<) 
               or maximize (>)).
   dim - number of dimensions in solution.
   nplants - number of initial plants.
   nplants-max - maximum number of plants in population.
   seed-min - minimum number of seeds generated per plant.
   seed-max - maximum number of seeds generated per plant.
   sigma-init - inital value for sigma (standard deviation).
   sigma-final - final value for sigma (standard deviation).
   modulation -  modulation index (usually 3). 
   max-iterations - maximum number of iterations.
   max-feat - maximum value for one feature."
  [fitness ftype dim nplants nplants-max seed-min seed-max
   sigma-init sigma-final modulation max-iterations max-feat]
  (let [population (init-population nplants dim seed-min seed-max sigma-init
                                    sigma-final max-feat)]
    ;; one time eval of initial plants
    (dorun (map #(send % eval-seed fitness) (:plantlist @population)))
    (apply await (:plantlist @population))
    (send population set-bestworst ftype)
    (await population)
    ; start IWO
    (dorun (map (fn [i]
                  (grow fitness ftype population max-iterations nplants-max
                        modulation i))
                (range max-iterations)))
    ; return best solution
    @(:gbest @population)))

;;;;PARTICLE SWARM OPTIMIZATION;;;
;;; http://en.wikipedia.org/wiki/Particle_swarm_optimization
;;;1. port my old .NET Version
;;;2. JSwarm version is in externals ns 
;;;   -- need to write macro that calls gen-class or something
;;;(defn hierarchical-particle-swarm ;alternatively, inertia could instead 
;;;decrease from 0.9 to 0.4 as algorithm progresses
;;;  ""
;;;  [f lowers uppers lazy-rnd iter & {:keys [guess inertia cognitive 
;;;var-socials level-particles] :or {inertia 0.729, cognitive 1.49445, 
;;;var-socials [1.49445], level-particles [10]}}]
;;;  ())

;;replace these with maps
(defstruct particle :position :velocity :pbestpos :pbestfit)
(defstruct swarm :particlelist :gbest :vmax :hfit)

(defn- init-particle
  [nDimensions max]
  (agent (struct particle
                 (double-array (for [i (range nDimensions)] (rand max)))
                 (double-array (for [i (range nDimensions)] (rand max)))
                 (double-array (replicate nDimensions 1.0))
                 m/max-dbl)))

(defn- init-swarm
  [nParticles nDimensions vmaxDelta nHistoryFitness max]
  (let [particles (map (fn [_] (init-particle nDimensions max))
                       (range nParticles))
        best (init-particle nDimensions max)
        vmax (agent (double-array (repeat nDimensions (* vmaxDelta max))))
        hfit (agent (double-array (replicate nHistoryFitness m/max-dbl)))]
    (struct swarm particles best vmax hfit)))

(defn- inertia-weight
  []
  (double (+ 0.5 (/ (rand) 2.0))))

(defn- update-particle
  [particle fitness ftype gbest vmaxVec iteration]
  (let [^doubles pos (:position particle)
        ^doubles vel (:velocity particle)
        ^doubles pbestpos (:pbestpos particle)
        pbestfit (double (:pbestfit particle))
        ^doubles gbestpos (:pbestpos @gbest)
        gbestfit (double (:pbestfit @gbest))

        ; acceleration constants
        c1 1.494
        c2 1.494
        ; random elements
        phi1 (double (rand))
        phi2 (double (rand))
        ; inertia weight
        w (inertia-weight)

        ; update velocity
        ^doubles newvel (amap vel idx ret
                              (+ (* c1 phi1 (- (aget gbestpos idx) (aget pos idx)))
                                 (* c2 phi2 (- (aget pbestpos idx) (aget pos idx)))
                                 (* (aget vel idx))))
        ;velocity clamping
        ^doubles vmax @vmaxVec
        ^doubles clampvel (amap newvel i ret
                                (if (< (m/abs (aget newvel i)) (aget vmax i))
                                  (aget newvel i)
                                  (* (/ (aget vmax i) (m/abs (aget newvel i)))
                                     (aget newvel i))))

        ; update position
        newpos (amap pos idx ret (+ (aget pos idx) (aget clampvel idx)))
        ; update pbest
        fit (fitness newpos)
        newpbestfit (if (ftype fit pbestfit) fit pbestfit)
        newpbestpos (if (ftype fit pbestfit) newpos pbestpos)]
    ; finally assign updates to agent
    (assoc particle :position newpos :velocity newvel
                    :pbestpos newpbestpos :pbestfit newpbestfit)))

(defn- update-gbest [gbest ftype particles]
  (let [newbest (first (sort-by :pbestfit ftype (map #(deref %) particles)))
        #^doubles newpbestpos (:pbestpos newbest)
        newpbestfit (double (:pbestfit newbest))]
    (assoc gbest :pbestpos newpbestpos
                 :pbestfit newpbestfit)))

(defn- update-hfit [^doubles hfit gbest]
  (let [gbestfit (double (:pbestfit @gbest))]
    (amap hfit i ret
          (if (= i 0)
            gbestfit
            (aget hfit (- i 1))))))

(defn- update-vmax [^doubles vmax swarm iteration]
  (let [#^doubles nlastfit @(:hfit swarm)
        gbestfit (double (:pbestfit @(:gbest swarm)))
        prebeta (- 1.0 (* iteration 0.0001))
        beta (if (<= prebeta 0.0001) 0.0001 prebeta)]
    (amap vmax i ret
          (if (every? (fn [x]
                        (>= gbestfit x))
                      nlastfit)
            (* beta (aget vmax i))
            (aget vmax i)))))

(defn- reset-vmax [vmax vmaxDelta max-feat]
  (double-array (repeat (count vmax) (* vmaxDelta max-feat))))

(defn- reset-particle [particle]
  (let [nDimensions (count (:position particle))
        newpos (double-array (for [i (range nDimensions)] (rand 255)))
        ;;newvel (double-array (repeat nDimensions 0.0))
        newvel (double-array (for [i (range nDimensions)] (rand 255)))
        #^doubles newpbest (:pbestpos particle)]
    (assoc particle :position newpos :velocity newvel :pbestpos newpbest)))

(defn- fly [fitness ftype swarm iteration max-feat]
  ;;reset x random particles every 4 iteration  
  ;; (if (= (rem iteration 25) 0)
  ;;   (do
  ;;     (send (:vmax swarm) reset-vmax 1.0 max-feat)
  ;;     (await (:vmax swarm))
  ;;     (let [plist (map #(nth (:particlelist swarm) %)
  ;;                        (take 15 
  ;;      (scramble (range 0 (count (:particlelist swarm))))))]
  ;;         (dorun (map #(send %1 reset-particle) plist))
  ;;         (apply await plist))))

  ; update particles
  (dorun (map #(send % update-particle fitness ftype (:gbest swarm)
                     (:vmax swarm) iteration)
              (:particlelist swarm)))
  (apply await (:particlelist swarm))

  ; update gbest
  (send (:gbest swarm) update-gbest ftype (:particlelist swarm))
  (await (:gbest swarm))
  (print (:pbestfit @(:gbest swarm)) " ")

  ;update  history fitness
  (send (:hfit swarm) update-hfit (:gbest swarm))
  (await (:hfit swarm))

  ;update vmax
  (send (:vmax swarm) update-vmax swarm iteration)
  (await (:vmax swarm))
  (println (aget ^doubles @(:vmax swarm) 0) " "))

(defn pso
  "Starts the PSO algorithm.
  Parameter: 
  fitness - defines the fitness function used in the PSO. 
              The only argument is the position of the 
              particle (a double-vector)
  ftype - defines the type of optimization 
             problem (minimize (<) or maximize (>)).
  dim - number of dimensions in solution.
  nparticles - number of particles in swarm.
  vc-init - initial value for velocity clamping (usually 1.0)
  vc-hist - number of history value for velocity 
                   clamping (usually 0.8 * max-iterations).
  max-iterations - maximum number of iterations.
  max-feat - maximum value for one feature."
  [fitness ftype dim nparticles vc-init vc-hist max-iterations max-feat]
  (let [swarm (init-swarm nparticles dim vc-init vc-hist max-feat)]
    (dorun (map (fn [i]
                  (fly fitness ftype swarm i max-feat))
                (range max-iterations)))
    @(:gbest swarm)))
;;;;;;;;GA and ES share a ton of code...

(defn- euclidean [^doubles v1 ^doubles v2]
  ;;use arrays ns here....
  (m/sqrt (areduce v1 i ret 0.0 (+ ret (m/sq (- (aget v1 i) (aget v2 i)))))))

(defstruct individual :tag :chromosome :steps :fitness)

(defstruct population :poplist)

;;;special GA
(defn- share
  [dist sigma alpha]
  (if (<= dist sigma)
    (- 1 (m/pow (/ dist sigma) alpha))
    0))

(defn- fitness-sharing
  [ind popu]
  (let [sigma 100                                           ;; 100 seems to be fine
        alpha 1
        dist-sum (reduce +
                         (for [other (:poplist popu)]
                           (share (euclidean
                                    (double-array (:chromosome ind))
                                    (double-array (:chromosome other)))
                                  sigma
                                  alpha)))]
    (assoc ind :fitness (/ (:fitness ind) dist-sum))))

;;;shared
(defn- rand-int-es
  [max number]
  (map (fn [_] (* max (rand)))
       (range number)))                                     ;need lazy-rnd

(defn- rand-int-ga
  [max number]
  (map (fn [_] (int (* max (rand))))
       (range number)))

(defn- generate-chromosome-es ^doubles
[^long n ^long bits]
  (let [max (int (- (m/pow 2 bits) 1))]
    (double-array (rand-int-es max n))))

(defn- generate-chromosome-ga ^ints
[^long n ^long bits]
  (let [max (int (- (m/pow 2 bits) 1))]
    (int-array (rand-int-ga max n))))

(defn- init-individual-es
  [^long n ^long bits]
  {:tag        0
   :chromosome (generate-chromosome-es n bits)
   :steps      (double-array (map #(* 10 %)
                                  (take 22 (repeatedly rand))))
   :fitness    0})

(defn- init-individual-ga
  [n bits]
  (struct individual 0 (generate-chromosome-ga n bits)
          (map #(* 10 %)
               (take 22 (repeatedly rand))) 0))

(defn- init-population-es [n dim bits]
  (let [poplist (map (fn [_]
                       (init-individual-es dim bits))
                     (range n))]
    {:poplist poplist}))

(defn- init-population-ga [n dim bits]
  (let [poplist (map (fn [_]
                       (init-individual-ga dim bits))
                     (range n))]
    {:poplist poplist}))

(defn chromo-to-phenotype-es
  [^doubles chromosome]
  (let [sum (areduce chromosome i ret 0.0 (+ ret (aget chromosome i)))]
    (amap chromosome i ret (/ (aget chromosome i) sum))))

(defn chromo-to-phenotype-ga
  [chromosome]
  (let [sum (reduce + chromosome)]
    (double-array (for [gene chromosome] (/ gene sum)))))

(defn- evaluate-individual-es
  [ind fitness]
  (let [w-pheno (chromo-to-phenotype-es (:chromosome ind))
        fitval (fitness w-pheno)]
    (assoc ind :fitness fitval)))

(defn- evaluate-individual-ga
  [ind fitness]
  (let [w-pheno (chromo-to-phenotype-ga (:chromosome ind))
        fitval (fitness w-pheno)]
    (assoc ind :fitness fitval)))

(defn- evaluate-all-es
  [popu fitness]
  (if DEBUG? (println "Evaluating population..."))
  (let [agentlist (for [ind (:poplist popu) :when (= (:tag ind) 1)]
                    (agent ind))]
    (dorun (map #(send %1 evaluate-individual-es fitness)
                agentlist))
    (apply await agentlist)
    (let [children (for [agent agentlist]
                     @agent)
          parents (for [p (:poplist popu) :when (= (:tag p) 0)]
                    p)]
      (assoc popu :poplist (concat parents children)))))

(defn- evaluate-all-ga
  [popu fitness fs?]
  (let [agentlist (for [ind (:poplist popu) :when (= (:tag ind) 1)]
                    (agent ind))]
    (dorun (map #(send %1 evaluate-individual-ga fitness)
                agentlist))
    (apply await agentlist)
    ;; fitness sharing
    (if fs? (do
              (dorun (map #(send %1 fitness-sharing popu)
                          agentlist))
              (apply await agentlist)))
    (let [children (for [agent agentlist]
                     @agent)
          parents (for [p (:poplist popu) :when (= (:tag p) 0)]
                    p)]
      (assoc popu :poplist (concat parents children)))))

(defn- evaluate-all-firstrun-es
  [popu fitness]
  (let [agentlist (for [ind (:poplist popu)]
                    (agent ind))]
    (dorun (map #(send %1 evaluate-individual-es fitness)
                agentlist))
    (apply await agentlist)
    (assoc popu :poplist (for [agent agentlist]
                           @agent))))
(defn- evaluate-all-firstrun-ga
  [popu fitness]
  (let [agentlist (for [ind (:poplist popu)]
                    (agent ind))]
    (dorun (map #(send %1 evaluate-individual-ga fitness)
                agentlist))
    (apply await agentlist)
    (assoc popu :poplist (for [agent agentlist]
                           @agent))))

;;;;;Evolution Strategies (ES)

(defn- crossover-inter
  "Intermediate recombination"
  [ind1 ind2]
  (let [dim (count (:chromosome ind1))
        ^doubles p1ch (:chromosome ind1)
        ^doubles p1st (:steps ind1)
        ^doubles p2ch (:chromosome ind2)
        ^doubles p2st (:steps ind2)
        childc (amap p1ch i ret (/ (+ (aget p1ch i) (aget p2ch i)) 2.0))
        childs (amap p1st i ret (/ (+ (aget p1st i) (aget p2st i)) 2.0))]
    (list childc childs)))

(defn- crossover-discrete
  "Discrete recombination"
  [ind1 ind2]
  (let [dim (count (:chromosome ind1)),
        ^doubles p1ch (:chromosome ind1),
        ^doubles p1st (:steps ind1),
        ^doubles p2ch (:chromosome ind2),
        ^doubles p2st (:steps ind2),
        childc (amap p1ch i ret (if (= (int (* 2 (rand))) 0) (aget p1ch i)
                                                             (aget p2ch i))),
        childs (amap p1st i ret (if (= (int (* 2 (rand))) 0) (aget p1st i)
                                                             (aget p2st i)))]
    (list childc childs)))

(defn- adapted-mutation-es
  "Uncorrelated Mutation with n Step Sizes"
  [^doubles chromosome ^doubles steps]
  (let [dim (count chromosome)
        bound 0.5
        tau (/ 1 (m/sqrt (* 2 (m/sqrt dim))))
        tauprime-rand (* (/ 1 (m/sqrt (* 2 dim)))
                         (normal 0 1))
        ^doubles new-steps-raw (amap steps i ret
                                     (* (aget steps i)
                                        (m/exp m/E (+ tauprime-rand
                                                      (* tau (normal 0 1))))))
        ^doubles new-steps (amap new-steps-raw i ret
                                 (if (< (aget new-steps-raw i) bound)
                                   bound
                                   (aget new-steps-raw i)))
        new-pos (amap chromosome i ret
                      (+ (aget chromosome i)
                         (* (aget new-steps i)
                            (normal 0 1))))]
    (list new-pos new-steps)))

(defn- do-offspring-es
  [acc parents-pair]
  (let [cross (crossover-inter (first parents-pair) (second parents-pair))
        mutated (adapted-mutation-es (first cross) (second cross))
        child {:tag        1
               :chromosome (first mutated)
               :steps      (second mutated)
               :fitness    0}
        old-pop-list (:poplist acc)]
    (assoc acc :poplist (conj old-pop-list child))))

(defn- generate-offspring-es
  [acc parents]
  (when DEBUG? (println "Generating Offspring..."))
  (reduce do-offspring-es acc parents))

(defn- survivor-selection-es
  [popu ftype popsize]
  (when DEBUG? (println "Survivor Selection..."))
  (let [survivors (take popsize (sort-by :fitness ftype (for [i (:poplist popu) :when (= (:tag i) 1)]
                                                          i)))]
    (assoc popu :poplist survivors)))

(defn- parent-selection-es
  [popu popsize factor]
  (if DEBUG? (println "Parent Selection..."))
  (let [parents (scramble (take (* factor popsize) (cycle (:poplist popu))))
        splitted (split-at (* (/ factor 2) popsize) parents)]
    (map #(list %1 %2)
         (nth splitted 0)
         (nth splitted 1))))

;;(def foo (init-population 10 22 8))
;;(count (parent-selection foo 0.4 10))
(defn es
  "Starts the ES algorithm.
  Parameter: 
  fitness - defines the fitness function used. 
              The only argument is the position of the 
               particle (a double-vector).
  ftype - defines the type of optimization 
               problem (minimize (<) or maximize (>)).
  dim - number of dimensions in solution.
  popsize - size of the population.
  offsp-factor - factor for the number of offspring created (usually 2).
  max-iterations - maximum number of iterations."
  [fitness ftype dim popsize offsp-factor max-iterations]
  (let [popu (init-population popsize dim 8),
        popu-evaluated (evaluate-all-firstrun-es popu fitness)]
    (loop [runs max-iterations popu popu-evaluated]
      (if DEBUG? (println "best in generation:"
                          (:fitness (first (sort-by :fitness ftype
                                                    (:poplist popu))))))
      (if (zero? runs)
        (first (sort-by :fitness ftype (:poplist popu)))
        (let [parents (parent-selection-es popu popsize offsp-factor)
              popu-with-children (generate-offspring-es popu parents)
              popu-eva (evaluate-all-es popu-with-children fitness)
              popu-surv (survivor-selection-es popu-eva ftype popsize)
              new-pop (assoc popu-surv :poplist (for [i (:poplist popu-surv)]
                                                  (assoc i :tag 0)))]
          (if DEBUG? (do (println "no. parent pairs:" (count parents))
                         (println "no. parents + children:"
                                  (count (:poplist popu-with-children)))
                         (println "no. survivors: "
                                  (count (:poplist popu-surv)))))
          (recur (dec runs) new-pop))))))

;;;GENETIC ALGORITHM
(defn- gene-crossover
  [^long gene1 ^long gene2 ^long len]
  (let [mask (int (- (m/pow 2 len) 1))
        last-g1 (bit-and gene1 mask)
        last-g2 (bit-and gene2 mask)
        new-g1 (bit-or (bit-shift-left
                         (bit-shift-right gene1 len) len) last-g2)
        new-g2 (bit-or (bit-shift-left
                         (bit-shift-right gene2 len) len) last-g1)]
    (list new-g1 new-g2)))

;;(crossover '(255 12 18 238 210 199 88) '(4 8 15 16 23 42 108) 8)
(defn- crossover
  [chromo1 chromo2 ^long bits]
  (let [split-pos (nth (range 1 bits) (rand-int (- bits 1)))
        split-length (- bits split-pos)
        result (map #(gene-crossover %1 %2 split-length)
                    chromo1
                    chromo2)
        c1-new (map #(first %1)
                    result)
        c2-new (map #(second %1)
                    result)]
    (list c1-new c2-new)))

;;(mutation '(34 21 57 56 11 10) 8 0.25) 
;;-- notice percentage is replaced by 0.25 below
(defn- mutation
  [chromosome ^long bits ^double percentage]
  (for [gene chromosome]
    (loop [g gene b 0]
      (if (= b bits)
        g
        (let [new-g (if (< (rand) 0.25)
                      (bit-flip g b)
                      g)]
          (recur new-g (inc b)))))))

(defn- do-offspring-ga
  [acc parents-pair]
  (let [percentage 0.3
        p1c (:chromosome (first parents-pair))
        p2c (:chromosome (second parents-pair))
        children (crossover p1c p2c 8)
        child1 {:tag        1
                :chromosome (mutation (first children) 8 percentage)
                :steps      (list)
                :fitness    0}
        child2 {:tag        1
                :chromosome (mutation (second children) 8 percentage)
                :steps      (list)
                :fitness    0}
        old-pop-list (:poplist acc)]
    (assoc acc :poplist (conj old-pop-list child1 child2))))

(defn- adapted-crossover
  [p1 p2]
  (let [dim (count (:chromosome p1))
        dim2 (/ dim 2)
        p1ch (:chromosome p1)
        p1st (:steps p1)
        p2ch (:chromosome p2)
        p2st (:steps p2)
        childc (map (fn [c1 c2]
                      (if (= (int (* 2 (rand))) 0) c1 c2))
                    p1ch
                    p2ch)
        childs (map (fn [s1 s2]
                      (if (= (int (* 2 (rand))) 0) s1 s2))
                    p1st
                    p2st)]
    (list childc childs)))

(defn- adapted-mutation-ga
  [chromosome steps]
  (let [dim 22
        bound 0.5
        tau (/ 1 (m/sqrt (* 2 (m/sqrt dim))))
        tauprime-rand (* (/ 1 (m/sqrt (* 2 dim))) (normal 0 1))
        new-steps-raw (for [s steps]
                        (* s (m/exp (+ tauprime-rand
                                       (* tau (normal 0 1))))))
        new-steps (map (fn [s]
                         (if (< s bound) bound s))
                       new-steps-raw)
        new-pos (map (fn [c s]
                       (+ c (* s (normal 0 1))))
                     chromosome new-steps)]
    (list new-pos new-steps)))

(defn- do-offspring-adapted
  [acc parents-pair]
  (let [p1 (first parents-pair)
        p2 (second parents-pair)
        crossed (adapted-crossover p1 p2)
        mutated (adapted-mutation-ga (first crossed) (second crossed))
        child {:tag        1
               :chromosome (first mutated)
               :steps      (second mutated)
               :fitness    0}
        old-pop-list (:poplist acc)]
    (assoc acc :poplist (conj old-pop-list child))))

(defn- generate-offspring-ga
  [acc parents adapted?]
  (if adapted?
    (reduce do-offspring-adapted acc parents)
    (reduce do-offspring-ga acc parents)))

(defn- survivor-selection-ga
  [popu ftype percent popsize]
  (let [sorted-all (sort-by :fitness ftype (:poplist popu))
        size-all (count sorted-all)
        size-top (* percent size-all)
        size-rest (- popsize size-top)
        splitted (split-at size-top sorted-all)
        top-percent-members (nth splitted 0)
        rest-offspring (take size-rest (sort-by :fitness ftype
                                                (for [i (nth splitted 1)
                                                      :when (if (= (:tag i) 1) true)] i)))]
    (assoc popu :poplist (concat top-percent-members rest-offspring))))

(defn- parent-selection-ga
  [popu ftype percent popsize adapted?]
  (let [sorted (sort-by :fitness ftype (:poplist popu))
        size-parents (int (* percent popsize))
        parents (take size-parents sorted)
        ext-parents (if adapted?
                      (scramble (take (* 2 popsize) (cycle parents)))
                      (scramble (take popsize (cycle parents))))
        splitted (if adapted?
                   (split-at popsize ext-parents)
                   (split-at (/ popsize 2) ext-parents))]
    (map #(list %1 %2)
         (nth splitted 0)
         (nth splitted 1))))

(defn ga
  "Starts the GA algorithm.
  Parameter: 
  fitness - defines the fitness function used. 
                 The only argument is the position of the 
                  particle (a double-vector).
  ftype - defines the type of optimization 
                problem (minimize (<) or maximize (>)).
  dim - number of dimensions in solution.
  popsize - size of the population.
  par-perc - percentage of parents chosen for parent selection.
  surv-perc - percentage of population chosen for next genertion.
  adpated? - use mutation step size adaption (boolean).
  max-iterations - maximum number of iterations."
  [fitness ftype dim popsize par-perc surv-perc adapted? max-iterations]
  (let [popu (init-population popsize 22 8)
        popu-evaluated (evaluate-all-firstrun-ga popu fitness)]
    (loop [runs max-iterations popu popu-evaluated]
      (when DEBUG?
        (println "best:"
                 (:fitness (first (sort-by :fitness ftype
                                           (:poplist popu))))))
      (if (zero? runs)
        (first (sort-by :fitness ftype (:poplist popu)))
        (let [parents (parent-selection-ga popu ftype par-perc popsize adapted?)
              popu-with-children (generate-offspring-ga popu parents adapted?)
              popu-eva (evaluate-all-ga popu-with-children fitness adapted?)
              popu-surv (survivor-selection-ga popu-eva ftype surv-perc
                                               popsize)
              new-pop (assoc popu-surv :poplist (for [i (:poplist popu-surv)]
                                                  (assoc i :tag 0)))]
          (when DEBUG?
            (do (println "no. parent pairs:" (count parents))
                (println "no. parents + children:"
                         (count (:poplist popu-with-children)))
                (println "no. survivors: "
                         (count (:poplist popu-surv)))))
          (recur (dec runs) new-pop))))))

