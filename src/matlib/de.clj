(ns matlib.de
  "Differential evolution heuristic gradient-free optimisation.
  See, for example, https://en.wikipedia.org/wiki/Differential_evolution ."
  (:require
      [matlib.core :as ml-core]
      [matlib.linalg :as ml-linalg]
      [uncomplicate.neanderthal.real :refer [entry entry!]]
      [uncomplicate.neanderthal.native :as native]
      [uncomplicate.neanderthal.linalg :as n-linalg]
      [uncomplicate.neanderthal.core :as n-core]
      [uncomplicate.neanderthal.vect-math :as vect-math]
      [uncomplicate.neanderthal.random :as random]))  

(defn perturb
  "Perturb a state vector `x`."
  [x]
  (let [x0 (native/dv x)]
    (while (= 0.0 (reduce * x0)) (n-core/axpy! (random/rand-normal! (native/dv x)) x0))
    (vect-math/mul x0 (random/rand-normal! (native/dv x)))))

(defn population
  "Initialise a population based on an example `x`.
  `x` cannot have zero entries.
  `NP` must be at least 5."
  ([x]
   (let [d (n-core/dim x)
         NP (cond (< d 10) (* 5 d)
                  (> d 100) (int (/ d 2))
                  :else 50)]
     (population x NP)))
  ([x NP]
   (for [n (range NP)] (perturb x))))

(defn update-individual
  "Returns updated `x`."
  [x a b c CR F]
  (let [R (rand-int (n-core/dim x))
        y (n-core/copy x)]
    (doseq
      [i (range (n-core/dim x))]
      (if (or (< (rand) CR) (= i R))
        (entry! y i (+ (entry a i) (* F (- (entry b i) (entry c i)))))))
    y))

(defn solve
  "Find the minimum of `f` given a population (vec) of `x`s.
  Constraints should be handled in `f`.
  If evaluations of `f` are expensive, consider memoizing `f` with lru."
  ([f xs CR F maxiter scores n]
   (let [x (first xs)
         ys (shuffle (rest xs))
         a (first ys)
         b (nth ys 2)
         c (nth ys 3)]
     ;(print n scores "\r")
     (if (and (> (- (reduce max scores) (reduce min scores)) ml-core/sq-eps) (< n maxiter))
       ; rotate xs with new one at end
       (let [y (update-individual x a b c CR F)
             fx (f x)
             fy (f y)
             new-xs (conj (vec (rest xs)) (if (or (> fy fx) (Double/isNaN fy)) x y))
             new-scores (conj (vec (rest scores)) (min fx fy))]
         (recur f new-xs CR F maxiter new-scores (inc n)))
       {:sol x
        :f (f x)
        :iterations n
        :maxiter maxiter
        :CR CR
        :F F
        :NP (count xs)
        :success (< n maxiter)}))))

