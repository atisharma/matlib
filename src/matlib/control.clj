(ns matlib.control
  "Linear algebra operations on matrices."
  (:require
    [matlib.core :refer :all]
    [matlib.linalg :refer :all]
    [uncomplicate.neanderthal.real :refer [entry entry!]]
    [uncomplicate.neanderthal.native :refer :all :exclude [sv]]
    [uncomplicate.neanderthal.linalg :refer :all]
    [uncomplicate.neanderthal.core :refer :all :exclude [entry entry!]]
    [uncomplicate.neanderthal.vect-math :as vect-math]
    [uncomplicate.neanderthal.random :as random]))

(defn obsv
  "Extended observability matrix up to order `i`. If `i` is not specified, up to
  order `n` (the usual controllability matrix)."
  ([A C]
   (obsv A C (mrows A)))
  ([A C i]
   (let [n (mrows A)
         m (mrows C)
         Gamma_i (dge (* m i) n)
         Apow (dge-eye n)]
     (doseq [block-row (range i)]
       (copy! (mm C Apow) (submatrix Gamma_i (* block-row m) 0 m n))
       (copy! (mm Apow A) Apow))  ; Apow = A^k, k=0...(i-1)
     Gamma_i)))

(defn ctrb
  "Controllability matrix."
  [A B]
  :not-implemented)

(defn sylv
  "Solution to the Sylvester equation  
  `AX + XB = C`."
  [A B C]
  :not-implemented)

(defn dgram
  "Discrete-time Gramian."
  [A B]
  :not-implemented)

(defn cgram
  "Continuous-time Gramian."
  [A B]
  :not-implemented)

(defn dlyap
  "Solution `X` to the discrete-time Lyapunov equation  
  `AXA' - X + Q = 0`."
  [A Q]
  :not-implemented)

(defn clyap
  "Solution `X` to the continuous-time Lyapunov equation  
  `AX + XA' + Q = 0`."
  [A Q]
  :not-implemented)

(defn dare
  "Solution to the discrete-time algebraic Riccati equation."
  [A B Q]
  :not-implemented)

(defn care
  "Solution to the continuous-time algebraic Riccati equation."
  [A B Q]
  :not-implemented)
