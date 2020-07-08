(ns matlib.state-space
  "Basic Linear Time-Invariant (LTI) state space system operations and patterns.
  
  A discrete-time state-space system is expressed as a set of matrices
  `A`, `B`, `C`, `D` which can be integrated as:  
  `x(k+1) = A x(k) + B u(k) + w(k)`  
  `  y(k) = C x(k) + D u(k) + v(k)`,  
  with integer `k` and  
  `u(k) ∈ ℝ^m,`  
  `y(k) ∈ ℝ^l,`  
  `x(k) ∈ ℝ^n`.  

  The noise terms `w(k)` and `v(k)` are unobserved, Gaussian distributed,
  zero-mean, non-zero white noise, with covariances  
  `ℰ ( w_k w_l' ) = Q  δ_kl >= 0`  
  `ℰ ( w_k v_l' ) = S  δ_kl >= 0`  
  `ℰ ( v_k w_l' ) = S' δ_kl >= 0`  
  `ℰ ( v_k v_l' ) = R  δ_kl >= 0`  

  A continuous-time state-space system is a set of matrices `A`, `B`, `C`, `D`
  which can be integrated as:
  
  `dx(t)/dt = A x(t) + B u(t) + w(t)`  
      `y(t) = C x(t) + D u(t) + v(t)`,  
  with `t ∈ ℝ` and  
  `u(t) ∈ ℝ^m,`  
  `y(t) ∈ ℝ^l,`  
  `x(t) ∈ ℝ^n`.

  Internally, discrete-time state-space systems are stored as maps,
  `{:A A :B B :C C :D D}`,  
  optionally with  
  `{:x x(k), :u u(k), :y y(k), :x+ x(k+1)}`  
  and optionally with `{:E E}`.  
  If `E` is supplied, `w(k)` and `v(k)` are generated with covariance matrix  
  `[ Q  S ] = E'E`  
  `[ S' R ]`  
  so that
  `( w(k) ) = E n(k)`  
  `( v(k) )`  
  with `n(k)` serially uncorrelated and drawn from `N(0,1)`.

  The snapshot matrix of inputs `{u(k)}` is stored under `:U`.
  Once the system is integrated using `make-snapshots`, snapshot matrix `Y`
  will be added to the map.
  "
  (:require
    [matlib.core :refer :all]
    [matlib.linalg :refer :all]
    [uncomplicate.neanderthal.real :refer [entry entry!]]
    [uncomplicate.neanderthal.native :refer :all :exclude [sv]]
    [uncomplicate.neanderthal.linalg :refer :all]
    [uncomplicate.neanderthal.core :refer :all :exclude [entry entry!]]
    [uncomplicate.neanderthal.vect-math :refer :all]
    [uncomplicate.neanderthal.random :as random]))

(defn step-discrete-time
  "Integrate a discrete-time system one step, by calculating  
  `x+ = Ax + Bu + w`  
  ` y = Cx + Du + v`  
  `x` is the state vector (a column vector).
  Argument `u` is u(k)."
  ([ss u]
   (let [x  (get ss :x+ (:x ss))
         y  (axpy
             1 (mm (:C ss) x)
             1 (mm (:D ss) (col-vector u))) 
         x+ (axpy
             1 (mm (:A ss) x)
             1 (mm (:B ss) (col-vector u)))]
    (merge ss
     {:k  (inc (get ss :k -1))
      :x  x
      :x+ x+
      :y  y
      :u  (col-vector u)}
     (if (:E ss)
       (let [n  (mm (:E ss) (random/rand-normal! (dge (ncols (:E ss)) 1)))
             w  (submatrix n (dim x) 1)
             v  (submatrix n (dim x) 0 (dim y) 1)]
         {:x+ (axpy x+ w)
          :y  (axpy y  v)
          :v  v
          :w  w})
      {})))))

(defn make-snapshots
  "Generate snapshot matrices of states `X`, outputs `Y` and noise `e` from a
  snapshot matrix of inputs, `U`, and a model, `ss`.
  `:k` of `ss` corresponds to the column of `U` and `Y`."
  ([ss]
   (cond
     (not (:U ss)) :input-snapshots-missing
     (not (:Y ss)) (recur (merge ss {:Y (dge (mrows (:C ss)) (ncols (:U ss)))}))
     (not (:X ss)) (recur (merge ss {:X (dge (mrows (:A ss)) (ncols (:U ss)))}))
     (not (:W ss)) (recur (merge ss {:W (dge (mrows (:A ss)) (ncols (:U ss)))}))
     (not (:V ss)) (recur (merge ss {:V (dge (mrows (:C ss)) (ncols (:U ss)))}))
     (not (:k ss)) (recur (step-discrete-time ss (col (:U ss) 0)))
     (< (:k ss) (ncols (:U ss))) (let [ss+ (step-discrete-time ss (col (:U ss) (:k ss)))]
                                   ; side effects on X, Y
                                   (copy! (:x ss) (submatrix (:X ss) 0 (:k ss) (mrows (:A ss)) 1))
                                   (copy! (:w ss) (submatrix (:W ss) 0 (:k ss) (mrows (:A ss)) 1))
                                   (copy! (:y ss) (submatrix (:Y ss) 0 (:k ss) (mrows (:C ss)) 1))
                                   (copy! (:v ss) (submatrix (:V ss) 0 (:k ss) (mrows (:C ss)) 1))
                                   (recur ss+))
     :else ss)))

(defn tf
  "Transfer function (input-output) of system,  
  `G(z) = D + C (zI - A)^-1 B`, or  
  `G(s) = D + C (sI - A)^-1 B`."
  ([ss z]
   (let [A (:A ss)
         B (:B ss)
         C (:C ss)
         D (:D ss)
         n (ncols A)
         I (transfer! (eye n) (dge n n))]
    (axpy 1 (mm C (minv (axpby! z I -1 (copy A))) B) D))))

(defn sigmas
  "Singular values of the transfer function over a range of `z` (or `s`)."
  [ss zs]
  (apply hcat (map #(col-vector (:sigma (svd (tf ss %)))) zs)))
  
;;; below is for testing

(def i 8000)

(def ss-model {:A (dge 2 2 [0.9 0.2 0 -0.995])
                ;:B (dge 3 2 [1 0 1, 0 1 1])
                :B (dge 2 2 [1 2, 0 1])
                :C (dge 1 2 [1 1, 0 1])
                :D (dge 1 2 [0 0, 0 0])
                :E (scal! 0.005 (dge 3 2 (range)))
                :x (dge 2 1 [0 0])
                ; need a better input:
                ; - persistently exciting
                ; - different frequencies
                :U (dge 2 i)})

; first input signal
(axpy! (sin (scal! 0.050 (dge 1 i (range)))) (submatrix (:U ss-model) 0 0 1 i))
(axpy! (sin (scal! 0.055 (dge 1 i (range)))) (submatrix (:U ss-model) 0 0 1 i))
(axpy! (sin (scal! 0.150 (dge 1 i (range)))) (submatrix (:U ss-model) 0 0 1 i))
(axpy! (sin (scal! 1.150 (dge 1 i (range)))) (submatrix (:U ss-model) 0 0 1 i))
(axpy! (sin (scal! 0.002 (dge 1 i (range)))) (submatrix (:U ss-model) 0 0 1 i))
                                                                                      
; second input signal                                                                 
(axpy! (sin (scal! 0.500 (dge 1 i (range)))) (submatrix (:U ss-model) 1 0 1 i))
(axpy! (sin (scal! 0.530 (dge 1 i (range)))) (submatrix (:U ss-model) 1 0 1 i))
(axpy! (sin (scal! 1.503 (dge 1 i (range)))) (submatrix (:U ss-model) 1 0 1 i))
(axpy! (sin (scal! 4.503 (dge 1 i (range)))) (submatrix (:U ss-model) 1 0 1 i))
(axpy! (sin (scal! 0.007 (dge 1 i (range)))) (submatrix (:U ss-model) 1 0 1 i))

; arbitrarily step some input
(scal! 5.0 (submatrix (:U ss-model) 0 200 2 1000))

; input noise
(axpy! (random/rand-normal! 0 0.1 (dge (mrows (:U ss-model)) (ncols (:U ss-model)))) (:U ss-model))
