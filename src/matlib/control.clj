(ns matlib.control
  "Linear algebra operations on matrices."
  (:require
    [matlib.core :refer :all]
    [matlib.linalg :refer :all]
    [matlib.optim :refer [l-bfgs]]
    [uncomplicate.commons.core :refer [release with-release let-release]]
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
  "Extended controllability matrix up to order `i`. If `i` is not specified, up to
  order `n` (the usual controllability matrix)."
  ([A B]
   (ctrb A B (mrows A)))
  ([A B i]
   (let [n (mrows A)
         m (ncols B)
         Ctrb (dge n (* m i))
         Apow (dge-eye n)]
     (doseq [block-col (range i)]
       (copy! (mm B Apow) (submatrix Ctrb 0 (* block-col n) n m))
       (copy! (mm Apow A) Apow))  ; Apow = A^k, k=0...(i-1)
     Ctrb)))

(defn sylv
  "Find a solution to the Sylvester equation  
  `AX + XB = C`."
  ([x A B C]
   (let [X (view-ge x (ncols A) (mrows B))]
     (nrm2 (axpy -1 C (mm A X) (mm X B)))))
  ([A B C]
   (let [x0 (view-vctr (dge (ncols A) (mrows B)))
         x (:sol (l-bfgs #(sylv % A B C) x0))
         X (view-ge x (ncols A) (mrows B))]
     X)))

(defn dlyap
  "Solution `X` to the discrete-time Lyapunov equation  
  `AXA' - X + Q = 0`.  
  Should be solved by Schur method, but uses L-BFGS for now."
  ([x A Q]
   (let [X (view-ge x (mrows Q) (ncols Q))]
     (nrm2 (axpy -1 X (mm A X (trans A)) Q))))
  ([A Q]
   (let [x0 (view-vctr Q)
         x (:sol (l-bfgs #(dlyap % A Q) x0))
         X (view-ge x (mrows Q) (ncols Q))]
     (axpy 0.5 X 0.5 (trans X)))))

(defn clyap
  "Solution `X` to the continuous-time Lyapunov equation  
  `AX + XA' + Q = 0`.  
  Should be solved by Bartels-Stewart method, but uses L-BFGS for now."
  ([x A Q]
   (let [X (view-ge x (mrows Q) (ncols Q))]
     (nrm2 (xpy (mm A X) (mm X (trans A)) Q))))
  ([A Q]
   (let [x0 (view-vctr Q)
         x (:sol (l-bfgs #(clyap % A Q) x0))
         X (view-ge x (mrows Q) (ncols Q))]
     (axpy 0.5 X 0.5 (trans X)))))

(defn dcgram
  "Discrete-time controllability Gramian."
  [A C]
  (dlyap A (mm (trans C) C)))

(defn dogram
  "Discrete-time observability Gramian."
  [A B]
  (scal! -1 (dlyap A (scal! -1 (mm B (trans B))))))

(defn ccgram
  "Solution `Wc` to the continuous-time Lyapunov equation  
  `A Wc + Wc A' + B B' = 0`."
  ([A B]
   (clyap A (mm B (trans B)))))

(defn cogram
  "Solution `Wo` to the continuous-time Lyapunov equation  
  `A' Wo + Wo A + C' C = 0`."
  ([A C]
   (clyap (trans A) (scal -1 (mm (trans C) C)))))

(defn dare
  "Solution to the discrete-time algebraic Riccati equation,  
  `P = A' P A - (A' P B) (R + B' P B)^-1 (B' P A) + Q`.  
  Also returns the stabilising controller  
  `K = (R + B' P B)^-1 B' P A`  
  and the closed-loop state transfer matrix `A-BK`.  
  The solutions is found via the eigenvalue decomposition of the symplectic matrix  
  `[A + B R^-1 B' (A^-1)' Q   |  -B R^-1 B' (A^-1)'  ] = Z  
   [       -(A^-1)' Q         |  (A^-1)'             ]`  
  TODO: extend to generalised DARE.  
  TODO: anti-stabilising solution.
  "
  ([A B Q R]
   (let [n (mrows A)
         Ainv' (trans (minv A))
         B' (trans B)
         Rinv (minv R)
         Z11  (xpy A (mm B Rinv B' Ainv' Q))
         Z12  (scal -1 (mm B Rinv B' Ainv'))
         Z21  (scal -1 (mm Ainv' Q))
         Z22  Ainv'
         Z (vcat (hcat Z11 Z12) (hcat Z21 Z22)) ; symplectic matrix
         U (dge (mrows Z) (ncols Z))  ; Schur vectors
         eigs (dge (mrows Z) 2)
         T (copy Z) ; real Schur form of Z
         _ (es! T eigs U) ; overwrite U, eigs, T
         ; the stable subspace of the Schur form
         U-stab (apply hcat (for [c (range (mrows eigs)) :when (< (entry (vect-math/abs eigs) c 0) 1.0)] (col-vector (col T c))))]
         ; split U into stable and anti-stable subspaces
         ;U1 (submatrix U-stab n n)
         ;U2 (submatrix U-stab 0 n n n)
         ;P (mm U2 (minv U1))]
         ;K (mm (minv (xpy R (mm B' P B))) B' P A)]
     {;:U1 U1
      ;:U2 U2
      :U-stable U-stab
      :U U
      :T T
      ;:K K
      ;:P P
      ;:e (trans (col eigs 0))
      ;:A_CL (axpy -1 (mm B K) A)
      :residual :not-calculated}))) 

(defn care
  "Solution to the continuous-time algebraic Riccati equation."
  [A B Q]
  :not-implemented)
