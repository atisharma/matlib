(ns matlib.ident
  "Identify a state-space model given input and output snapshot matrices,
  using deterministic-stochastic subspace methods (both MOESP and N4SID).
  
  Notation and assumptions follow [vOdM-95].

  Notation:

  `xₖ₊₁ = A xₖ + B uₖ + wₖ`  
  `  yₖ = C xₖ + D uₖ + vₖ`

  `uₖ ∈ ℝᵐ,`  
  `yₖ ∈ ℝˡ,`  
  `xₖ ∈ ℝⁿ,`  
  and `wₖ, vₖ` are unobserved, Gaussian distributed, zero-mean, non-zero white noise.

  They have covariances  
  `E ( wₖ wₗ' ) = Q  δₖₗ>= 0`  
  `E ( wₖ vₗ' ) = S  δₖₗ>= 0`  
  `E ( vₖ wₗ' ) = S' δₖₗ>= 0`  
  `E ( vₖ vₗ' ) = R  δₖₗ>= 0`  
  
  Snapshot matrices `U` and `Y` of `uₖ` and `yₖ` are taken with `k=0...(t-1)`.

  References:

  [V-94]  
  Identification of the deterministic part of MIMO state space models given in
  innovations form from input-output data,  
  M Verhaegen,  
  Automatica, Vol. 30, No. 1, pp. 61-74 (1994)

  [vOdM-94]  
  'N4SID: Subspace Algorithms for the Identification of Combined
  Deterministic-Stochastic Systems'  
  P van Overschee & B de Moor  
  Automatica Vol. 30 No. 1 pp. 75-93 (1993)

  [vOdM-95]  
  A Unifying Theorem for Three Subspace System Identification Algorithms'  
  P van Overschee & B de Moor  
  Automatica, Vol. 31, No. 12, pp. 1853-1861 (1995)

  [vOdM-96]  
  Subspace Identification for Linear Systems  
  Theory — Implementation — Applications  
  P van Overschee, B de Moor  
  Edition 1, Springer US, (1996)  
  ISBN 978-1-4613-8061-0

  [SSvH-04]  
  'High-Performance Numerical Algorithms and Software for Subspace-Based Linear
  Multivariable System Identification'  
  V Sima, DM Sima and S van Huffel  
  J. Comp. Appl. Math., Vol. 170, pp. 371-397 (2004)

  [DSC-06]
  'A New Insight to the Matrices Extraction in a MOESP Type Subspace
  Identification Algorithm'  
  CJM Delgado, P Lopes dos Santos and J L Martins de Carvalho  
  Int. J. Systems Science, Vol. 37, No. 8, pp. 565-574 (2006)

  "
  (:require
   [matlib.core :refer [dge-eye vcat hcat take-rows take-cols]]
   [matlib.control :refer [obsv]]
   [matlib.linalg :refer [rsp rsp-perp oblique-rsp pinv rsvd condition]]
   [uncomplicate.neanderthal
    [core :refer [mm mm! transfer! copy! scal! axpy ncols mrows trans dia view-tr view-vctr view-ge subvector submatrix]]
    [vect-math :refer [sqrt inv]]
    [linalg :as linalg]
    [native :refer [dge dgd]]]))

;;; TODO: follow through algebra for W2 according to [vOdM-96]
;;; TODO: find BDQRS in N4SID-1
;;; TODO: find ABCDQRS in N4SID-robust
;;; TODO: find ABCDQRS in MOESP
;;; TODO: standardise methods to [vOdM-96] presentation

(defn block-hankel
  "Construct a block-Hankel matrix `U_a,b` from snapshot matrix `U`.
  Notation follows [SSvH-04] which corresponds to [vOdM-*]
  as `U_a,b,c` <-> `U_a|b` (the `c` being implicit in the latter).
;
  Arguments:
;
  `a, b, c`   determine size of output
  `U`         `m x j` snapshot matrix `U = [u_0 ... u_j-1]`.
              Row `q` of `U` contains the history of the `q`th input, and each
              column of `U` is a snapshot of inputs at sample time `k`.
;
  Returns:
;
  `U_a,b,c`   block-Hankel matrix of size `(1+b-a)m x (1+c-a)`, where
  `U_a,b,c = [ u_a      u_a+1   u_a+2   ...   u_c     ]
             [ u_a+1    u_a+2   u_a+3   ...   u_c+1   ]
             [ ...                ...                 ]
             [ u_b                ...         u_c+b-a ]`
  "
  ([U a b c]
   (let [m (mrows U)
         bhm-rows (* m (- b a -1))
         bhm-cols (- c a -1)
         bhm (dge bhm-rows bhm-cols)]
     (doseq [k0 (range a (+ b 1))]
       (copy!
        (submatrix U 0 k0 m bhm-cols)
        (submatrix bhm (* (- k0 a) m) 0 m bhm-cols)))
     bhm)))

(defn block-hankel-matrices
  "Construct a block-Hankel matrices from input, output snapshot matrices
  `U` and `Y`,
  `U ∈ ℝ^(m x t)`, `Y ∈ ℝ^(l x t)`.

  `W_p := [ U_0|i-1 ]` is past inputs and outputs,
  `       [ Y_0|i-1 ]`

  `Y_f := Y_i|2i-1` is future outputs,

  `U_f := U_i|2i-1` is future inputs.

  `W_p+`, `U_f-` and `Y_f-` are defined similarly but with the boundary between
  past and future shifted one down (later).
  
  Notation follows [vOdM-95/96]"
  ([ss i]
   (block-hankel-matrices (:U ss) (:Y ss) i))
  ([U Y i]
   (let [t (ncols U)
         N (- t i i)
         i2 (* i 2)]
     {:N N
      :W_p  (vcat (block-hankel U 0 (- i 1) N)
                  (block-hankel Y 0 (- i 1) N))
      :W_p+ (vcat (block-hankel U 0 i       N)
                  (block-hankel Y 0 i       N))
      :Y_f  (block-hankel Y i       (- i2 1) (+ N i))
      :Y_f- (block-hankel Y (+ i 1) (- i2 1) (+ N i 1))
      :U_f  (block-hankel U i       (- i2 1) (+ N i))
      :U_f- (block-hankel U (+ i 1) (- i2 1) (+ N i 1))})))
      
(defn W_2
  "Weighting matrix `W_2` is based on the method being used.
  `W_1` is not calculated because it's only used in CVA, which is not
  implemented here.
  `method` should be one of :MOESP or :N4SID.
  Follows notation in [vOdM-96]"
  ([W_p Y_f U_f method]
   (let [M (rsp-perp W_p U_f)]
     (case method
       :MOESP (mm (pinv M) M) ; = Π_U_f-perp
       :N4SID nil
       :CVA :not-implemented))))

(defn model-spectrum
  "Return the singular values of `(Y_f /_U_f W_p) Π_U_f-perp`.
  The user should determine model order `n` from the largest logarithmic gap in
  the spectrum. See (4.24-4.25) or Fig 4.8 of [vOdM-96]."
  ([ss i]
   (model-spectrum (:U ss) (:Y ss) i))
  ([U Y i]
   (let [{N :N W_p :W_p Y_f :Y_f U_f :U_f} (block-hankel-matrices U Y i)
         W (W_2 W_p Y_f U_f :MOESP)
         O (mm (rsp-perp Y_f U_f) W)] ; weighted oblique projection
     (scal! (/ 1.0 N) (:sigma (linalg/svd (mm O W)))))))

(defn intermediates
  "Intermediate calculations for all methods.
  `W2` should be `nil` if not used."
  ([W_p Y_f U_f W2 l n]
   (let [O_i (oblique-rsp Y_f U_f W_p)      ; eq (4.24)
         {S1 :sigma U1 :u V1' :vt
          S2 :sigma_perp U2 :u_perp V2' :vt_perp} (if W2
                                                    (rsvd (mm O_i W2) :rank n)
                                                    (rsvd O_i :rank n))
         Gamma_i (mm U1 (sqrt S1))
         Gamma_up (submatrix Gamma_i l 0 (- (mrows Gamma_i) l) n)
         Gamma_down (submatrix Gamma_i (- (mrows Gamma_i) l) n)
         Z_i (rsp Y_f (vcat W_p U_f))       ; eq (4.18)
         X_i (mm (sqrt S1) V1')]
     {:O_i O_i
      :Gamma_i Gamma_i
      :Gamma_up Gamma_up
      :Gamma_down Gamma_down
      :order n
      :X_i X_i
      :Z_i Z_i
      :S1 S1 :U1 U1 :V1' V1'})))

(defn Hd_i
  "Block Toeplitz matrix of system matrices."
  ([B D Gamma_i l m n i]
   (let [li (* l i)
         mi (* m i)
         H (dge li mi)
         DGamma_iB (vcat D (mm Gamma_i B))]
      (doseq [stripe (range (- i 1))]
        (let [block-col (take-rows DGamma_iB (- (mrows DGamma_iB) (* l stripe)))]
          (copy! block-col (submatrix H (* stripe l) (* stripe m) (mrows block-col) (ncols block-col)))))
      H)))
  
(defn BD_target
  "Target function that is minimised when `B` and `D` are arguments."
  [BD Gamma_i Xhat_i l m n]
  (let [B (submatrix BD n l)
        D (submatrix BD 0 n m l)
        Z_i (axpy (mm Gamma_i Xhat_i) ((Hd_i A B C D) U_f))]
    (norm2 (mm M_R Z_i))))

(defn find-BD
  "Find `B` and `D` by optimisation method of [DSC-06] for use when `U_f U_f'`
  is poorly conditioned."
  [A C Xhat_i]
  (let [Gamma_i (obsv A C (- i 1))]
    (argmin BD_target)))

(defn n4sid
  "Algorithm 1 (Figure 4.6) of [vOdM-96]."
  ([ss i n]
   (n4sid-1 (:U ss) (:Y ss) i n))
  ([U Y i n]
   (let [l (mrows Y)
         m (mrows U)
         {:keys [N W_p Y_f U_f W_p+ Y_f- U_f-]} (block-hankel-matrices U Y i)
         ;W2 (W_2 W_p Y_f U_f :MOESP) ; untested, check algebra
         W2 nil
         {Gamma_i :Gamma_i, Gamma_down :Gamma_down, Gamma_up :Gamma_up,
          O_i :O_i, Z_i :Z_i} (intermediates W_p Y_f U_f W2 l n)
         {O_i+1 :O_i, Z_i+1 :Z_i} (intermediates W_p+ Y_f- U_f- W2 l n)
         Gamma_i-1 Gamma_down
         Jlu (mm (pinv Gamma_i-1) Z_i+1)
         Jru (mm (pinv Gamma_i) Z_i)
         Y_i|i (block-hankel Y i i (- (ncols Jlu) i))
         U_i|i (block-hankel U i i (- (ncols Jru) i))
         ; solve linear eqns for A, C, K
         ACK (mm (vcat Jlu Y_i|i) (pinv (vcat Jru U_f)))
         A (submatrix ACK 0 0 n n)
         C (submatrix ACK n 0 l n)
         K (submatrix ACK 0 n (+ n l) (- (ncols ACK) n))]
         ; solve for B, D using (4.53), (4.61), (4.62)
         ; determine covariance matrices
     {:A A
      :C C
      :K K
      :i i
      :order n
      :method :N4SID-1})))

(defn n4sid-biased
  "Algorithm 2 (Figure 4.7) of [vOdM-96].
  This algorithm gives asymptotically biased solutions."
  ([ss i n]
   (n4sid-2 (:U ss) (:Y ss) i n))
  ([U Y i n]
   (let [l (mrows Y)
         m (mrows U)
         ; block Hankel matrices, time-shifted and not
         {:keys [N W_p Y_f U_f W_p+ Y_f- U_f-]} (block-hankel-matrices U Y i)
         W2 nil
         ;W2 (W_2 W_p Y_f U_f :MOESP) ; untested, check algebra
         {Gamma_i :Gamma_i, Gamma_down :Gamma_down, Gamma_up :Gamma_up,
          O_i :O_i, Z_i :Z_i} (intermediates W_p Y_f U_f W2 l n)
         {O_i+1 :O_i, Z_i+1 :Z_i} (intermediates W_p+ Y_f- U_f- W2 l n)
         Gamma_i-1 Gamma_down
         X_i+1 (mm (pinv Gamma_i-1) O_i+1)
         X_i (mm (pinv Gamma_i) O_i)
         U_i|i (block-hankel U i i (- (ncols X_i) i))
         Y_i|i (block-hankel Y i i (- (ncols X_i) i))
         ; solve LSQ for A, B, C, D
         Jl (vcat X_i+1 Y_i|i)
         Jr (vcat X_i U_i|i)
         ABCD (mm Jl (pinv Jr))
         ; determine covariance matrices
         samples (ncols X_i)
         residual (axpy -1 (mm ABCD Jr) Jl)
         covariance (scal! (/ 1.0 samples) (mm residual (trans residual)))]
     {:A (submatrix ABCD 0 0 n n)
      :B (submatrix ABCD 0 n n m)
      :C (submatrix ABCD n 0 l n)
      :D (submatrix ABCD n n l m)
      :Q (submatrix covariance 0 0 n n)
      :S (submatrix covariance 0 n n l)
      :R (submatrix covariance n n l l)
      :order n
      :samples samples
      :i i
      :method :N4SID-2})))

(defn robust
  "Explicit calculation, per Figure 4.8 of [vOdM-96].
  This is the so-called 'robust' algorithm.
  The expressions in [DSC-08] are also used.
  `i` is the model delay order and `n` is the model order."
  ([ss i n]
   (robust (:U ss) (:Y ss) i n))
  ([U Y i n]
   (let [l (mrows Y)
         m (mrows U)
         t (ncols U)
         ; block Hankel matrices, time-shifted and not:
         {:keys [N W_p Y_f U_f W_p+ Y_f- U_f-]} (block-hankel-matrices U Y i)
         ; intermediate calculations
         Pi (W_2 W_p Y_f U_f :MOESP)
         {O_i :O_i, Gamma_i :Gamma_i, Gamma_i-1 :Gamma_down,
          Z_i :Z_i
          S1 :S1, U1 :U1, V1' :V1'} (intermediates W_p Y_f U_f Pi l n)
          Z_i+1 (rsp Y_f- (vcat W_p+ U_f-))
          Jlu (mm (pinv Gamma_i-1) Z_i+1)
          Jru (mm (pinv Gamma_i) Z_i)
          Y_i|i (block-hankel Y i i (- (ncols Jlu) i))
          U_i|i (block-hankel U i i (- (ncols Jru) i))
          ; solve linear eqns for A, C, K
          ACK (mm (vcat Jlu Y_i|i) (pinv (vcat Jru U_f)))
          A (submatrix ACK 0 0 n n)
          C (submatrix ACK n 0 l n)
          K (submatrix ACK 0 n (+ n l) (- (ncols ACK) n))
          ;; recompute Gamma_i, Gamma_i-1 from A, C
          ; here we depart a bit to follow [DSC08]
          {U2 :u_perp} (rsvd Gamma_i)]
          ;M_R (mm (hcat (dge-eye l) (dge l (nrows U2))) U2 (trans U2))
          ;P (axpy -1 (mm (vcat A C) Jru) (vcat Jlu Y_i|i))
          ; If condition number of future input covariances > 1e4, it's poorly
          ; conditioned for our purposes and we use this approach:
          ;U_f_con (condition (mm U_f (trans U_f)))
          ; solve for B, D
          ;BD (argmin (norm2 (mm M_R Z_i)))]
          ; determine covariance matrices
     {:A A
      :C C
      ;:U_f_con U_f_con
      :spectrum (seq (dia S1))
      :order n
      :i i
      :method :N4SID-robust})))

(defn hankel-rq
  "(Lower) RQ factorisation of block-Hankel matrices, given snapshot matrices
  `U` and `Y`."
  ([ss i]
   (hankel-rq (:U ss) (:Y ss) i))
  ([U Y i]
   (let [p (- (* 2 i) 1)
         t (ncols U) 
         N (- t i i)
         H (scal! (/ 1.0 (Math/sqrt N)) (vcat (block-hankel U 0 p N) (block-hankel Y 0 p N)))
         ; default no pivoting
         qr (linalg/qrf (trans H))
         Q (trans (linalg/org qr))
         R (trans (view-tr (:or qr) {:uplo :upper}))]
     {:R R :Q Q})))

(defn partition-R
  "Given the lower-triangular `L` of the LQ-factorisation of `H`,
  Follows Chapter 6 / p.164 [vOdM-96]."
  ([ss i]
   (let [L (:R (hankel-rq ss i))
         l (ncols (:B ss))
         m (mrows (:C ss))]
     (partition-R L l m)))
  ([L l m]
   (let [R (view-ge L)
         ; infer i from size of R, R is 2i(m+l) x 2i(m+l)
         i (/ (mrows R) (* 2 (+ m l)))
         ; block row/column sizes
         x1 (* m i)
         x2 m
         x3 (- (* m i) m)
         x4 (* l i)
         x5 l
         x6 (- (* l i) l)
         ; block row/column starting indices
         ix1 0,
         ix2 (+ ix1 x1),
         ix3 (+ ix2 x2),
         ix4 (+ ix3 x3),
         ix5 (+ ix4 x4),
         ix6 (+ ix5 x5)
         ; spans
         x1-2 (+ x1 x2)
         x1-3 (+ x1 x2 x3)
         x1-4 (+ x1 x2 x3 x4)
         x1-5 (+ x1 x2 x3 x4 x5)
         x1-6 (+ x1 x2 x3 x4 x5 x6)
         x2-3 (+ x2 x3)
         x5-6 (+ x5 x6)]
    ; first block column
     {:i i
      :l l
      :m m
      :R11 (submatrix R ix1 ix1 x1 x1)
      :R21 (submatrix R ix2 ix1 x2 x1)
      :R31 (submatrix R ix3 ix1 x3 x1)
      :R41 (submatrix R ix4 ix1 x4 x1) 
      :R51 (submatrix R ix5 ix1 x5 x1)
      :R61 (submatrix R ix6 ix1 x6 x1)
     ; second block column
      :R22 (submatrix R ix2 ix2 x2 x2)
      :R32 (submatrix R ix3 ix2 x3 x2)
      :R42 (submatrix R ix4 ix2 x4 x2) 
      :R52 (submatrix R ix5 ix2 x5 x2)
      :R62 (submatrix R ix6 ix2 x6 x2)
     ; third block column
      :R33 (submatrix R ix3 ix3 x3 x3)
      :R43 (submatrix R ix4 ix3 x4 x3) 
      :R53 (submatrix R ix5 ix3 x5 x3)
      :R63 (submatrix R ix6 ix3 x6 x3)
     ; fourth block column
      :R44 (submatrix R ix4 ix4 x4 x4) 
      :R54 (submatrix R ix5 ix4 x5 x4)
      :R64 (submatrix R ix6 ix4 x6 x4)
     ; fifth block column
      :R55 (submatrix R ix5 ix5 x5 x5)
      :R65 (submatrix R ix6 ix5 x6 x5)
     ; sixth block column
      :R66 (submatrix R ix6 ix6 x6 x6)
     ; useful partitions
      :R5515 (submatrix R ix5 ix1 x5   x1-5)
      :R1515 (submatrix R ix1 ix1 x1-5 x1-5)
      :R5614 (submatrix R ix5 ix1 x5-6 x1-4)
      :R1414 (submatrix R ix1 ix1 x1-4 x1-4)
      :R6614 (submatrix R ix6 ix1 x6   x1-4)
      :R5514 (submatrix R ix5 ix1 x5   x1-4)
      :R2314 (submatrix R ix2 ix1 x2-3 x1-4)
      :R2313 (submatrix R ix2 ix1 x2-3 x1-3)
      :R1113 (submatrix R ix1 ix1 x1   x1-3)
      :R4413 (submatrix R ix4 ix1 x4   x1-3)
      :R6615 (submatrix R ix6 ix1 x6   x1-5)
      :R5615 (submatrix R ix5 ix1 x5-6 x1-5)
      :R2315 (submatrix R ix2 ix1 x2-3 x1-5)})))
      
(defn robust-rq
  "Robust identification algorithm of Chapter 6 of [vOdM-96], using the
  QR decomposition. See Figure 6.1."
  ; TODO: validate R decomposition. Gamma -> AC is correct, but the method overall
  ; gives wrong results.
  ([ss i n] (robust-rq (partition-R ss i) n))
  ([Rd n]
   (let [l (:l Rd),     m (:m Rd),    i (:i Rd)
         li (* l i),    mi (* m i),   mi2 (* m i 2)
         R44 (:R44 Rd)
         R1113 (:R1113 Rd)
         R1414 (:R1414 Rd)
         R2313 (:R2313 Rd)
         R2315 (:R2315 Rd)
         R4413 (:R4413 Rd)
         R5614 (:R5614 Rd)
         R5515 (:R5515 Rd)
         R5615 (:R5615 Rd)
         R6615 (:R6615 Rd)
         R2313' (trans R2313)
         I2mi (dge-eye mi2)
         ; step 1
         Ls (mm R5614 (pinv R1414))
         L_Up (submatrix Ls li mi)
         L_Uf (submatrix Ls 0 mi li mi)
         L_Yp (submatrix Ls 0 mi2 li li)
         ; step 2 (step 3 is to determine n, which is given
         ;Pi (axpy -1 (mm R2313' (pinv (mm R2313 R2313')) R2313) I2mi)
         Pi (axpy -1 (mm (pinv R2313) R2313) I2mi)
         F (hcat
             (axpy (mm L_Up R1113 Pi) (mm L_Yp R4413 Pi))
             (mm L_Yp R44))
         {S1 :sigma U1 :u V1' :vt} (rsvd F :rank n)
         ; step 4
         Gamma_i (mm U1 (sqrt S1))
         Gamma_up (submatrix Gamma_i l 0 (- (mrows Gamma_i) l) n)
         Gamma_down (submatrix Gamma_i (- (mrows Gamma_i) l) n)
         Gamma_i-1 Gamma_down
         Tl (vcat (mm (pinv Gamma_i-1) R6615)
                  R5515)
         Tr (vcat (mm (pinv Gamma_i) R5615)
                  R2315)
         S (mm Tl (pinv Tr))
         A (submatrix S n n)
         C (submatrix S n 0 m n)]
     {:S1 S1
      :A A
      :C C
      :A-old (mm (pinv Gamma_down) Gamma_up)
      :C-old (submatrix Gamma_i m n)
      :F F
      :S S
      :Gamma_i Gamma_i
      :Gamma_i-1 Gamma_i-1})))
  
(defn moesp
  "Explicit calculation, per (6-8) of [vOdM-95],
  following notation of [SSvH-04]."
  ([ss i n]
   (moesp (:U ss) (:Y ss) i n))
  ([U Y i n]
   (let [l (mrows Y)
         m (mrows U)
         t (ncols U)
         {:keys [N W_p Y_f U_f W_p+ Y_f- U_f-]} (block-hankel-matrices U Y i)
         H (scal! (/ 1.0 (Math/sqrt N)) (vcat U_f W_p Y_f))
         x1 (mrows U_f)
         x2 (mrows W_p)
         x3 (mrows Y_f)
         ; default no pivoting
         qr (linalg/qrf (trans H))
         ;Q (trans (linalg/org qr))
         L (view-ge (trans (view-ge (view-tr (:or qr) {:uplo :upper}))))
         L11 (submatrix L 0         0         x1 x1)
         L21 (submatrix L x1        0         x2 x1)
         L31 (submatrix L (+ x1 x2) 0         x3 x1)
         L22 (submatrix L x1        x1        x2 x2)
         L32 (submatrix L (+ x1 x2) x1        x3 x2)
         L33 (submatrix L (+ x1 x2) (+ x1 x2) x3 x3)
         {sigma_1 :sigma U_1 :u V_1' :vt} (rsvd L32 :rank n)
         Gamma_i (mm U_1 (sqrt sigma_1))
         Gamma_up (submatrix Gamma_i l 0 (- (mrows Gamma_i) l) n)
         Gamma_down (submatrix Gamma_i (- (mrows Gamma_i) l) n)
         A (mm (pinv Gamma_down) Gamma_up)
         C (submatrix Gamma_i l n)]
     {:A A
      :C C
      :order n
      :i i
      :method :MOESP})))

