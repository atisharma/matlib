(ns matlib.linalg
  "Linear algebra operations on matrices.
  
  Useful references:
  
  [S-06]  
  'Linear Algebra and Its Applications (4th Ed.)'  
  G Strang  
  Wellesley-Cambridge Press (2016)  
  ISBN 10: 0980232775  
  ISBN 13: 9780980232776

  [TB-97]  
  'Numerical Linear Algebra'  
  Lloyd N. Trefethen, David Bau III  
  SIAM: Society for Industrial and Applied Mathematics (1997)
  ISBN 10: 0898713617  
  ISBN 13: 9780898713619

  "
  (:require
    [matlib.core :refer :all]
    [uncomplicate.neanderthal.real :refer [entry entry!]]
    [uncomplicate.neanderthal.native :refer :all :exclude [sv]]
    [uncomplicate.neanderthal.linalg :refer :all]
    [uncomplicate.neanderthal.core :refer :all :exclude [entry entry!]]
    [uncomplicate.neanderthal.vect-math :as vect-math]
    [uncomplicate.neanderthal.random :as random]))

(defn minv
  "Matrix inverse of `M` based on LU decomposition."
  [M]
  (tri! (trf M)))

(defn rsvd
  "The reduced SVD of `M`,  
  `M = [ U₁ | U₂ ] [ Σ₁ | 0  ] [ V₁' ]`  
  `                [ 0  | Σ₂ ] [ V₂' ]`  
  where the leading subspace is given with `u, vt, sigma` and similarly with 
  subscript `2` as `u_perp` etc."
  ([M & options]
   (let [{:keys [tol rank] :or {tol sq-eps, rank nil}} options
         {:keys [u vt sigma]} (svd M true true)
         r (if rank
             rank
             (count (take-while #(> % tol) (lazy-seq (dia sigma)))))
         p (- (min (mrows M) (ncols M)) r)
         sigma_r (transfer!
                   (subvector (view-vctr sigma) 0 r)
                   (dgd r))
         sigma_perp (transfer!
                      (subvector (view-vctr sigma) r p)
                      (dgd p))
         vt_r (submatrix vt r (ncols vt))
         u_r (submatrix u (mrows u) r)
         vt_perp (submatrix vt r 0 (- (mrows vt) r) (ncols vt))
         u_perp (submatrix u 0 r (mrows u) (- (ncols u) r))]
     {:master true
      :sigma sigma_r
      :u u_r
      :vt vt_r
      :u_perp u_perp
      :vt_perp vt_perp
      :sigma_perp sigma_perp})))

(defn pinv
  "`M⁺`, the Moore-Penrose pseudoinverse of real matrix `M`, such that
  `M M⁺ M = M`,  
  `M⁺ M M⁺ = M⁺`  
  and `M M⁺` and `M⁺ M` are symmetric.
  The subspace associated with singular values less than or equal
  to `:tol` (or, improperly, after rank `:rank`) will be projected away."
  ([M & options]
   (let [{:keys [tol rank] :or {tol sq-eps, rank nil}} options
         {:keys [u vt sigma]} (rsvd M :tol tol :rank rank)
         inv-sigma (vect-math/inv sigma)]
     (mm (trans vt) inv-sigma (trans u)))))

(defn rank
  "Estimate rank of a matrix based on the ratio of singular values."
  ([M & options]
   (let [{:keys [tol rank] :or {tol sq-eps, rank nil}} options
         sigma (seq (dia (:sigma (svd M))))
         s1 (first sigma)]
     (count (filter #(> % tol) sigma)))))
       
(defn cokern
  "(Left) null-space or kernel `X` of a matrix `M`, such that `x'M=0`, where 
  `x` is any linear combination of columns of `X`."
  ([M & options]
   (let [{:keys [tol] :or {tol sq-eps}} options
         m (mrows M)
         n (ncols M)
         d (max m n)
         Z (dge m (- d n))
         {:keys [u sigma]} (svd (hcat M Z) true false)
         r (count (filter #(> % tol) (dia sigma)))
         u (submatrix u 0 r (mrows u) (- (ncols u) r))]
       u)))

(defn kern
  "(Right) null-space or kernel `X` of a matrix `M`, such that `Mx=0`, where 
  `x` is any linear combination of columns of `X`."
  ([M & options]
   (let [{:keys [tol] :or {tol sq-eps}} options
         m (mrows M)
         n (ncols M)
         d (max m n)
         Z (dge (- d m) n)
         {:keys [vt sigma]} (svd (vcat M Z) false true)
         r (count (filter #(> % tol) (dia sigma)))
         vt_perp (submatrix vt r 0 (- (mrows vt) r) (ncols vt))]
       (trans vt_perp))))

(defn span
  "Span (range) `R` of matrix `M`, such that for `r=Mx` and any vector `x`,
  `r` is a linear combination of columns of `R`."
  ([M & options]
   (let [{:keys [tol rank] :or {tol sq-eps, rank nil}} options]
     (:u (rsvd M :tol tol :rank rank)))))

(defn rsp
  "Projection of row-space of `E` onto row space of `F`, `E/F = EF⁺F`."
  [E F]
  (mm E (pinv F) F))

(defn rsp-perp
  "Orthogonal complement of the row-space projection of `E` onto `F`,
  `E/F^⟂ = E-E/F`."
  [E F]
  (axpy 1 E -1 (rsp E F)))

(defn oblique-rsp
  "Oblique projection of row space of `E` onto row space of `F` along row space
  of `G`, `E/_GF = (E/G^⟂)(F/G^⟂)⁺F`."
  [E G F]
  (mm (rsp-perp E G) (pinv (rsp-perp F G)) F))

(defn hada
  "Hadamard (element-wise) product of matrices."
  ([A]
   A)
  ([A B]
   (dge (mrows A) (ncols B) (vect-math/mul (view-vctr A) (view-vctr B))))
  ([A B & rst]
   (reduce hada (hada A B) rst)))
  
(defn kron
  "Kroneker product of matrices."
  [A B]
  :not-implemented)

(defn hess
  "Hessenberg form of matrix `M`"
  [M]
  :not-implemented)

(defn trace
  "Trace of `M`."
  [M]
  (sum (dia M)))

(defn schur-ordered
  "Sorted Schur decomposition."
  [M]
  :not-implemented)
