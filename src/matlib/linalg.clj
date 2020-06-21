(ns matlib.linalg
  "Linear algebra operations on matrices."
  (:require
   [matlib.core :refer [eps sq-eps]]
   [uncomplicate.neanderthal
    [core :refer [transfer! copy! scal! axpy entry nrm2 sum mm dia dim ncols mrows trans view-vctr subvector submatrix]]
    [vect-math :refer [sqrt inv mul]]
    [native :refer [dgd dge]]
    [linalg :refer :all]]))

(defn minv
  "Matrix inverse of `M` based on LU decomposition."
  [M]
  (tri! (trf M)))

(defn rsvd
  "The reduced SVD."
  ([M & options]
   (let [{:keys [tol rank] :or {tol sq-eps, rank nil}} options
         {:keys [u vt sigma]} (svd M true true)
         r (if rank
             rank
             (count (take-while #(> % tol) (lazy-seq (dia sigma)))))
         sigma_r (transfer!
                   (subvector (view-vctr sigma) 0 r)
                   (dgd r))
         vt_r (submatrix vt r (ncols vt))
         u_r (submatrix u (mrows u) r)]
     {:master true
      :sigma sigma_r
      :u u_r
      :vt vt_r})))

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
         inv-sigma (inv sigma)]
     (mm (trans vt) inv-sigma (trans u)))))

(defn rank
  "Estimate rank of a matrix based on the ratio of singular values."
  ([M & options]
   (let [{:keys [tol rank] :or {tol sq-eps, rank nil}} options
         sigma (seq (dia (:sigma (svd M))))
         s1 (first sigma)]
     (count (filter #(> % tol) sigma)))))
       
(defn condition
  "Condition number of a matrix."
  [M]
  (let [sigma (seq (dia (:sigma (svd M))))]
    (/ (first sigma) (last sigma))))

(defn kern
  "(Right) null-space or kernel `X` of a matrix `M`, such that `Mx=0`, where 
  `x` is any linear combination of columns of `X`."
  ([M & options]
   (let [{:keys [tol] :or {tol sq-eps}} options
         {:keys [u vt sigma]} (svd M true true)
         r (count (take-while #(> % tol) (dia sigma)))
         V-perp (submatrix (trans vt) 0 r (ncols vt) (- (mrows vt) r))]
     V-perp)))

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
   (dge (mrows A) (ncols B) (mul (view-vctr A) (view-vctr B))))
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

(defn schur
  "Schur form of matrix `M`."
  [M]
  :not-implemented)

(defn trace
  "Trace of `M`."
  [M]
  (sum (dia M)))
