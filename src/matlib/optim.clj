(ns matlib.optim
  "Various optimisation algorithms.

  [NW-06]
  Numerical Optimization (second Ed.)
  J Nocedal, S Wright
  Springer-Verlag New York (2006)
  ISBN: 978-0-387-40065-5
  DOI: 10.1007/978-0-387-40065-5
  ISBN: 978-0-387-30303-1, 978-1-4939-3711-0

  [L-BFGS]
  On the limited memory BFGS method for large scale optimization methods  
  DC Liu and J Nocedal  
  Mathematical Programming Vol. 45, pp. 503-528 (1989).
  "
  (:require
    [matlib.core :refer [eps sq-eps shift-update dge-eye last-cols col-vector take-cols drop-cols]]
    [data.plot :as plot]
    [uncomplicate.neanderthal
     [core :refer [transfer! copy! copy scal! scal xpy axpy axpy! dot nrm1 nrm2 nrmi sum mm dia dim col ncols mrows trans view-vctr subvector submatrix vctr?]]
     [vect-math :refer [sqrt inv mul div]]
     [native :refer [dgd dge dv]]
     [real :refer [entry entry!]]
     [random :refer :all]
     [linalg :refer :all]]))

(def ip (/ (- (Math/sqrt 5.0) 1.0) 2.0))
(def ip2 (/ (- 3.0 (Math/sqrt 5.0)) 2.0))

(defn vctr-grad
  "Gradient of `f` at `x`, approximated by central finite-difference,
  where `x` is a Neanderthal vctr.
  `f: ℝⁿ -> ℝ`."
  ([f x]
   (vctr-grad f x (* sq-eps (max 1.0 (nrm2 x)))))
  ([f x h]
   (let [df_dx (dv (repeat (dim x) 0))]
     (doseq [i (range (dim x))]
       (let [dx (dv (repeat (dim x) 0))]
         (entry! dx i h)  ; perturb in one direction
         (let [x+ (axpy  1.0 dx x)
               x- (axpy -1.0 dx x)
               df_dx_i (/ (- (f x+) (f x-)) (+ h h))]
           (entry! df_dx i df_dx_i))))
     df_dx)))

(defn scalar-grad
  "Approximate gradient of a function `phi` by forward finite-difference, where
  `phi: ℝ -> ℝ."
  ([phi ^double x]
   (scalar-grad phi x (* sq-eps (max 1.0 (Math/abs x)))))
  ([phi ^double x ^double h]
   (let [x+ (+ x h)
         x- (- x h)]
     (/ (- (phi x+) (phi x-)) (* 2 h)))))

(defn grad
  "Approximate vector or scalar gradient."
  [f x]
  (if (vctr? x)
    (vctr-grad f x)
    (scalar-grad f x)))

(defn golden-section
  "Return double `x` that minimises differentiable `phi(x)` using golden mean search.
  If bounds `x-` and `x+` are not given, they are assumed to be +-1e10."
  ([phi]
   (golden-section phi -1e10 1e10 {}))
  ([phi x- x+ args]
   (let [a (min x- x+)
         b (max x- x+)
         {:keys [h c d fc fd tol] :or {h (- b a)
                                       c (+ a (* ip2 h))
                                       d (+ a (* ip h))
                                       fc (phi c)
                                       fd (phi d)
                                       tol sq-eps}} args]
     (cond (<= h tol) a
           (< fc fd) (recur phi a d {:h (* h ip) :d c :fd fc :tol tol})
           (>= fc fd) (recur phi c b {:h (* h ip) :c d :fc fd :tol tol})))))

(defn optimal?
  "Check for local optimality, using infinity norm of the gradient."
  [f x tol]
  (< (nrmi (grad f x)) tol))

(defn- zoom
  "Algorithm 3.6 of [NW-06] using bisection.  
  `phi` phi: ℝ -> ℝ, usually `phi(a) = f(x + ap)`  
  `a-`  lower bound on `a`  
  `a+`  upper bound on `a`  
  `0 < c1 < c2 < 1`."
  [phi a- a+ c1 c2 i]
  (let [max-iter 20
        aj (/ (+ a+ a-) 2.0)
        phiaj (phi aj)
        phi'0 (scalar-grad phi 0.0)
        phi'aj (scalar-grad phi aj)]
    (cond (> i max-iter)                          aj
          (> phiaj (+ (phi 0.0) (* c1 aj phi'0))) (recur phi a- aj c1 c2 (inc i))
          (>= phiaj (phi a-))                     (recur phi a- aj c1 c2 (inc i))
          (> (Math/abs phi'aj) (- (* c2 phi'0)))  aj
          (>= (* phi'aj (- a+ a-)) 0.0)           (recur phi aj a- c1 c2 (inc i))
          :else                                   (recur phi aj a+ c1 c2 (inc i)))))
               
(defn- wolfe
  "Perform a line search to find a step length `a` satisfying the strong Wolfe
  conditions.  
  `phi` function of scalar a  
  `a` starting value.  
  See p.60 [NW-06]."
  ([phi]
   (wolfe phi 1.0 0.0 1))
  ([phi a_i a_i-1 i]
   (let [c1 1e-4  ; from [NW-06]
         c2 0.9   ; from [NW-06]
         a_max 1.0
         maxiter 20
         phi_a_i (phi a_i)
         phi'a_i (scalar-grad phi a_i)
         phi'0 (scalar-grad phi 0)
         phi_a_i-1 (phi a_i-1)]
     (cond (> i maxiter)                              a_i
           (> phi_a_i (+ (phi 0) (* c1 a_i phi'0)))   (zoom phi a_i-1 a_i c1 c2 0)
           (and (> phi_a_i phi_a_i-1) (> i 1))        (zoom phi a_i-1 a_i c1 c2 0)
           (<= (Math/abs phi'a_i) (* -1 c2 phi'0))    a_i
           (>= phi'a_i 0)                             (zoom phi a_i a_i-1 c1 c2 0)
           :else                 (recur phi (/ (+ a_i a_max) 2.0) a_i (inc i))))))

(defn backtrack
  "Backtracking line search."
  ([phi] (backtrack phi 5.0 0))
  ([phi a j]
   (let [tau 0.5
         c 0.5
         m (scalar-grad phi a)
         t (* -1 c m)]
     (if (or (< (- (phi 0) (phi a)) (* a t)) (< a 1e-8))
       a
       (recur phi (* tau a) j)))))

(defn- alg-7-4
  "Two-loop recursion in L-BFGS to calculate product of Hessian with gradient.  
  `S` is a matrix whose columns are `sᵢ`  
  `Y` is a matrix whose columns are `yᵢ`  
  `q` is the supplied gradient of `f` at `xₖ`.
  See Algorithm 7.4 of [NW-06]."
  [S Y q k]
  (assert (= (ncols S) (ncols Y)))
  (assert (= (mrows S) (mrows Y) (dim q)))
  ; i = k-1, ..., k-m
  (let [c (ncols S)
        m (min c k)
        is (range (- c m) c)  ; k-m, ..., k-1
        rhos (dv (repeat c 0.0))
        as (copy rhos)
        sl (last-cols S)
        yl (last-cols Y)
        yl2 (dot yl yl)
        gamma (/ (dot yl sl) (max yl2 1e-2))
        r (dv (repeat (dim q) 0.0))] 
    ; i = k-1, ... k-m
    (doseq [i (reverse is)]
      (let [s (col S i)
            y (col Y i)
            rho (/ 1.0 (dot y s))
            a (* rho (dot s q))]
        ; update as, rhos, q in-place
        (entry! rhos i rho)
        (entry! as i a)
        (axpy! (- a) y q)))
    (copy! (scal gamma q) r)  ; r0 = H^0_k, H^0_k = gamma I
    ; i = k-m, ..., k-1
    (doseq [i is]
      (let [rho (entry rhos i)
            a (entry as i)
            s (col S i)
            y (col Y i)
            b (* rho (dot y r))]
        (axpy! (- b a) s r)))
    r))

(defn- alg-7-5
  "Main L-BFGS loop. See Algorithm 7.5 of [NW-06]."
  ([f x S Y k tol maxiter]
   (let [q (vctr-grad f x)
         p (alg-7-4 S Y q k)
         ;a (backtrack #(f (axpy % p x)))
         ;a (wolfe #(f (axpy % p x)))
         a (golden-section #(f (axpy % p x)))
         x+ (axpy (- a) p x)
         s (axpy -1 x x+)
         y (axpy -1 (grad f x+) (grad f x))
         X (dge (dim x) maxiter)]
     (shift-update S (col-vector s))
     (shift-update Y (col-vector y))
     (shift-update X (col-vector x))
     (printf "k: %05d\ta: %8.5f\t|s|: %10.3f\t|y|: %10.3f\t|q| %10.3f\t|p|: %10.3f\t|f(x+)|: %10.4f\n"
             k a (nrm2 s) (nrm2 y) (nrmi q) (nrm2 p) (f x+))
     ;(print "p, s: " p s "\n")
     (cond (< (nrmi q) tol) {:sol x+ :k k :X X :S S :Y Y :success true}
           (> k maxiter) {:sol x+ :k k :X X :S S :Y Y :success false}
           :else (recur f x+ S Y (inc k) tol maxiter)))))
        
(defn l-bfgs
  "A pure Neanderthal implementation of L-BFGS, using finite-difference gradient
  approximation.
  `f`   function to be solved, `f: ℝⁿ -> ℝ`.  
  `x`   initial solution guess  
  `m`   number of last iterations to store."
  ([f x m & options]
   (let [{:keys [tol maxiter] :or {tol (* sq-eps (nrm2 x))}, maxiter 100} options
         X (rand-normal! (dge (dim x) (+ m 1)))  ; some random states
         S (axpy -1.0 (submatrix X (dim x) m) (submatrix X 0 1 (dim x) m))
         Y (dge (dim x) m)]
     (printf "tol %f\tmaxiter: %d\n" tol maxiter)
     (print "\t\ta: steplength\ts: step\t\ty: ddf/dx\tq: df/dx\tp: search dirn\n")
     (doseq [i (range m)]
       (copy!
         (axpy -1
               (vctr-grad f (col X i))
               (vctr-grad f (col X (+ i 1))))
         (col Y i)))
     (alg-7-5 f (copy x) S Y 0 tol maxiter))))
  
(defn f
  "A test function f: ℝⁿ -> ℝ."
  [x]
  (+ 10.0 (nrm2 x)))

(defn make-phi
  "A test function phi(a) = f(x + ap)."
  [f x p]
  #(f (axpy % p x)))
