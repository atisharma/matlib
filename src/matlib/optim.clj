(ns matlib.optim
  "Various optimisation algorithms.
  L-BFGS and gradient descent are implemented.

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
    [matlib.core :refer [eps sq-eps shift-update dge-eye last-cols col-vector ones take-cols drop-cols]]
    [uncomplicate.neanderthal
     [core :refer [transfer! copy! copy scal! scal xpy axpy axpy! dot nrm1 nrm2 nrmi sum mm dia dim col ncols mrows trans view-vctr subvector submatrix vctr?]]
     [vect-math :refer [sqrt inv mul div]]
     [native :refer [dgd dge dv]]
     [real :refer [entry entry!]]
     [random :refer :all]
     [linalg :refer :all]]))

(def ^:private ip (/ (- (Math/sqrt 5.0) 1.0) 2.0))
(def ^:private ip2 (/ (- 3.0 (Math/sqrt 5.0)) 2.0))

(defn- tolerance
  "Default tolerance for iterative functions, `sqrt eps * |x|₂`."
  [x]
  (if (vctr? x)
    (* sq-eps (max 1.0 (nrm2 x)))
    (* sq-eps (max 1.0 (Math/abs x)))))

(defn vctr-grad
  "Gradient of `f` at `x`, approximated by central finite-difference,
  where `x` is a Neanderthal vctr.
  `f: ℝⁿ -> ℝ`."
  ([f x]
   (vctr-grad f x (tolerance x)))
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
  "Approximate gradient of a function `phi` by central finite-difference, where
  `phi: ℝ -> ℝ."
  ([phi ^double x]
   (scalar-grad phi x (tolerance x)))
  ([phi ^double x ^double h]
   (let [x+ (+ x h)
         x- (- x h)]
     (/ (- (phi x+) (phi x-)) (* 2 h)))))

(defn grad
  "Approximate vector or scalar gradient of function `f` using central
  finite-difference."
  [f x]
  (if (vctr? x)
    (vctr-grad f x)
    (scalar-grad f x)))

(defn golden-section
  "Return double `x` that minimises differentiable `phi(x)` using golden mean search.
  If bounds `x-` and `x+` are not given, they are assumed to be +-1e20."
  ([phi]
   (golden-section phi -1e20 1e20))
  ([phi x- x+ & args]
   (let [a (min x- x+)
         b (max x- x+)
         {:keys [h c d fc fd k] :or {h (- b a)
                                     c (+ a (* ip2 h))
                                     d (+ a (* ip h))
                                     fc (phi c)
                                     fd (phi d)
                                     k 0}} args]
     (cond (<= h (tolerance x+))  {:sol (/ (+ a b) 2.0) :iterations k}
           (< fc fd)              (recur phi a d {:h (* h ip) :d c :fd fc :k (inc k)})
           (>= fc fd)             (recur phi c b {:h (* h ip) :c d :fc fd :k (inc k)})))))

(defn- zoom
  "Algorithm 3.6 of [NW-06] using bisection.  
  `phi` phi: ℝ -> ℝ, usually `phi(a) = f(x + ap)`  
  `a-`  lower bound on `a`  
  `a+`  upper bound on `a`  
  `0 < c1 < c2 < 1`."
  [phi a- a+ phi0 phi'0 c1 c2 i]
  (let [maxiter 20
        a (/ (+ a+ a-) 2.0)
        phia (phi a)
        phi'a (scalar-grad phi a)]
    (cond (> i maxiter)                           {:sol a :success false}
          (> phia (+ phi0 (* c1 a phi'0)))        (recur phi a- a phi0 phi'0 c1 c2 (inc i))
          (>= phia (phi a-))                      (recur phi a- a phi0 phi'0 c1 c2 (inc i))
          (<= (Math/abs phi'a) (* -1 c2 phi'0))   {:sol a :success true}
          (>= (* phi'a (- a+ a-)) 0.0)            (recur phi a a- phi0 phi'0 c1 c2 (inc i))
          :else                                   (recur phi a a+ phi0 phi'0 c1 c2 (inc i)))))
               
(defn wolfe
  "Perform a line search to find a step length `a` satisfying the strong Wolfe
  conditions.  
  `phi` function of scalar a  
  `a` starting value.  
  See p.60 [NW-06]."
  ([phi & options]
   (let [{:keys [a a- amax maxiter k phi0 phi'0 c1 c2]
          :or {a 1.0
               a- 0.0
               amax 2.0
               maxiter 20
               k 0
               phi0 (phi 0)
               phi'0 (scalar-grad phi 0)
               c1 1e-4
               c2 0.9}} options
         phi_a (phi a)
         phi'a (scalar-grad phi a)
         phi_a- (phi a-)]
     (cond (> k maxiter) {:sol a :success false}
           (> phi_a (+ phi0 (* c1 a phi'0))) (zoom phi a- a phi0 phi'0 c1 c2 0)
           (and (> phi_a phi_a-) (> k 1)) (zoom phi a- a phi0 phi'0 c1 c2 0)
           (<= (Math/abs phi'a) (* -1 c2 phi'0)) {:sol a :success true}
           (>= phi'a 0.0) (zoom phi a a- phi0 phi'0 c1 c2 0)
           :else (recur phi {:a (/ (+ a amax) 2.0)
                             :a- a
                             :amax amax
                             :maxiter maxiter
                             :k (inc k)
                             :phi0 phi0
                             :phi'0 phi'0
                             :c1 c1
                             :c2 c2})))))

(defn backtrack
  "Perform a backtracking line search to find a step length `a` satisfying the
  Armijo-Goldstein conditions.  
  `phi` function of scalar a  
  `a` starting value.  
  See p.60 [NW-06]."
  ([phi & options]
   (let [{:keys [a maxiter k tau c phi0 t]
          :or {a 10.0
               maxiter 1000
               k 0
               tau 0.5
               c 0.5
               phi0 (phi 0.0)
               t (* -1 c (scalar-grad phi 0))}} options
         armijo (>= (- phi0 (phi a)) (* a t))]
     (cond armijo {:sol a :success true :iterations k}
           (> k maxiter) {:sol a :success false :iterations k}
           :else (recur phi {:a (* tau a)
                             :maxiter maxiter
                             :k (inc k)
                             :tau tau
                             :c c
                             :phi0 phi0
                             :t t})))))

(defn gradient-descent
  "Gradient descent with line search.  
  `f`   function to be solved, `f: ℝⁿ -> ℝ`.  
  `x`   initial solution guess  
  options (default):
  `:tol`      solution converges when `(nrm1 (grad f x)) < tol)`, (`sqrt eps * |x₀|²`)  
  `:maxiter`  maximum iterations (1000)  
  `:output`   print progress every iteration (`false`)  
  `:lsmethod` line-search method for step length, one of `:wolfe`, `:gs`, `:backtrack` (`wolfe`)  
  "
  ([f x & options]
   (let [{:keys [tol maxiter output k lsmethod]
          :or {tol (tolerance x)
               maxiter 1000
               output false
               k 0
               lsmethod :wolfe}} options
         linesearch (get {:wolfe wolfe :gs golden-section :backtrack backtrack} lsmethod wolfe)
         g (vctr-grad f x)
         q (scal -1 g)
         a (:sol (wolfe #(f (axpy % q x))))
         x+ (axpy a q x)
         success (< (nrm1 g) tol)]
     (when output
       (when (= k 0) (print options "\n"))
       (printf "k: %05d\ta: %8.5f\t|grad| %10.3f\t|dx|: %10.3f\t|f(x+)|: %10.4f\n"
               k a (nrmi g) (nrm2 (axpy -1 x x+)) (f x+)))
     (cond success        (merge options {:sol x+ :f (f x+) :grad g :iterations k :success true :lsmethod lsmethod})
           (> k maxiter)  (merge options {:sol x+ :f (f x+) :grad g :iterations k :success false :lsmethod lsmethod})
           :else          (recur f x+ {:tol tol :maxiter maxiter :output output :k (inc k) :lsmethod lsmethod})))))

(defn- alg-7-4
  "Two-loop recursion in L-BFGS to calculate descent direction as the
  product of Hessian with the gradient.  
  `S` is a matrix whose columns are `sᵢ`  
  `Y` is a matrix whose columns are `yᵢ`  
  `q` is the supplied gradient of `f` at `xₖ`.
  See Algorithm 7.4 of [NW-06]."
  [S Y q k]
  ; i = k-1, ..., k-m
  (let [c (ncols S)
        m (min c k)
        is (range (- c m) c)  ; k-m, ..., k-1
        irhos (dv (repeat c 0.0))
        as (dv (repeat c 0.0))
        sl (col S (- c 1))
        yl (col Y (- c 1))
        yl2 (dot yl yl)
        gamma (/ (dot yl sl) yl2)
        r (dv (repeat (dim q) 0.0))] 
    ; i = k-1, ... k-m
    (doseq [i (reverse is)]
      (let [s (col S i)
            y (col Y i)
            irho (dot y s)
            a (/ (dot s q) irho)]
        ; update as, 1/rhos, q in-place
        (entry! irhos i irho)
        (entry! as i a)
        (axpy! (- a) y q)))
    (copy! (scal gamma q) r)  ; r0 = H^0_k, H^0_k = gamma I
    ; i = k-m, ..., k-1
    (doseq [i is]
      (let [irho (entry irhos i)
            a (entry as i)
            s (col S i)
            y (col Y i)
            b (/ (dot y r) irho)]
        (axpy! (- a b) s r)))
    (scal! 1 r)))

(defn- alg-7-5
  "Internal L-BFGS solver. See Algorithm 7.5 of [NW-06]."
  ([f x S Y k linesearch tol maxiter history output]
   (let [q (vctr-grad f x)
         ; search direction, start off downhill
         p (if (= 0 k) (scal -1 q) (alg-7-4 S Y q k))
         a (:sol (linesearch #(f (axpy % p x))))
         x+ (axpy a p x)
         s (axpy -1 x x+)
         y (axpy -1 (grad f x+) (grad f x))
         X (dge (dim x) maxiter)
         success (< (nrmi q) tol)]
     (when output
       (printf "k: %5d\ta: %8.5f\t|q| %10.3f\tp.s: %10.3f\t|f(x+)|: %10.4f\n"
               k a (nrmi q) (dot p s) (f x+)))
     (shift-update S (col-vector s))
     (shift-update Y (col-vector y))
     (when history
       (shift-update X (col-vector x)))
     (cond success        (merge {:sol x+ :f (f x+) :iterations k :tol tol :maxiter maxiter :output output :success true} (if history {:X X :S S :Y Y} {}))
           (> k maxiter)  (merge {:sol x+ :f (f x+) :iterations k :tol tol :maxiter maxiter :output output :success false} (if history {:X X :S S :Y Y} {}))
           :else          (recur f x+ S Y (inc k) linesearch tol maxiter history output)))))
        
(defn l-bfgs
  "L-BFGS, using finite-difference gradient approximation.
  `f`   function to be solved, `f: ℝⁿ -> ℝ`.  
  `x`   initial solution guess  
  options (default):
  `:m`        number of last iterations to store for approximate Hessian (20)  
  `:tol`      solution converges when `(nrm1 (grad f x)) < tol)`, (`sqrt eps * |x₀|²`)  
  `:maxiter`  maximum iterations ( 1000)  
  `:history`  return search history (`false`)  
  `:output`   print progress every iteration (`false`)  
  `:lsmethod` line-search method for step length, one of `:wolfe`, `:gs`, `:backtrack` (`wolfe`)  
  "
  ([f x & options]
   (let [{:keys [tol maxiter output m history lsmethod]
          :or {tol (tolerance x)
               maxiter 1000
               output false
               m 20
               history false
               lsmethod :wolfe}} options
         linesearch (get {:wolfe wolfe :gs golden-section :backtrack backtrack} lsmethod wolfe)
         S (dge (dim x) m)
         Y (dge (dim x) m)]
     (when output
       (print options "\n")
       (print "\t\ta: steplength\ts: step\t\ty: ddf/dx\tq: df/dx\tp: search dirn\n"))
     (merge (alg-7-5 f (copy x) S Y 0 linesearch tol maxiter history output) {:lsmethod lsmethod}))))

(defn- booth
  "Booth function f: ℝ² -> ℝ, minimum at f(1, 3) = 0."
  [v]
  (let [x (entry v 0)
        y (entry v 1)]
    (+ (Math/pow (+ x y y -7) 2)
       (Math/pow (+ x x y -5) 2))))

(defn- beale
  "Beale function f: ℝ² -> ℝ, minimum at f(3, 0.5) = 0."
  [v]
  (let [x (entry v 0)
        y (entry v 1)
        mx (- x)]
    (+ (Math/pow (+ 1.5 (* x y) mx) 2)
       (Math/pow (+ 2.25 (* x y y) mx) 2)
       (Math/pow (+ 2.625 (* x y y y) mx) 2))))

(defn- rosenbrock
  "Rosenbrock function f: ℝ² -> ℝ.
  Since there is a large flat valley, requires a low tolerance."
  ([v]
   (rosenbrock v 1 100))
  ([v a b]
   (let [x (entry v 0)
         y (entry v 1)]
     (+ (Math/pow (- a x) 2)
        (* b (Math/pow (- y (* x x)) 2))))))
