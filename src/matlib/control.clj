(ns matlib.control
  "Linear algebra operations on matrices."
  (:require
   [matlib.core :refer [eps sq-eps dge-eye eye]]
   [uncomplicate.neanderthal
    [core :refer [transfer! copy! scal! axpy entry nrm2 sum mm dia dim ncols mrows trans view-vctr subvector submatrix]]
    [vect-math :refer [sqrt inv mul]]
    [native :refer [dgd dge]]
    [linalg :refer :all]]))

(defn obsv
  "Extended observability matrix up to order `i`. If `i` is not specified, up to
  order `n`."
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

