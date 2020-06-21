(ns matlib.optim
  "Optimisation algorithms."
  (:require
   [uncomplicate.neanderthal
    [core :refer [transfer! copy! scal! axpy entry nrm2 sum mm dia dim ncols mrows trans view-vctr subvector submatrix]]
    [vect-math :refer [sqrt inv mul]]
    [native :refer [dgd dge]]
    [linalg :refer :all]]))


(defn bfgs
  []
  :not-implemented)

(defn l-bfgs
  []
  :not-implemented)


(defn l-obfgs
  []
  :not-implemented)

(defn gmres
  []
  :not-implemented)
