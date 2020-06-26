(ns matlib.stats
  "Basic statistical matrix operations."
  (:require
    [matlib.core :refer :all]
    [matlib.linalg :refer :all]
    [uncomplicate.neanderthal.real :refer [entry entry!]]
    [uncomplicate.neanderthal.native :refer :all :exclude [sv]]
    [uncomplicate.neanderthal.linalg :refer :all]
    [uncomplicate.neanderthal.core :refer :all :exclude [entry entry!]]
    [uncomplicate.neanderthal.vect-math :refer :all]
    [uncomplicate.neanderthal.random :refer :all]))

(defn row-mean
  "Mean across rows."
  ([M]
   (scal! (/ 1.0 (mrows M)) (dge (map sum (cols M))))))

(defn col-mean
  "Mean across columns."
  ([M]
   (scal! (/ 1.0 (ncols M)) (dge (map sum (rows M))))))

(defn subtract-col-mean
  "Subtract the mean over columns from each row."
  [M]
  (let [ave (col-vector (col-mean M))
        m-ones (dge 1 (ncols M) (repeat (double -1)))
        result (mm ave m-ones)]
    (axpy! M result)))

(defn scale-col-mean
  "Divide each row by its mean."
  [M]
  (let [ave (col-mean M)
        ave-inv (div! (ones ave) ave)
        d-ave-inv (transfer! ave-inv (dgd (dim ave)))]
    (mm d-ave-inv M)))

(defn covar
  "Covariance matrix from a matrix of observations `M`.
  Columns of `M` correspond to each sample."
  [M]
  (let [R (subtract-col-mean M)]
    (scal! (/ 1.0 (ncols R)) (mm R (trans R)))))

(defn sd
  "Unscaled sqrt of the diagonal of the covariance matrix."
  [M]
  (sqrt! (copy (dia (covar M)))))

(defn snr
  "Unscaled mean / std dev."
  [M]
  (div (col-mean M) (sd M)))

(defn fuzz
  "Add normally distributed random numbers to `M`. The distribution has mean of 
  `mu` (default 0) and standard deviation `sigma`."
  ([M mu sigma]
   (axpy! M (rand-normal! mu sigma (copy M))))
  ([M sigma]
   (fuzz M 0 sigma)))
