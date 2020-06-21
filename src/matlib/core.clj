(ns matlib.core
  "Basic matrix constructions, operations and patterns."
  (:require
   [uncomplicate.neanderthal
    [native :refer [dge dgd native-double]]
    [core :refer [alter! entry entry! copy copy! scal! dim vctr vctr? view-vctr mm axpy ncols mrows trans submatrix]]]))

(defn machine-epsilon
  "Approximate machine epsilon."
  ; Machine precision is about 1.11e-16 on 64bit dp according to LAPACK user guide
  ([]
   (machine-epsilon 1e-5))
  ([e]
   (if (= (+ 1.0 e) 1.0)
     e
     (recur (* 0.99 e)))))

(def eps (machine-epsilon 1.0))

(def sq-eps (Math/sqrt eps))

(defn nan?
  "True iff `##NaN`."
  [^double x]
  (and (number? x) (Double/isNaN x)))

(defn numeric?
  "True if not `##NaN` and is a number. True for `##Inf`."
  [x]
  (and (number? x) (Double/isFinite x) (not (Double/isNaN x))))

(defn ones
  "Matrix (vector) of ones with same dimension as `M`."
  ([M]
   (if (vctr? M)
     (vctr native-double (repeat (dim M) 1))
     (dge (mrows M) (ncols M) (repeat (double 1)))))
  ([m n]
   (dge m n (repeat (double 1)))))

(defn eye
  "Identity matrix."
  [n]
  (dgd n (repeat 1)))

(defn sign
  "Sign of each element in a matrix or vector."
  [M]
  (alter! (copy M) (fn ^double [^double x] (double (Math/signum x)))))

(defn negatives
  "Elements are `1` where matrix/vector `M` has negative elements."
  [M]
  (alter! (copy M) (fn ^double [^double x] (if (> 0.0 x) 1.0 0.0))))

(defn positives
  "Elements are `1` where `M` has positive elements."
  [M]
  (alter! (copy M) (fn ^double [^long i ^long j ^double x] (if (< 0 x) 1.0 0.0))))

(defn replace-nan!
  "Replace `##NaN` in-place (with default 0)."
  ([M x]
   (let [d (double x)]
     (alter! M (fn ^double [^double y] (if (nan? y) d y)))))
  ([M]
   (replace-nan! M 0)))

(defn replace-non-numeric!
  "Replace `##NaN` and `+-##Inf` in-place (with default 0)."
  ([M x]
   (let [d (double x)]
     (alter! M (fn ^double [^double y] (if (numeric? y) y d)))))
  ([M]
   (replace-nan! M 0)))

(defn get-zoh
  "Return the first non-`##NaN` numeric value found by decrementing column `j`.
  If one is not found, return 0."
  [M i j]
  (if (< j 0)
    0
    (let [x (entry M i j)]
      (if (Double/isNaN x)
        (get-zoh M i (dec j))
        x))))

(defn zoh!
  "Set elements to the first non-`##NaN` numeric value found by decrementing column `j`.
  If one is not found, use 0."
  ;TODO: performance here needs improvement
  ([M r c]
   (entry! M r c (get-zoh M r c)))
  ([M r]
   (for [c (range (ncols M))]
     (last
      (zoh! M r c))))
  ([M]
   (last
    (for [r (range (mrows M))
          c (range (ncols M))]
      (zoh! M r c)))))

(defn hcat
  "Concatenate A and B (and D...) to a new matrix [A B]."
  ([A B]
   (let [C (dge (mrows A) (+ (ncols A) (ncols B)))]
     (copy! A (submatrix C (mrows A) (ncols A)))
     (copy! B (submatrix C 0 (ncols A) (mrows B) (ncols B)))
     C))
  ([A B & D]
   (reduce hcat (hcat A B) D)))

(defn vcat
  "Concatenate A and B (and D...) to a new matrix [A' B']'."
  ([A B]
   (let [C (dge (+ (mrows A) (mrows B)) (ncols A))]
     (copy! A (submatrix C (mrows A) (ncols A)))
     (copy! B (submatrix C (mrows A) 0 (mrows B) (ncols B)))
     C))
  ([A B & D]
   (reduce vcat (vcat A B) D)))

(defn col-vector
  "Copy as a column vector."
  [v]
  (dge (dim v) 1 v))

(defn row-vector
  "Copy as a row vector."
  [v]
  (dge 1 (dim v) v))

(defn last-cols
  "Last `n` (1) columns of matrix `M`."
  ([M n]
   (submatrix M 0 (- (ncols M) n) (mrows M) n))
  ([M]
   (last-cols M 1)))

(defn last-rows
  "Last `n` rows of matrix `M`."
  ([M n]
   (submatrix M (- (mrows M) n) 0 n (ncols M)))
  ([M]
   (last-rows M 1)))

(defn take-cols
  "Take `q` columns starting from `p` or first `n` columns of matrix `M`."
  ([M p q]
   (submatrix M 0 p (mrows M) q))
  ([M n]
   (take-cols M 0 n)))

(defn take-rows
  "Take `q` rows starting from `p` or first `n` rows of matrix `M`."
  ([M p q]
   (submatrix M p 0 q (ncols M)))
  ([M n]
   (take-rows M 0 n)))

(defn drop-cols
  "All columns but the first (1 or `n`)."
  ([M]
   (drop-cols M 1))
  ([M n]
   (submatrix M 0 n (mrows M) (- (ncols M) n))))

(defn drop-rows
  "All rows but the first (1 or `n`)."
  ([M]
   (drop-rows M 1))
  ([M n]
   (submatrix M n 0 (- (mrows M) n) (ncols M))))

(defn drop-last-cols
  "All columns but the last (1 or `n`)."
  ([M]
   (drop-last-cols M 1))
  ([M n]
   (submatrix M 0 0 (mrows M) (- (ncols M) n))))

(defn drop-last-rows
  "All rows but the last (1 or `n`)."
  ([M]
   (drop-last-rows M 1))
  ([M n]
   (submatrix M 0 0 (- (mrows M) n) (ncols M))))

(defn col-diff
  "Elements are the diff of M from its previous column."
  [M]
  (axpy -1.0 (hcat (dge (mrows M) 1) (take-cols M (- (ncols M) 1))) M))
