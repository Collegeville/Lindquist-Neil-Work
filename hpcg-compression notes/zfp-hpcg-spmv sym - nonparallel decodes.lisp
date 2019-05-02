#!/bin/sbcl --script

(defparameter *l1-time* 5)
(defparameter *l2-time* 12)
(defparameter *main-mem-time* 1656/10)
(defparameter *bytes-per-mat-ind* 4)
(defparameter *bytes-per-mat-val* 8)
(defparameter *bytes-per-vect* 8)
(defparameter *inds-decode-time* 0)
(defparameter *vals-decode-time* 0)
(defparameter *vect-decode-time* 0)
(defparameter *vect-encode-time* 0)

(defmethod fetch ((obj (eql :mat-inds)) (i integer))
  (if (= (floor (* (1- i) 64) *bytes-per-mat-ind*)
         (floor (* i 64) *bytes-per-mat-ind*))
    *main-mem-time*
    *l1-time*))
(defmethod fetch ((obj (eql :mat-vals)) (i integer))
  (if (= (floor (* (1- i) 64) *bytes-per-mat-val*)
         (floor (* i 64) *bytes-per-mat-val*))
    *main-mem-time*
    *l1-time*))
(defmethod fetch ((obj (eql :vect)) (i integer))
  (cond
    ((> (mod i 3) 0) *l1-time*)
    ((< (/ i 3) 6) *l2-time*)
    (t (+ (* (/ *bytes-per-vect* 64) *main-mem-time*)
          (* (/ (- 64 *bytes-per-vect*) 64) *l1-time*)))))

(defun 1-row ()
  (/ (+ (loop :for i :from 0 :below (* 5 27)
               :for inds-fetch-time = (fetch :mat-inds i)
               :for vals-fetch-time = (fetch :mat-vals i)
               :for vect-fetch-time = (fetch :vect i)
               :for vals-decode-start = (max vals-fetch-time (+ inds-fetch-time *inds-decode-time*))
               :for total-vals-time = (+ vals-decode-start *vals-decode-time*)
               :for vect-decode-start = (max total-vals-time (+ inds-fetch-time *inds-decode-time* vect-fetch-time))
               :for total-vect-time = (+ vect-decode-start *vect-decode-time*)
          :do (assert (<= total-vect-time (+ inds-fetch-time *inds-decode-time* vals-fetch-time *vals-decode-time* vect-fetch-time *vect-decode-time*)))
          :summing total-vect-time)
        (* 5 *vect-encode-time*))
    5))

(defun 1-row-with-props (ind-size val-size vect-size ind-decode val-decode vect-decode vect-encode)
  (let ((*bytes-per-mat-ind* ind-size)
        (*bytes-per-mat-val* val-size)
        (*bytes-per-vect* vect-size)
        (*inds-decode-time* ind-decode)
        (*vals-decode-time* val-decode)
        (*vect-decode-time* vect-decode)
        (*vect-encode-time* vect-encode))
    (1-row)))

(defparameter options
  (loop :for vect-size :from 1 :to 7
    :append (loop :for vect-decode :from 0 :to 50
            :append (loop :for vect-encode :from 0 :to 500
                    :if (> 86997/200 (1-row-with-props 4 8 vect-size 0 0 vect-decode vect-encode))
                      :collect (list vect-size vect-decode vect-encode)))
    :do (format t "~&; vect size ~S complete~%" vect-size)))

(defparameter filtered-options
  (remove-if (lambda (lst)
               (some (lambda (lst2)
                       (and (not (eql lst lst2))
                            (>= (the fixnum (first lst2)) (the fixnum (first lst)))
                            (>= (the fixnum (second lst2)) (the fixnum (second lst)))
                            (>= (the fixnum (third lst2)) (the fixnum (third lst)))))
                     options))
             options))


; To get the important information by running this as a script
(format t "Combinations of (bytesPerVectorValue vectDecodeTime vectEncodeTime) that \"perform\" faster than the baseline:
Note that any combination that has all three values less than or equal to the respective values of another successful combination is ommitted~%")
(print filtered-options)
(format t "~%")
