#!/bin/sbcl --script

(defparameter *l1-time* 5)
(defparameter *l2-time* 12)
(defparameter *main-mem-time* 1656/10)
(defparameter *rows-to-check* 128)
(defparameter *bytes-per-mat-ind* 4)
(defparameter *bytes-per-mat-val* 8)
(defparameter *bytes-per-vect* 8)
(defparameter *inds-decode-time* 0)
(defparameter *vals-decode-time* 0)
(defparameter *vect-decode-time* 0)
(defparameter *vect-encode-time* 0)

(defmethod fetch ((obj (eql :mat-inds)) (i integer))
  (if (/= (floor (* (1- i) *bytes-per-mat-ind*) 64)
          (floor (* i *bytes-per-mat-ind*) 64))
    *main-mem-time*
    *l1-time*))
(defmethod fetch ((obj (eql :mat-vals)) (i integer))
  (if (/= (floor (* (1- i) *bytes-per-mat-val*) 64)
          (floor (* i *bytes-per-mat-val*) 64))
    *main-mem-time*
    *l1-time*))
(defmethod fetch ((obj (eql :vect)) (i integer))
  (cond
    ; 2/3rds of values were used by the previous index
    ((/= (mod i 3) 2) *l1-time*)
    ; 2/9ths of values were used by y-1
    ((< (/ i 3) 6) *l2-time*)
    ; 1/9th of values are being used for the first time
    (t (if (/= (floor (* (1- (/ i 27)) *bytes-per-vect*) 64)
               (floor (* (/ i 27) *bytes-per-vect*) 64))
         *main-mem-time*
         *l1-time*))))

(defun 1-row ()
  (+
    (/ (loop :for i :from 0 :below (* *rows-to-check* 27)
               :for inds-fetch-time = (fetch :mat-inds i)
               :for vals-fetch-time = (fetch :mat-vals i)
               :for vect-fetch-time = (fetch :vect i)
               :for total-vals-time = (+ vals-fetch-time *vals-decode-time*)
               :for total-vect-time = (+ inds-fetch-time *inds-decode-time* vect-fetch-time *vect-decode-time*)
          :summing (min total-vals-time total-vect-time))
      *rows-to-check*)
    *vect-encode-time*))

(defun 1-row-with-props (ind-size val-size vect-size ind-decode val-decode vect-decode vect-encode)
  (let ((*bytes-per-mat-ind* ind-size)
        (*bytes-per-mat-val* val-size)
        (*bytes-per-vect* vect-size)
        (*inds-decode-time* ind-decode)
        (*vals-decode-time* val-decode)
        (*vect-decode-time* vect-decode)
        (*vect-encode-time* vect-encode))
    (1-row)))

(defparameter *baseline-time* (1-row))

(defparameter options
  (loop :for ind-size :from 1 :to 4
    :append (loop :for val-size :from 1 :to 8
              :append (loop :for vect-size :from 8 :to 8
                        :append (loop :for vect-decode :downfrom 0 :to 0 :by 1
                                  :append (loop :for vect-encode :downfrom 0 :to 0 :by 1
                                            :append (loop :for ind-decode :downfrom 200 :to 0 :by 2
                                                      :append (loop :for val-decode :downfrom 100 :to 0 :by 2
                                                                :if (>= *baseline-time* (1-row-with-props ind-size val-size vect-size ind-decode val-decode vect-decode vect-encode))
                                                                  :collect (list ind-size val-size vect-size vect-decode vect-encode ind-decode val-decode)))))
                        :do (format t "~&;   vect size ~S complete~%" vect-size))
              :do (format t "~&;  val size ~S complete~%" val-size))
    :do (format t "~&; ind size ~S complete ~%" ind-size)))


(defparameter filtered-options
  (remove-if (lambda (lst)
               (some (lambda (lst2)
                       (and (not (eql lst lst2))
                            (= (first lst2) (first lst))
                            (= (second lst2) (second lst))
                            (= (third lst2) (third lst))
                            (>= (nth 3 lst2) (nth 3 lst))
                            (>= (nth 4 lst2) (nth 4 lst))
                            (>= (nth 5 lst2) (nth 5 lst))
                            (>= (nth 6 lst2) (nth 6 lst))))
                    options))
             options))

(defparameter filtered-options2
  (loop :with result = nil
        :for lst :in options
    :unless (or (some (lambda (lst2)
                        (and (not (eql lst lst2))
                             (= (first lst2) (first lst))
                             (= (second lst2) (second lst))
                             (= (third lst2) (third lst))
                             (>= (nth 3 lst2) (nth 3 lst))
                             (>= (nth 4 lst2) (nth 4 lst))
                             (>= (nth 5 lst2) (nth 5 lst))
                             (>= (nth 6 lst2) (nth 6 lst))))
                      result)
                (some (lambda (lst2)
                        (and (not (eql lst lst2))
                             (= (first lst2) (first lst))
                             (= (second lst2) (second lst))
                             (= (third lst2) (third lst))
                             (>= (nth 3 lst2) (nth 3 lst))
                             (>= (nth 4 lst2) (nth 4 lst))
                             (>= (nth 5 lst2) (nth 5 lst))
                             (>= (nth 6 lst2) (nth 6 lst))))
                      options))
      :do (push lst result)
    :finally (return (print result))))

; To get the important information by running this as a script
(format t "Combinations of (bytesPerVectorValue vectDecodeTime vectEncodeTime) that \"perform\" faster than the baseline:
Note that any combination that has all three values less than or equal to the respective values of another successful combination is ommitted~%")
(print filtered-options)
(format t "~%")


(defun inds-performance (ind-size)
    (loop :for ind-decode :downfrom 20 :to 0
      :collect (1-row-with-props ind-size 8 8 ind-decode 0 0 0)))

(loop :for val-decode :from 0 :to 20
  :do (format t "~A: " val-decode)
  :do (loop :for val-size :from 1 :to 8
        :do (format t "  ~A" (float (1-row-with-props 4 val-size 8 0 val-decode 0 0))))
  :do (format t "~%"))
