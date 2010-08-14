#lang scheme

#| matrix-test.ss: Test suite for the matrix.ss module.
Copyright (C) 2007 Will M. Farr

This program is free software; you can redistribute it and/or modify it
under the terms of the GNU General Public License as published by the
Free Software Foundation; either version 2 of the License, or (at your
option) any later version.

This program is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
General Public License for more details.

You should have received a copy of the GNU General Public License along
with this program; if not, write to the Free Software Foundation, Inc.,
51 Franklin St, Fifth Floor, Boston, MA 02110-1301, USA.

|#


(require (planet "test.ss" ("schematics" "schemeunit.plt" 2))
         (planet "text-ui.ss" ("schematics" "schemeunit.plt" 2))
         (lib "42ref.ss" "srfi")
         (lib "4.ss" "srfi")
         "matrix.ss")

(define-simple-check (check-close? eps a b)
  (< (abs (- a b)) (abs eps)))

(define-simple-check (check-mnorm-close? eps m1 m2)
  (< (abs (- (matrix-norm m1)
             (matrix-norm m2)))
     (abs eps)))

(provide matrix-test-suite)

(define matrix-test-suite
  (test-suite
   "matrix.ss test suite"
   (test-case
    "basic matrix operations"
    (let ((m (matrix-ec 3 3 (:range i 9) (random)))
          (i (random 3))
          (j (random 3))
          (elt (random)))
      (matrix-set! m i j elt)
      (check-equal? (matrix-ref m i j) elt)))
   (test-case
    "column-major order"
    (let ((m (matrix-ec 3 3 (:range i 9) i)))
      (check-equal? (matrix-ref m 0 0) 0.0)
      (check-equal? (matrix-ref m 1 0) 1.0)
      (check-equal? (matrix-ref m 0 1) 3.0))
    (let ((m (matrix 2 2 1 2 3 4)))
      (check-equal? (matrix-ref m 1 0) 2.0)
      (check-equal? (matrix-ref m 0 1) 3.0)))
   (test-case
    "add, subtract and scale matrix"
    (let ((m1 (matrix-ec 3 3 (:range i 9) (random)))
          (m2 (matrix-ec 3 3 (:range i 9) (random))))
      (check-mnorm-close? 1e-10 (matrix-add m1 m1) (matrix-scale m1 2))
      (check-mnorm-close? 1e-10 (matrix-scale m2 2) (matrix-scale m2 2))
      (check-mnorm-close? 1e-10 (matrix-sub m1 m2) (matrix-add m1 (matrix-scale m2 -1)))))
   (test-case
    "norm, inverse, and matrix-mul"
    (let ((m (matrix-ec 11 11 (:range i 121) (random)))
          (m2 (matrix-ec 11 11 (:range i 121) (random))))
      (check-close? 1e-10
                    (matrix-norm (matrix-sub (matrix-identity 11)
                                             (matrix-mul (matrix-inverse m) m)))
                    0)
      (check-close? 1e-10
                    (matrix-norm (matrix-sub (matrix-identity 11)
                                             (matrix-mul (matrix-inverse m2) m2)))
                    0)))
   (test-case
    "transpose and vector-matrix multiplication"
    (let ((m (matrix-ec 10 10 (:range i 100) (random)))
          (v (list->f64vector (list-ec (:range i 10) (random)))))
      (let ((mv (matrix-f64vector-mul m v))
            (vm (f64vector-matrix-mul v (matrix-transpose m))))
        (do-ec (:range i 10)
          (check-close? 1e-10 (f64vector-ref mv i) (f64vector-ref vm i))))))
   (test-case
    "matrix-solve"
    (let* ((m (matrix-ec 10 10 (:range i 100) (random)))
           (b (list->f64vector (list-ec (:range i 10) (random))))
           (x (matrix-solve m b))
           (mx (matrix-f64vector-mul m x)))
      (do-ec (:range i 10)
        (check-close? 1e-10 (f64vector-ref b i) (f64vector-ref mx i)))))
   (test-case
    "matrix-solve-least-squares"
    (let* ((m (matrix-ec 10 11 (:range i 110) (random)))
           (b (list->f64vector (list-ec (:range i 10) (random)))))
      (do-ec (:range i 10)
        (matrix-set! m i 10 0.0))
      (let ((x (matrix-solve-least-squares m b)))
        (let ((mx (matrix-f64vector-mul m x)))
          (do-ec (:range i 10)
            (check-close? 1e-8 (f64vector-ref mx i) (f64vector-ref b i)))))))
   (test-case
    "Noel Welsh's seg fault"
    (f64vector-matrix-mul (f64vector 1. .0 .0 .0 .0 .0 .0 .0 .0 .0) (make-matrix 10 10 1.)))
   (test-case
    ":matrix test"
    (let ((m (matrix-ec 3 3 (:range i 3) (:range j 3) (+ i j))))
      (let ((elts (list-ec (:matrix x m) x)))
        (check-equal? elts (list 0.0 1.0 2.0 1.0 2.0 3.0 2.0 3.0 4.0)))))))