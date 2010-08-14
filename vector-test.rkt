#lang scheme

#|  vector-test.ss: Test suite for vector.ss
    Copyright (C) 2007 Will M. Farr <farr@mit.edu>

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License along
    with this program; if not, write to the Free Software Foundation, Inc.,
    51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
|#

(require (planet "test.ss" ("schematics" "schemeunit.plt" 2))
         "vector.ss"
         (lib "math.ss")
         (planet "srfi-4-comprehensions.ss" ("wmfarr" "srfi-4-comprehensions.plt" 1)))

(provide vector-test-suite)

(define-simple-check (check-close? eps a b)
  (< (abs (- a b)) (abs eps)))

(define vector-test-suite
  (test-suite
   "vector.ss test suite"
   (test-case
    "f64vector-norm"
    (check-close? 1e-10
                  (f64vector-norm (f64vector 1.0 2.0 3.0 4.0))
                  (sqrt 30)))
   (test-case
    "f64vector-add, f64vector-sub, f64vector-scale"
    (let ((v1 (f64vector-of-length-ec 10 (:range i 10) (random)))
          (v2 (f64vector-of-length-ec 10 (:range i 10) (random))))
      (let ((x (f64vector-sub v1 v2))
            (y (f64vector-add v1 (f64vector-scale v2 -1))))
        (do ((i 0 (add1 i)))
          ((= i (f64vector-length x)) #t)
          (check-close? 1e-10 (f64vector-ref x i) (f64vector-ref y i))))))
   (test-case
    "f64vector-dot"
    (let ((v (f64vector-of-length-ec 10 (:range i 10) (random))))
      (check-close? 1e-10 (f64vector-dot v v) (sqr (f64vector-norm v)))))))