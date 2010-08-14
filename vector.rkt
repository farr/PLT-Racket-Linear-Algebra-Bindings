#lang scheme

#|  vector.ss: f64vectors in linear algebra.
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

(require (lib "foreign.ss")
         (except-in (lib "contract.ss") ->)
         "blas-lapack.ss"
         (rename-in (lib "contract.ss") (-> ->/c))
         (lib "4.ss" "srfi"))

(unsafe!)

(provide (all-from-out (lib "4.ss" "srfi")) f64vector-same-length/c)

(define (f64vector-same-length/c v1)
  (let ((n (f64vector-length v1)))
    (flat-named-contract
     (format "length ~a f64vector" (f64vector-length v1))
     (lambda (v2) (= (f64vector-length v2) n)))))

(provide/contract
 (f64vector-norm (->/c f64vector? (>=/c 0)))
 (f64vector-copy (->/c f64vector? f64vector?))
 (f64vector-scale (->/c f64vector? real? f64vector?))
 (f64vector-add (->r ((v1 f64vector?)
                      (v2 (and/c f64vector?
                                 (f64vector-same-length/c v1))))
                     f64vector?))
 (f64vector-sub (->r ((v1 f64vector?)
                      (v2 (and/c f64vector?
                                 (f64vector-same-length/c v1))))
                     f64vector?))
 (f64vector-dot (->/c f64vector? f64vector? number?)))

(define f64vector-norm
  (get-ffi-obj 'cblas_dnrm2 *blas*
               (_fun (v) ::
                     (_int = (f64vector-length v))
                     (_f64vector = v)
                     (_int = 1) ->
                     _double)))

(define f64vector-copy
  (get-ffi-obj 'cblas_dcopy *blas*
               (_fun (v) ::
                     (_int = (f64vector-length v))
                     (_f64vector = v)
                     (_int = 1)
                     (v-out : (_f64vector o (f64vector-length v)))
                     (_int = 1) ->
                     _void ->
                     v-out)))

(define f64vector-scale
  (get-ffi-obj 'cblas_dscal *blas*
               (_fun (v s) ::
                     (_int = (f64vector-length v))
                     (_double* = s)
                     (v-out : _f64vector = (f64vector-copy v))
                     (_int = 1) ->
                     _void ->
                     v-out)))

(define f64vector-add
  (get-ffi-obj 'cblas_daxpy *blas*
               (_fun (v1 v2) ::
                     (_int = (f64vector-length v1))
                     (_double = 1.0)
                     (_f64vector = v1)
                     (_int = 1)
                     (v-out : _f64vector = (f64vector-copy v2))
                     (_int = 1) ->
                     _void ->
                     v-out)))

(define f64vector-sub
  (get-ffi-obj 'cblas_daxpy *blas*
               (_fun (v1 v2) ::
                     (_int = (f64vector-length v1))
                     (_double = -1.0)
                     (_f64vector = v2)
                     (_int = 1)
                     (v-out : _f64vector = (f64vector-copy v1))
                     (_int = 1) ->
                     _void ->
                     v-out)))

(define f64vector-dot
  (get-ffi-obj 'cblas_ddot *blas*
               (_fun (v1 v2) ::
                     (_int = (f64vector-length v1))
                     (_f64vector = v1)
                     (_int = 1)
                     (_f64vector = v2)
                     (_int = 1) ->
                     _double)))