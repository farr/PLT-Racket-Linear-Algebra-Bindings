#lang racket

#|  matrix.rkt: Matrices and matrix operations using BLAS and LAPACK
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
         (lib "etc.ss")
         (except-in (lib "contract.ss") ->)
         (rename-in (lib "contract.ss") (-> ->/c))
         (except-in (lib "42ref.ss" "srfi") :)
         (except-in (planet "srfi-4-comprehensions.ss" ("wmfarr" "srfi-4-comprehensions.plt")) :)
         "blas-lapack.ss"
         "vector.ss"
         racket/fixnum
         racket/flonum)

(define (list/length/c n)
  (flat-named-contract
   (format "list of length ~a" n)
   (lambda (l) (= (length l) n))))

(define (matrix-multiplication-compatible/c m)
  (flat-named-contract
   (format "compatible for multiplication by a ~a by ~a matrix" (matrix-rows m) (matrix-cols m))
   (lambda (m2)
     (= (matrix-cols m) (matrix-rows m2)))))

(define (matrix-same-dimensions/c m)
  (flat-named-contract
   (format "~a by ~a matrix" (matrix-rows m) (matrix-cols m))
   (lambda (m2)
     (and (= (matrix-rows m) (matrix-rows m2))
          (= (matrix-cols m) (matrix-cols m2))))))

(define (matrix-valid-row-index/c m)
  (let ((r (matrix-rows m)))
    (flat-named-contract
     (format "valid row index for a ~a by ~a matrix" r (matrix-cols m))
     (lambda (i) (and (>= i 0)
                      (< i r))))))

(define (matrix-valid-col-index/c m)
  (let ((c (matrix-cols m)))
    (flat-named-contract
     (format "valid column index for ~a by ~a matrix" (matrix-rows m) c)
     (lambda (j) (and (>= j 0)
                      (< j c))))))

(define (matrix-col-vector-compatible/c m)
  (let ((c (matrix-cols m)))
    (flat-named-contract
     (format "column vector of length ~a" c)
     (lambda (v) (= (f64vector-length v) c)))))

(define (matrix-row-vector-compatible/c m)
  (let ((r (matrix-rows m)))
    (flat-named-contract
     (format "row vector of length ~a" r)
     (lambda (v) (= (f64vector-length v) r)))))

(define matrix-square/c
  (flat-named-contract
   "square matrix"
   (lambda (m) (= (matrix-rows m) (matrix-cols m)))))

(define-struct matrix
  (ptr rows cols) #:transparent)

(provide matrix? matrix-multiplication-compatible/c matrix-valid-row-index/c
         matrix-valid-col-index/c matrix-square/c matrix-same-dimensions/c
         matrix-col-vector-compatible/c matrix-row-vector-compatible/c
         matrix-ec :matrix in-matrix
         _matrix
         struct:matrix
         (rename-out (matrix s:matrix))
         exn:singular-matrix)

(provide/contract
 (rename my-make-matrix make-matrix
         (->/c natural-number/c natural-number/c number? matrix?))
 (rename my-matrix matrix
         (->r ((i natural-number/c)
               (j natural-number/c))
              elts (and/c (listof number?)
                          (list/length/c (* i j)))
              matrix?))
 (matrix-rows (->/c matrix? natural-number/c))
 (matrix-cols (->/c matrix? natural-number/c))
 (matrix-ref (->r ((m matrix?)
                   (i (and/c natural-number/c
                             (matrix-valid-row-index/c m)))
                   (j (and/c natural-number/c
                             (matrix-valid-col-index/c m))))
                  number?))
 (matrix-set! (->r ((m matrix?)
                    (i (and/c natural-number/c
                              (matrix-valid-row-index/c m)))
                    (j (and/c natural-number/c
                              (matrix-valid-col-index/c m)))
                    (x number?))
                   any))
 (matrix-add (->r ((m1 matrix?)
                   (m2 (and/c matrix?
                              (matrix-same-dimensions/c m1))))
                  matrix?))
 (matrix-sub (->r ((m1 matrix?)
                   (m2 (and/c matrix?
                              (matrix-same-dimensions/c m1))))
                  matrix?))
 (matrix-scale (->/c matrix? number? matrix?))
 (matrix-mul (->r ((m1 matrix?)
                   (m2 (and/c matrix?
                              (matrix-multiplication-compatible/c m1))))
                  matrix?))
 (matrix-f64vector-mul (->r ((m matrix?)
                             (v (and/c f64vector?
                                       (matrix-col-vector-compatible/c m))))
                            f64vector?))
 (f64vector-matrix-mul (->r ((v (and/c f64vector?
                                       (matrix-row-vector-compatible/c m)))
                             (m matrix?))
                            f64vector?))
 (matrix-inverse (->/c (and/c matrix? matrix-square/c) matrix?))
 (matrix-norm (->/c matrix? number?))
 (matrix-identity (->/c natural-number/c matrix?))
 (matrix-transpose (->/c matrix? matrix?))
 (matrix-solve (->r ((m matrix-square/c)
                     (v (and/c f64vector?
                               (matrix-row-vector-compatible/c m))))
                    (and/c f64vector? (matrix-col-vector-compatible/c m))))
 (matrix-solve-many (->r ((m1 matrix-square/c)
                          (m2 (and/c matrix?
                                     (matrix-multiplication-compatible/c m1))))
                         matrix?))
 (matrix-solve-least-squares
  (->r ((m matrix?)
        (b (and/c f64vector?
                  (matrix-row-vector-compatible/c m))))
       (and/c f64vector? (matrix-col-vector-compatible/c m))))
 (eigensystem 
  (->/c matrix-square/c (values f64vector? f64vector? matrix? matrix?)))
 (eigenvalues->vector
  (->/c f64vector? f64vector? (vectorof number?))))

(unsafe!)

(define-struct (exn:singular-matrix exn)
  (elt) #:transparent)

(define my-make-matrix
  (case-lambda
    ((rows cols)
     (let* ((n (* rows cols))
            (p (malloc n _double 'atomic)))
       (memset p 0 n _double)
       (make-matrix p rows cols)))
    ((rows cols elt)
     (let* ((m (my-make-matrix rows cols))
            (p (matrix-ptr m)))
       (do-ec (:range i (* rows cols))
         (ptr-set! p _double* i elt))
       m))))

(define (my-matrix i j . elts)
  (let* ((m (my-make-matrix i j))
         (p (matrix-ptr m)))
    (do-ec (:parallel (:range k (* i j))
                      (:list elt elts))
      (ptr-set! p _double* k elt))
    m))

(define _matrix*
  (make-ctype _pointer matrix-ptr
              (lambda (x)
                (error '_matrix
                       "cannot convert C output to _matrix"))))

(define-fun-syntax _matrix
  (syntax-id-rules (i o io)
    ((_matrix i)
     _matrix*)
    ((_matrix o rows cols)
     (type: _pointer
            pre: (let* ((n (* rows cols))
                        (p (malloc n _double 'atomic)))
                   (memset p 0 n _double)
                   p)
            post: (p => (make-matrix p rows cols))))
    ((_matrix io)
     (type: _pointer
            bind: m
            pre: (m => (matrix-ptr m))
            post: m))
    (_matrix _matrix*)))

(define (matrix-ptr-index m i j)
  (+ i (* j (matrix-rows m))))

(define (matrix-ref m i j)
  (ptr-ref (matrix-ptr m) _double (matrix-ptr-index m i j)))

(define (matrix-set! m i j elt)
  (ptr-set! (matrix-ptr m) _double* (matrix-ptr-index m i j) elt))

(define (matrix-length m)
  (* (matrix-cols m)
     (matrix-rows m)))

(define matrix-copy
  (get-ffi-obj 'cblas_dcopy *blas*
               (_fun (m) ::
                     (_int = (matrix-length m))
                     (m : _matrix) (_int = 1)
                     (m-out : (_matrix o (matrix-rows m) (matrix-cols m))) (_int = 1) ->
                     _void ->
                     m-out)))

(define matrix-add
  (get-ffi-obj 'cblas_daxpy *blas*
               (_fun (m1 m2) ::
                     (_int = (matrix-length m1))
                     (_double = 1.0)
                     (_matrix = m1) (_int = 1)
                     (m-out : _matrix = (matrix-copy m2)) (_int = 1) ->
                     _void ->
                     m-out)))

(define matrix-sub
  (get-ffi-obj 'cblas_daxpy *blas*
               (_fun (m1 m2) ::
                     (_int = (matrix-length m1))
                     (_double = -1.0)
                     (_matrix = m2) (_int = 1)
                     (m-out : _matrix = (matrix-copy m1)) (_int = 1) ->
                     _void ->
                     m-out)))

(define matrix-scale
  (get-ffi-obj 'cblas_dscal *blas*
               (_fun (m s) ::
                     (_int = (matrix-length m))
                     (_double* = s)
                     (m-out : _matrix = (matrix-copy m)) (_int = 1) ->
                     _void ->
                     m-out)))

(define matrix-mul
  (get-ffi-obj 'cblas_dgemm *blas*
               (_fun (m1 m2) ::
                     (_cblas-order = 'col-major)
                     (_cblas-transpose = 'no-trans)
                     (_cblas-transpose = 'no-trans)
                     (_int = (matrix-rows m1))
                     (_int = (matrix-cols m2))
                     (_int = (matrix-rows m2))
                     (_double = 1.0)
                     (_matrix = m1)
                     (_int = (matrix-rows m1))
                     (_matrix = m2)
                     (_int = (matrix-rows m2))
                     (_double = 0.0)
                     (m-out : (_matrix o (matrix-rows m1) (matrix-cols m2)))
                     (_int = (matrix-rows m1)) ->
                     _void ->
                     m-out)))

(define matrix-f64vector-mul
  (get-ffi-obj 'cblas_dgemv *blas*
               (_fun (m v) ::
                     (_cblas-order = 'col-major)
                     (_cblas-transpose = 'no-trans)
                     (_int = (matrix-rows m))
                     (_int = (matrix-cols m))
                     (_double = 1.0)
                     (_matrix = m)
                     (_int = (matrix-rows m))
                     (_f64vector = v)
                     (_int = 1)
                     (_double = 0.0)
                     (v-out : (_f64vector o (matrix-rows m)))
                     (_int = 1) ->
                     _void ->
                     v-out)))

(define f64vector-matrix-mul
  (get-ffi-obj 'cblas_dgemv *blas*
               (_fun (v m) ::
                     (_cblas-order = 'col-major)
                     (_cblas-transpose = 'trans)
                     (_int = (matrix-rows m))
                     (_int = (matrix-cols m))
                     (_double = 1.0)
                     (_matrix = m)
                     (_int = (matrix-rows m))
                     (_f64vector = v)
                     (_int = 1)
                     (_double = 0.0)
                     (v-out : (_f64vector o (matrix-cols m)))
                     (_int = 1) ->
                     _void ->
                     v-out)))

(define matrix-lu-decomp
  (get-ffi-obj 'dgetrf_ *lapack*
               (_fun (m) ::
                     ((_ptr i _int) = (matrix-rows m))
                     ((_ptr i _int) = (matrix-cols m))
                     (m-out : _matrix = (matrix-copy m))
                     ((_ptr i _int) = (matrix-rows m))
                     (ipiv : (_u32vector o (matrix-rows m)))
                     (_ptr o _int) ->
                     _void ->
                     (values m-out ipiv))))

(define dgetri-lwork
  (get-ffi-obj 'dgetri_ *lapack*
               (_fun (m ipiv) ::
                     (n : (_ptr i _int) = (matrix-rows m))
                     (_matrix = m)
                     ((_ptr i _int) = n)
                     (_u32vector = ipiv)
                     (lwork : (_ptr o _double))
                     ((_ptr i _int) = -1)
                     (res : (_ptr o _int)) ->
                     _void ->
                     (values (inexact->exact (round lwork)) res))))

(define dgetri/lwork
  (get-ffi-obj 'dgetri_ *lapack*
               (_fun (m ipiv lwork) ::
                     (n : (_ptr i _int) = (matrix-rows m))
                     (m-out : _matrix = (matrix-copy m))
                     ((_ptr i _int) = n)
                     (_u32vector = ipiv)
                     (_f64vector o lwork)
                     ((_ptr i _int) = lwork)
                     (_ptr o _int) ->
                     _void ->
                     m-out)))

(define (matrix-inverse m)
  (let-values (((m-lu ipiv)
                (matrix-lu-decomp m)))
    (let-values (((lwork res)
                  (dgetri-lwork m-lu ipiv)))
      (dgetri/lwork m-lu ipiv lwork))))

(define matrix-norm
  (get-ffi-obj 'cblas_dnrm2 *blas*
               (_fun (m) ::
                     (_int = (matrix-length m))
                     (_matrix = m)
                     (_int = 1) ->
                     _double)))

(define (matrix-transpose m)
  (let ((r (matrix-rows m))
        (c (matrix-cols m)))
    (matrix-ec c r (:range j r) (:range i c) (matrix-ref m j i))))

(define matrix-solve-many
  (get-ffi-obj 'dgesv_ *lapack*
               (_fun (m b) ::
                     ((_ptr i _int) = (matrix-rows m))
                     ((_ptr i _int) = (matrix-cols b))
                     (_matrix = (matrix-copy m))
                     ((_ptr i _int) = (matrix-rows m))
                     (_u32vector o (matrix-rows m))
                     (x : _matrix = (matrix-copy b))
                     ((_ptr i _int) = (matrix-rows b))
                     (info : (_ptr o _int)) ->
                     _void ->
                     (if (> info 0)
                         (raise (make-exn:singular-matrix (format "singular matrix (elt ~a in U = 0) in matrix-solve-many"
                                                                  info)
                                                          (current-continuation-marks)
                                                          info))
                         x))))

(define matrix-solve
  (get-ffi-obj 'dgesv_ *lapack*
               (_fun (m v) ::
                     ((_ptr i _int) = (matrix-rows m))
                     ((_ptr i _int) = 1)
                     (_matrix = (matrix-copy m))
                     ((_ptr i _int) = (matrix-rows m))
                     (_u32vector o (matrix-rows m))
                     (x : _f64vector = (f64vector-copy v))
                     ((_ptr i _int) = (f64vector-length v))
                     (info : (_ptr o _int)) ->
                     _void ->
                     (if (> info 0)
                         (raise (make-exn:singular-matrix (format "singular matrix (elt ~a in U = 0) in matrix-solve"
                                                                  info)
                                                          (current-continuation-marks)
                                                          info))
                         x))))

(define dgelsd-lwork
  (get-ffi-obj 'dgelsd_ *lapack*
               (_fun (m v) ::
                     ((_ptr i _int) = (matrix-rows m))
                     ((_ptr i _int) = (matrix-cols m))
                     ((_ptr i _int) = 1)
                     (_matrix = m)
                     ((_ptr i _int) = (matrix-rows m))
                     (_f64vector = v)
                     ((_ptr i _int) = (max (matrix-rows m) (matrix-cols m)))
                     (_ptr o _double)
                     ((_ptr i _double) = 1e-8)
                     (_ptr o _int)
                     (work : (_ptr o _double))
                     ((_ptr i _int) = -1)
                     (iwork : (_ptr o _int))
                     (_ptr o _int) ->
                     _void ->
                     (values (inexact->exact (round work))
                             iwork))))

(define dgelsd/lwork
  (get-ffi-obj 'dgelsd_ *lapack*
               (_fun (m v lwork iwork) ::
                     ((_ptr i _int) = (matrix-rows m))
                     ((_ptr i _int) = (matrix-cols m))
                     ((_ptr i _int) = 1)
                     (_matrix = (matrix-copy m))
                     ((_ptr i _int) = (matrix-rows m))
                     (b : _f64vector = (let ((nb (max (matrix-rows m) (matrix-cols m)))
                                             (nv (f64vector-length v)))
                                         (f64vector-of-length-ec nb (:range i nb)
                                                                 (if (< i nv)
                                                                     (f64vector-ref v i)
                                                                     0.0))))
                     ((_ptr i _int) = (max (matrix-rows m) (matrix-cols m)))
                     (_f64vector o (min (matrix-rows m) (matrix-cols m)))
                     ((_ptr i _double) = -1.0)
                     (_ptr o _int)
                     (_f64vector o lwork)
                     ((_ptr i _int) = lwork)
                     (_u64vector o iwork)
                     (_ptr o _int) ->
                     _void ->
                     (f64vector-of-length-ec (matrix-cols m) (:range i (matrix-cols m))
                                             (f64vector-ref b i)))))
;; Some sort of segfault---removed from the provided functions until fixed. 
(define (matrix-solve-least-squares m b)
  (error 'matrix-solve-least-squares "there is a hard-to-find segfault in matrix-solve-least-squares; contact w-farr@northwestern.edu and ask him to track it down if you *really* need access to least-squares solvers")

  #;(let-values (((lwork iwork) (dgelsd-lwork m b)))
    (printf "Found lwork, iwork: ~a, ~a~%" lwork iwork)
    (let ((answer (dgelsd/lwork m b lwork lwork))) ; The final lwork should be iwork
      (printf "Returned from dgelsd/lwork.~%")
      answer)))

(define (matrix-identity n)
  (let ((m (my-make-matrix n n 0.0)))
    (do-ec (:range i n) (matrix-set! m i i 1.0))
    m))

(define-syntax matrix-ec
  (syntax-rules ()
    ((matrix-ec rrows ccols etc ...)
     (let ((rows rrows)
           (cols ccols))
       (apply my-matrix rows cols (list-ec etc ...))))))

(define-syntax :matrix
  (syntax-rules (index)
    ((:matrix cc var arg)
     (:matrix cc var (index i j) arg))
    ((:matrix cc var (index i j) arg)
     (:do cc
          (let ((m arg)
                (rows #f)
                (cols #f))
            (set! rows (matrix-rows m))
            (set! cols (matrix-cols m)))
          ((i 0) (j 0))
          (< j cols)
          (let ((i+1 (+ i 1))
                (j+1 (+ j 1))
                (wrapping? #f)
                (var (matrix-ref m i j)))
            (set! wrapping? (>= i+1 rows)))
          #t
          ((if wrapping? 0 i+1)
           (if wrapping? j+1 j))))))

(define ptr->matrix make-matrix)
(provide* (unsafe ptr->matrix))

(define (in-matrix* m)
  (let ((nrows (matrix-rows m))
        (ncols (matrix-cols m)))
    (make-do-sequence
     (lambda ()
       (values (lambda (ij) (matrix-ref m (car ij) (cdr ij)))
               (lambda (ij)
                 (match ij
                   ((cons i j)
                    (let ((ii (add1 i)))
                      (if (< ii nrows)
                          (cons ii j)
                          (cons 0 (add1 j)))))))
               (cons 0 0)
               (lambda (ij) 
                 (match ij
                   ((cons i j)
                    (not (>= j ncols)))))
               (lambda args #t)
               (lambda args #t))))))

(define-sequence-syntax in-matrix
  (lambda () (syntax in-matrix*))
  (lambda (stx)
    (syntax-case stx ()
      (((x-id) (in-matrix mexpr))
       (syntax/loc stx
         ((m) (:do-in (((m) mexpr))
                      #t
                      ((nrows (matrix-rows m))
                       (ncols (matrix-cols m))
                       (i 0)
                       (j 0))
                      (fx< j ncols)
                      (((x-id) (matrix-ref m i j))
                       ((ii) (add1 i)))
                      #t
                      #t
                      (nrows ncols (if (fx>= ii nrows) 0 ii) (if (fx>= ii nrows) (add1 j) j)))))))))

(define dgeev-lwork
  (get-ffi-obj 
   'dgeev_ *lapack*
   (_fun (m) ::
         (_string = "N")
         (_string = "N")
         ((_ptr i _int) = (matrix-rows m))
         (_matrix = m)
         ((_ptr i _int) = (matrix-rows m))
         (_ptr o _double)
         (_ptr o _double)
         (_ptr o _double)
         ((_ptr i _int) = 1)
         (_ptr o _double)
         ((_ptr i _int) = 1)
         (work : (_ptr o _double))
         ((_ptr i _int) = -1)
         (_ptr o _int) -> 
         _void -> 
         (inexact->exact (round work)))))

(define dgeev
  (get-ffi-obj
   'dgeev_ *lapack*
   (_fun (m lwork) :: 
         (_string = "V")
         (_string = "V")
         ((_ptr i _int) = (matrix-rows m))
         (_matrix = (matrix-copy m))
         ((_ptr i _int) = (matrix-rows m))
         (wr : (_f64vector o (matrix-rows m)))
         (wi : (_f64vector o (matrix-rows m)))
         (vl : _matrix = (matrix-copy m))
         ((_ptr i _int) = (matrix-rows m))
         (vr : _matrix = (matrix-copy m))
         ((_ptr i _int) = (matrix-rows m))
         (_f64vector o lwork)
         ((_ptr i _int) = lwork)
         (_ptr o _int) -> 
         _void -> 
         (values wr wi vl vr))))

(define (eigensystem m)
  (let ((lwork (dgeev-lwork m)))
    (dgeev m lwork)))

(define (eigenvalues->vector er ei)
  (list->vector
   (for/list ((i (in-range (f64vector-length er))))
     (if (= 0.0 (f64vector-ref ei i))
         (f64vector-ref er i)
         (make-rectangular (f64vector-ref er i) (f64vector-ref ei i))))))

(define-unsafer matrix-unsafe!)