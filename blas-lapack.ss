#lang scheme

#| blas-lapack.ss: Library locations for BLAS/LAPACK.
Copyright (C) 2007 Will M. Farr <farr@mit.edu>

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License along
with this program; if not, write to the Free Software Foundation, Inc.,
51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
|# 

#| Special thanks to Noel Welsh, who contributed the library-searching
code below. |#


(require (lib "foreign.ss")
         srfi/1)

(unsafe!)

(define-unsafer blas-lapack-unsafe!)

(provide *blas* *lapack*
         _cblas-order _cblas-transpose _cblas-uplo
         _cblas-diag _cblas-side)

;; search-paths : (listof string)
(define search-paths
  (case (system-type)
    [(macosx)
     (list
      "/System/Library/Frameworks/Accelerate.framework/Versions/Current/Frameworks/vecLib.framework/Versions/Current"
      "/System/Library/Frameworks/vecLib.framework/Versions/Current")]
    [(unix)
     (list
      "/usr/lib"
      "/home/pg/nhw/data/lib" ;;;; NB: NHW specific
      "/usr/libsse2" ;; Bug report #196 on PLaneT Issue Tracking System
      )]
    [(windows)
     (list
      "c:\\windows\\system32")]))

(define default-path "")

;; lib-blas : (listof string)
;;
;; Possible names for the BLAS library
;; 
;; libcblas.so.3gf is apparently required for Debian.  Thanks to
;; Viktor Winschel for pointing this out.
(define lib-blas
  (case (system-type)
    [(macosx) '("libBLAS")]
    [(unix)  '("libcblas" "libgslcblas" "libcblas.so.3gf")]
    (else '("libblas" "libcblas"))))

;; lib-lapack : (listof string)
;;
;; Possible names for the LAPACK library
;;
;; liblapack.so.3gf is apparently required for Debian.  Thanks to
;; Viktor Winschel for pointing this out.
(define lib-lapack
  (case (system-type)
    [(macosx) '("libLAPACK")]
    [(unix)  '("liblapack" "liblapack.so.3gf")]
    (else '("liblapack"))))

(define (string-empty? s)
  (= (string-length s) 0))

;; base-paths : (listof (U path string))
(define base-paths
  (filter (lambda (path)
            (and (not (string-empty? path))
                 (directory-exists? path)))
          (append search-paths (list default-path))))

(define (build-path* . paths-or-empty)
  (apply build-path (filter (lambda (p-or-e) (not (string-empty?
                                                   p-or-e))) paths-or-empty)))

(define (find-libraries-that-exist names search-paths)
  (let ([found
         (remove not
                 (append-map
                  (lambda (name)
                    (map (lambda (search-path)
                           (with-handlers
                               ([exn? (lambda (exn) #f)])
                             (ffi-lib (build-path* search-path name))))
                         search-paths))
                  names))])
    (if (null? found)
        (error
         (format "Could not find any of ~a under paths ~a~n" names
                 search-paths))
        found)))

(define *blas* (car (find-libraries-that-exist lib-blas base-paths)))
(define *lapack* (car (find-libraries-that-exist lib-lapack base-paths)))

(define _cblas-order (_enum '(row-major = 101 col-major = 102)))
(define _cblas-transpose (_enum '(no-trans = 111 trans = 112 conj-trans = 113 atlas-conj = 114)))
(define _cblas-uplo (_enum '(upper = 121 lower = 122)))
(define _cblas-diag (_enum '(non-unit = 131 unit = 132)))
(define _cblas-side (_enum '(left = 141 right = 142)))

(define-for-syntax (append-to-syntax stx . objs)
  (let ((strings (map (lambda (obj)
                        (cond
                          ((symbol? obj) (symbol->string obj))
                          ((syntax? obj) (symbol->string (syntax->datum obj)))
                          (else obj)))
                      objs)))
    (datum->syntax stx
                          (string->symbol
                           (apply string-append strings)))))

(define-syntax (define-blas stx)
  (syntax-case stx ()
    ((define-blas name args ...)
     (with-syntax ((_TAGvector (datum->syntax stx '_TAGvector))
                   (TAGvector-length (datum->syntax stx 'TAGvector-length))
                   (_type (datum->syntax stx '_type))
                   (sname (append-to-syntax stx 's (syntax name)))
                   (cblas_sname (append-to-syntax stx 'cblas_s (syntax name)))
                   (dname (append-to-syntax stx 'd (syntax name)))
                   (cblas_dname (append-to-syntax stx 'cblas_d (syntax name))))
       (syntax/loc stx
         (begin
           (provide* (unsafe sname))
           (define sname
             (let ((_TAGvector _f32vector)
                   (TAGvector-length f32vector-length)
                   (_type _float))
               (get-ffi-obj 'cblas_sname *blas*
                            (_fun args ...)
                            (lambda ()
                              (lambda x
                                (error 'blas (string-append "function "
                                                            (symbol->string 'cblas_sname)
                                                            " not found in blas library.")))))))
           (provide* (unsafe dname))
           (define dname
             (let ((_TAGvector _f64vector)
                   (TAGvector-length f64vector-length)
                   (_type _double*))
               (get-ffi-obj 'cblas_dname *blas*
                            (_fun args ...)
                            (lambda ()
                              (lambda x
                                (error 'blas (string-append "function "
                                                            (symbol->string 'cblas_dname)
                                                            " not found in blas library.")))))))))))))

(define-blas dot _int _TAGvector _int _TAGvector _int -> _type)
(define-blas nrm2 _int _TAGvector _int -> _type)
(define-blas asum _int _TAGvector _int -> _type)
(define-blas swap _int _TAGvector _int _TAGvector _int -> _void)
(define-blas copy _int _TAGvector _int _TAGvector _int -> _void)
(define-blas axpy _int _type _TAGvector _int _TAGvector _int -> _void)
(define-blas scal _int _type _TAGvector _int -> _void)
(define-blas gemv _cblas-order _cblas-transpose _int _int _type _TAGvector _int _TAGvector _int _type _TAGvector _int -> _void)
(define-blas gbmv _cblas-order _cblas-transpose _int _int _int _int _type _TAGvector _int _TAGvector _int _type _TAGvector _int -> _void)
(define-blas trmv _cblas-order _cblas-uplo _cblas-transpose _cblas-diag _int _TAGvector _int _TAGvector _int -> _void)
(define-blas tbmv _cblas-order _cblas-uplo _cblas-transpose _cblas-diag _int _int _TAGvector _int _TAGvector _int -> _void)
(define-blas tpmv _cblas-order _cblas-uplo _cblas-transpose _cblas-diag _int _TAGvector _TAGvector _int -> _void)
(define-blas trsv _cblas-order _cblas-uplo _cblas-transpose _cblas-diag _int _TAGvector _int _TAGvector _int -> _void)
(define-blas tbsv _cblas-order _cblas-uplo _cblas-transpose _cblas-diag _int _int _TAGvector _int _TAGvector _int -> _void)
(define-blas tpsv _cblas-order _cblas-uplo _cblas-transpose _cblas-diag _int _TAGvector _TAGvector _int -> _void)
(define-blas symv _cblas-order _cblas-uplo _int _type _TAGvector _int _TAGvector _int _type _TAGvector _int -> _void)
(define-blas sbmv _cblas-order _cblas-uplo _int _int _type _TAGvector _int _TAGvector _int _type _TAGvector _int -> _void)
(define-blas spmv _cblas-order _cblas-uplo _int _type _TAGvector _TAGvector _int _type _TAGvector _int -> _void)
(define-blas ger _cblas-order _int _int _type _TAGvector _int _TAGvector _int _TAGvector _int -> _void)
(define-blas syr _cblas-order _cblas-uplo _int _type _TAGvector _int _TAGvector _int -> _void)
(define-blas spr _cblas-order _cblas-uplo _int _type _TAGvector _int _TAGvector -> _void)
(define-blas syr2 _cblas-order _cblas-uplo _int _type _TAGvector _int _TAGvector _int _TAGvector _int -> _void)
(define-blas spr2 _cblas-order _cblas-uplo _int _type _TAGvector _int _TAGvector _int _TAGvector -> _void)
(define-blas gemm _cblas-order _cblas-transpose _cblas-transpose _int _int _type _TAGvector _int
  _TAGvector _int _type _TAGvector _int -> _void)
(define-blas symm _cblas-order _cblas-side _cblas-uplo _int _int _type _TAGvector _int
  _TAGvector _int _type _TAGvector _int -> _void)
(define-blas syrk _cblas-order _cblas-uplo _cblas-transpose _int _int _type
  _TAGvector _int _type _TAGvector _int -> _void)
(define-blas syr2k _cblas-order _cblas-uplo _cblas-transpose _int _int _type _TAGvector
  _int _TAGvector _int _type _TAGvector _int -> _void)
(define-blas trmm _cblas-order _cblas-side _cblas-uplo _cblas-transpose _cblas-diag _int _int
  _type _TAGvector _int _TAGvector _int -> _void)
(define-blas trsm _cblas-order _cblas-side _cblas-uplo _cblas-transpose _cblas-diag _int _int
  _type _TAGvector _int _TAGvector _int -> _void)