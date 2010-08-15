#lang scribble/doc

@(require scribble/manual
          (for-label scheme)
          (for-label "matrix.ss")
          (for-label srfi/4)
          (for-label (only-in scheme/foreign _fun)))

@title{A Linear Algebra Library for PLT Scheme}

This document describes the @filepath{plt-linalg.plt} @PLaneT package.
This package provides linear-algebra procedures based on underlying
native BLAS and LAPACK libraries.  

The library makes a reasonable effort to find the BLAS and LAPACK
libraries that are installed on your system, but this can fail if they
are in a non-standard place.  If so, you should modify the search
paths in @filepath{blas-lapack.ss}; if you do, please contact
@link["mailto:w-farr@northwestern.edu"]{Will M. Farr} with the updated
search paths so they can be added to the library and benefit other
users.

This library uses the older, reference implementation of
@link["http://srfi.schemers.org/srfi-42/"]{SRFI-42}, which you can
access using @scheme[(require srfi/42ref)] instead of @scheme[(require
srfi/42)].  (This is required because of some syntax certificate
errors that occur in the new library.)

The @filepath{plt-linalg.plt} package is released under the
@link["http://www.gnu.org/licenses/old-licenses/gpl-2.0.html"]{GPL
Version 2} (see the file @filepath{COPYING} in the collection
directory for more information).  If you find any bugs, or have
feature requests, please contact @link["w-farr@northwestern.edu"]{Will
M. Farr}.

@section{Vectors}

@defmodule[(planet wmfarr/plt-linalg:1:13/vector)]

The @filepath{plt-linalg.plt} package represents vectors as
f64vectors.

@subsection{Contracts for Vectors}

@defproc[(f64vector-same-length/c (v f64vector?)) flat-contract?]{

Produces a contract which matches f64vectors of the same length as
@scheme[v].}

@subsection{Procedures for Vectors}

@defproc[(f64vector-norm (v f64vector?)) (>=/c 0)]{

Computes the two-norm of @scheme[v].}

@defproc[(f64vector-copy (v f64vector?)) f64vector?]{

Produces a copy of @scheme[v].}

@defproc[(f64vector-scale (v f64vector?) (x real?)) f64vector?]{

Produces a new f64vector whose elements are those of @scheme[v] scaled
(multiplied) by @scheme[x].}

@defproc[(f64vector-add (v1 f64vector?) 
                        (v2 (f64vector-same-length/c v1))) 
         f64vector?]{

Produces a new f64vector whose elements are the sum of corresponding
elements of @scheme[v1] and @scheme[v2].}

@defproc[(f64vector-sub (v1 f64vector?)
                        (v2 (f64vector-same-length/c v1)))
          f64vector?]{

Produces a new f64vector whose elements are the difference of
corresponding elements of @scheme[v1] and @scheme[v2].}

@defproc[(f64vector-dot (v1 f64vector?)
                        (v2 (f64vector-same-length/c v1)))
          real?]{

Computes the dot product of @scheme[v1] and @scheme[v2].}

@section{Matrices}

@defmodule[(planet wmfarr/plt-linalg:1:13/matrix)]

A matrix is represented in the @filepath{plt-linalg.plt} package by a
special datastructure: @scheme[(make-matrix rows cols init)].  It
contains a contiguous block of memory to hold @scheme[(* rows cols)]
double-precision floating-point numbers.  Operations on matrices (and
vectors) are implemented by interfacing to native BLAS and LAPACK
libraries.

@subsection{The Matrix Structrue}

@defproc[(make-matrix (rows natural-number/c)
                      (cols natural-number/c)
                      (init real?))
         matrix?]{

Constructs a matrix with @scheme[rows] rows and @scheme[cols] columns.
The elements of this matrix are all @scheme[init].}

@defproc[(matrix (rows natural-number/c)
                 (cols natural-number/c)
                 (elt real?)
                 ...)
          matrix?]{

Constructs a matrix with the given @scheme[elt ...] as elements.
There must be exactly @scheme[(* rows cols)] elements given; the
matrix is constructed in column-major order (note @emph{this is not
left-to-right order}).  Column-major order is sometimes also called
"FORTRAN" order, as opposed to "C" order.}

@defproc[(matrix-identity (n natural-number/c)) matrix?]{Constructs an
@scheme[n] by @scheme[n] identity matrix.}

@defproc[(matrix? (obj any/c)) boolean?]{

Type predicate for matrices.}

@defproc*[(((matrix-rows (m matrix?)) natural-number/c)
           ((matrix-cols (m matrix?)) natural-number/c))]{

Selectors.}

@defproc*[(((matrix-ref (m matrix?) 
                        (i (matrix-valid-row-index/c m))
                        (j (matrix-valid-col-index/c m)))
            real?)
           ((matrix-set! (m matrix?)
                         (i (matrix-valid-row-index/c m))
                         (j (matrix-valid-col-index/c m))
                         (x real?))
            any))]{

Returns or sets the @scheme[i]-@scheme[j]-th element of @scheme[m]. }

@defthing[struct:matrix any/c]{Struture type descriptor for matrix
structs.}
        
@defthing[s:matrix any/c]{Transformer binding for matrix struct.}

@defthing[_matrix any/c]{Foreign type for matrices.  Translates to a
pointer to a column-major (i.e. "FORTRAN", not "C" order) block of
memory containing the elements of the matrix.  Can be used inside
@scheme[_fun] syntax in input, output, or input-output forms.  In
output form, the syntax is @scheme[(_matrix o row-expr col-expr)].}

@subsection{Matrix Contracts}

Some contracts.

@defproc[(matrix-multiplication-compatible/c (m matrix?))
         flat-contract?]{Contract matching matrices which can be
         left-multiplied by @scheme[m].}

@defproc*[(((matrix-valid-row-index/c (m matrix?)) flat-contract?)
           ((matrix-valid-col-index/c (m matrix?)) flat-contract?))]{
           Contracts which matches natural numbers which are valid
           row/col indices for @scheme[m].  (That is, @scheme[x]s such
           that @scheme[(<= 0 x (sub1 (matrix-rows m)))] or
           @scheme[(<= 0 x (sub1 (matrix-cols m)))].)}

@defthing[matrix-square/c flat-contract?]{Matches square matrices.}

@defproc[(matrix-same-dimensions/c (m matrix?)) flat-contract?]{
Matches matrices with the same dimensions as @scheme[m].}

@defproc[(matrix-col-vector-compatible/c (m matrix?)) flat-contract?]{
Matches f64vectors which are compatible for multiplication on the left
by @scheme[m].}

@defproc[(matrix-row-vector-compatible/c (m matrix?)) flat-contract?]{
Matches f64vectors which are compatible for multiplication on the
right by @scheme[m].}

@defstruct[(exn:singular-matrix exn) ((elt natural-number/c))
#:inspector #f]{Thrown by LU-decomposition routines.  @scheme[elt] is
the diagonal element of U which is zero.}

@subsection{Single-Matrix Operations}

@defproc[(matrix-norm (m matrix?)) (>=/c 0)]{Returns the two-norm of a
given matrix.}

@defproc[(matrix-transpose (m matrix?)) matrix?]{Returns a new matrix
which is the transpose of @scheme[m].}

@defproc[(matrix-inverse (m square-matrix/c)) matrix?]{Returns the
matrix inverse of @scheme[m].  @emph{If you are trying to solve linear
equations, it is much more stable (and efficient) to use the
@scheme[matrix-solve] or @scheme[matrix-solve-many] procedures.}}

@subsection{Matrix-Matrix and Matrix-Vector Operations}

@defproc*[(((matrix-add (m1 matrix?) (m2 (matrix-same-dimensions/c
m1)))
            matrix?)
           ((matrix-sub (m1 matrix?) (m2 (matrix-same-dimensions/c
           m1)))
            matrix?))]{ Matrix addition and subtraction.  Does not
            modify arguments.}

@defproc[(matrix-scale (m matrix?) (x real?)) matrix?]{Returns a
matrix whose elements are those of @scheme[m] scaled by @scheme[x].}

@defproc[(matrix-mul (m1 matrix?) (m2
(matrix-multiplication-compatible/c m1)))
         matrix?]{Matrix multiplication; does not modify its arguments.}

@subsection{Solving Linear Systems}

@defproc[(matrix-solve (m matrix-square/c)
                       (b (matrix-row-vector-compatible/c m)))
         (matrix-col-vector-compatible/c m)]{Solves the system
                       @scheme[m]*x = @scheme[b] for x.}

@defproc[(matrix-solve-many (m matrix-square/c) (b
                            (matrix-multiplication-compatible/c m)))
         (matrix-multiplication-compatible/c m)]{ 

Solves simultaneously many linear systems of the form @scheme[m]*x =
@scheme[b] for x.}