#| info.ss: Information file for plt-linalg package.
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

(module info (lib "infotab.ss" "setup")
  
  (define name "PLT-Linalg")
  
  (define blurb
    (list "Simple linear algebra operations in double-precision for PLT Scheme.  Uses LAPACK and BLAS.  Contains bindings for all single- and double-precision BLAS operations."))
  
  (define primary-file "all.ss")
  
  (define categories '(scientific datastructures))
  
  (define can-be-loaded-with 'all)
  
  (define release-notes
    '("Updated linux paths to add /usr/libsse2 per bug report #196."))
  
  (define repositories '("4.x"))
  
  (define scribblings '(("plt-linalg.scrbl" (multi-page))))
  
  (define compile-omit-paths
    '("_darcs")))