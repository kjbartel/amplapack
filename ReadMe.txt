AMP LAPACK: C++ AMP LAPACK Library

This library contains a subset of the LAPACK library for C++ AMP.  This library
depends on the C++ AMP BLAS library (http://ampblas.codeplex.com), which must be
present on your computer in order to build this LAPACK library.

The library is in the form of C++ header files that are available in "inc" directory.

To use the library:
1) add "inc" to INCLUDE, or use compiler switch "/I inc", 
  or add "inc" to "Include Directories" in visual studio project properties.
2) Add "include <amp_lapack.h>" in cpp source file. 
  (There may be new headers that will be introduced in future)

 