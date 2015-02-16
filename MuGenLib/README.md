Introduction       {#mainpage}
=============

MuGen is a C++ library that implements a comprehensive approach to Bayesian inference of multi-trait models for quantitative genetics.  It enables genome-wide association 
studies and genome-enabled prediction, allows for complicated and unbalanced experimental designs, outlier observations, and missing data.  However, while motivation for the construction of this library was to model genetic data, other types of data that could benefit from hierarchical treatment (see \cite gelman04 \cite gelman07 ) can also be analyzed.

The basic idea in hierarchical modeling is that so-called location parameters (i.e., group means and regression coefficients) can be arranged so that levels are nested within each other.  This arrangement is particularly suited for Bayesian inference because model hierarchy fits a hierarchy of prior distributions \cite gelman04 \cite gelman07 \cite greenberg10 .
In addition to the location parameters, scale parameters (in our case, covariance matrices) are also modeled.  Furthermore, some location parameters may not be hierarchically nested within other location parameters.  Most parameters of interest can be fit into one of these categories, and given an implementation that can form a part of a model. These individual parameters are represented by objects of various classes in MuGen.
Having provided implementations for Gibbs sampling of a variety of types of parameters, we enable the user to quickly build a flexible model that accomodates most applications without having to worry about computational complexity.

The structure of the hierarchy is represented by index classes (RanIndex and its derivatives) that store information about which location parameter is related to which upper or lower member of the model hierarchy.  These index parameters can themselves stochastically change to accomodate mixture models where membership of a group is uncertain.  However, the latter feature is still experimental.

Describing statistical hierarchical models unfortunately creates a collision of terms, since our classes are also arranged in a hierarchy.  The latter is the C++ class derivation hierarchy and is not related to the statistical term.  In the documentation we use the term "hierarchy" we to describe the statistical model relationships among variables unless otherwise specified.

Internal implementation of MuGen is multithreaded (using OpenMP, <http://openmp.org/wp/>) and uses GNU Scientific Library and BLAS for computation.  These computational details are hidden behind the interface which is designed for users familiar 
with basic C++ programming.  A publication describing and validating the library is in preparation.

File formats
-------------

Unless explicitly specified otherwise, all files are saved and read as binary files created by GSL functions of the class `gsl_matrix_fwrite`.  To read and write these files from R, we created a set of functions found in the file GSLRIO.cpp.  To compile for R, use

	R CMD SHLIB GSLRIO.cpp -lgsl -O3 -DHAVE_INLINE -DGSL_RANGE_CHECK_OFF

and put the resulting GSLRIO.o where R can find it.  It can then be loaded using `dyn.load()`.

Dependencies
-------------

Aside from the C++ standard library, the implementation depends only on the GNU Scienetific Library (GSL) (<http://www.gnu.org/software/gsl/>) and Basic Linear Algebra Subprograms (BLAS, <http://www.netlib.org/blas/>).
We are currently working with GSL version 1.16, but any version that supports linear algebra routines like Cholesky decomposition will work.  Choosing the right BLAS implementation is crucial for performance.
Generic BLAS is not recommended.  We have successfully compiled MuGen with Apple's Accelerate framework, Intel's MKL (<https://software.intel.com/en-us/intel-mkl>) and AMD's ACML (<http://developer.amd.com/tools-and-sdks/cpu-development/amd-core-math-library-acml/>), all with good results.
MKL seems to perform similarly to ACML on AMD Opteron processors we use to perform most computationally intensive tasks.

Compilation instructions
-------------------------

We use the GNU compiler collection (both version 4.7 and 4.9) for compilation.  We did not use any C++11 features, so older compilers should work, too.  We also successfully compiled the library with Intel icc and Apple's LLVM version 6.0.  However, the latter does not support OpenMP, so the resulting code is not multithreaded.

To compile, make sure GSL is where the compiler can see it and is pointing to your preferred implemetation of BLAS.  Then, simply run

	g++ -c MuGen.cpp -lz -lgslcblas -lgsl -lpthread -lm -fopenmp \\
		-DHAVE_INLINE -DGSL_RANGE_CHECK_OFF -O3 -march=native
    
	ar crv libMuGen.a MuGen.o

and move libMuGen.a to, say, /usr/local/lib/ and MuGen.h to, say, /usr/local/include/.

MuGen has been successfully compiled under RedHat, Amazon Linux and Mac OS X Mavericks and Yosemite.  On Mac OS X, at least when using gcc-4.7 and above, using -march=native does not always work.  In that case, try including -m64 -msse4.2 -mfancy-math-387.  
Windows compilation has not been attempted.
