/*
*  libMuGen.h
*  libMuGen
*
*  Created by ajg67 on 10/30/12.
*  Copyright (c) 2012 SEELE. All rights reserved.
*
*  Classes for Hierarchical Bayesian quantitative-genetic models.  Using gsl_vector and gsl_matrix as internal storage types.
*
*
*/

/// C++ classes for Hierarchical Bayesian Multi-trait quantitative-genetic models.
/** \file
 * \author Anthony J. Greenberg
 * \copyright GNU public license
 * \version 0.9.0
 * 
 * MuGen is a library that implements a comprehensive approach to Bayesian inference of multi-trait models for quantitative genetics.  It enables genome-wide association 
 * studies and genome-enabled prediction, allows for complicated and unbalanced experimental designs, outlier observations, and missing data.  Internal implementation is
 * multithreaded and uses GNU Scientific Library and BLAS for computation.  These computational details are hidden behind the interface which is designed for users familiar 
 * with basic C++ programming.
 *
 */

#ifndef LIBMUGEN_H
#define LIBMUGEN_H
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_permutation.h>
#include <vector>
#include <string>
#include <list>

using std::vector;
using std::list;
using std::string;

class MVnorm;
class MVnormMu;
class MVnormMuPEX;
class MVnormBetaPEX;
class MVnormMuMiss;
class MVnormBeta;
class MVnormBetaMiss;
class MVnormBetaFt;
class MVnormBetaSc;
class MVnormMuBlk;
class MVnormBetaBlk;
class MVnormBetaFtBlk;
class RanIndex;
class RanIndexVS;
class Apex;
class Grp;
class MuGrp;
class MuGrpPEX;
class BetaGrpPEX;
class MuGrpMiss;
class MuGrpEE;
class MuGrpEEmiss;
class BetaGrpSnp;
class BetaGrpSnpCV;
class BetaGrpPSR;
class BetaGrpSnpMiss;
class BetaGrpSnpMissCV;
class BetaGrpPSRmiss;
class BetaGrpFt;
class BetaGrpSc;
class BetaGrpPC;
class BetaGrpPCpex;
class BetaGrpBVSR;
class MuBlk;
class BetaBlk;
class SigmaI;
class SigmaIblk;
class SigmaIpex;
class StTq;
class Qgrp;
class QgrpPEX;
class MixP;

/** \defgroup distributions Various distribution functions
 *
 * Sampling functions for distributions not already in GSL.
 *
 * \warning These are internal functions that are invoked from within classes. Because they need to be as efficient as possible, no range-checking is performed (especially if GSL range checking is turned off at compilation).  They rely on the user to keep track of dimensions.
 *
 *
 * @{
 */

/** \brief Multivariate Gaussian sampling function
 * 
 * This function gives a sample from the MV Gaussian distribution with a given mean and covariance matrix.
 *
 * \param[in] gsl_vector* pointer to a vector of means
 * \param[in] gsl_matrix* pointer to a matrix with the Cholesky-decomposed covariance matrix
 * \param[in] gsl_rng* pointer to a pseudo-random number generator (PNG)
 * \param[out] gsl_vector* pointer to a destination vector for the sample
 *
 */
void MVgauss(const gsl_vector*, const gsl_matrix*, const gsl_rng*, gsl_vector*);

/** \defgroup Overloaded Wishart sampling functions
 *
 * These functions give samples from the Wishart distribution, the only difference between them is the type of the degrees of freedom paramter.
 *
 * \param[in] gsl_matrix* pointer to a matrix with the Cholesky-transformed covariance matrix parameter
 * \param[in] int/size_t& degrees of freedom parameter
 * \param[in] gsl_rng* pointer to a PNG
 * \param[out] gsl_matrix* pointer to a matrix with the sample
 * 
 * @{
 */
void Wishart(const gsl_matrix*, const int&, const gsl_rng*, gsl_matrix*);
void Wishart(const gsl_matrix*, const size_t&, const gsl_rng*, gsl_matrix*);
/** @} */

/** \brief Truncated geometric distribution
 *
 * Sampling from the truncated geometric distribution, returning a size_t type index value.
 *
 * \param[in] double& probability of success
 * \param[in] size_t& maximum value
 * \param[in] gsl_rng* pointer to a PNG
 * \return a value of the type size_t that is a sample from the distribution.
 *
 */
size_t rtgeom(const double&, const size_t&, const gsl_rng*);
/** @} */

/** \defgroup auxFun Auxiliary functions
 *	 
 * Various functions that are needed internally.
 *
 * @{
 */
/** \brief Mahalanobis distance
 *
 * Calculates the Mahalanobis distance of a vector from zero.
 *
 * \param[in] gsl_vector* the vector
 * \param[in] gsl_matrix* the corresponding inverse-covariance matrix
 * \return a scalar value of type double.
 *
 */
double mhl(const gsl_vector *beta, const gsl_matrix *SigI);
/** \defgroup centerFun Centering functions
 *
 * Functions that center matrix columns and vectors.
 *
 * @{
 */
/** \brief Matrix centering in-place
 *
 * \param[in,out] gsl_matrix* the matrix to be modified in-place
 *
 */
void colCenter(gsl_matrix *inplace);
/** \brief Matrix centering with copy
 *
 * \param[in] gsl_matrix* the matrix to be centered, left unmodified
 * \param[out] gsl_matrix* the modified matrix
 *
 */
void colCenter(const gsl_matrix *source, gsl_matrix *res);
/** \brief Matrix centering in-place with missing values
 *
 * \param[in,out] gsl_matrix* matrix to be modified in-place
 * \param[in] double& label for missing values
 *
 */
void colCenter(gsl_matrix *inplace, const double &absLab);
/** \brief Matrix centering with copy and missing values
 *
 *	In addition to centering does mean-imputation of missing values.
 *
 * \param[in] gsl_matrix* matrix to be centered, but not changed
 * \param[out] gsl_matrix* centered matrix
 * \param[in] double& label for missing values
 *
 */
void colCenter(const gsl_matrix *source, gsl_matrix *res, const double &absLab);
/** \brief Vector centering in-place
 *
 * \param[in,out] gsl_vector* vector to be modified
 * \param[in] double& label for missing values
 *
 */
void vecCenter(gsl_vector *inplace);
/** \brief Vector centering with copy
 *
 * \param[in] gsl_vector* vector to be centered, but not changed
 * \param[out] gsl_vector* modified vector
 *
 */
void vecCenter(const gsl_vector *source, gsl_vector *res);
/** \brief Mean imputation without centering
 *
 * \param[in,out] gsl_matrix* matrix to be modified in-place
 * \param[in] double& label for missing values
 *
 */
void meanImpute(gsl_matrix *inplace, const double &absLab);
/** @} */
/** \brief Print matrix to screen
 *
 * Printing a GSL matrix to screen, with each row on a line, values seprated by spaces.
 *
 * \param[in] gsl_matrix* matrix to be displayed
 * \param[in] double& label for missing data
 *
 */
void printMat(const gsl_matrix *m);
/** \brief Accessing the processor RTDSC instruction
 *
 * \return a value of type unsigned long long.
 *
 * This function outputs the processor RTDSC instruction for use in seeding random number generators using assembly code.
 * \warning this function is unlikely to work under Windows.
 *
 */
unsigned long long rdtsc();
/** @} */

/** \defgroup rowLocParam Individual location parameters
 *
 * A hierarchy of classes that refer to rows of location parameter matrices. These classes are internal to Grp classes and are not directly declared by the user of the library.  
 * They are typically updated with samples from the multivariate normal distribution.  Because the methods are used for much of the computation, they are implemented for speed and almost no chacking for errors is done.
 * The assumption is that the encapsualting classes will do all of that correctly. Constructors are used to set initial values.
 *
 * @{
 */
/** \brief The abstract base class for location parameter rows
 *
 * This is the generic location parameter class and cannot be envoked directly.  All the constructors are protected.
 */
class MVnorm {
protected:
	/** \brief Data vector
	 * 
	 * This is actually a pointer (implemented as a GSL *vector_view*) to a row in the corresponding matrix of location parameters.
	 *
	 */
	gsl_vector_view _vec;
	
	/** \brief Length of the data vector */
	size_t _d;
	
	/** \brief Default constructor
	 *
	 * This is a deterministic constructor that results in a pointer to nowhere and the dimension member set to zero.
	 */
	MVnorm() : _d(0) {};
	/** \brief Dimension-only constructor
	 *
	 * Sets the dimension to a value, but the *vector_view* is not assigned a target.
	 * 
	 * \param[in] size_t& dimension value
	 */
	MVnorm(const size_t &d) : _d(d) {};
	/** \brief Dimension and vector value constructor
	 *
	 * A deterministic constructor that sets the _vec member to point to the provided *gsl_vector* and the dimension member to the size of the provided *gsl_vector*.
	 *
	 * \param[in] gsl_vector* target vector
	 */
	MVnorm(gsl_vector *mn);
	/** \brief Univariate Gaussian constructor
	 *
	 * Sets the _vec member to point to the provided *gsl_vector* and ads univariate Gaussian noise independtly to each element, modifying the target.
	 *
	 * \param[in,out] gsl_vector* vector of means
	 * \param[in] gsl_vector* vector of standard deviations for the univariate Gaussian, has to be no shorter than the mean vector.  If it is longer, the extra values are ignored.
	 * \param[in] gsl_rng* pointer to a PNG
	 *
	 */
	MVnorm(gsl_vector *mn, const gsl_vector *sd, const gsl_rng *r);
	/** \brief Multivariate Gaussian constructor
	 *
	 * Sets the _vec member to point to the provided *gsl_vector* and ads multivariate Gaussian noise, modifying the target.
	 *
	 * \param[in,out] gsl_vector* vector of means
	 * \param[in] gsl_matrix* covariance matrix for the multivariate Gaussian, has to match the mean vector in dimensions and has to be symmetric positive-definite
	 * \param[in] gsl_rng* pointer to a PNG
	 *
	 */
	MVnorm(gsl_vector *mn, const gsl_matrix *Sig, const gsl_rng *r);
	/** \brief Dimension and vector value constructor
	 *
	 * A deterministic constructor that sets the _vec member to point to a row (indicated by the row index) of the provided *gsl_matrix* and the dimension member to the number of columns in the provided *gsl_matrix*.
	 *
	 * \param[in] gsl_matrix* target matrix
	 * \param[in] size_t& row index of the matrix
	 */
	MVnorm(gsl_matrix *mn, const size_t &iRw);
	/** \brief Univariate Gaussian constructor with a matrix
	 *
	 * Sets the _vec member to point to a row (indicated by the row index) of the provided *gsl_matrix* and ads univariate Gaussian noise independtly to each element, modifying the target.
	 *
	 * \param[in,out] gsl_matrix* target matrix
	 * \param[in] size_t& row index of the target matrix
	 * \param[in] gsl_vector* vector of standard deviations for the univariate Gaussian, has to be no shorter than a row of the mean matrix.  If it is longer, the extra values are ignored.
	 * \param[in] gsl_rng* pointer to a PNG
	 *
	 */
	MVnorm(gsl_matrix *mn, const size_t &iRw, const gsl_vector *sd, const gsl_rng *r);
	/** \brief Multivariate Gaussian constructor with a matrix
	 *
	 * Sets the _vec member to point to a row (indicated by the row index) of the provided *gsl_matrix* and ads multivariate Gaussian noise, modifying the target.
	 *
	 * \param[in,out] gsl_matrix* target matrix
	 * \param[in] size_t& row index of the target matrix
	 * \param[in] gsl_matrix* covariance matrix for the multivariate Gaussian, has to be symmetric positive-definite, with diemsions equal to the number of columns in the target matrix
	 * \param[in] gsl_rng* pointer to a PNG
	 *
	 */
	MVnorm(gsl_matrix *mn, const size_t &iRw, const gsl_matrix *Sig, const gsl_rng *r);
	
public:
	/** \brief Copy constructor
	 *
	 * Provided for completeness and has not been extensively tested.  The constructor is deterministic, i.e.\ the values are copied exactly, with no stochastic perturbation.
	 *
	 * \param[in] MVnorm& object to be copied
	 * \return object of class _MVnorm_.
	 */
	MVnorm(const MVnorm &);
	/** \brief Assignement operator
	 *
	 * Provided for completeness and has not been extensively tested.  The operator is deterministic, i.e.\ the values are copied exactly, with no stochastic perturbation.
	 *
	 * \param[in] MVnorm& object to be assigned
	 * \return reference to an object of class _MVnorm_.
	 */
	MVnorm & operator=(const MVnorm &);
	
	/** \brief Virtual destructor
	 *
	 * The destructor has nothing to do since members of this class cannot be de-allocated.
	 */
	virtual ~MVnorm() {};
	
	/** \defgroup updateFun Overloaded update functions
	 *
	 * These functions generate updated values of the corresponding lines in the location matrix via Gibbs sampling, given the current states of the parameters listed as arguments. 
	 * In this class, these are all pure virtual functions.
	 *
	 * @{
	 */
	/** \defgroup locIP Improper prior methods
	 *
	 * @{
	 */
	/** \brief Gaussian likelihood
	 *
	 * \param[in] Grp& data
	 * \param[in] SigmaI& inverse-covariance matrix for the likelihood
	 * \param[in] gsl_rng* pointer to a PNG
	 */
	virtual void update(const Grp &, const SigmaI &, const gsl_rng *) = 0;
	/** \brief Sudent-\f$t\f$ likelihood
	 *
	 * \param[in] Grp& data
	 * \param[in] Qgrp& vector of Student-\f$t\f$ covariance scale parameters for the likelihood covariance
	 * \param[in] SigmaI& inverse-covariance matrix for the likelihood
	 * \param[in] gsl_rng* pointer to a PNG
	 */
	virtual void update(const Grp &, const Qgrp &, const SigmaI &, const gsl_rng *) = 0;
	/** @} */
	/** \addtogroup locZMn 0-mean prior methods
	 *
	 * @{
	 */
	/** \brief Gaussian likelihood, Gaussian prior
	 *
	 * \param[in] Grp& data for the likelihood
	 * \param[in] SigmaI& inverse-covariance matrix for the likelihood
	 * \param[in] SigmaI& prior inverse-covariance matrix
	 * \param[in] gsl_rng* pointer to a PNG
	 */
	virtual void update(const Grp &, const SigmaI &, const SigmaI &, const gsl_rng *) = 0;
	/** \brief Gaussian likelihood, Student-\f$t\f$ prior
	 *
	 * \param[in] Grp& data for the likelihood
	 * \param[in] SigmaI& inverse-covariance matrix for the likelihood
	 * \param[in] double& Student-\f$t\f$ scale parameter for the prior covariance
	 * \param[in] SigmaI& prior inverse-covariance matrix
	 * \param[in] gsl_rng* pointer to a PNG
	 */
	virtual void update(const Grp &, const SigmaI &, const double &, const SigmaI &, const gsl_rng *) = 0;
	/** \brief Student-\f$t\f$ likelihood, Gaussian prior
	 *
	 * \param[in] Grp& data for the likelihood
	 * \param[in] Qgrp& vector of Student-\f$t\f$ covariance scale parameters for the likelihood covariance
	 * \param[in] SigmaI& inverse-covariance matrix for the likelihood
	 * \param[in] SigmaI& prior inverse-covariance matrix
	 * \param[in] gsl_rng* pointer to a PNG
	 */
	virtual void update(const Grp &, const Qgrp &, const SigmaI &, const SigmaI &, const gsl_rng *) = 0;
	/** \brief Student-\f$t\f$ likelihood, Student-\f$t\f$ prior
	 *
	 * \param[in] Grp& data for the likelihood
	 * \param[in] Qgrp& vector Student-\f$t\f$ covariance scale parameter for the likelihood covariance
	 * \param[in] SigmaI& inverse-covariance matrix for the likelihood
	 * \param[in] double& Student-\f$t\f$ scale parameter for the prior covariance
	 * \param[in] SigmaI& prior inverse-covariance matrix
	 * \param[in] gsl_rng* pointer to a PNG
	 */
	virtual void update(const Grp &, const Qgrp &, const SigmaI &, const double &, const SigmaI &, const gsl_rng *) = 0;
	/** @} */
	/** \addtogroup locNZMn non-0-mean prior methods
	 *
	 * @{
	 */
	/** \brief Gaussian likelihood, Gaussian prior
	 *
	 * \param[in] Grp& data for the likelihood
	 * \param[in] SigmaI& inverse-covariance matrix for the likelihood
	 * \param[in] Grp& prior mean
	 * \param[in] SigmaI& prior inverse-covariance matrix
	 * \param[in] gsl_rng* pointer to a PNG
	 */
	virtual void update(const Grp &, const SigmaI &, const Grp &, const SigmaI &, const gsl_rng *) = 0;
	/** \brief Gaussian likelihood, Student-\f$t\f$ prior
	 *
	 * \param[in] Grp& data for the likelihood
	 * \param[in] SigmaI& inverse-covariance matrix for the likelihood
	 * \param[in] Grp& prior mean
	 * \param[in] double& Student-\f$t\f$ scale parameter for the prior covariance
	 * \param[in] SigmaI& prior inverse-covariance matrix
	 * \param[in] gsl_rng* pointer to a PNG
	 */
	virtual void update(const Grp &, const SigmaI &, const Grp &, const double &, const SigmaI &, const gsl_rng *) = 0;
	/** \brief Student-\f$t\f$ likelihood, Gaussian prior
	 *
	 * \param[in] Grp& data for the likelihood
	 * \param[in] Qgrp& vector Student-\f$t\f$ covariance scale parameter for the likelihood covariance
	 * \param[in] SigmaI& inverse-covariance matrix for the likelihood
	 * \param[in] Grp& prior mean
	 * \param[in] SigmaI& prior inverse-covariance matrix
	 * \param[in] gsl_rng* pointer to a PNG
	 */
	virtual void update(const Grp &, const Qgrp &, const SigmaI &, const Grp &, const SigmaI &, const gsl_rng *) = 0;
	/** \brief Student-\f$t\f$ likelihood, Student-\f$t\f$ prior
	 *
	 * \param[in] Grp& data for the likelihood
	 * \param[in] Qgrp& vector Student-\f$t\f$ covariance scale parameter for the likelihood covariance
	 * \param[in] SigmaI& inverse-covariance matrix for the likelihood
	 * \param[in] Grp& prior mean
	 * \param[in] double& Student-\f$t\f$ scale parameter for the prior covariance
	 * \param[in] SigmaI& prior inverse-covariance matrix
	 * \param[in] gsl_rng* pointer to a PNG
	 */
	virtual void update(const Grp &, const Qgrp &, const SigmaI &, const Grp &, const double &, const SigmaI &, const gsl_rng *) = 0;
	/** @} */
	/** @} */
	// end of the update group
	
	/** \defgroup mahal Mahalanobis distance functions
	 *
	 * Overloaded functions to calculate the Mahalanobis distance of the vector of location parameters stored in the class to another vector or to zero.
	 *
	 * @{
	 */
	/** \brief Mahalanobis distance to a vector
	 *
	 * Calculates the Mahalanobis distance between the vector stored in the class and the provided vector of class MVnorm.
	 *
	 * \param[in] MVnorm& vector to calculate the distance from
	 * \param[in] SigmaI& inverse covariance matrix for scaling
	 * \return the distance value of type _double_.
	 *
	 */
	virtual double mhl(const MVnorm *x, const SigmaI &SigI);
	/** \brief Mahalanobis distance to a vector
	 *
	 * Calculates the Mahalanobis distance between the vector stored in the class and the provided vector of class MVnorm.  The const version of the above function.
	 *
	 * \param[in] MVnorm& vector to calculate the distance from
	 * \param[in] SigmaI& inverse covariance matrix for scaling
	 * \return the distance value of type _double_.
	 *
	 */
	virtual double mhl(const MVnorm *x, const SigmaI &SigI) const;
	/** \brief Mahalanobis distance to a vector
	 *
	 * Calculates the Mahalanobis distance between the vector stored in the class and the provided GSL vector.
	 *
	 * \param[in] MVnorm& vector to calculate the distance from
	 * \param[in] SigmaI& inverse covariance matrix for scaling
	 * \return the distance value of type _double_.
	 *
	 */
	virtual double mhl(const gsl_vector *x, const SigmaI &SigI);
	/** \brief Mahalanobis distance to a vector
	 *
	 * Calculates the Mahalanobis distance between the vector stored in the class and the provided GSL vector. The const version of the above function.
	 *
	 * \param[in] MVnorm& vector to calculate the distance from
	 * \param[in] SigmaI& inverse covariance matrix for scaling
	 * \return the distance value of type _double_.
	 *
	 */
	virtual double mhl(const gsl_vector *x, const SigmaI &SigI) const;
	/** \brief Mahalanobis distance to zero
	 *
	 * Calculates the Mahalanobis distance between the vector stored in the class zero.
	 *
	 * \param[in] SigmaI& inverse covariance matrix for scaling
	 * \return the distance value of type _double_.
	 *
	 */
	virtual double mhl(const SigmaI &SigI); // distance from 0
	/** \brief Mahalanobis distance to zero
	 *
	 * Calculates the Mahalanobis distance between the vector stored in the class zero.  The const version of the above function.
	 *
	 * \param[in] SigmaI& inverse covariance matrix for scaling
	 * \return the distance value of type _double_.
	 *
	 */
	virtual double mhl(const SigmaI &SigI) const; // distance from 0
	/** @} */
	/** \defgroup MVnormDens Multivariate Gaussian density functions
	 *
	 * These functions calculate multivariate Gaussian density given the mean and inverse-covariance provided and the vector of values stored in the class as data.
	 *
	 * @{
	 */
	/** \brief Multivariate Gaussian density
	 *
	 * \param[in] gsl_vector* vector of means
	 * \param[in] SigmaI& inverse-covariance matrix
	 * \return density value of type _double_.
	 *
	 */
	double density(const gsl_vector *theta, const SigmaI &SigI);
	/** \brief Multivariate Gaussian density
	 *
	 * \note A const version of the above function.
	 *
	 * \param[in] gsl_vector* vector of means
	 * \param[in] SigmaI& inverse-covariance matrix
	 * \return density value of type _double_.
	 *
	 */
	double density(const gsl_vector *theta, const SigmaI &SigI) const;
	/** \brief Multivariate Gaussian density
	 *
	 * \param[in] MVnorm* vector of means
	 * \param[in] SigmaI& inverse-covariance matrix
	 * \return density value of type _double_
	 *
	 */
	double density(const MVnorm *theta, const SigmaI &SigI);
	/** \brief Multivariate Gaussian density
	 *
	 * \note A const version of the above function.
	 *
	 * \param[in] MVnorm* vector of means
	 * \param[in] SigmaI& inverse-covariance matrix
	 * \return density value of type _double_.
	 *
	 */
	double density(const MVnorm *theta, const SigmaI &SigI) const;
	/** @} */
	
	/** \brief Save function
	 *
	 * Saves the current value of the location parameter to a file, appending by default.
	 *
	 * \param[in] string& file name
	 * \param[in] char* save mode, "a" for append by default. The allowable values are as for a standard C file stream
	 *
	 */
	void save(const string &fileNam, const char *how = "a");
	/** \brief Save function
	 *
	 * Saves the current value of the location parameter to a file.
	 *
	 * \param[in] FILE* C file stream
	 *
	 */
	void save(FILE *fileStr);
	
	/** \brief Subscript operator
	 *
	 * Overloaded subscript operator for the vector of location parameter values.
	 *
	 * \param[in] size_t index value
	 * \return value of the member vector corresponding to the index value, of type _double_.
	 *
	 */
	double operator[](const size_t i) const{return gsl_vector_get(&_vec.vector, i); };
	/** \brief Setting an element to a value
	 *
	 * Sets the _i_-th element of the location vector to a value deterministically
	 *
	 * \param[in] size_t index
	 * \param[in] double new value of the element
	 *
	 */
	void valSet(const size_t i, const double x){gsl_vector_set(&_vec.vector, i, x); };
	/** \brief Access the location vector
	 *
	 * Provides access to the internal parameter vector.
	 *
	 * \return GSL vector pointing to the corresponding location parameter.
	 */
	const gsl_vector *getVec() const{return &_vec.vector;};
	/** \brief Length of the location vector
	 *
	 * \return Length of the parameter vector, of type _size_t_. Corresponds to the number of traits.
	 *
	 */
	size_t len() const{return _d; };
	/** \brief Number of missing values
	 *
	 * Accesses the  number of missing phenotype values.  Non-zero only for classes that implement treatment of missing values.
	 *
	 * \return Number of missing phenotypes, of type _size_t_.
	 */
	virtual size_t nMissP() const{return 0; };
	/** \brief Indexes of missing values
	 *
	 * Accesses a vector with indexes of missing phenotypes.
	 *
	 * \return Vector of _size_t_ that contains indexes that correspond to missing values. For classes that do not allow missing values it is empty.
	 *
	 */
	virtual const vector<size_t> getMisPhen() const{ return vector<size_t>(0); };
	/** \brief Points to the corresponding data
	 *
	 * Access to the pointer to a vector that indexes the data used to update an instance of the class.  Is non-zero only for classes which can be part of a hierarchical model.
	 *
	 * \return Pointer to a vector of _size_t_.
	 *
	 */
	virtual const vector<size_t> *down() const{return 0; };
	/** \brief Points to the prior
	 *
	 * Access to the pointer to the correspoding vector of priors.  Is non-zero only for classes where a prior is implemented.
	 *
	 * \return Pointer to _size_t_.
	 *
	 */
	virtual const size_t *up() const{return 0; };
	/** \brief Scale parameter
	 *
	 * Access to a scale parameter member, defined only for some derived regression-type classes.  If not defined, returns 1.0.
	 *
	 * \return Scale paramter value of type _double_.
	 *
	 */
	virtual double scalePar() const {return 1.0; };
};

/** \brief Individual vector of means
 *
 * Refers to a row of a matrix of location parameters that are updated with a likelihood that is an inverse-covariance weighted mean of the data. Variables of this class are typically imbedded in a model hierarchy, and in that case are updated with non-0-mean priors. 
 * This class allows us to keep track of which rows in data and prior matrices to use through pointers to index vectors (for data, where there may be more than one row) and to a single index (for the prior mean, since there is only one).  
 * We use pointers because they allow us to modify the indexes in other parts of the model, thus allowing for things like mixture models and model selection schemes.
 *
 */
class MVnormMu : public MVnorm {
protected:
	/** \brief Data indexes
	 *
	 * A pointer to a vector of _size_t_. The elements are row indexes of the data matrix.  Means among these rows are used in the likelihood of the _udpate_ functions.
	 * The pointer is _const_ to make sure we are not modifying it from this class.
	 *
	 */
	const vector<size_t> *_lowLevel;
	
	/** \brief Prior index
	 *
	 * Pointer to a value of type _size_t_ that is a row index of the prior location matrix. The pointer is _const_ to make sure we are not modifying it from this class.
	 *
	 */
	const size_t *_upLevel;
	
public:
	/** \brief Default constructor
	 *
	 * The _vec member is unassigned and the low-level and upper-level pointers are set to zero.
	 *
	 */
	MVnormMu();
	/** \brief Zero vector constructor
	 *
	 * Sets _vec to point to a vector of zeros, of length *d*.
	 *
	 * \param[in] size_t& size of the vector
	 */
	MVnormMu(const size_t &d) : MVnorm(d) {};
	/** \brief Zero vector with pointers
	 *
	 * Sets _vec to point to a vector of zeros, of length *d*.  Pointers to the lower and upper level of the model hierarchy are set.
	 *
	 * \param[in] size_t& size of the vector
	 * \param[in] vector<size_t>& vector of indexes of rows in the data matrix
	 * \param[in] size_t& index in the matrix of priors
	 *
	 */
	MVnormMu(const size_t &d, const vector<size_t> &low, const size_t &up);
	/** \brief Deterministic constructor
	 *
	 * Sets _vec to point to a GSL vector of values, without modifying the target during the invocation of the constructor. The pointers to the data and the prior indexes are set.
	 *
	 * \param[in] gsl_vector* vector of values
	 * \param[in] vector<size_t>& vector of indexes of rows in the data matrix
	 * \param[in] size_t& index in the matrix of priors
	 *
	 */
	MVnormMu(gsl_vector *mn, const vector<size_t> &low, const size_t &up);
	/** \brief Deterministic constructor, prior index only
	 *
	 * Sets _vec to point to a GSL vector of values, without modifying the target during the invocation of the constructor. Only the pointer to the prior index is set. This constructor can be used for the lowest hierarchy level.
	 *
	 * \param[in] gsl_vector* vector of values
	 * \param[in] size_t& index in the matrix of priors
	 *
	 */
	MVnormMu(gsl_vector *mn, const size_t &up);
	/** \brief Deterministic constructor with a matrix
	 *
	 * Sets _vec to point to a row of the GSL matrix of values, without modifying the target during the invocation of the constructor.
	 *
	 * \param[in] gsl_matrix* matrix of values
	 * \param[in] size_t& row index of the matrix of values, has to be no smaller than 0 and strictly smaller than the number of rows in _mn_
	 *
	 */
	MVnormMu(gsl_matrix *mn, const size_t &iRw);
	/** \brief Deterministic constructor with a matrix and an index to data
	 *
	 * Sets _vec to point to a row of the GSL matrix of values, without modifying the target during the invocation of the constructor.  None of the index pointers are set.
	 *
	 * \param[in] gsl_matrix* matrix of values
	 * \param[in] size_t& row index of the matrix of values, has to be no smaller than 0 and strictly smaller than the number of rows in _mn_
	 * \param[in] vector<size_t>& vector of indexes of rows in the data matrix
	 *
	 */
	MVnormMu(gsl_matrix *mn, const size_t &iRw, const vector<size_t> &low);
	/** \brief Deterministic constructor with a matrix and indexes to data and a prior
	 *
	 * Sets _vec to point to a row of the GSL matrix of values, without modifying the target during the invocation of the constructor.
	 *
	 * \param[in] gsl_matrix* vector of values
	 * \param[in] size_t& row index of the matrix of values, has to be no smaller than 0 and strictly smaller than the number of rows in _mn_
	 * \param[in] vector<size_t>& vector of indexes of rows in the data matrix
	 * \param[in] size_t& index in the matrix of priors
	 *
	 */
	MVnormMu(gsl_matrix *mn, const size_t &iRw, const vector<size_t> &low, const size_t &up);
	/** \brief Deterministic constructor with a matrix and an index to a prior
	 *
	 * Sets _vec to point to a row of the GSL matrix of values, without modifying the target during the invocation of the constructor.
	 *
	 * \param[in] gsl_matrix* vector of values
	 * \param[in] size_t& row index of the matrix of values, has to be no smaller than 0 and strictly smaller than the number of rows in _mn_
	 * \param[in] size_t& index in the matrix of priors
	 *
	 */
	MVnormMu(gsl_matrix *mn, const size_t &iRw, const size_t &up);
	/** \brief Univariate random constructor with a vector and indexes to data and a prior
	 *
	 * Sets _vec to point to a GSL vector of values, modifying the target during the invocation of the constructor. Modification is by adding independent samples from a univariate normal distribution with the provided standard deviations.
	 *
	 * \param[in] gsl_vector* vector of values
	 * \param[in] gsl_vector* vector of standard deviations
	 * \param[in] gsl_rng* pointer to a PNG
	 * \param[in] vector<size_t>& vector of indexes of rows in the data matrix
	 * \param[in] size_t& row index of the matrix of values, has to be no smaller than 0 and strictly smaller than the number of rows in _mn_
	 *
	 */
	MVnormMu(gsl_vector *mn, const gsl_vector *sd, const gsl_rng *r, const vector<size_t> &low, const size_t &up);
	/** \brief Multivariate random constructor with a vector and indexes to data and a prior
	 *
	 * Sets _vec to point to a GSL vector of values, modifying the target during the invocation of the constructor. Modification is by adding samples from a multivariate normal distribution with the provided covariance matrix.
	 *
	 * \param[in] gsl_vector* vector of values
	 * \param[in] gsl_matrix* covariance matrix
	 * \param[in] gsl_rng* pointer to a PNG
	 * \param[in] vector<size_t>& vector of indexes of rows in the data matrix
	 * \param[in] size_t& index in the matrix of priors
	 *
	 */
	MVnormMu(gsl_vector *mn, const gsl_matrix *Sig, const gsl_rng *r, const vector<size_t> &low, const size_t &up);
	/** \brief Univariate random constructor with a matrix and indexes to data and a prior
	 *
	 * Sets _vec to point to a row of the GSL matrix of values, modifying the target during the invocation of the constructor. Modification is by adding independent samples from a univariate normal distribution with the provided standard deviations.
	 *
	 * \param[in] gsl_matrix* matrix of values
	 * \param[in] size_t& row index of the matrix of values, has to be no smaller than 0 and strictly smaller than the number of rows in _mn_
	 * \param[in] gsl_vector* vector of standard deviations
	 * \param[in] gsl_rng* pointer to a PNG
	 * \param[in] vector<size_t>& vector of indexes of rows in the data matrix
	 * \param[in] size_t& index in the matrix of priors
	 *
	 */
	MVnormMu(gsl_matrix *mn, const size_t &iRw, const gsl_vector *sd, const gsl_rng *r, const vector<size_t> &low, const size_t &up);
	/** \brief Multivariate random constructor with a matrix and indexes to data and a prior
	 *
	 * Sets _vec to point to a row of the GSL matrix of values, modifying the target during the invocation of the constructor. Modification is by adding samples from a multivariate normal distribution with the provided covariance matrix.
	 *
	 * \param[in] gsl_matrix* matrix of values
	 * \param[in] size_t& row index of the matrix of values, has to be no smaller than 0 and strictly smaller than the number of rows in _mn_
	 * \param[in] gsl_matrix* covariance matrix
	 * \param[in] gsl_rng* pointer to a PNG
	 * \param[in] vector<size_t>& vector of indexes of rows in the data matrix
	 * \param[in] size_t& index in the matrix of priors
	 *
	 */
	MVnormMu(gsl_matrix *mn, const size_t &iRw, const gsl_matrix *Sig, const gsl_rng *r, const vector<size_t> &low, const size_t &up);
	
	/** \brief Copy constructor
	 *
	 * Deterministic copy constructor
	 *
	 * \param[in] MVnormMu& object of type MVnormMu
	 *
	 */
	MVnormMu(const MVnormMu &);
	/** \brief Assignment operator
	 *
	 * Deterministic assignement operator
	 *
	 * \param[in] MVnormMu& object of type MVnormMu
	 *
	 * \return A reference to an object of type MVnormMu
	 *
	 */
	MVnormMu & operator=(const MVnormMu &);
	
	/** \brief Destructor
	 *
	 */
	virtual ~MVnormMu();
	
	// improper prior methods
	virtual void update(const Grp &dat, const SigmaI &SigIm, const gsl_rng *r);
	virtual void update(const Grp &dat, const Qgrp &q, const SigmaI &SigIm, const gsl_rng *r);
	// 0-mean prior methods
	virtual void update(const Grp &dat, const SigmaI &SigIm, const SigmaI &SigIp, const gsl_rng *r);
	virtual void update(const Grp &dat, const SigmaI &SigIm, const double &qPr, const SigmaI &SigIp, const gsl_rng *r);
	virtual void update(const Grp &dat, const Qgrp &q, const SigmaI &SigIm, const SigmaI &SigIp, const gsl_rng *r);
	void update(const Grp &dat, const Qgrp &q, const SigmaI &SigIm, const double &qPr, const SigmaI &SigIp, const gsl_rng *r);
	// non-zero mean prior methods
	virtual void update(const Grp &dat, const SigmaI &SigIm, const Grp &muPr, const SigmaI &SigIp, const gsl_rng *r);
	virtual void update(const Grp &dat, const SigmaI &SigIm, const Grp &muPr, const double &qPr, const SigmaI &SigIp, const gsl_rng *r);
	void update(const Grp &dat, const Qgrp &q, const SigmaI &SigIm, const Grp &muPr, const SigmaI &SigIp, const gsl_rng *r);
	void update(const Grp &dat, const Qgrp &q, const SigmaI &SigIm, const Grp &muPr, const double &qPr, const SigmaI &SigIp, const gsl_rng *r);
	
	virtual size_t nMissP() const{return 0; };
	virtual const vector<size_t> getMisPhen() const{ return vector<size_t>(0); };
	/** \brief Points to the corresponding data
	 *
	 * Access to the pointer to a vector that indexes the data used to update an instance of the class.  Is implemented for this class.
	 *
	 * \return Pointer to a vector of _size_t_.
	 *
	 */
	const vector<size_t> *down() const{return _lowLevel; };
	/** \brief Points to the prior
	 *
	 * Access to the pointer to the correspoding vector of priors.  Is implemented for this class.
	 *
	 * \return Pointer to _size_t_.
	 *
	 */
	const size_t *up() const{return _upLevel; };
};

/** \brief Individual vector of means with parameter expansion
 *
 * This class is the same as MVnormMu, but implements multivariate parameter expansion for updating (Greenberg, in prep.)
 *
 */
class MVnormMuPEX : public MVnormMu {
protected:
	/** \brief Pointer to the redundant parameter
	 *
	 * A pointer to an object of class Apex, that stores the matrix of redundant parameters.
	 *
	 */
	Apex *_A;
	/** \brief Pre-computed auxiliary matrix
	 *
	 * Vector view of a matrix that stores a pre-computed \f$ (SA)^{\textsc{t}} \f$ matrix that is common for all members of the encapsulating Grp class.
	 *
	 */
	gsl_matrix_view _tSAprod;
	
public:
	/** \brief Default constructor
	 *
	 */
	MVnormMuPEX() : MVnormMu() {_A = 0; };
	/** \brief Deterministic constructor
	 *
	 * \param[in] gsl_matrix* matrix of values
	 * \param[in] size_t& row index of the matrix of values
	 * \param[in] vector<size_t>& vector of indexes of rows in the data matrix
	 * \param[in] size_t& index in the matrix of priors
	 * \param[in] Apex& object storing the redundant parameter matrix
	 * \param[in] gsl_matrix* pre-computed auxiliary matrix
	 *
	 */
	MVnormMuPEX(gsl_matrix *mn, const size_t &iRw, const vector<size_t> &low, const size_t &up, Apex &A, gsl_matrix *tSigIAt);
	/** \brief Univariate random constructor
	 *
	 * \param[in] gsl_matrix* matrix of values
	 * \param[in] size_t& row index of the matrix of values
	 * \param[in] gsl_vector* vector of standard deviations for the univariate Gaussian
	 * \param[in] gsl_rng* pointer to a PNG
	 * \param[in] vector<size_t>& vector of indexes of rows in the data matrix
	 * \param[in] size_t& index in the matrix of priors
	 * \param[in] Apex& object storing the redundant parameter matrix
	 * \param[in] gsl_matrix* pre-computed auxiliary matrix
	 *
	 */
	MVnormMuPEX(gsl_matrix *mn, const size_t &iRw, const gsl_vector *sd, const gsl_rng *r, const vector<size_t> &low, const size_t &up, Apex &A, gsl_matrix *tSigIAt);
	
	/** \brief Deterministic copy constructor
	 *
	 * \param[in] MVnormMuPEX& object of type MVnormMuPEX
	 *
	 */
	MVnormMuPEX(const MVnormMuPEX &mu);
	/** \brief Assignment operator
	 *
	 * Deterministic assignement operator
	 *
	 * \param[in] MVnormMuPEX& object of type MVnormMuPEX
	 *
	 * \return A reference to an object of type MVnormMuPEX
	 *
	 */
	MVnormMuPEX & operator=(const MVnormMuPEX &mu);
	
	/** \brief Destructor */
	virtual ~MVnormMuPEX() {};
	
	virtual void update(const Grp &dat, const SigmaI &SigIm, const SigmaI &SigIp, const gsl_rng *r);
	virtual void update(const Grp &dat, const SigmaI &SigIm, const double &qPr, const SigmaI &SigIp, const gsl_rng *r);
	
	virtual void update(const Grp &dat, const SigmaI &SigIm, const Grp &muPr, const SigmaI &SigIp, const gsl_rng *r);
	virtual void update(const Grp &dat, const SigmaI &SigIm, const Grp &muPr, const double &qPr, const SigmaI &SigIp, const gsl_rng *r);
};

/** \brief Individual regression with parameter expansion
 *
 * Implements handling of one-at-a-time regressions, as in MVnormBetaFt but with parameter expansion.
 *
 */
class MVnormBetaPEX : public MVnormMuPEX {
protected:
	/** \brief Predictor
	 *
	 * Points to a vector of predictor values, typically a column of the predictor matrix that is in the encapsulating Grp class.
	 *
	 */
	gsl_vector_view _X;
	/** \brief Scale parameter
	 *
	 * Typically \f$ x^{\textsc{t}}x \f$, but some exceptions exist in special cases.
	 *
	 */
	double _scale;
	/** \brief Length of the predictor */
	size_t _N;
	/** \brief Fitted values
	 *
	 * Points to a matrix of partial fitted values that exclude the current regression (\f$ \boldsymbol{X}_{\cdot -i}\boldsymbol{B}_{-i \cdot} \f$).
	 *
	 */
	gsl_matrix_view _fitted;
	
public:
	/** \brief Default constructor */
	MVnormBetaPEX() : MVnormMuPEX() {};
	/** \brief Random constructor
	 *
	 * Calculates initial values for the regression using the provided responce and covariance matrix.  The initial values modify the corresponding row of the value matrix, stored in the encapsulating class.
	 *
	 * \param[in] gsl_matrix* response matrix
	 * \param[in] gsl_matrix* predictor matrix
	 * \param[in] size_t& column index of the current predictor in the predictor matrix
	 * \param[in] vector<double>& vectorized partial fitted matrix
	 * \param[in] gsl_matrix* covariance matrix
	 * \param[in] gsl_rng* pointer to a PNG
	 * \param[in] size_t& row index for the prior matrix
	 * \param[in] gsl_matrix* value matrix
	 * \param[in] size_t& row index of the value matrix
	 * \param[in] Apex& object storing the redundant parameter matrix
	 * \param[in] gsl_matrix* pre-computed auxiliary matrix
	 *
	 */
	MVnormBetaPEX(const gsl_matrix *resp, gsl_matrix *pred, const size_t &iCl, vector<double> &eaFt, const gsl_matrix *Sig, const gsl_rng *r, const size_t &up, gsl_matrix *bet, const size_t &iRw, Apex &A, gsl_matrix *tSigIAt);
	/** \brief Deterministic constructor
	 *
	 * Sets up the proper pointers to parameters of the regression that have already been calculated in the encapsulating class.  The corresponding value matrix is not modified.
	 *
	 * \param[in] gsl_matrix* response matrix
	 * \param[in] gsl_matrix* predictor matrix
	 * \param[in] size_t& column index of the current predictor in the predictor matrix
	 * \param[in] vector<double>& vectorized partial fitted matrix
	 * \param[in] size_t& row index for the prior matrix
	 * \param[in] gsl_matrix* value matrix
	 * \param[in] size_t& row index of the value matrix
	 * \param[in] Apex& object storing the redundant parameter matrix
	 * \param[in] gsl_matrix* pre-computed auxiliary matrix
	 *
	 */
	MVnormBetaPEX(const size_t &d, gsl_matrix *pred, const size_t &iCl, vector<double> &eaFt, const size_t &up, gsl_matrix *bet, const size_t &iRw, Apex &A, gsl_matrix *tSigIAt);
	
	/** \brief Deterministic copy constructor
	 *
	 * \param[in] MVnormBetaPEX& object of type MVnormBetaPEX
	 *
	 */
	MVnormBetaPEX(const MVnormBetaPEX &bet);
	/** \brief Assignment operator
	 *
	 * Deterministic assignement operator
	 *
	 * \param[in] MVnormBetaPEX& object of type MVnormBetaPEX
	 *
	 * \return A reference to an object of type MVnormBetaPEX
	 *
	 */
	MVnormBetaPEX & operator=(const MVnormBetaPEX &bet);
	
	/** \brief Destructor */
	~MVnormBetaPEX() {};
	
	void update(const Grp &dat, const SigmaI &SigIm, const SigmaI &SigIp, const gsl_rng *r);
	void update(const Grp &dat, const SigmaI &SigIm, const double &qPr, const SigmaI &SigIp, const gsl_rng *r);
	
	void update(const Grp &dat, const SigmaI &SigIm, const Grp &muPr, const SigmaI &SigIp, const gsl_rng *r);
	void update(const Grp &dat, const SigmaI &SigIm, const Grp &muPr, const double &qPr, const SigmaI &SigIp, const gsl_rng *r);
	
};

/** \brief Individual vector of means with missing data
 *
 * Implements missing phenotype data imputation.  Some parameters for update methods have a different meaning than for other MVnorm classes.
 *
 */
class MVnormMuMiss : public MVnormMu {
protected:
	/** \brief Missing data index
	 *
	 * Vector of indexes tagging missing data.
	 *
	 */
	vector<size_t> _misPhenInd;
	/** \brief own index
	 *
	 * Index of the row the current instance of this class is pointing to.  Only defined with constructors that take _iRw_.
	 *
	 */
	size_t _myInd;
	
public:
	/** \brief Default constructor */
	MVnormMuMiss() : MVnormMu(){};
	/** \brief 0-mean deterministic constructor
	 *
	 * Sets up the instance to point to a vector of zeros
	 *
	 * \param[in] size_t& dimension of the vector
	 *
	 */
	MVnormMuMiss(const size_t &d) : MVnormMu(d) {};
	/** \brief Deterministic constructor
	 *
	 * \param[in] size_t& dimension of the vector
	 * \param[in] vector<size_t>& vector of row indexes for data
	 * \param[in] size_t& row index for the prior matrix
	 * \param[in] vector<size_t>& indexes of missing data
	 *
	 */
	MVnormMuMiss(const size_t &d, const vector<size_t> &low, const size_t &up, const vector<size_t> &mis);
	/** \brief Deterministic constructor with a vector
	 *
	 * \param[in] gsl_vector* vector of values
	 * \param[in] size_t& dimension of the vector
	 * \param[in] vector<size_t>& vector of row indexes for data
	 * \param[in] size_t& row index for the prior matrix
	 * \param[in] vector<size_t>& indexes of missing data
	 *
	 */
	MVnormMuMiss(gsl_vector *mn, const vector<size_t> &low, const size_t &up, const vector<size_t> &mis);
	/** \brief Deterministic constructor with a vector
	 *
	 * \param[in] gsl_vector* vector of values
	 * \param[in] size_t& row index for the prior matrix
	 * \param[in] vector<size_t>& indexes of missing data
	 *
	 */
	MVnormMuMiss(gsl_vector *mn, const size_t &up, const vector<size_t> &mis);
	/** \brief Univariate random constructor with a vector
	 *
	 * \param[in] gsl_vector* vector of values
	 * \param[in] gsl_vector* vector of standard deviations
	 * \param[in] gsl_rng* pointer to a PNG
	 * \param[in] size_t& dimension of the vector
	 * \param[in] vector<size_t>& vector of row indexes for data
	 * \param[in] size_t& row index for the prior matrix
	 * \param[in] vector<size_t>& indexes of missing data
	 *
	 */
	MVnormMuMiss(gsl_vector *mn, const gsl_vector *sd, const gsl_rng *r, const vector<size_t> &low, const size_t &up, const vector<size_t> &mis);
	/** \brief Multivariate random constructor with a vector
	 *
	 * \param[in] gsl_vector* vector of values
	 * \param[in] gsl_matrix* covariance matrix
	 * \param[in] gsl_rng* pointer to a PNG
	 * \param[in] size_t& dimension of the vector
	 * \param[in] vector<size_t>& vector of row indexes for data
	 * \param[in] size_t& row index for the prior matrix
	 * \param[in] vector<size_t>& indexes of missing data
	 *
	 */
	MVnormMuMiss(gsl_vector *mn, const gsl_matrix *Sig, const gsl_rng *r, const vector<size_t> &low, const size_t &up, const vector<size_t> &mis);
	/** \brief Deterministic constructor with a matrix
	 *
	 * \param[in] gsl_matrix* vector of values
	 * \param[in] size_t& row index of the matrix of values
	 * \param[in] vector<size_t>& vector of row indexes for data
	 * \param[in] size_t& row index for the prior matrix
	 * \param[in] vector<size_t>& indexes of missing data
	 *
	 */
	MVnormMuMiss(gsl_matrix *mn, const size_t &iRw, const vector<size_t> &low, const size_t &up, const vector<size_t> &mis);
	/** \brief Deterministic constructor with a matrix
	 *
	 * \param[in] gsl_matrix* vector of values
	 * \param[in] size_t& row index of the matrix of values
	 * \param[in] size_t& row index for the prior matrix
	 * \param[in] vector<size_t>& indexes of missing data
	 *
	 */
	MVnormMuMiss(gsl_matrix *mn, const size_t &iRw, const size_t &up, const vector<size_t> &mis);
	/** \brief Univariate random constructor with a matrix
	 *
	 * \param[in] gsl_matrix* matrix of values
	 * \param[in] size_t& row index of the matrix of values
	 * \param[in] gsl_vector* vector of standard deviations
	 * \param[in] gsl_rng* pointer to a PNG
	 * \param[in] size_t& dimension of the vector
	 * \param[in] vector<size_t>& vector of row indexes for data
	 * \param[in] size_t& row index for the prior matrix
	 * \param[in] vector<size_t>& indexes of missing data
	 *
	 */
	MVnormMuMiss(gsl_matrix *mn, const size_t &iRw, const gsl_vector *sd, const gsl_rng *r, const vector<size_t> &low, const size_t &up, const vector<size_t> &mis);
	/** \brief Multivariate random constructor with a matrix
	 *
	 * \param[in] gsl_matrix* matrix of values
	 * \param[in] size_t& row index of the matrix of values
	 * \param[in] gsl_matrix* covariance matrix
	 * \param[in] gsl_rng* pointer to a PNG
	 * \param[in] size_t& dimension of the vector
	 * \param[in] vector<size_t>& vector of row indexes for data
	 * \param[in] size_t& row index for the prior matrix
	 * \param[in] vector<size_t>& indexes of missing data
	 *
	 */
	MVnormMuMiss(gsl_matrix *mn, const size_t &iRw, const gsl_matrix *Sig, const gsl_rng *r, const vector<size_t> &low, const size_t &up, const vector<size_t> &mis);
	
	/** \brief Deterministic copy constructor
	 *
	 * \param[in] MVnormMuMiss& object of type MVnormMuMiss
	 *
	 */
	MVnormMuMiss(const MVnormMuMiss &); // copy constructor
	/** \brief Assignment operator
	 *
	 * Deterministic assignement operator
	 *
	 * \param[in] MVnormMuMiss& object of type MVnormMuMiss
	 *
	 * \return A reference to an object of type MVnormMuMiss
	 *
	 */
	MVnormMuMiss & operator=(const MVnormMuMiss &);
	
	/** \brief Destructor */
	~MVnormMuMiss() {};
	
	/** \brief Gaussian likelihood
	 *
	 * Here, the Grp variable contains the means for imputation rather than the data
	 *
	 * \param[in] Grp& data
	 * \param[in] SigmaI& inverse-covariance matrix for the likelihood
	 * \param[in] gsl_rng* pointer to a PNG
	 */
	void update(const Grp &mu, const SigmaI &SigIm, const gsl_rng *r);
	/** \brief Gaussian likelihood, Gaussian 0-mean prior
	 *
	 * Here, the Grp variable contains the means for imputation rather than the data
	 *
	 * \param[in] Grp& data
	 * \param[in] SigmaI& inverse-covariance matrix for the likelihood
	 * \param[in] SigmaI& inverse-covariance matrix for the prior
	 * \param[in] gsl_rng* pointer to a PNG
	 */
	void update(const Grp &mu, const SigmaI &SigIm, const SigmaI &SigIp, const gsl_rng *r);
	
	/** \brief Number of missing values
	 *
	 * Accesses the  number of missing phenotype values.  Implemented in this class.
	 *
	 * \return Number of missing phenotypes, of type _size_t_.
	 */
	size_t nMissP() const{return _misPhenInd.size(); };
	/** \brief Indexes of missing values
	 *
	 * Accesses a vector with indexes of missing phenotypes.  Implemented for this class.
	 *
	 * \return Vector of _size_t_ that contains indexes that correspond to missing values.
	 *
	 */
	const vector<size_t> getMisPhen() const{ return _misPhenInd; };
};

/** \brief Generic regression
 *
 * Generic implementation of a regression. By itself to be used strictly for regressions with a single predictor because it ignores the effects of other predictors.
 *
 */
class MVnormBeta : public MVnorm {
protected:
	/** \brief Predictor
	 *
	 * Points to a vector of predictor values, typically a column of the predictor matrix that is in the encapsulating Grp class.
	 *
	 */
	gsl_vector_view _X;
	/** \brief Scale parameter
	 *
	 * Typically \f$ x^{\textsc{t}}x \f$, but some exceptions exist in special cases.
	 *
	 */
	double _scale;
	/** \brief Length of the predictor */
	size_t _N;
	
	/** \brief Row index of the prior
	 *
	 * Points to the index of a row in the prior matrix
	 *
	 */
	const size_t *_upLevel; // index into the vector of priors
	
public:
	/** \brief Default constructor */
	MVnormBeta();
	/** \brief Dimension-only constructor
	 *
	 * Sets the dimension to a value, but the _vec is not assigned a target.
	 *
	 */
	MVnormBeta(const size_t d);
	
	/** \brief Univariate random constructor with vector
	 *
	 * Sets up _vec to point to the provided vector of values, and the predictor to a column in the provided matrix of predictors.
	 *
	 * \param[in] gsl_vector* vector of values
	 * \param[in] gsl_vector* vector of standard deviations
	 * \param[in] gsl_matrix* matrix of predictors
	 * \param[in] size_t& column index for the predictor matrix
	 * \param[in] gsl_rng* pointer to a PNG
	 *
	 */
	MVnormBeta(gsl_vector *b, const gsl_vector *sd, gsl_matrix *pred, const size_t &iCl, const gsl_rng *r);
	/** \brief Univariate random constructor with vector and pointer to prior
	 *
	 * Sets up _vec to point to the provided vector of values, the predictor to a column in the provided matrix of predictors, and the pointer to the corresponding row in the prior matrix.
	 *
	 * \param[in] gsl_vector* vector of values
	 * \param[in] gsl_vector* vector of standard deviations
	 * \param[in] gsl_matrix* matrix of predictors
	 * \param[in] size_t& column index for the predictor matrix
	 * \param[in] gsl_rng* pointer to a PNG
	 * \param[in] size_t& index of the appropriate row in the prior matrix
	 *
	 */
	MVnormBeta(gsl_vector *b, const gsl_vector *sd, gsl_matrix *pred, const size_t &iCl, const gsl_rng *r, const size_t &up);
	
	/** \brief Multivariate random constructor with vector
	 *
	 * Sets up _vec to point to the provided vector of values, and the predictor to a column in the provided matrix of predictors.
	 *
	 * \param[in] gsl_vector* vector of values
	 * \param[in] gsl_matrix* covariance matrix
	 * \param[in] gsl_matrix* matrix of predictors
	 * \param[in] size_t& column index for the predictor matrix
	 * \param[in] gsl_rng* pointer to a PNG
	 *
	 */
	MVnormBeta(gsl_vector *b, const gsl_matrix *Sig, gsl_matrix *pred, const size_t &iCl, const gsl_rng *r);
	/** \brief Multivariate random constructor with vector and pointer to prior
	 *
	 * Sets up _vec to point to the provided vector of values, the predictor to a column in the provided matrix of predictors, and the pointer to the corresponding row in the prior matrix.
	 *
	 * \param[in] gsl_vector* vector of values
	 * \param[in] gsl_matrix* covariance matrix
	 * \param[in] gsl_matrix* matrix of predictors
	 * \param[in] size_t& column index for the predictor matrix
	 * \param[in] gsl_rng* pointer to a PNG
	 * \param[in] size_t& index of the appropriate row in the prior matrix
	 *
	 */
	MVnormBeta(gsl_vector *b, const gsl_matrix *Sig, gsl_matrix *pred, const size_t &iCl, const gsl_rng *r, const size_t &up);
	
	/** \brief Multivariate constructor with response matrix
	 *
	 * Sets up _vec to point to a row in the provided matrix of regression coefficients, and the predictor to a column of the matrix of predictors.  Initializes regression coefficients using the provided response matrix.
	 *
	 * \param[in] gsl_matrix* response matrix
	 * \param[in] gsl_matrix* predictor martix
	 * \param[in] size_t& column index of the predictor matrix
	 * \param[in] gsl_matrix* covariance matrix
	 * \param[in] gsl_rng* pointer to a PNG
	 * \param[in] gsl_matrix* matrix of regression coefficients
	 * \param[in] size_t& row index of the regression coefficient matrix
	 *
	 */
	MVnormBeta(const gsl_matrix *resp, gsl_matrix *pred, const size_t &iCl, const gsl_matrix *Sig, const gsl_rng *r, gsl_matrix *bet, const size_t &iRw);
	/** \brief Multivariate constructor with response matrix and index for a prior
	 *
	 * Sets up _vec to point to a row in the provided matrix of regression coefficients, and the predictor to a column of the matrix of predictors.  Initializes regression coefficients using the provided response matrix.  
	 * Initializes _upLevel to point to the appropriate row in the matrix of priors.
	 *
	 * \param[in] gsl_matrix* response matrix
	 * \param[in] gsl_matrix* predictor martix
	 * \param[in] size_t& column index of the predictor matrix
	 * \param[in] gsl_matrix* covariance matrix
	 * \param[in] gsl_rng* pointer to a PNG
	 * \param[in] size_t& row index of the prior matrix
	 * \param[in] gsl_matrix* matrix of regression coefficients
	 * \param[in] size_t& row index of the regression coefficient matrix
	 *
	 */
	MVnormBeta(const gsl_matrix *resp, gsl_matrix *pred, const size_t &iCl, const gsl_matrix *Sig, const gsl_rng *r, const size_t &up, gsl_matrix *bet, const size_t &iRw);
	
	/** \brief Multivariate constructor with a Grp type response
	 *
	 * Sets up _vec to point to a row in the provided matrix of regression coefficients, and the predictor to a column of the matrix of predictors.  Initializes regression coefficients using the provided response matrix.
	 *
	 * \param[in] Grp& response
	 * \param[in] gsl_matrix* predictor martix
	 * \param[in] size_t& column index of the predictor matrix
	 * \param[in] gsl_matrix* covariance matrix
	 * \param[in] gsl_rng* pointer to a PNG
	 * \param[in] gsl_matrix* matrix of regression coefficients
	 * \param[in] size_t& row index of the regression coefficient matrix
	 *
	 */
	MVnormBeta(const Grp &resp, gsl_matrix *pred, const size_t &iCl, const gsl_matrix *Sig, const gsl_rng *r, gsl_matrix *bet, const size_t &iRw);
	/** \brief Multivariate constructor with a Grp type response and index for a prior
	 *
	 * Sets up _vec to point to a row in the provided matrix of regression coefficients, and the predictor to a column of the matrix of predictors.  Initializes regression coefficients using the provided response matrix.
	 * Initializes _upLevel to point to the appropriate row in the matrix of priors.
	 *
	 * \param[in] Grp& response
	 * \param[in] gsl_matrix* predictor martix
	 * \param[in] size_t& column index of the predictor matrix
	 * \param[in] gsl_matrix* covariance matrix
	 * \param[in] gsl_rng* pointer to a PNG
	 * \param[in] size_t& row index of the prior matrix
	 * \param[in] gsl_matrix* matrix of regression coefficients
	 * \param[in] size_t& row index of the regression coefficient matrix
	 *
	 */
	MVnormBeta(const Grp &resp, gsl_matrix *pred, const size_t &iCl, const gsl_matrix *Sig, const gsl_rng *r, const size_t &up, gsl_matrix *bet, const size_t &iRw);
	
	/** \brief Deterministic constructor
	 *
	 * Does not initialize the regression coefficients, but simply points to the already-initialized matrix of values.
	 *
	 * \param[in] gsl_matrix* predictor matrix
	 * \param[in] size_t& column index of the predictor matrix
	 * \param[in] gsl_matrix* regression coefficient matrix
	 * \param[in] size_t& row index of the regression coefficient matrix
	 *
	 */
	MVnormBeta(gsl_matrix *pred, const size_t &iCl, gsl_matrix *bet, const size_t &iRw);
	/** \brief Deterministic constructor with prior index
	 *
	 * Does not initialize the regression coefficients, but simply points to the already-initialized matrix of values.
	 * Sets the _upLevel to point to the corresponding row of the prior matrix
	 *
	 * \param[in] gsl_matrix* predictor matrix
	 * \param[in] size_t& column index of the predictor matrix
	 * \param[in] size_t& row index of the prior matrix
	 * \param[in] gsl_matrix* regression coefficient matrix
	 * \param[in] size_t& row index of the regression coefficient matrix
	 *
	 */
	MVnormBeta(gsl_matrix *pred, const size_t &iCl, const size_t &up, gsl_matrix *bet, const size_t &iRw);
	
	/** \brief Deterministic copy constructor
	 *
	 * \param[in] MVnormBeta& object to be copied
	 *
	 */
	MVnormBeta(const MVnormBeta &);
	/** \brief Deterministic assignement operator
	 *
	 * \param[in] MVnormBeta& object to be assigned
	 *
	 * \return A reference to a MVnormBeta object
	 *
	 */
	MVnormBeta & operator=(const MVnormBeta &);
	
	/** \brief Destructor */
	virtual ~MVnormBeta(); 	// destructor
	
	// improper prior methods
	virtual void update(const Grp &dat, const SigmaI &SigIb, const gsl_rng *r);
	virtual void update(const Grp &dat, const Qgrp &q, const SigmaI &SigIb, const gsl_rng *r);
	// 0-mean prior methods
	virtual void update(const Grp &dat, const SigmaI &SigIb, const SigmaI &SigIp, const gsl_rng *r); // G data, G prior
	virtual void update(const Grp &dat, const SigmaI &SigIb, const double &qPr, const SigmaI &SigIp, const gsl_rng *r); // G data, t prior
	virtual void update(const Grp &dat, const Qgrp &q, const SigmaI &SigIb, const SigmaI &SigIp, const gsl_rng *r); // t data, G prior
	virtual void update(const Grp &dat, const Qgrp &q, const SigmaI &SigIb, const double &qPr, const SigmaI &SigIp, const gsl_rng *r); // t data, t prior
	// non-zero mean prrior methods
	virtual void update(const Grp &dat, const SigmaI &SigIb, const Grp &muPr, const SigmaI &SigIp, const gsl_rng *r); // G data, G prior
	virtual void update(const Grp &dat, const SigmaI &SigIb, const Grp &muPr, const double &qPr, const SigmaI &SigIp, const gsl_rng *r); // G data, t prior
	virtual void update(const Grp &dat, const Qgrp &q, const SigmaI &SigIb, const Grp &muPr, const SigmaI &SigIp, const gsl_rng *r); // t data, G prior
	virtual void update(const Grp &dat, const Qgrp &q, const SigmaI &SigIb, const Grp &muPr, const double &qPr, const SigmaI &SigIp, const gsl_rng *r); // t data, t prior
	
	const size_t *up() const{return _upLevel; };
	double scalePar() const {return _scale; };
};

/** \brief Regression with multiple predictors
 *
 * Full implementation of a multiple regression using one-at-a-time updating. Updates regression coefficients for the current predictor, while accounting for the effects of other predictors
 *
 */
class MVnormBetaFt: public MVnormBeta {
protected:
	/** \brief Matrix of already fitted values
	 *
	 * \f$ \boldsymbol{XB} \f$ matrix for all predictors but the current one (i.e., \f$ \boldsymbol{X}_{\cdot -j}\boldsymbol{B}_{-j\cdot} \f$).  It is subtracted from the response and the regression is performed on the residual.
	 *
	 */
	gsl_matrix_view _fitted;
	
public:
	/** \brief Default constructor */
	MVnormBetaFt(): MVnormBeta() {};
	
	/** \brief Univariate random constructor with vector
	 *
	 * Sets up _vec to point to the provided vector of values, and the predictor to a column in the provided matrix of predictors.  
	 * The _fitted member points to the indicated submatrix of a matrix that contains all the partial fitted matrices.
	 *
	 * \param[in] gsl_vector* vector of values
	 * \param[in] gsl_vector* vector of standard deviations
	 * \param[in] gsl_matrix* matrix of predictors
	 * \param[in] size_t& column index for the predictor matrix
	 * \param[in] gsl_matrix* matrix with all the partial fitted matrices stacked
	 * \param[in] size_t& index of the row where the relevant fitted matrix begins
	 * \param[in] gsl_rng* pointer to a PNG
	 *
	 */
	MVnormBetaFt(gsl_vector *b, const gsl_vector *sd, gsl_matrix *pred, const size_t &iCl, gsl_matrix *allFt, const size_t &begRw, const gsl_rng *r) : MVnormBeta(b, sd, pred, iCl, r) {_fitted = gsl_matrix_submatrix(allFt, begRw, 0, pred->size1, _d); };
	/** \brief Univariate random constructor with vector and prior index
	 *
	 * Sets up _vec to point to the provided vector of values, and the predictor to a column in the provided matrix of predictors.
	 * The _fitted member points to the indicated submatrix of a matrix that contains all the partial fitted matrices. The _upLevel member is set to point to a row in the prior matrix.
	 *
	 * \param[in] gsl_vector* vector of values
	 * \param[in] gsl_vector* vector of standard deviations
	 * \param[in] gsl_matrix* matrix of predictors
	 * \param[in] size_t& column index for the predictor matrix
	 * \param[in] gsl_matrix* matrix with all the partial fitted matrices stacked
	 * \param[in] size_t& index of the row where the relevant fitted matrix begins
	 * \param[in] gsl_rng* pointer to a PNG
	 * \param[in] size_t& row index of the prior matrix
	 *
	 */
	MVnormBetaFt(gsl_vector *b, const gsl_vector *sd, gsl_matrix *pred, const size_t &iCl, gsl_matrix *allFt, const size_t &begRw, const gsl_rng *r, const size_t &up) : MVnormBeta(b, sd, pred, iCl, r, up) {_fitted = gsl_matrix_submatrix(allFt, begRw, 0, pred->size1, _d); };
	
	/** \brief Multivariate random constructor with vector
	 *
	 * Sets up _vec to point to the provided vector of values, and the predictor to a column in the provided matrix of predictors.
	 * The _fitted member points to the indicated submatrix of a matrix that contains all the partial fitted matrices.
	 *
	 * \param[in] gsl_vector* vector of values
	 * \param[in] gsl_matrix* covariance matrix
	 * \param[in] gsl_matrix* matrix of predictors
	 * \param[in] size_t& column index for the predictor matrix
	 * \param[in] gsl_matrix* matrix with all the partial fitted matrices stacked
	 * \param[in] size_t& index of the row where the relevant fitted matrix begins
	 * \param[in] gsl_rng* pointer to a PNG
	 *
	 */
	MVnormBetaFt(gsl_vector *b, const gsl_matrix *Sig, gsl_matrix *pred, const size_t &iCl, gsl_matrix *allFt, const size_t &begRw, const gsl_rng *r) : MVnormBeta(b, Sig, pred, iCl, r) {_fitted = gsl_matrix_submatrix(allFt, begRw, 0, pred->size1, _d); };
	/** \brief Multivariate random constructor with vector and prior index
	 *
	 * Sets up _vec to point to the provided vector of values, and the predictor to a column in the provided matrix of predictors.
	 * The _fitted member points to the indicated submatrix of a matrix that contains all the partial fitted matrices. The _upLevel member is set to point to a row in the prior matrix.
	 *
	 * \param[in] gsl_vector* vector of values
	 * \param[in] gsl_matrix* covariance matrix
	 * \param[in] gsl_matrix* matrix of predictors
	 * \param[in] size_t& column index for the predictor matrix
	 * \param[in] gsl_matrix* matrix with all the partial fitted matrices stacked
	 * \param[in] size_t& index of the row where the relevant fitted matrix begins
	 * \param[in] gsl_rng* pointer to a PNG
	 * \param[in] size_t& row index of the prior matrix
	 *
	 */
	MVnormBetaFt(gsl_vector *b, const gsl_matrix *Sig, gsl_matrix *pred, const size_t &iCl, gsl_matrix *allFt, const size_t &begRw, const gsl_rng *r, const size_t &up) : MVnormBeta(b, Sig, pred, iCl, r, up) {_fitted = gsl_matrix_submatrix(allFt, begRw, 0, pred->size1, _d); };
	
	/** \brief Multivariate random constructor with response matrix
	 *
	 * Sets up _vec to point to the provided vector of values, and the predictor to a column in the provided matrix of predictors.   Initializes regression coefficients using the provided response matrix.
	 * The _fitted member points to the indicated submatrix of a matrix that contains all the partial fitted matrices.
	 *
	 * \param[in] gsl_matrix* response matrix
	 * \param[in] gsl_matrix* predictor matrix
	 * \param[in] size_t& column index for the predictor matrix
	 * \param[in] gsl_matrix* matrix with all the partial fitted matrices stacked
	 * \param[in] size_t& index of the row where the relevant fitted matrix begins
	 * \param[in] gsl_matrix* covariance matrix
	 * \param[in] gsl_rng* pointer to a PNG
	 * \param[in] gsl_matrix* matrix of regression coefficients
	 * \param[in] size_t& row index of the regression coefficient matrix
	 *
	 */
	MVnormBetaFt(const gsl_matrix *resp, gsl_matrix *pred, const size_t &iCl, gsl_matrix *allFt, const size_t &begRw, const gsl_matrix *Sig, const gsl_rng *r, gsl_matrix *bet, const size_t &iRw) : MVnormBeta(resp, pred, iCl, Sig, r, bet, iRw) {_fitted = gsl_matrix_submatrix(allFt, begRw, 0, pred->size1, _d); };
	/** \brief Multivariate random constructor with response matrix and prior index
	 *
	 * Sets up _vec to point to the provided vector of values, and the predictor to a column in the provided matrix of predictors. Initializes regression coefficients using the provided response matrix.
	 * The _fitted member points to the indicated submatrix of a matrix that contains all the partial fitted matrices.  The _upLevel member is set to point to a row in the prior matrix.
	 *
	 * \param[in] gsl_matrix* response matrix
	 * \param[in] gsl_matrix* predictor matrix
	 * \param[in] size_t& column index for the predictor matrix
	 * \param[in] gsl_matrix* matrix with all the partial fitted matrices stacked
	 * \param[in] size_t& index of the row where the relevant fitted matrix begins
	 * \param[in] gsl_matrix* covariance matrix
	 * \param[in] gsl_rng* pointer to a PNG
	 * \param[in] size_t& row index of the prior matrix
	 * \param[in] gsl_matrix* matrix of regression coefficients
	 * \param[in] size_t& row index of the regression coefficient matrix
	 *
	 */
	MVnormBetaFt(const gsl_matrix *resp, gsl_matrix *pred, const size_t &iCl, gsl_matrix *allFt, const size_t &begRw, const gsl_matrix *Sig, const gsl_rng *r, const size_t &up, gsl_matrix *bet, const size_t &iRw) : MVnormBeta(resp, pred, iCl, Sig, r, up, bet, iRw) {_fitted = gsl_matrix_submatrix(allFt, begRw, 0, pred->size1, _d); };
	
	/** \brief Multivariate random constructor with response matrix
	 *
	 * Sets up _vec to point to the provided vector of values, and the predictor to a column in the provided matrix of predictors. Initializes regression coefficients using the provided response matrix.
	 * The _fitted member points to the indicated vector that has the appropriate partial fitted values.
	 *
	 * \param[in] gsl_matrix* response matrix
	 * \param[in] gsl_matrix* predictor matrix
	 * \param[in] size_t& column index for the predictor matrix
	 * \param[in] vector<double>& vectorized partial fitted matrix
	 * \param[in] gsl_matrix* covariance matrix
	 * \param[in] gsl_rng* pointer to a PNG
	 * \param[in] gsl_matrix* matrix of regression coefficients
	 * \param[in] size_t& row index of the regression coefficient matrix
	 *
	 */
	MVnormBetaFt(const gsl_matrix *resp, gsl_matrix *pred, const size_t &iCl, vector<double> &eaFt, const gsl_matrix *Sig, const gsl_rng *r, gsl_matrix *bet, const size_t &iRw) : MVnormBeta(resp, pred, iCl, Sig, r, bet, iRw) {_fitted = gsl_matrix_view_array(eaFt.data(), pred->size1, _d); };
	/** \brief Multivariate random constructor with response matrix and prior index
	 *
	 * Sets up _vec to point to the provided vector of values, and the predictor to a column in the provided matrix of predictors. Initializes regression coefficients using the provided response matrix.
	 * The _fitted member points to the indicated vector that has the appropriate partial fitted values. The _upLevel member is set to point to a row in the prior matrix.
	 *
	 * \param[in] gsl_matrix* response matrix
	 * \param[in] gsl_matrix* predictor matrix
	 * \param[in] size_t& column index for the predictor matrix
	 * \param[in] vector<double>& vectorized partial fitted matrix
	 * \param[in] gsl_matrix* covariance matrix
	 * \param[in] gsl_rng* pointer to a PNG
	 * \param[in] size_t& row index of the prior matrix
	 * \param[in] gsl_matrix* matrix of regression coefficients
	 * \param[in] size_t& row index of the regression coefficient matrix
	 *
	 */
	MVnormBetaFt(const gsl_matrix *resp, gsl_matrix *pred, const size_t &iCl, vector<double> &eaFt, const gsl_matrix *Sig, const gsl_rng *r, const size_t &up, gsl_matrix *bet, const size_t &iRw) : MVnormBeta(resp, pred, iCl, Sig, r, up, bet, iRw) {_fitted = gsl_matrix_view_array(eaFt.data(), pred->size1, _d); };
	
	/** \brief Multivariate random constructor with Grp type response
	 *
	 * Sets up _vec to point to the provided vector of values, and the predictor to a column in the provided matrix of predictors.   Initializes regression coefficients using the provided response matrix.
	 * The _fitted member points to the indicated submatrix of a matrix that contains all the partial fitted matrices.
	 *
	 * \param[in] Grp& response
	 * \param[in] gsl_matrix* predictor matrix
	 * \param[in] size_t& column index for the predictor matrix
	 * \param[in] gsl_matrix* matrix with all the partial fitted matrices stacked
	 * \param[in] size_t& index of the row where the relevant fitted matrix begins
	 * \param[in] gsl_matrix* covariance matrix
	 * \param[in] gsl_rng* pointer to a PNG
	 * \param[in] gsl_matrix* matrix of regression coefficients
	 * \param[in] size_t& row index of the regression coefficient matrix
	 *
	 */
	MVnormBetaFt(const Grp &resp, gsl_matrix *pred, const size_t &iCl, gsl_matrix *allFt, const size_t &begRw, const gsl_matrix *Sig, const gsl_rng *r, gsl_matrix *bet, const size_t &iRw) : MVnormBeta(resp, pred, iCl, Sig, r, bet, iRw) {_fitted = gsl_matrix_submatrix(allFt, begRw, 0, pred->size1, _d); }; // deprecated
	/** \brief Multivariate random constructor with Grp type response and prior index
	 *
	 * Sets up _vec to point to the provided vector of values, and the predictor to a column in the provided matrix of predictors. Initializes regression coefficients using the provided response matrix.
	 * The _fitted member points to the indicated submatrix of a matrix that contains all the partial fitted matrices.  The _upLevel member is set to point to a row in the prior matrix.
	 *
	 * \param[in] Grp& response
	 * \param[in] gsl_matrix* predictor matrix
	 * \param[in] size_t& column index for the predictor matrix
	 * \param[in] gsl_matrix* matrix with all the partial fitted matrices stacked
	 * \param[in] size_t& index of the row where the relevant fitted matrix begins
	 * \param[in] gsl_matrix* covariance matrix
	 * \param[in] gsl_rng* pointer to a PNG
	 * \param[in] size_t& row index of the prior matrix
	 * \param[in] gsl_matrix* matrix of regression coefficients
	 * \param[in] size_t& row index of the regression coefficient matrix
	 *
	 */
	MVnormBetaFt(const Grp &resp, gsl_matrix *pred, const size_t &iCl, gsl_matrix *allFt, const size_t &begRw, const gsl_matrix *Sig, const gsl_rng *r, const size_t &up, gsl_matrix *bet, const size_t &iRw) : MVnormBeta(resp, pred, iCl, Sig, r, up, bet, iRw) {_fitted = gsl_matrix_submatrix(allFt, begRw, 0, pred->size1, _d); };
	
	/** \brief Multivariate random constructor with Grp type response
	 *
	 * Sets up _vec to point to the provided vector of values, and the predictor to a column in the provided matrix of predictors. Initializes regression coefficients using the provided response matrix.
	 * The _fitted member points to the indicated vector that has the appropriate partial fitted values.
	 *
	 * \param[in] Grp& response
	 * \param[in] gsl_matrix* predictor matrix
	 * \param[in] size_t& column index for the predictor matrix
	 * \param[in] vector<double>& vectorized partial fitted matrix
	 * \param[in] gsl_matrix* covariance matrix
	 * \param[in] gsl_rng* pointer to a PNG
	 * \param[in] gsl_matrix* matrix of regression coefficients
	 * \param[in] size_t& row index of the regression coefficient matrix
	 *
	 */
	MVnormBetaFt(const Grp &resp, gsl_matrix *pred, const size_t &iCl, vector<double> &eaFt, const gsl_matrix *Sig, const gsl_rng *r, gsl_matrix *bet, const size_t &iRw) : MVnormBeta(resp, pred, iCl, Sig, r, bet, iRw) {_fitted = gsl_matrix_view_array(eaFt.data(), pred->size1, _d); };
	/** \brief Multivariate random constructor with Grp type response and prior index
	 *
	 * Sets up _vec to point to the provided vector of values, and the predictor to a column in the provided matrix of predictors. Initializes regression coefficients using the provided response matrix.
	 * The _fitted member points to the indicated vector that has the appropriate partial fitted values. The _upLevel member is set to point to a row in the prior matrix.
	 *
	 * \param[in] Grp& response
	 * \param[in] gsl_matrix* predictor matrix
	 * \param[in] size_t& column index for the predictor matrix
	 * \param[in] vector<double>& vectorized partial fitted matrix
	 * \param[in] gsl_matrix* covariance matrix
	 * \param[in] gsl_rng* pointer to a PNG
	 * \param[in] size_t& row index of the prior matrix
	 * \param[in] gsl_matrix* matrix of regression coefficients
	 * \param[in] size_t& row index of the regression coefficient matrix
	 *
	 */
	MVnormBetaFt(const Grp &resp, gsl_matrix *pred, const size_t &iCl, vector<double> &eaFt, const gsl_matrix *Sig, const gsl_rng *r, const size_t &up, gsl_matrix *bet, const size_t &iRw) : MVnormBeta(resp, pred, iCl, Sig, r, up, bet, iRw) {_fitted = gsl_matrix_view_array(eaFt.data(), pred->size1, _d); };
	
	/** \brief Deterministic constructor
	 *
	 * Does not initialize the regression coefficients, but simply points to the already-initialized matrix of values.
	 * The _fitted member points to the indicated vector that has the appropriate partial fitted values.
	 *
	 * \param[in] gsl_matrix* predictor matrix
	 * \param[in] size_t& column index of the predictor matrix
	 * \param[in] vector<double>& vectorized partial fitted matrix
	 * \param[in] gsl_matrix* regression coefficient matrix
	 * \param[in] size_t& row index of the regression coefficient matrix
	 *
	 */
	MVnormBetaFt(gsl_matrix *pred, const size_t &iCl, vector<double> &eaFt, gsl_matrix *bet, const size_t &iRw) : MVnormBeta(pred, iCl, bet, iRw) {_fitted = gsl_matrix_view_array(eaFt.data(), pred->size1, _d); };
	/** \brief Deterministic constructor with a prior index
	 *
	 * Does not initialize the regression coefficients, but simply points to the already-initialized matrix of values.
	 * The _fitted member points to the indicated vector that has the appropriate partial fitted values.  The _upLevel member is set to point to a row in the prior matrix.
	 *
	 * \param[in] gsl_matrix* predictor matrix
	 * \param[in] size_t& column index of the predictor matrix
	 * \param[in] vector<double>& vectorized partial fitted matrix
	 * \param[in] size_t& row index of the prior matrix
	 * \param[in] gsl_matrix* regression coefficient matrix
	 * \param[in] size_t& row index of the regression coefficient matrix
	 *
	 */
	MVnormBetaFt(gsl_matrix *pred, const size_t &iCl, vector<double> &eaFt, const size_t &up, gsl_matrix *bet, const size_t &iRw) : MVnormBeta(pred, iCl, up, bet, iRw) {_fitted = gsl_matrix_view_array(eaFt.data(), pred->size1, _d); };
	
	/** \brief Deterministic copy constructor
	 *
	 * \param[in] MVnormBetaFt& object to be copied
	 *
	 */
	MVnormBetaFt(const MVnormBetaFt &);
	/** \brief Deterministic assignment operator
	 *
	 * \param[in] MVnormBetaFt& object to be copied
	 *
	 * \return Refernce to an object of class MVnormBetaFt
	 *
	 */
	MVnormBetaFt & operator=(const MVnormBetaFt &);
	
	/** \brief Destructor */
	virtual ~MVnormBetaFt() {};
	
	// update() methods
	// Gaussian model, improper prior (for "fixed effect" regression)
	void update(const Grp &resp, const SigmaI &SigIb, const gsl_rng *r);
	// Student-t model, improper prior (for "fixed effect" regression).  To use for replicate-level covariate regression
	void update(const Grp &resp, const Qgrp &q, const SigmaI &SigIb, const gsl_rng *r);
	
	// Gaussian model, 0-mean Gaussian prior
	void update(const Grp &resp, const SigmaI &SigIb, const SigmaI &SigIp, const gsl_rng *r);
	// Gaussian model, 0-mean Student-t prior
	void update(const Grp &resp, const SigmaI &SigIb, const double &qPr, const SigmaI &SigIp, const gsl_rng *r);
	// Student-t model, 0-mean Gaussian prior
	void update(const Grp &resp, const Qgrp &q, const SigmaI &SigIb, const SigmaI &SigIp, const gsl_rng *r);
	// Student-t model, 0-mean Student-t prior
	void update(const Grp &resp, const Qgrp &q, const SigmaI &SigIb, const double &qPr, const SigmaI &SigIp, const gsl_rng *r);
	
	// Gaussian model, Gaussian prior
	void update(const Grp &resp, const SigmaI &SigIb, const Grp &muPr, const SigmaI &SigIp, const gsl_rng *r);
	// Gaussian model, Student-t prior
	void update(const Grp &resp, const SigmaI &SigIb, const Grp &muPr, const double &qPr, const SigmaI &SigIp, const gsl_rng *r);
	// Student-t model, Gaussian prior
	void update(const Grp &resp, const Qgrp &q, const SigmaI &SigIb, const Grp &muPr, const SigmaI &SigIp, const gsl_rng *r);
	// Student-t model, Student-t prior
	void update(const Grp &resp, const Qgrp &q, const SigmaI &SigIb, const Grp &muPr, const double &qPr, const SigmaI &SigIp, const gsl_rng *r);
	
};
/** \brief Individual vector of means with blocks of traits
 *
 * Implements separate models for blocks of traits.  The likelihood covariance matrix is block-diagonal.  Traits of the same block have to be contiguous within the vector.
 *
 */
class MVnormMuBlk : public MVnorm {
protected:
	/** \brief Pointer to a row index of the prior
	 *
	 * This has to be the same for all the blocks.
	 *
	 */
	const size_t *_upLevel;
	/** \brief Start positions of blocks
	 *
	 * Pointer to the vector that's in the corresponding Grp class.
	 *
	 */
	const vector< size_t > *_blkStart;
	/** \brief Trait blocks
	 *
	 * Vector views of vector pieces, each corresponding to a block of variables. Pointer to the vector that's in the corresponding Grp class.
	 *
	 */
	vector< gsl_vector_view > _eachVec;
	/** \brief Data matrix row indexes
	 *
	 * Each block of traits can have a different number and identity of rows in the data matrix it refers to.  Each vector of *size_t* in this vector corresponds to a block, so the size of this should be the same as the size of *_eachVec*.
	 * This member is a pointer to the vector that's in the corresponding Grp class.
	 *
	 */
	const vector< vector<size_t> > *_eachLL;
	
public:
	/** \brief Default constructor */
	MVnormMuBlk() : MVnorm() {};
	/** \brief Deterministic constructor
	 *
	 * Sets up the member variables to point to blocks of traits in the matrix _mn_, but does not perform stochastic intitialization.
	 *
	 * \param[in] gsl_matrix* mean values matrix
	 * \param[in] size_t& row index of the mean values matrix
	 * \param[in] vector<size_t>& vector of start indexes for each block
	 * \param[in] size_t& row index of the prior matrix
	 *
	 */
	MVnormMuBlk(gsl_matrix *mn, const size_t &iRw, const vector<size_t> &blkStart, const size_t &up);
	/** \brief Univariate stochastic constructor
	 *
	 * Sets initial values independently for each trait, modifying the matrix row that corresponds to this vector.
	 *
	 * \param[in] gsl_matrix* mean values matrix
	 * \param[in] size_t& row index of the mean values matrix
	 * \param[in] gsl_vector* vector of standard deviations
	 * \param[in] gsl_rng* pointer to a PNG
	 * \param[in] vector<size_t>& vector of start indexes for each block
	 * \param[in] vector< vector<size_t> >& vector of lower-level index vectors
	 * \param[in] size_t& row index of the prior matrix
	 *
	 */
	MVnormMuBlk(gsl_matrix *mn, const size_t &iRw, const gsl_vector *sd, const gsl_rng *r, const vector<size_t> &blkStart, const vector< vector<size_t> > &eachLL, const size_t &up);
	
	/** \brief Destructor */
	virtual ~MVnormMuBlk() {};
	
	const size_t *up() const{return _upLevel; };
	
	void update(const Grp &dat, const SigmaI &SigIm, const gsl_rng *r);
	void update(const Grp &dat, const Qgrp &q, const SigmaI &SigIm, const gsl_rng *r);
	// 0-mean prior methods
	void update(const Grp &dat, const SigmaI &SigIm, const SigmaI &SigIp, const gsl_rng *r);
	void update(const Grp &dat, const SigmaI &SigIm, const double &qPr, const SigmaI &SigIp, const gsl_rng *r);
	void update(const Grp &dat, const Qgrp &q, const SigmaI &SigIm, const SigmaI &SigIp, const gsl_rng *r);
	void update(const Grp &dat, const Qgrp &q, const SigmaI &SigIm, const double &qPr, const SigmaI &SigIp, const gsl_rng *r);
	// non-zero mean prior methods
	void update(const Grp &dat, const SigmaI &SigIm, const Grp &muPr, const SigmaI &SigIp, const gsl_rng *r);
	void update(const Grp &dat, const SigmaI &SigIm, const Grp &muPr, const double &qPr, const SigmaI &SigIp, const gsl_rng *r);
	void update(const Grp &dat, const Qgrp &q, const SigmaI &SigIm, const Grp &muPr, const SigmaI &SigIp, const gsl_rng *r);
	void update(const Grp &dat, const Qgrp &q, const SigmaI &SigIm, const Grp &muPr, const double &qPr, const SigmaI &SigIp, const gsl_rng *r);
};

class MVnormBetaBlk : public MVnorm {
protected:
	const vector< size_t > *_blkStart;  // starts of vector pieces; pointer to the vector that's in the corresponding Grp class
	vector< gsl_vector_view > _eachVec; // views of vector pieces, each corresponding to a block of variables
	vector< gsl_vector_view > _X;       // vector view of predictor values (typically column of a predictor matrix that is in the encapsulating group class).  The matrix this points to can be modified through this, as when I scale the predictor
	vector<double> _scale;              // Scale -- typically XtX, but not always (special cases dealt with in derived classes)
	
	const size_t *_upLevel;             // index into the vector of priors

	
public:
	MVnormBetaBlk() : MVnorm () {};
	MVnormBetaBlk(const Grp &resp, gsl_matrix *pred, const vector<size_t> &iCl, const gsl_matrix *Sig, const vector<size_t> &blkStart, const gsl_rng *r, gsl_matrix *bet, const size_t &iRw);
	MVnormBetaBlk(const Grp &resp, gsl_matrix *pred, const vector<size_t> &iCl, const gsl_matrix *Sig, const vector<size_t> &blkStart, const gsl_rng *r, const size_t &up, gsl_matrix *bet, const size_t &iRw);
	MVnormBetaBlk(const gsl_matrix *resp, gsl_matrix *pred, const vector<size_t> &iCl, const gsl_matrix *Sig, const vector<size_t> &blkStart, const gsl_rng *r, gsl_matrix *bet, const size_t &iRw);
	MVnormBetaBlk(const gsl_matrix *resp, gsl_matrix *pred, const vector<size_t> &iCl, const gsl_matrix *Sig, const vector<size_t> &blkStart, const gsl_rng *r, const size_t &up, gsl_matrix *bet, const size_t &iRw);
	
	virtual ~MVnormBetaBlk() {};
	
	const size_t *up() const{return _upLevel; };
	
	virtual void update(const Grp &dat, const SigmaI &SigIb, const gsl_rng *r);
	virtual void update(const Grp &dat, const Qgrp &q, const SigmaI &SigIb, const gsl_rng *r);
	// 0-mean prior methods
	virtual void update(const Grp &dat, const SigmaI &SigIb, const SigmaI &SigIp, const gsl_rng *r);
	virtual void update(const Grp &dat, const SigmaI &SigIb, const double &qPr, const SigmaI &SigIp, const gsl_rng *r);
	virtual void update(const Grp &dat, const Qgrp &q, const SigmaI &SigIb, const SigmaI &SigIp, const gsl_rng *r);
	virtual void update(const Grp &dat, const Qgrp &q, const SigmaI &SigIb, const double &qPr, const SigmaI &SigIp, const gsl_rng *r);
	// non-zero mean prior methods
	virtual void update(const Grp &dat, const SigmaI &SigIb, const Grp &muPr, const SigmaI &SigIp, const gsl_rng *r);
	virtual void update(const Grp &dat, const SigmaI &SigIb, const Grp &muPr, const double &qPr, const SigmaI &SigIp, const gsl_rng *r);
	virtual void update(const Grp &dat, const Qgrp &q, const SigmaI &SigIb, const Grp &muPr, const SigmaI &SigIp, const gsl_rng *r);
	virtual void update(const Grp &dat, const Qgrp &q, const SigmaI &SigIb, const Grp &muPr, const double &qPr, const SigmaI &SigIp, const gsl_rng *r);
};

class MVnormBetaFtBlk : public MVnormBetaBlk {
protected:
	gsl_matrix_view _fitted;
	
public:
	MVnormBetaFtBlk() : MVnormBetaBlk() {};
	MVnormBetaFtBlk(const Grp &resp, gsl_matrix *pred, vector<double> &eaFt, const vector<size_t> &iCl, const gsl_matrix *Sig, const vector<size_t> &blkStart, const gsl_rng *r, gsl_matrix *bet, const size_t &iRw) : MVnormBetaBlk(resp, pred, iCl, Sig, blkStart, r, bet, iRw) {_fitted = gsl_matrix_view_array(eaFt.data(), pred->size1, _d);};
	MVnormBetaFtBlk(const Grp &resp, gsl_matrix *pred, vector<double> &eaFt, const vector<size_t> &iCl, const gsl_matrix *Sig, const vector<size_t> &blkStart, const gsl_rng *r, const size_t &up, gsl_matrix *bet, const size_t &iRw) : MVnormBetaBlk(resp, pred, iCl, Sig, blkStart, r, up, bet, iRw) {_fitted = gsl_matrix_view_array(eaFt.data(), pred->size1, _d);};
	MVnormBetaFtBlk(const gsl_matrix *resp, gsl_matrix *pred, vector<double> &eaFt, const vector<size_t> &iCl, const gsl_matrix *Sig, const vector<size_t> &blkStart, const gsl_rng *r, gsl_matrix *bet, const size_t &iRw) : MVnormBetaBlk(resp, pred, iCl, Sig, blkStart, r, bet, iRw) {_fitted = gsl_matrix_view_array(eaFt.data(), pred->size1, _d);};
	MVnormBetaFtBlk(const gsl_matrix *resp, gsl_matrix *pred, vector<double> &eaFt, const vector<size_t> &iCl, const gsl_matrix *Sig, const vector<size_t> &blkStart, const gsl_rng *r, const size_t &up, gsl_matrix *bet, const size_t &iRw) : MVnormBetaBlk(resp, pred, iCl, Sig, blkStart, r, up, bet, iRw) {_fitted = gsl_matrix_view_array(eaFt.data(), pred->size1, _d);};
	
	~MVnormBetaFtBlk() {};
	
	void update(const Grp &dat, const SigmaI &SigIb, const gsl_rng *r);
	void update(const Grp &dat, const Qgrp &q, const SigmaI &SigIb, const gsl_rng *r);
	// 0-mean prior methods
	void update(const Grp &dat, const SigmaI &SigIb, const SigmaI &SigIp, const gsl_rng *r);
	void update(const Grp &dat, const SigmaI &SigIb, const double &qPr, const SigmaI &SigIp, const gsl_rng *r);
	void update(const Grp &dat, const Qgrp &q, const SigmaI &SigIb, const SigmaI &SigIp, const gsl_rng *r);
	void update(const Grp &dat, const Qgrp &q, const SigmaI &SigIb, const double &qPr, const SigmaI &SigIp, const gsl_rng *r);
	// non-zero mean prior methods
	void update(const Grp &dat, const SigmaI &SigIb, const Grp &muPr, const SigmaI &SigIp, const gsl_rng *r);
	void update(const Grp &dat, const SigmaI &SigIb, const Grp &muPr, const double &qPr, const SigmaI &SigIp, const gsl_rng *r);
	void update(const Grp &dat, const Qgrp &q, const SigmaI &SigIb, const Grp &muPr, const SigmaI &SigIp, const gsl_rng *r);
	void update(const Grp &dat, const Qgrp &q, const SigmaI &SigIb, const Grp &muPr, const double &qPr, const SigmaI &SigIp, const gsl_rng *r);

};
/** @} */
// end of lineLoc group

/** \brief Random index class that allows for mixture models.
 *	Can be used as a deterministic index if there is deterministic initiation and no updating
 */
class RanIndex {
protected:
	vector< vector<size_t> > _idx; // rugged 2D array that relates the upper-level elements to the indexes of the lower levels
	vector<size_t> _vecInd;        // vector that has the upper-level element IDs for each lower-level instance
	gsl_rng *_r;
	
public:
	// constructors
	RanIndex();
	RanIndex(const gsl_vector_int *lInd, const size_t &Ntot, const size_t &Nup);
	RanIndex(const size_t &Ntot);
	RanIndex(const size_t &Ntot, const size_t &Nup); // one-to-one for all the Ntot that are =< Nup
	RanIndex(const size_t &Ntot, const size_t &Nup, const string &fileNam); // reading it directly from a file
	RanIndex(const size_t &Ntot, const size_t &Nup, FILE *fileStr);
	
	virtual ~RanIndex() {gsl_rng_free(_r);};
	
	const vector<size_t> &operator[](const size_t i) const{return _idx[i];};
	const vector<size_t> &getIndVec() const{return _vecInd; };
	const size_t &priorInd(const size_t i) const{return _vecInd[i]; };
	vector<size_t> &operator[](const size_t i) {return _idx[i];};
	vector<size_t> &getIndVec() {return _vecInd; };
	size_t &priorInd(const size_t i) {return _vecInd[i]; };
	
	size_t getNtot(){return _vecInd.size(); };
	size_t getNgrp(){return _idx.size(); };
	size_t getNtot() const{return _vecInd.size(); };
	size_t getNgrp() const{return _idx.size(); };
	virtual const vector<double> props() const;
	
	virtual void init(const vector< vector<size_t> > &idx, const vector< vector<size_t> > &rLD); // determnistic initiation; rLD is ignored in this class
	virtual void save(const string &outFlNam) {};
	virtual void save(const Grp &y, const BetaGrpBVSR *theta, const SigmaI &SigIe) {};
	virtual void save(const Grp &y, const Grp &theta, const SigmaI &SigIe) {};
	virtual void dump(){};
	
	void update(const Grp &theta, const Grp &mu, const vector<SigmaI> &SigI, const MixP &p);
	void update(const Grp &theta, const Grp &mu, const SigmaI &SigI, const MixP &p);
	virtual void update(const Grp &y, const SigmaI &SigIe, BetaGrpBVSR *theta, const SigmaI &SigIp) {}; // not supported at present
};

/*
 *	Index class that keeps track of the variable selection process.
 */
class RanIndexVS : public RanIndex {  //  is a friend of BetaGrpBVSR
	
protected:
	vector<double> _pep;              // posterior EXCLUSION probabilities, in the order of the picked X; keep adding to this at each save, then dump into a file in the end
	vector<bool> _mcmcTrack;          // accept/reject and which move tracking.  Each _mcmcTrack[i] has 4 elements: accepted? adding proposal? dropping proposal? swapping proposal?
	vector< vector<size_t> > _relLD;  // relates the new positions of selected predictors to the original positions; _relLD[i][0] is the original position, i.e. _rel has to be as long as the selected # of predictors; the rest of _relLD (if any) are positions of linked SNPs (or, in general, correlated predictors)
	size_t _Ntp;			          // original number of predictors
	
	MixP *_prior;               // the prior
	
	vector<bool> _acceptDrop;	// these vectors store the possible values of _mcmcTrak, so I can just copy over the correct one instead of doing element-wise assignment
	vector<bool> _rejectDrop;
	vector<bool> _acceptAdd;
	vector<bool> _rejectAdd;
	vector<bool> _acceptSwap;
	vector<bool> _rejectSwap;
	
	double _numSaves;     // tracking the number of saves; make it double so that no static_cast is necessary for the division at the end
	string _pepOutFlnam;  // out file for pep
	string _mcmcPutFlNam; // file name to save the accept/reject statistics
	
	size_t _proposal(const size_t &Ntot, const int &Nmn, const gsl_rng *r); // rank proposal-generating function for Metropolis updating; not clear if it is symmetric so don't use for now

public:
	RanIndexVS();
	RanIndexVS(const size_t &Ntot, const string &outFlNam, MixP &pr); // preliminary set-up; initiation to be finished with an init() function
	RanIndexVS(const size_t &Ntot, const string &outFlNam);
	
	~RanIndexVS(){};
	
	void props(vector<double> &prp) const;
	
	void save(const Grp &y, const BetaGrpBVSR *theta, const SigmaI &SigIe); // not saving to a file (mostly)
	void save(const Grp &y, const Grp &theta, const SigmaI &SigIe);
	void dump();
	void init(const vector< vector<size_t> > &idx, const vector< vector<size_t> > &rLD); // idx is _idx, rLD is _relLD
	
	void update(const Grp &y, const SigmaI &SigIe, BetaGrpBVSR *theta, const SigmaI &SigIp);
};

/*
 *	Parameter-expansion variable for 0-mean MVnormMu
 */

class Apex {
private:
	gsl_matrix *_Amat;
	gsl_matrix *_SigIpr;
	vector<vector<double> > *_fitEach;
	gsl_rng *_r;
	
public:
	Apex();
	Apex(const size_t &d, vector<vector <double> > &fE);
	Apex(const double &dV, const size_t &d, vector<vector <double> > &fE);
	Apex(const SigmaI &SigI, vector<vector <double> > &fE);
	
	~Apex();
	
	Apex(const Apex &A);
	Apex &operator=(const Apex &A);
	
	gsl_matrix *getMat(){return _Amat; };
	const gsl_matrix *getMat() const{return _Amat; };
	void setPr(const double &pr);
	
	void update(const Grp &y, const gsl_matrix *xi, const SigmaI &SigIm);
	void update(const Grp &y, const gsl_matrix *xi, const SigmaI &SigIm, const RanIndex &ind);
};

MuGrp operator+(const Grp &m1, const Grp &m2);
MuGrp operator-(const Grp &m1, const Grp &m2);

class Grp {
	friend MuGrp operator+(const Grp &m1, const Grp &m2);
	friend MuGrp operator-(const Grp &m1, const Grp &m2);
	friend class MuGrp;
	
protected:
	vector<MVnorm *> _theta;
	gsl_matrix *_valueMat; // the individual _theta members point to rows of this and modify it on initialization and during updates
	RanIndex *_lowLevel;   // points to the level below; if unset will be 0
	RanIndex *_upLevel;    // points to the level above; if unset will be 0
	vector<gsl_rng *> _rV; // vector of RNGs; mostly will be length 1, but occasionally length == num_thread()
	
	string _outFlNam;
	
	Grp();
	
public:
	virtual ~Grp();
	
	// improper prior
	virtual void update(const Grp &, const SigmaI &) = 0;
	virtual void update(const Grp &, const Qgrp &, const SigmaI &) = 0;
	// 0-mean prior
	virtual void update(const Grp &, const SigmaI &, const SigmaI &) = 0;
	virtual void update(const Grp &, const Qgrp &, const SigmaI &, const SigmaI &) = 0;
	virtual void update(const Grp &, const SigmaI &, const Qgrp &, const SigmaI &) = 0;
	virtual void update(const Grp &, const Qgrp &, const SigmaI &, const Qgrp &, const SigmaI &) = 0;
	// non-0-mean prior
	virtual void update(const Grp &, const SigmaI &, const Grp &, const SigmaI &) = 0;
	virtual void update(const Grp &, const Qgrp &, const SigmaI &, const Grp &, const SigmaI &) = 0;
	virtual void update(const Grp &, const SigmaI &, const Grp &, const Qgrp &, const SigmaI &) = 0;
	virtual void update(const Grp &, const Qgrp &, const SigmaI &, const Grp &, const Qgrp &, const SigmaI &) = 0;
	
	virtual void save();
	virtual void save(const string &outFlNam);
	virtual void save(const string &outMuFlNam, const string &outSigFlNam, const SigmaI &SigI);
	virtual void save(const SigmaI &SigI) {};                // these two are for situations where we want to dump to a file only at the end of the run
	virtual void save(const Grp &y, const SigmaI &SigI) {};  // they have an effect only in the child classes that implement them
	void mhlSave(const string &outFlNam, const SigmaI SigI);
	
	virtual void dump(){}; // dump to a file in cases where save() saves to a matrix
	
	const vector<MVnorm *> &dataVec() const{return _theta; };
	virtual const gsl_matrix *dMat() const{return _valueMat; };
	virtual const gsl_matrix *fMat() const{return _valueMat; };
	size_t Ndata() const{return _theta.size(); };
	size_t phenD() const{return _valueMat->size2; };
	
	virtual double lnOddsRat(const Grp &y, const SigmaI &SigI, const size_t i) const {return 0.0; }; // makes sense only for regression classes
	
	const MVnorm * operator[](const size_t i) const{return _theta[i]; };
	MVnorm * operator[](const size_t i) {return _theta[i]; };
	
	virtual MuGrp mean(RanIndex &grp);
	virtual const MuGrp mean(RanIndex &grp) const;
	
	virtual MuGrp mean(RanIndex &grp, const Qgrp &q);
	virtual const MuGrp mean(RanIndex &grp, const Qgrp &q) const;
	
	void center(){ colCenter(_valueMat); };
	
};

class MuGrp : public Grp {
	friend MuGrp operator+(const Grp &m1, const Grp &m2);
	friend MuGrp operator-(const Grp &m1, const Grp &m2);
protected:
	
public:
	MuGrp() : Grp(){};
	MuGrp(RanIndex &low, const size_t &d);  // 0-mean prior constructor
	MuGrp(const string &datFlNam, RanIndex &low, RanIndex &up, const size_t &d); // RanIndex's are not const because I have a pointer to it as a member
	MuGrp(const string &datFlNam, RanIndex &up, const size_t d); // for setting up the lowest level deterministically
	MuGrp(const vector<MVnorm *> &dat, RanIndex &low, RanIndex &up);
	MuGrp(const Grp &dat, RanIndex &low, RanIndex &up);
	MuGrp(const vector<MVnorm *> &dat, RanIndex &low, RanIndex &up, const string &outFlNam);
	MuGrp(const Grp &dat, RanIndex &low, RanIndex &up, const string &outFlNam);
	MuGrp(const Grp &dat, RanIndex &low); // deterministic constructor to calculate means of dat according to groups in low
	MuGrp(const Grp &dat, const Qgrp &q, RanIndex &low); // deterministic constructor to calculate weighted means of dat according to groups in low
	MuGrp(const gsl_matrix *dat); // simple deterministic constructor
	MuGrp(const gsl_matrix *dat, RanIndex &low); // deterministic constructor to calculate means of dat according to groups in low
	MuGrp(const gsl_matrix *dat, const Qgrp &q, RanIndex &low); // deterministic constructor to calculate weighted means of dat according to groups in low
	
	virtual ~MuGrp() {};
	
	MuGrp(const MuGrp &mG); // copy constructor
	MuGrp(const Grp &g); // copy constructor
	MuGrp & operator=(const MuGrp &mG);
	
	// improper prior
	virtual void update(const Grp &dat, const SigmaI &SigIm);
	virtual void update(const Grp &dat, const Qgrp &q, const SigmaI &SigIm);
	// 0-mean prior
	virtual void update(const Grp &dat, const SigmaI &SigIm, const SigmaI &SigIp);
	virtual void update(const Grp &dat, const Qgrp &q, const SigmaI &SigIm, const SigmaI &SigIp);
	virtual void update(const Grp &dat, const SigmaI &SigIm, const Qgrp &qPr, const SigmaI &SigIp);
	virtual void update(const Grp &dat, const Qgrp &q, const SigmaI &SigIm, const Qgrp &qPr, const SigmaI &SigIp);
	// non-0-mean prior
	virtual void update(const Grp &dat, const SigmaI &SigIm, const Grp &muPr, const SigmaI &SigIp);
	virtual void update(const Grp &dat, const Qgrp &q, const SigmaI &SigIm, const Grp &muPr, const SigmaI &SigIp);
	virtual void update(const Grp &dat, const SigmaI &SigIm, const Grp &muPr, const Qgrp &qPr, const SigmaI &SigIp);
	virtual void update(const Grp &dat, const Qgrp &q, const SigmaI &SigIm, const Grp &muPr, const Qgrp &qPr, const SigmaI &SigIp);
	
};

class MuGrpPEX : public MuGrp {
protected:
	gsl_matrix *_adjValMat;
	gsl_matrix *_tSigIAt;   // t(SigIm%*%t(_A)) that is in common among all individual MVnormMuPEX's
	Apex _A;
	vector<vector<double> > _ftA;
	int _nThr;
	
	string _outSigFlNam; // file name to save the Sigma that goes with this Mu
	
	void _updateFitted();
	
public:
	MuGrpPEX();
	MuGrpPEX(const Grp &dat, RanIndex &low, RanIndex &up, const double &Spr, const int &nThr);
	MuGrpPEX(const Grp &dat, RanIndex &low, RanIndex &up, const string outMuFlNam, const double &Spr, const int &nThr);
	MuGrpPEX(const Grp &dat, RanIndex &low, RanIndex &up, const string outMuFlNam, const string outSigFlNam, const double &Spr, const int &nThr);
	
	~MuGrpPEX();
	
	MuGrpPEX(const MuGrpPEX &mGp);
	MuGrpPEX &operator=(const MuGrpPEX &mGp);
	
	const gsl_matrix *fMat() const{return _adjValMat; }; // _adjValMat returns the interesting values muPEX%*%A, but the dMat() still returns _valueMat that stores the "raw" (per Gelman and Hill) values
	void save(); // interested in saving the adjusted values, so the fMat()
	void save(const SigmaI &SigI); // saving the location parameters and the cognate Sigma
	Apex &getA(){return _A; };
	void setApr(const double &pr){_A.setPr(pr); };
	
	MuGrp mean(RanIndex &grp);
	const MuGrp mean(RanIndex &grp) const;
	MuGrp mean(RanIndex &grp, const Qgrp &q);
	const MuGrp mean(RanIndex &grp, const Qgrp &q) const;
	
	void update(const Grp &dat, const SigmaI &SigIm, const SigmaI &SigIp);
	void update(const Grp &dat, const SigmaI &SigIm, const Qgrp &qPr, const SigmaI &SigIp);
	
	void update(const Grp &dat, const SigmaI &SigIm, const Grp &muPr, const SigmaI &SigIp);
	void update(const Grp &dat, const SigmaI &SigIm, const Grp &muPr, const Qgrp &qPr, const SigmaI &SigIp);
};

class MuGrpMiss : public MuGrp {
protected:
	vector<size_t> _misInd;
	
public:
	MuGrpMiss() : MuGrp() {};
	MuGrpMiss(const string &datFlNam, const string &misMatFlNam, const string &misVecFlNam, RanIndex &up, const size_t &d);
	
	~MuGrpMiss() {};
	
	MuGrpMiss(const MuGrpMiss &mG); // copy constructor
	MuGrpMiss & operator=(const MuGrpMiss &mG);
	
	// only Gaussian updates work for imputation (in general)
	virtual void update(const Grp &mu, const SigmaI &SigIm);
	virtual void update(const Grp &mu, const SigmaI &SigIm, const SigmaI &SigIp);
	
	size_t nMis() const{return _misInd.size();};
	size_t nMis(){return _misInd.size();};
};

/** \brief Data with measurement error
 *
 */

class MuGrpEE : public MuGrp {
protected:
	gsl_matrix *_errorVar;
	vector<size_t> _errInd; // each row is the same
	
public:
	MuGrpEE() : MuGrp() {_errorVar = gsl_matrix_alloc(1, 1); };
	MuGrpEE(const string &datFlNam, const string &varFlNam, const string &indFlNam, RanIndex &up, const size_t &d);
	MuGrpEE(const string &datFlNam, const string &varFlNam, const vector<size_t> &varInd, RanIndex &up, const size_t &d);
	
	~MuGrpEE() {gsl_matrix_free(_errorVar); };
	
	MuGrpEE(const MuGrpEE &mG); // copy constructor
	MuGrpEE & operator=(const MuGrpEE &mG);
	
	// all of these are actually priors
	virtual void update(const Grp &muPr, const SigmaI &SigIm);
	virtual void update(const Grp &muPr, const Qgrp &q, const SigmaI &SigIm);
	
};

/** \brief Data with measurement error and missing phenotypes
 *
 */

class MuGrpEEmiss : public MuGrpMiss {
protected:
	vector< list<double> > _errorVar;
	vector< list<size_t> > _missErrMat;
	
public:
	MuGrpEEmiss() : MuGrpMiss() {};
	MuGrpEEmiss(const string &datFlNam, const string &varFlNam, const string &indFlNam, const string &misMatFlNam, const string &misVecFlNam, RanIndex &up, const size_t &d);
	MuGrpEEmiss(const string &datFlNam, const string &varFlNam, const vector<size_t> &varInd, const string &misMatFlNam, const string &misVecFlNam, RanIndex &up, const size_t &d);
	
	~MuGrpEEmiss() { };
	
	MuGrpEEmiss(const MuGrpEEmiss &mG); // copy constructor
	MuGrpEEmiss & operator=(const MuGrpEEmiss &mG);
	
	// all of these are actually priors
	void update(const Grp &muPr, const SigmaI &SigIm);
	void update(const Grp &muPr, const Qgrp &q, const SigmaI &SigIm);
	
};

class BetaGrpFt : public Grp {
protected:
	vector< vector<double> > _fittedEach; // each member of the outer vector stores the element-specific fitted matrix as vector in the row-major format, to be accessed as a matrix_view of an array
	gsl_matrix *_fittedAll;
	gsl_matrix *_valueSum; // sum of all the saved estimates; allocated only if needed (under the condition that _numSaves != 0.0)
	gsl_matrix *_Xmat;
	int _nThr;
	
	double _numSaves;      // number of saves made, to calculate the mean at the end; double b/c will need to divide a matrix of doubles by it
	
	virtual void _updateFitted();
	void _rankPred(const gsl_matrix *y, const SigmaI &SigI, gsl_vector *XtX, gsl_permutation *prm);  // ranking the predictors in preparation for throwing out the ones below a certain rank
	void _rankPred(const gsl_matrix *y, const SigmaI &SigI, const double &absLab, gsl_vector *XtX, gsl_permutation *prm);
	void _ldToss(const gsl_vector *var, const gsl_permutation *prm, const double &rSqMax, const size_t &Npck, vector< vector<size_t> > &idx, vector< vector<size_t> > &rLd, gsl_matrix *Xpck);  // tossing the predictors in LD with top hits, noting their identity and topping off the selection with "unlinked" SNPs
	virtual double _MGkernel(const Grp &dat, const SigmaI &SigI) const;   // MV Gaussian kernel calculation for all Xbeta
	virtual double _MGkernel(const Grp &dat, const SigmaI &SigI, const size_t &prInd) const;  // MV Gaussian kernel, dropping predictor prInd
public:
	BetaGrpFt();
	BetaGrpFt(const Grp &rsp, const string &predFlNam, const size_t &Npred, const int &nThr);
	BetaGrpFt(const Grp &rsp, const string &predFlNam, const size_t &Npred, RanIndex &up, const int &nThr);
	BetaGrpFt(const Grp &rsp, const string &predFlNam, const size_t &Npred, RanIndex &low, RanIndex &up, const int &nThr);
	BetaGrpFt(const Grp &rsp, const string &predFlNam, const size_t &Npred, const string &outFlNam, const int &nThr);
	BetaGrpFt(const Grp &rsp, const string &predFlNam, const size_t &Npred, RanIndex &up, const string &outFlNam, const int &nThr);
	BetaGrpFt(const Grp &rsp, const string &predFlNam, const size_t &Npred, RanIndex &low, RanIndex &up, const string &outFlNam, const int &nThr);
	
	// with missing predictor values (labeled by absLab)
	BetaGrpFt(const Grp &rsp, const string &predFlNam, const size_t &Npred, const double &absLab, const int &nThr);
	BetaGrpFt(const Grp &rsp, const string &predFlNam, const size_t &Npred, const double &absLab, RanIndex &up, const int &nThr);
	BetaGrpFt(const Grp &rsp, const string &predFlNam, const size_t &Npred, const double &absLab, RanIndex &low, RanIndex &up, const int &nThr);
	BetaGrpFt(const Grp &rsp, const string &predFlNam, const size_t &Npred, const double &absLab, const string &outFlNam, const int &nThr);
	BetaGrpFt(const Grp &rsp, const string &predFlNam, const size_t &Npred, const double &absLab, RanIndex &up, const string &outFlNam, const int &nThr);
	BetaGrpFt(const Grp &rsp, const string &predFlNam, const size_t &Npred, const double &absLab, RanIndex &low, RanIndex &up, const string &outFlNam, const int &nThr);
	
	// pre-screening of predictors based on initial rank.  Unlike BVSR, there is no update of the selected set
	BetaGrpFt(const Grp &rsp, const SigmaI &SigI, const string &predFlNam, const size_t &Npred, const double &Nmul, const double &rSqMax, RanIndex &up, const string &outFlNam, const int &nThr);
	BetaGrpFt(const Grp &rsp, const SigmaI &SigI, const string &predFlNam, const size_t &Npred, const double &Nmul, const double &rSqMax, RanIndex &low, RanIndex &up, const string &outFlNam, const int &nThr);
	BetaGrpFt(const Grp &rsp, const SigmaI &SigI, const string &predFlNam, const size_t &Npred, const double &Nmul, const double &rSqMax, const double &absLab, RanIndex &up, const string &outFlNam, const int &nThr);
	BetaGrpFt(const Grp &rsp, const SigmaI &SigI, const string &predFlNam, const size_t &Npred, const double &Nmul, const double &rSqMax, const double &absLab, RanIndex &low, RanIndex &up, const string &outFlNam, const int &nThr);
	
	virtual ~BetaGrpFt();
	
	BetaGrpFt(const BetaGrpFt &mG); // copy constructor
	BetaGrpFt & operator=(const BetaGrpFt &mG);
	
	virtual const gsl_matrix *fMat() const{return _fittedAll; };
	void save(const SigmaI &SigI);
	void dump();
	
	double lnOddsRat(const Grp &y, const SigmaI &SigI, const size_t i) const;
	
	// improper prior
	void update(const Grp &dat, const SigmaI &SigIm);
	void update(const Grp &dat, const Qgrp &q, const SigmaI &SigIm);
	// 0-mean prior
	virtual void update(const Grp &dat, const SigmaI &SigIm, const SigmaI &SigIp);
	void update(const Grp &dat, const Qgrp &q, const SigmaI &SigIm, const SigmaI &SigIp);
	virtual void update(const Grp &dat, const SigmaI &SigIm, const Qgrp &qPr, const SigmaI &SigIp);
	void update(const Grp &dat, const Qgrp &q, const SigmaI &SigIm, const Qgrp &qPr, const SigmaI &SigIp);
	// non-0-mean prior
	virtual void update(const Grp &dat, const SigmaI &SigIm, const Grp &muPr, const SigmaI &SigIp);
	void update(const Grp &dat, const Qgrp &q, const SigmaI &SigIm, const Grp &muPr, const SigmaI &SigIp);
	virtual void update(const Grp &dat, const SigmaI &SigIm, const Grp &muPr, const Qgrp &qPr, const SigmaI &SigIp);
	void update(const Grp &dat, const Qgrp &q, const SigmaI &SigIm, const Grp &muPr, const Qgrp &qPr, const SigmaI &SigIp);
	
};

class BetaGrpPEX : virtual public BetaGrpFt {
protected:
	gsl_matrix *_tSigIAt;   // t(SigIm%*%t(_A)) that is in common among all individual MVnormMuPEX's
	Apex _A;
	vector<vector<double> > _ftA;
	gsl_matrix *_fittedAllAdj;  // this is the X%*%Zeta%*%A, the regular scale matrix.  The _fittedAll is the "raw" matrix, as is the _valueMat (we don't need the adjusted _valueMat)
	
	void _finishConstruct(const double &Spr);
	void _finishFitted();
	void _updateAfitted();
	
	BetaGrpPEX(const double &Spr) {_finishConstruct(Spr); }; // only used by BetaGrpPCpex
public:
	BetaGrpPEX() : BetaGrpFt() {};
	BetaGrpPEX(const Grp &rsp, const string &predFlNam, const size_t &Npred, const double &Spr, const int &nThr) : BetaGrpFt(rsp, predFlNam, Npred, nThr) {_finishConstruct(Spr); };
	BetaGrpPEX(const Grp &rsp, const string &predFlNam, const size_t &Npred, const double &Spr, RanIndex &up, const int &nThr) : BetaGrpFt(rsp, predFlNam, Npred, up, nThr) {_finishConstruct(Spr); };
	BetaGrpPEX(const Grp &rsp, const string &predFlNam, const size_t &Npred, const double &Spr, RanIndex &low, RanIndex &up, const int &nThr) : BetaGrpFt(rsp, predFlNam, Npred, low, up, nThr) {_finishConstruct(Spr); };
	BetaGrpPEX(const Grp &rsp, const string &predFlNam, const size_t &Npred, const double &Spr, const string &outFlNam, const int &nThr) : BetaGrpFt(rsp, predFlNam, Npred, outFlNam, nThr) {_finishConstruct(Spr); };
	BetaGrpPEX(const Grp &rsp, const string &predFlNam, const size_t &Npred, const double &Spr, RanIndex &up, const string &outFlNam, const int &nThr) : BetaGrpFt(rsp, predFlNam, Npred, up, outFlNam, nThr) {_finishConstruct(Spr); };
	BetaGrpPEX(const Grp &rsp, const string &predFlNam, const size_t &Npred, const double &Spr, RanIndex &low, RanIndex &up, const string &outFlNam, const int &nThr) : BetaGrpFt(rsp, predFlNam, Npred, low, up, outFlNam, nThr) {_finishConstruct(Spr); };
	
	// with missing predictor values (labeled by absLab)
	BetaGrpPEX(const Grp &rsp, const string &predFlNam, const size_t &Npred, const double &Spr, const double &absLab, const int &nThr) : BetaGrpFt(rsp, predFlNam, Npred, absLab, nThr) {_finishConstruct(Spr); };
	BetaGrpPEX(const Grp &rsp, const string &predFlNam, const size_t &Npred, const double &Spr, const double &absLab, RanIndex &up, const int &nThr) : BetaGrpFt(rsp, predFlNam, Npred, absLab, up, nThr) {_finishConstruct(Spr); };
	BetaGrpPEX(const Grp &rsp, const string &predFlNam, const size_t &Npred, const double &Spr, const double &absLab, RanIndex &low, RanIndex &up, const int &nThr) : BetaGrpFt(rsp, predFlNam, Npred, absLab, low, up, nThr) {_finishConstruct(Spr); };
	BetaGrpPEX(const Grp &rsp, const string &predFlNam, const size_t &Npred, const double &Spr, const double &absLab, const string &outFlNam, const int &nThr) : BetaGrpFt(rsp, predFlNam, Npred, absLab, outFlNam, nThr) {_finishConstruct(Spr); };
	BetaGrpPEX(const Grp &rsp, const string &predFlNam, const size_t &Npred, const double &Spr, const double &absLab, RanIndex &up, const string &outFlNam, const int &nThr) : BetaGrpFt(rsp, predFlNam, Npred, absLab, up, outFlNam, nThr) {_finishConstruct(Spr); };
	BetaGrpPEX(const Grp &rsp, const string &predFlNam, const size_t &Npred, const double &Spr, const double &absLab, RanIndex &low, RanIndex &up, const string &outFlNam, const int &nThr) : BetaGrpFt(rsp, predFlNam, Npred, absLab, low, up, outFlNam, nThr) {_finishConstruct(Spr); };
	
	// pre-screening of predictors based on initial rank.  Unlike BVSR, there is no update of the selected set
	BetaGrpPEX(const Grp &rsp, const SigmaI &SigI, const string &predFlNam, const size_t &Npred, const double &Spr, const double &Nmul, const double &rSqMax, RanIndex &up, const string &outFlNam, const int &nThr) : BetaGrpFt(rsp, SigI, predFlNam, Npred, Nmul, rSqMax, up, outFlNam, nThr) {_finishConstruct(Spr); };
	BetaGrpPEX(const Grp &rsp, const SigmaI &SigI, const string &predFlNam, const size_t &Npred, const double &Spr, const double &Nmul, const double &rSqMax, RanIndex &low, RanIndex &up, const string &outFlNam, const int &nThr) : BetaGrpFt(rsp, SigI, predFlNam, Npred, Nmul, rSqMax, low, up, outFlNam, nThr) {_finishConstruct(Spr); };
	BetaGrpPEX(const Grp &rsp, const SigmaI &SigI, const string &predFlNam, const size_t &Npred, const double &Spr, const double &Nmul, const double &rSqMax, const double &absLab, RanIndex &up, const string &outFlNam, const int &nThr) : BetaGrpFt(rsp, SigI, predFlNam, Npred, Nmul, rSqMax, absLab, up, outFlNam, nThr) {_finishConstruct(Spr); };
	BetaGrpPEX(const Grp &rsp, const SigmaI &SigI, const string &predFlNam, const size_t &Npred, const double &Spr, const double &Nmul, const double &rSqMax, const double &absLab, RanIndex &low, RanIndex &up, const string &outFlNam, const int &nThr) : BetaGrpFt(rsp, SigI, predFlNam, Npred, Nmul, rSqMax, absLab, low, up, outFlNam, nThr) {_finishConstruct(Spr); };
	
	virtual ~BetaGrpPEX();
	
	virtual const gsl_matrix *fMat() const{return _fittedAllAdj; };
	
	void save();
	void save(const string &outFlNam);
	
	virtual void update(const Grp &dat, const SigmaI &SigIm, const SigmaI &SigIp);
	virtual void update(const Grp &dat, const SigmaI &SigIm, const Qgrp &qPr, const SigmaI &SigIp);
	
	virtual void update(const Grp &dat, const SigmaI &SigIm, const Grp &muPr, const SigmaI &SigIp);
	virtual void update(const Grp &dat, const SigmaI &SigIm, const Grp &muPr, const Qgrp &qPr, const SigmaI &SigIp);
};


class BetaGrpPC : virtual public BetaGrpFt {
protected:
	
public:
	BetaGrpPC() : BetaGrpFt() {};
	BetaGrpPC(const Grp &rsp, const string &predFlNam, const string &evFlNam, const size_t &Npred, RanIndex &up, const int &nThr);
	BetaGrpPC(const Grp &rsp, const string &predFlNam, const string &evFlNam, const size_t &Npred, RanIndex &low, RanIndex &up, const int &nThr);
	BetaGrpPC(const Grp &rsp, const string &predFlNam, const string &evFlNam, const size_t &Npred, RanIndex &up, const string &outFlNam, const int &nThr);
	BetaGrpPC(const Grp &rsp, const string &predFlNam, const string &evFlNam, const size_t &Npred, RanIndex &low, RanIndex &up, const string &outFlNam, const int &nThr);
	
	virtual ~BetaGrpPC() {};
	
	BetaGrpPC(const BetaGrpPC &mG); // copy constructor
	BetaGrpPC & operator=(const BetaGrpPC &mG);
	
	//update methods taken care of by BetaGrpFt
};

class BetaGrpPCpex : public BetaGrpPC, public BetaGrpPEX {
protected:
	
public:
	BetaGrpPCpex(){};
	BetaGrpPCpex(const Grp &rsp, const string &predFlNam, const string &evFlNam, const size_t &Npred, const double &Spr, RanIndex &up, const int &nThr) : BetaGrpFt(), BetaGrpPC(rsp, predFlNam, evFlNam, Npred, up, nThr), BetaGrpPEX(Spr) {};
	BetaGrpPCpex(const Grp &rsp, const string &predFlNam, const string &evFlNam, const size_t &Npred, const double &Spr, RanIndex &low, RanIndex &up, const int &nThr) : BetaGrpFt(), BetaGrpPC(rsp, predFlNam, evFlNam, Npred, low, up, nThr), BetaGrpPEX(Spr) {};
	BetaGrpPCpex(const Grp &rsp, const string &predFlNam, const string &evFlNam, const size_t &Npred, const double &Spr, RanIndex &up, const string &outFlNam, const int &nThr) : BetaGrpFt(), BetaGrpPC(rsp, predFlNam, evFlNam, Npred, up, outFlNam, nThr), BetaGrpPEX(Spr) {};
	BetaGrpPCpex(const Grp &rsp, const string &predFlNam, const string &evFlNam, const size_t &Npred, const double &Spr, RanIndex &low, RanIndex &up, const string &outFlNam, const int &nThr) : BetaGrpFt(), BetaGrpPC(rsp, predFlNam, evFlNam, Npred, low, up, outFlNam, nThr), BetaGrpPEX(Spr) {};
	
	~BetaGrpPCpex() {};
	
	const gsl_matrix *fMat() const{return _fittedAllAdj; };
	
	void update(const Grp &dat, const SigmaI &SigIm, const SigmaI &SigIp) { BetaGrpPEX::update(dat, SigIm, SigIp); };
	void update(const Grp &dat, const SigmaI &SigIm, const Qgrp &qPr, const SigmaI &SigIp) { BetaGrpPEX::update(dat, SigIm, qPr, SigIp); };
	
	void update(const Grp &dat, const SigmaI &SigIm, const Grp &muPr, const SigmaI &SigIp) { BetaGrpPEX::update(dat, SigIm, muPr, SigIp); };
	void update(const Grp &dat, const SigmaI &SigIm, const Grp &muPr, const Qgrp &qPr, const SigmaI &SigIp) { BetaGrpPEX::update(dat, SigIm, muPr, qPr, SigIp); };
};

class BetaGrpSnp : public MuGrp {
protected:
	gsl_matrix *_Xmat;
	gsl_matrix *_fakeFmat; // is all 0, in case someone uses this class in an addition/subtraction operator
	
	gsl_matrix *_Ystore;   // storing the values for predictor until the end, when the mean is used for the SNP regression
	
	size_t _Npred;
	int _nThr;
	double _numSaves;
	double _priorVar;  // W in Wakefield's (2007) ABF paper
	
	string _inPredFl;
	
public:
	BetaGrpSnp();
	BetaGrpSnp(const string &predFlNam, const string &outFlNam, const size_t &Ndat, const size_t &Npred, const size_t &d, const int &Nthr);
	BetaGrpSnp(const string &predFlNam, const string &outFlNam, RanIndex &low, const size_t &Npred, const size_t &d, const int &Nthr);
	BetaGrpSnp(const string &predFlNam, const string &outFlNam, const size_t &Ndat, const size_t &Npred, const size_t &d, const int &Nthr, const double &prVar);
	BetaGrpSnp(const string &predFlNam, const string &outFlNam, RanIndex &low, const size_t &Npred, const size_t &d, const int &Nthr, const double &prVar);
	
	~BetaGrpSnp();
	
	virtual void dump();
	
	BetaGrpSnp(const BetaGrpSnp &mG); // copy constructor
	BetaGrpSnp & operator=(const BetaGrpSnp &mG);
	
	const gsl_matrix *fMat() const{return _fakeFmat; };

	void update(const Grp &dat, const SigmaI &SigIm);
};

/** \brief Regression with conditional variance
 *
 */
class BetaGrpSnpCV : public BetaGrpSnp {
protected:
	
public:
	BetaGrpSnpCV() : BetaGrpSnp() {};
	BetaGrpSnpCV(const string &predFlNam, const string &outFlNam, const size_t &Ndat, const size_t &Npred, const size_t &d, const int &Nthr) : BetaGrpSnp(predFlNam, outFlNam, Ndat, Npred, d, Nthr) {};
	BetaGrpSnpCV(const string &predFlNam, const string &outFlNam, RanIndex &low, const size_t &Npred, const size_t &d, const int &Nthr) : BetaGrpSnp(predFlNam, outFlNam, low, Npred, d, Nthr) {};
	BetaGrpSnpCV(const string &predFlNam, const string &outFlNam, const size_t &Ndat, const size_t &Npred, const size_t &d, const int &Nthr, const double &prVar) : BetaGrpSnp(predFlNam, outFlNam, Ndat, Npred, d, Nthr, prVar) {};
	BetaGrpSnpCV(const string &predFlNam, const string &outFlNam, RanIndex &low, const size_t &Npred, const size_t &d, const int &Nthr, const double &prVar) : BetaGrpSnp(predFlNam, outFlNam, low, Npred, d, Nthr, prVar) {};
	
	~BetaGrpSnpCV() {};
	
	BetaGrpSnpCV(const BetaGrpSnpCV &mG); // copy constructor
	BetaGrpSnpCV & operator=(const BetaGrpSnpCV &mG);
	
	void dump();
	
};

/** \brief regressions on a single SNP controlling for other traits
 *
 */
class BetaGrpPSR : public BetaGrpSnp {
protected:
	
public:
	BetaGrpPSR() : BetaGrpSnp() {};
	BetaGrpPSR(const string &predFlNam, const string &outFlNam, const size_t &Ndat, const size_t &Npred, const size_t &d, const int &Nthr) : BetaGrpSnp(predFlNam, outFlNam, Ndat, Npred, d, Nthr) {};
	BetaGrpPSR(const string &predFlNam, const string &outFlNam, RanIndex &low, const size_t &Npred, const size_t &d, const int &Nthr) : BetaGrpSnp(predFlNam, outFlNam, low, Npred, d, Nthr) {};
	BetaGrpPSR(const string &predFlNam, const string &outFlNam, const size_t &Ndat, const size_t &Npred, const size_t &d, const int &Nthr, const double &prVar) : BetaGrpSnp(predFlNam, outFlNam, Ndat, Npred, d, Nthr, prVar) {};
	BetaGrpPSR(const string &predFlNam, const string &outFlNam, RanIndex &low, const size_t &Npred, const size_t &d, const int &Nthr, const double &prVar) : BetaGrpSnp(predFlNam, outFlNam, low, Npred, d, Nthr, prVar) {};
	
	~BetaGrpPSR() {};
	
	void dump();
	
	BetaGrpPSR(const BetaGrpPSR &mG); // copy constructor
	BetaGrpPSR & operator=(const BetaGrpPSR &mG);
	
};

class BetaGrpSnpMiss : public MuGrp {
protected:
	vector<vector<double> > _Xmat; // a rugged array that only stores present genotypes
	gsl_matrix *_fakeFmat;
	
	gsl_matrix *_Ystore;   // storing the values for predictor until the end, when the mean is used for the SNP regression
	
	size_t _Npred;
	int _nThr;
	double _numSaves;
	double _absLab;
	double _priorVar;  // W in Wakefield's (2007) ABF paper
	
	string _inPredFl;
	
public:
	BetaGrpSnpMiss();
	BetaGrpSnpMiss(const string &predFlNam, const string &outFlNam, const size_t &Ndat, const size_t &Npred, const size_t &d, const int &Nthr, const double &absLab);
	BetaGrpSnpMiss(const string &predFlNam, const string &outFlNam, RanIndex &low, const size_t &Npred, const size_t &d, const int &Nthr, const double &absLab);
	BetaGrpSnpMiss(const string &predFlNam, const string &outFlNam, const size_t &Ndat, const size_t &Npred, const size_t &d, const int &Nthr, const double &prVar, const double &absLab);
	BetaGrpSnpMiss(const string &predFlNam, const string &outFlNam, RanIndex &low, const size_t &Npred, const size_t &d, const int &Nthr, const double &prVar, const double &absLab);
	
	~BetaGrpSnpMiss();
	
	BetaGrpSnpMiss(const BetaGrpSnpMiss &mG); // copy constructor
	BetaGrpSnpMiss & operator=(const BetaGrpSnpMiss &mG);
	
	virtual void dump();
	
	const gsl_matrix *fMat() const{return _fakeFmat; };
	
	void update(const Grp &dat, const SigmaI &SigIm);
	
};

/** \brief Missing-data regression with conditional variance
 *
 */
class BetaGrpSnpMissCV : public BetaGrpSnpMiss {
protected:
	
public:
	BetaGrpSnpMissCV() : BetaGrpSnpMiss() {};
	BetaGrpSnpMissCV(const string &predFlNam, const string &outFlNam, const size_t &Ndat, const size_t &Npred, const size_t &d, const int &Nthr, const double &absLab) : BetaGrpSnpMiss(predFlNam, outFlNam, Ndat, Npred, d, Nthr, absLab) {};
	BetaGrpSnpMissCV(const string &predFlNam, const string &outFlNam, RanIndex &low, const size_t &Npred, const size_t &d, const int &Nthr, const double &absLab) : BetaGrpSnpMiss(predFlNam, outFlNam, low, Npred, d, Nthr, absLab) {};
	BetaGrpSnpMissCV(const string &predFlNam, const string &outFlNam, const size_t &Ndat, const size_t &Npred, const size_t &d, const int &Nthr, const double &prVar, const double &absLab) : BetaGrpSnpMiss(predFlNam, outFlNam, Ndat, Npred, d, Nthr, prVar, absLab) {};
	BetaGrpSnpMissCV(const string &predFlNam, const string &outFlNam, RanIndex &low, const size_t &Npred, const size_t &d, const int &Nthr, const double &prVar, const double &absLab) : BetaGrpSnpMiss(predFlNam, outFlNam, low, Npred, d, Nthr, prVar, absLab) {};
	
	~BetaGrpSnpMissCV(){};
	
	BetaGrpSnpMissCV(const BetaGrpSnpMissCV &mG); // copy constructor
	BetaGrpSnpMissCV & operator=(const BetaGrpSnpMissCV &mG);
	
	void dump();
	
};

class BetaGrpPSRmiss : public BetaGrpSnpMiss {
protected:
	
public:
	BetaGrpPSRmiss() : BetaGrpSnpMiss() {};
	BetaGrpPSRmiss(const string &predFlNam, const string &outFlNam, const size_t &Ndat, const size_t &Npred, const size_t &d, const int &Nthr, const double &absLab) : BetaGrpSnpMiss(predFlNam, outFlNam, Ndat, Npred, d, Nthr, absLab) {};
	BetaGrpPSRmiss(const string &predFlNam, const string &outFlNam, RanIndex &low, const size_t &Npred, const size_t &d, const int &Nthr, const double &absLab) : BetaGrpSnpMiss(predFlNam, outFlNam, low, Npred, d, Nthr, absLab) {};
	BetaGrpPSRmiss(const string &predFlNam, const string &outFlNam, const size_t &Ndat, const size_t &Npred, const size_t &d, const int &Nthr, const double &prVar, const double &absLab) : BetaGrpSnpMiss(predFlNam, outFlNam, Ndat, Npred, d, Nthr, prVar, absLab) {};
	BetaGrpPSRmiss(const string &predFlNam, const string &outFlNam, RanIndex &low, const size_t &Npred, const size_t &d, const int &Nthr, const double &prVar, const double &absLab) : BetaGrpSnpMiss(predFlNam, outFlNam, low, Npred, d, Nthr, prVar, absLab) {};
	
	~BetaGrpPSRmiss() {};
	
	BetaGrpPSRmiss(const BetaGrpPSRmiss &mG); // copy constructor
	BetaGrpPSRmiss & operator=(const BetaGrpPSRmiss &mG);
	
	void dump();

};

class BetaGrpBVSR : public BetaGrpFt {
	friend void RanIndexVS::update(const Grp &, const SigmaI &, BetaGrpBVSR *, const SigmaI &);
protected:
	gsl_matrix *_selB;     // matrix of selected predictors
	gsl_matrix *_tmpXb;    // Xb matrix for the candidate dropped/added element
	
	// protected functions
	void _updateFitted();  // local version of the _updateFitted function used in all BetaGrpFt
	void _ldToss(const gsl_vector *var, const gsl_permutation *prm, const double &rSqMax, const size_t &Nsel, const size_t &Npck, vector< vector<size_t> > &idx, vector< vector<size_t> > &rLd, gsl_matrix *Xpck);  // tossing the predictors in LD with top hits, noting their identity and topping off the selection with "unlinked" SNPs
	double _MGkernel(const Grp &dat, const SigmaI &SigI) const;   // MV Gaussian kernel calculation for all Xbeta
	double _MGkernel(const Grp &dat, const SigmaI &SigI, const size_t &prInd) const;  // MV Gaussian kernel, dropping predictor prInd
	double _MGkernel(const Grp &dat, const SigmaI &SigIe, const SigmaI &SigIpr, const size_t &prInd);  // MV Gaussian kernel, adding predictor prInd with the prior covariance matrix SigIpr
public:
	BetaGrpBVSR() : BetaGrpFt() {};
	BetaGrpBVSR(const Grp &y, const SigmaI &SigI, const string &predFlNam, const double &Nmul, const double &rSqMax, RanIndex &up, const string &outFlNam, const int &nThr);
	BetaGrpBVSR(const Grp &y, const SigmaI &SigI, const string &predFlNam, const double &Nmul, const double &rSqMax, RanIndex &low, RanIndex &up, const string &outFlNam, const int &nThr);
	BetaGrpBVSR(const Grp &y, const SigmaI &SigI, const string &predFlNam, const double &Nmul, const double &rSqMax, const double &absLab, RanIndex &up, const string &outFlNam, const int &nThr); // absLab is the token that signifies missing data
	BetaGrpBVSR(const Grp &y, const SigmaI &SigI, const string &predFlNam, const double &Nmul, const double &rSqMax, const double &absLab, RanIndex &low, RanIndex &up, const string &outFlNam, const int &nThr);
	
	~BetaGrpBVSR();
	
	void save(const SigmaI &SigI);
	void save(const Grp &y, const SigmaI &SigI);
	void dump();
	const gsl_matrix *dMat() const{return _selB; };
	
	void update(const Grp &dat, const SigmaI &SigIm, const SigmaI &SigIp);
};

/** \brief Independent blocks of traits
 *
 */
class MuBlk : public MuGrp {
protected:
	vector< size_t > _blkStart;                  // vector of indexes into beginnings of each block
	vector< vector< vector<size_t> > > _blkLow;  // lowLevel indexes, separate for each block, not necessarily the same number of levels
	gsl_matrix *_expandedVM;                     // the _valuMat expanded according to _blkLow
	vector<size_t> _shortLevels;                 // for cases where the number of levels of low differ among blocks, store here index of the first element past the available one
	
	void _updateExp();
	void _fillIn();
	void _fillInUp();
	
public:
	MuBlk() : MuGrp() {};
	// Nval is the number of rows in _valueMat; if there are diffent numbers of levels of the low factor in each block, it's the maximum number
	MuBlk(const Grp &dat, const string &lowIndFlName, const size_t &Nval, RanIndex &up, const string &blkIndFileNam);
	MuBlk(const Grp &dat, const string &lowIndFlName, const size_t &Nval, RanIndex &up, const string &outFlNam, const string &blkIndFileNam);
	
	~MuBlk() {gsl_matrix_free(_expandedVM); };
	
	const gsl_matrix *fMat() const{ return _expandedVM; };
	
	// improper prior
	void update(const Grp &dat, const SigmaI &SigIm);
	void update(const Grp &dat, const Qgrp &q, const SigmaI &SigIm);
	// 0-mean prior
	void update(const Grp &dat, const SigmaI &SigIm, const SigmaI &SigIp);
	void update(const Grp &dat, const Qgrp &q, const SigmaI &SigIm, const SigmaI &SigIp);
	void update(const Grp &dat, const SigmaI &SigIm, const Qgrp &qPr, const SigmaI &SigIp);
	void update(const Grp &dat, const Qgrp &q, const SigmaI &SigIm, const Qgrp &qPr, const SigmaI &SigIp);
	// non-0-mean prior
	void update(const Grp &dat, const SigmaI &SigIm, const Grp &muPr, const SigmaI &SigIp);
	void update(const Grp &dat, const Qgrp &q, const SigmaI &SigIm, const Grp &muPr, const SigmaI &SigIp);
	void update(const Grp &dat, const SigmaI &SigIm, const Grp &muPr, const Qgrp &qPr, const SigmaI &SigIp);
	void update(const Grp &dat, const Qgrp &q, const SigmaI &SigIm, const Grp &muPr, const Qgrp &qPr, const SigmaI &SigIp);
};


class BetaBlk : public BetaGrpFt {

protected:
	vector< size_t > _blkStart;        // vector of indexes into beginnings of each block
	vector< gsl_matrix_view > _eachX;  // matrix-view submatrices of the predictors
	vector< gsl_matrix_view > _eachB;  // matrix-view submatrices of coefficients
	
	void _updateFitted();
public:
	BetaBlk() : BetaGrpFt() {};
	// total Ncol of the predictor matrix is Npred*Nblocks
	BetaBlk(const Grp &dat, const string &predFlName, const size_t &Npred, const string &blkIndFileNam, const int &nThr);
	BetaBlk(const Grp &dat, const string &predFlName, const size_t &Npred, const string &outFlNam, const string &blkIndFileNam, const int &nThr);
	BetaBlk(const Grp &dat, const string &predFlName, const size_t &Npred, RanIndex &up, const string &blkIndFileNam, const int &nThr);
	BetaBlk(const Grp &dat, const string &predFlName, const size_t &Npred, RanIndex &up, const string &outFlNam, const string &blkIndFileNam, const int &nThr);
	
	BetaBlk(const Grp &dat, const string &predFlName, const size_t &Npred, const double &absLab, const string &blkIndFileNam, const int &nThr);
	BetaBlk(const Grp &dat, const string &predFlName, const size_t &Npred, const double &absLab, const string &outFlNam, const string &blkIndFileNam, const int &nThr);
	BetaBlk(const Grp &dat, const string &predFlName, const size_t &Npred, const double &absLab, RanIndex &up, const string &blkIndFileNam, const int &nThr);
	BetaBlk(const Grp &dat, const string &predFlName, const size_t &Npred, const double &absLab, RanIndex &up, const string &outFlNam, const string &blkIndFileNam, const int &nThr);
	
	~BetaBlk() {};
	
	// updates handled by BetaGrpFt
};


class SigmaI {
protected:
	gsl_matrix *_mat;
	size_t _d;
	double _srDet;      // square root of the determinant; only calculate if necessary (e.g., for MV normal density calculation)
	
	gsl_matrix *_LamSc; // prior covariance matrix scaled by prior df
	double _n0;         // prior degrees of freedom
	
	string _outFlNam;   // file to save output
	
	gsl_rng *_r;
	
public:
	// overloaded constructors
	SigmaI(); // necessary default constructor
	SigmaI(const size_t &d, const double &invVar);  // constructor to give a vague high-variance prior to a mean
	SigmaI(const size_t &d, const double &invVar, const double &df);  // constructor to give a diagonal prior with a set df
	SigmaI(const size_t &d, const double &invVar, const double &df, const string &outFlNam);
	SigmaI(const gsl_matrix *mat);                  // deterministic constructor (use for the SNP regression, for example)
	// legacy C-based constructors
	SigmaI(const gsl_matrix *S, const size_t d, const size_t &df, const gsl_matrix *LamPr, const double &nu0);
	SigmaI(const gsl_matrix *S, const size_t d, const size_t &df, const double &diagPr, const double &nu0);
	
	// C++ abstract class-based constructor
	SigmaI(const Grp &dat, const double &prDiag, const double &nu0);
	SigmaI(const Grp &dat, const string &outFlNam, const double &prDiag, const double &nu0);
	
	SigmaI(const SigmaI &); // copy constructor
	SigmaI & operator=(const SigmaI &);
	
	virtual ~SigmaI(); // destructor
	
	// C++ abstract class-based updates
	// Gaussian
	virtual void update(const Grp &dat);
	virtual void update(const Grp &dat, const Grp &mu);
	// Student-t
	virtual void update(const Grp &dat, const Qgrp &q);
	virtual void update(const Grp &dat, const Grp &mu, const Qgrp &q);
	
	// save function for the stored file name
	virtual void save(const char *how = "a");
	// save function, taking file name, appending by default
	virtual void save(const string &fileNam, const char *how = "a");
	// PEX save
	void save(const string &fileNam, const Apex &A, const char *how = "a");
	// save function, taking file stream name
	void save(FILE *fileStr) {gsl_matrix_fwrite(fileStr, _mat); };
	
	const gsl_matrix *getMat() const{return _mat; };
	
	//dealing with the determinant
	
	void   srDetUpdate();
	double getSrDet() const {return _srDet; };
	
};

class SigmaIblk : public SigmaI {  // block-diagonal inverse-covariance matrix; tests show that it pays to do block-wise inversions
protected:
	vector< size_t > _blkStart;          // vector of indexes into beginnings of each block
	vector< gsl_matrix_view > _eachBlk;  // matrix-view submatrices of the whole thing
	vector< gsl_matrix_view > _eachLS;   // matrix-view submatrices of the prior
	
public:
	SigmaIblk();
	SigmaIblk(const size_t &d, const double &invVar, const double &df, const string &blkIndFileNam);
	SigmaIblk(const size_t &d, const double &invVar, const double &df, const string &blkIndFileNam, const string &outFlNam);
	SigmaIblk(const Grp &dat, const string &blkIndFileNam, const double &prDiag, const double &nu0);
	SigmaIblk(const Grp &dat, const string &blkIndFileNam, const string &outFlNam, const double &prDiag, const double &nu0);

	~SigmaIblk() {};
	
	// Gaussian
	void update(const Grp &dat);
	void update(const Grp &dat, const Grp &mu);
	
};

class SigmaIpex : public SigmaI {  // class for van Dyk and Meng's PEX for MV Student-t
protected:
	double _alpha;  // the PEX parameter
	
public:
	SigmaIpex() : _alpha(1.0), SigmaI() {};
	SigmaIpex(const size_t &d, const double &invVar, const double &df) : _alpha(1.0), SigmaI(d, invVar, df) {};
	SigmaIpex(const size_t &d, const double &invVar, const double &df, const string &outFlNam) : _alpha(1.0), SigmaI(d, invVar, df, outFlNam) {};
	SigmaIpex(const Grp &dat, const double &prDiag, const double &nu0) : _alpha(1.0), SigmaI(dat, prDiag, nu0) {};
	SigmaIpex(const Grp &dat, const string &outFlNam, const double &prDiag, const double &nu0) : _alpha(1.0), SigmaI(dat, outFlNam, prDiag, nu0) {};
	
	~SigmaIpex() {};
	
	void update(const Grp &dat, const Qgrp &q);
	void update(const Grp &dat, const Grp &mu, const Qgrp &q);
	
	void save(const char *how = "a");
	void save(const string &fileNam, const char *how = "a");
};

class StTq {	// the weigting parameter for Student-t sampling
private:
	double _q;
	const double *_nu; // Student-t degrees of freedom; pointer to a value kept in Qgrp, since all members of the group should have the same
	size_t _locInd;
public:
	// constructors
	StTq() : _nu(0), _q(1.0), _locInd(0) {};
	StTq(const double &nu) : _q(1.0), _locInd(0) { _nu = &nu; };
	StTq(const double &nu, const int &d, gsl_rng *r);
	StTq(const double &q, const double &nu, const int &d, gsl_rng *r);
	StTq(const double &nu, const size_t &ind) : _q(1.0), _locInd(ind) { _nu = &nu; };
	StTq(const double &q, const double &nu, const size_t &ind, const int &d, gsl_rng *r);
	
	~StTq() {};
	
	double getVal() const {return _q; };
	double getNu() const {return *_nu; };
	
	void update(const Grp &dat, const Grp &mu, const SigmaI &SigI, const gsl_rng *r);
	void update(const Grp &dat, const SigmaI &SigI, const gsl_rng *r);
	// for the PEX scheme
	void update(const Grp &dat, const Grp &mu, const SigmaI &SigI, const double &alpha, const gsl_rng *r);
	void update(const Grp &dat, const SigmaI &SigI, const double &alpha, const gsl_rng *r);
	
	// save function, taking file name, appending by default
	void save(const string &fileNam, const char *how = "a");
	// save function, taking file stream name
	void save(FILE *fileStr);
};

class Qgrp { // group of StTq
protected:
	vector<StTq> _qVec;
	vector<size_t> _presInd; // vector of indexes of samples with no missing phenotypes.  These are the ones we can update q for
	double _nu;  // degrees of freedom
	
	gsl_rng *_r;
	
public:
	Qgrp() : _nu(3.0) {};
	Qgrp(const size_t &N);
	Qgrp(const size_t &N, const double &nu);
	Qgrp(const size_t &N, const double &nu, const string &misVecFlNam);
	
	virtual ~Qgrp() {};
	
	double operator[](const size_t i) const{return _qVec[i].getVal(); };
	double operator[](const size_t i) {return _qVec[i].getVal(); };
	
	virtual double alpha() const {return 1.0; };
	virtual double alpha() {return 1.0; };
	
	virtual void update(const Grp &dat, const Grp &mu, const SigmaI &SigI);
	virtual void update(const Grp &dat, const SigmaI &SigI);
	
};

class QgrpPEX : public Qgrp {  // group of StTq for van Dyk and Meng's PEX
protected:
	double _alpha;  // the extra parameter
	
public:
	QgrpPEX() : _alpha(1.0), Qgrp() {};
	QgrpPEX(const size_t &N) : _alpha(1.0), Qgrp(N) {};
	QgrpPEX(const size_t &N, const double &nu) : _alpha(1.0), Qgrp(N, nu) {};
	QgrpPEX(const size_t &N, const double &nu, const string &misVecFlNam) : _alpha(1.0), Qgrp(N, nu, misVecFlNam) {};
	
	double alpha() const {return _alpha; };
	double alpha() {return _alpha; };
	
	void update(const Grp &dat, const Grp &mu, const SigmaI &SigI);
	void update(const Grp &dat, const SigmaI &SigI);
};

class MixP {
private:
	vector<double> _p;
	vector<double> _alpha; // prior
	gsl_rng *_r;
	
public:
	MixP();
	MixP(const vector<double> initP) : _p(initP), _alpha(initP.size(), 1.0) {};
	MixP(const gsl_vector *initP, const size_t &len);
	MixP(const vector<double> initP, const double &a);
	MixP(const gsl_vector *initP, const size_t &len, const double &a);
	MixP(const vector<double> initP, const double *a, const size_t &len);
	MixP(const gsl_vector *initP, const size_t &len, const double *a);
	MixP(const vector<double> initP, const vector<double> &a);
	MixP(const gsl_vector *initP, const size_t &len, const vector<double> &a);
	MixP(const RanIndex &ind, const double &a);
	
	~MixP() {};
	
	MixP(const MixP &);
	MixP & operator=(const MixP &);
	const double &operator[](const size_t i) const{return _p[i]; };
	
	void update(const RanIndex &Nvec);
	
};





#endif
