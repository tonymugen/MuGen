/*
*  HBqgen.h
*  HBqgen
*
*  Created by ajg67 on 10/30/12.
*  Copyright (c) 2012 SEELE. All rights reserved.
*
*  Classes for Hierarchical Bayesian quantitative-genetic models.  Using gsl_vector and gsl_matrix as internal storage types.
*  MVnormMu and MVnormBeta have overloaded update() methods to accomodate MV Student-t  and Gaussian models
*
*/

#ifndef HBQGEN_H
#define HBQGEN_H
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_permutation.h>
#include <vector>
#include <string>

using std::vector;
using std::string;

class MVnorm;
class MVnormMu;
class MVnormMuPEX;
class MVnormMuMiss;
class MVnormBeta;
class MVnormBetaMiss;
class MVnormBetaFt;
class MVnormBetaPC;
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
class MuGrpMiss;
class BetaGrpSnp;
class BetaGrpSnpMiss;
class BetaGrpFt;
class BetaGrpSc;
class BetaGrpPC;
class BetaGrpBVSR;
class MuBlk;
class BetaBlk;
class SigmaI;
class SigmaIblk;
class StTq;
class Qgrp;
class MixP;

/*
	various distribution functions
*/
void MVgauss(const gsl_vector*, const gsl_matrix*, const int, const gsl_rng*, gsl_vector*);
void Wishart(const gsl_matrix*, const int&, const int&, const gsl_rng*, gsl_matrix*);
void MVgauss(const gsl_vector*, const gsl_matrix*, const size_t, const gsl_rng*, gsl_vector*);
void Wishart(const gsl_matrix*, const size_t&, const size_t&, const gsl_rng*, gsl_matrix*);
size_t rtgeom(const double &, const size_t &, const gsl_rng *);

/*
	auxiliary functions for the sampler
*/
void rspRsd(const MVnormMu *rsp, const MVnormMu *ftd, const int &N, const int &d, gsl_matrix *rsd);
double plsOne(const gsl_matrix *resp, const gsl_vector *pred, const vector<int> &pres, const int &d, const gsl_rng *r);
void colCenter(gsl_matrix *inplace);
void colCenter(const gsl_matrix *source, gsl_matrix *res);
void colCenter(gsl_matrix *inplace, const double &absLab); // when there are missing data labeled by absLab
void colCenter(const gsl_matrix *source, gsl_matrix *res, const double &absLab);
void printMat(const gsl_matrix *);
unsigned long long rdtsc();


class MVnorm {
protected:
	gsl_vector_view _vec; // the vector of values
	size_t _d;	          // dimension, needs to be the same for all objects in a given program
	
	// constructors
	MVnorm() : _d(0) {};
	MVnorm(const size_t &d) : _d(d) {};
	MVnorm(gsl_vector *mn);  // deterministic constructor
	MVnorm(gsl_vector *mn, const gsl_vector *sd, const gsl_rng *r);
	MVnorm(gsl_vector *mn, const gsl_matrix *Sig, const gsl_rng *r);
	MVnorm(gsl_matrix *mn, const size_t &iRw);  // deterministic constructor
	MVnorm(gsl_matrix *mn, const size_t &iRw, const gsl_vector *sd, const gsl_rng *r);
	MVnorm(gsl_matrix *mn, const size_t &iRw, const gsl_matrix *Sig, const gsl_rng *r);
	
public:
	MVnorm(const MVnorm &); // copy constructor
	MVnorm & operator=(const MVnorm &);
	
	virtual ~MVnorm() {}; // necessary virtual destructor
	
	//	overloaded update() methods; prior q values are the const double &
	
	// improper prior methods
	virtual void update(const Grp &, const SigmaI &, const gsl_rng *) = 0;
	virtual void update(const Grp &, const Qgrp &, const SigmaI &, const gsl_rng *) = 0;
	// 0-mean prior methods
	virtual void update(const Grp &, const SigmaI &, const SigmaI &, const gsl_rng *) = 0; // G data, G prior
	virtual void update(const Grp &, const SigmaI &, const double &, const SigmaI &, const gsl_rng *) = 0; // G data, t prior
	virtual void update(const Grp &, const Qgrp &, const SigmaI &, const SigmaI &, const gsl_rng *) = 0; // t data, G prior
	virtual void update(const Grp &, const Qgrp &, const SigmaI &, const double &, const SigmaI &, const gsl_rng *) = 0; // t data, t prior
	// non-zero mean prrior methods
	virtual void update(const Grp &, const SigmaI &, const Grp &, const SigmaI &, const gsl_rng *) = 0; // G data, G prior
	virtual void update(const Grp &, const SigmaI &, const Grp &, const double &, const SigmaI &, const gsl_rng *) = 0; // G data, t prior
	virtual void update(const Grp &, const Qgrp &, const SigmaI &, const Grp &, const SigmaI &, const gsl_rng *) = 0; // t data, G prior
	virtual void update(const Grp &, const Qgrp &, const SigmaI &, const Grp &, const double &, const SigmaI &, const gsl_rng *) = 0; // t data, t prior
	
	virtual double mhl(const MVnorm *x, const SigmaI &SigI);
	virtual double mhl(const SigmaI &SigI); // distance from 0
	virtual double mhl(const gsl_vector *x, const SigmaI &SigI);
	virtual double mhl(const MVnorm *x, const SigmaI &SigI) const;
	virtual double mhl(const SigmaI &SigI) const; // distance from 0
	virtual double mhl(const gsl_vector *x, const SigmaI &SigI) const;
	
	double density(const gsl_vector *theta, const SigmaI &SigI);
	double density(const MVnorm *theta, const SigmaI &SigI);
	double density(const gsl_vector *theta, const SigmaI &SigI) const;
	double density(const MVnorm *theta, const SigmaI &SigI) const;
	
	void save(const string &fileNam, const char *how = "a");
	void save(FILE *fileStr);
	
	double operator[](const size_t i) const{return gsl_vector_get(&_vec.vector, i); };
	void valSet(const size_t i, const double x){gsl_vector_set(&_vec.vector, i, x); };
	const gsl_vector *getVec() const{return &_vec.vector;};
	size_t len() const{return _d; };
	virtual size_t nMissP() const{return 0; };
	virtual const vector<size_t> getMisPhen() const{ return vector<size_t>(0); };
	virtual const vector<size_t> *down() const{return 0; };
	virtual const size_t *up() const{return 0; };
	virtual double scalePar() const {return 1.0; };
};

class MVnormMu : public MVnorm {
protected:
	const vector<size_t> *_lowLevel; // pointer to a vector of lower level (data) indexes; pointers are const to make sure they are not used to modify the indexes, but the indexes themselves may change (as in the mixture model)
	const size_t *_upLevel;          // pointer to the upper-level (prior) index
	
public:
	// constructors
	MVnormMu();								   // default constructor, necessary if a user plans to allocate arrays of MVnormMu on the heap
	MVnormMu(const size_t &d) : MVnorm(d) {};  // constructor for a vague 0-mean prior
	MVnormMu(const size_t &d, const vector<size_t> &low, const size_t &up);
	MVnormMu(gsl_vector *mn, const vector<size_t> &low, const size_t &up); // deterministic constructor
	MVnormMu(gsl_vector *mn, const size_t &up); // deterministic constructor for the lowest hierarchy level: *_lowLevel will be set to 0
	MVnormMu(gsl_matrix *mn, const size_t &iRw); // deterministic constructor for, e.g., a sum. Both *_lowLevel and  *_upLevel will be set to 0
	MVnormMu(gsl_matrix *mn, const size_t &iRw, const vector<size_t> &low); // deterministic constructor for the highest hierarchy level: *_upLevel will be set to 0
	MVnormMu(gsl_vector *mn, const gsl_vector *sd, const gsl_rng *r, const vector<size_t> &low, const size_t &up);
	MVnormMu(gsl_vector *mn, const gsl_matrix *Sig, const gsl_rng *r, const vector<size_t> &low, const size_t &up);
	MVnormMu(gsl_matrix *mn, const size_t &iRw, const vector<size_t> &low, const size_t &up); // deterministic constructor
	MVnormMu(gsl_matrix *mn, const size_t &iRw, const size_t &up); // deterministic constructor for the lowest hierarchy level: *_lowLevel will be set to 0
	MVnormMu(gsl_matrix *mn, const size_t &iRw, const gsl_vector *sd, const gsl_rng *r, const vector<size_t> &low, const size_t &up);
	MVnormMu(gsl_matrix *mn, const size_t &iRw, const gsl_matrix *Sig, const gsl_rng *r, const vector<size_t> &low, const size_t &up);
	
	MVnormMu(const MVnormMu &); // copy constructor
	MVnormMu & operator=(const MVnormMu &);
	
	virtual ~MVnormMu(); // destructor
	
	//	overloaded update() methods
	// improper prior methods
	virtual void update(const Grp &dat, const SigmaI &SigIm, const gsl_rng *r);
	virtual void update(const Grp &dat, const Qgrp &q, const SigmaI &SigIm, const gsl_rng *r);
	// 0-mean prior methods
	virtual void update(const Grp &dat, const SigmaI &SigIm, const SigmaI &SigIp, const gsl_rng *r);
	void update(const Grp &dat, const SigmaI &SigIm, const double &qPr, const SigmaI &SigIp, const gsl_rng *r);
	virtual void update(const Grp &dat, const Qgrp &q, const SigmaI &SigIm, const SigmaI &SigIp, const gsl_rng *r);
	void update(const Grp &dat, const Qgrp &q, const SigmaI &SigIm, const double &qPr, const SigmaI &SigIp, const gsl_rng *r);
	// non-zero mean prior methods
	virtual void update(const Grp &dat, const SigmaI &SigIm, const Grp &muPr, const SigmaI &SigIp, const gsl_rng *r);
	void update(const Grp &dat, const SigmaI &SigIm, const Grp &muPr, const double &qPr, const SigmaI &SigIp, const gsl_rng *r);
	void update(const Grp &dat, const Qgrp &q, const SigmaI &SigIm, const Grp &muPr, const SigmaI &SigIp, const gsl_rng *r);
	void update(const Grp &dat, const Qgrp &q, const SigmaI &SigIm, const Grp &muPr, const double &qPr, const SigmaI &SigIp, const gsl_rng *r);
	
	virtual size_t nMissP() const{return 0; };
	virtual const vector<size_t> getMisPhen() const{ return vector<size_t>(0); };
	const vector<size_t> *down() const{return _lowLevel; };
	const size_t *up() const{return _upLevel; };
};

class MVnormMuPEX : public MVnormMu {
protected:
	Apex *_A;
	gsl_matrix_view _tSAprod;
	
public:
	MVnormMuPEX() : MVnormMu() {_A = 0; };
	MVnormMuPEX(gsl_matrix *mn, const size_t &iRw, const vector<size_t> &low, const size_t &up, Apex &A, gsl_matrix *tSigIAt); // deterministic constructor
	MVnormMuPEX(gsl_matrix *mn, const size_t &iRw, const gsl_vector *sd, const gsl_rng *r, const vector<size_t> &low, const size_t &up, Apex &A, gsl_matrix *tSigIAt);
	
	MVnormMuPEX(const MVnormMuPEX &mu); // copy constructor
	MVnormMuPEX & operator=(const MVnormMuPEX &mu);
	
	~MVnormMuPEX(){};
	
	void update(const Grp &dat, const SigmaI &SigIm, const SigmaI &SigIp, const gsl_rng *r);
	void update(const Grp &dat, const SigmaI &SigIm, const Grp &muPr, const SigmaI &SigIp, const gsl_rng *r);
};

class MVnormMuMiss : public MVnormMu { // further modifies MVnormMu to account for missing phenotypes
protected:
	vector<size_t> _misPhenInd; // vector of indexes tagging missing data
	size_t _myInd;              // it's own index; defined only for instances that take iRw
	
public:
	// constructors
	MVnormMuMiss() : MVnormMu(){};                   // default constructor, necessary if a user plans to allocate arrays of MVnormMuMiss on the heap
	MVnormMuMiss(const size_t &d) : MVnormMu(d) {};  // constructor for a vague 0-mean prior
	MVnormMuMiss(const size_t &d, const vector<size_t> &low, const size_t &up, const vector<size_t> &mis);
	MVnormMuMiss(gsl_vector *mn, const vector<size_t> &low, const size_t &up, const vector<size_t> &mis); // deterministic constructor
	MVnormMuMiss(gsl_vector *mn, const size_t &up, const vector<size_t> &mis); // deterministic constructor
	MVnormMuMiss(gsl_vector *mn, const gsl_vector *sd, const gsl_rng *r, const vector<size_t> &low, const size_t &up, const vector<size_t> &mis);
	MVnormMuMiss(gsl_vector *mn, const gsl_matrix *Sig, const gsl_rng *r, const vector<size_t> &low, const size_t &up, const vector<size_t> &mis);
	MVnormMuMiss(gsl_matrix *mn, const size_t &iRw, const vector<size_t> &low, const size_t &up, const vector<size_t> &mis); // deterministic constructor
	MVnormMuMiss(gsl_matrix *mn, const size_t &iRw, const size_t &up, const vector<size_t> &mis); // deterministic constructor
	MVnormMuMiss(gsl_matrix *mn, const size_t &iRw, const gsl_vector *sd, const gsl_rng *r, const vector<size_t> &low, const size_t &up, const vector<size_t> &mis);
	MVnormMuMiss(gsl_matrix *mn, const size_t &iRw, const gsl_matrix *Sig, const gsl_rng *r, const vector<size_t> &low, const size_t &up, const vector<size_t> &mis);
	
	MVnormMuMiss(const MVnormMuMiss &); // copy constructor
	MVnormMuMiss & operator=(const MVnormMuMiss &);
	
	~MVnormMuMiss() {}; 	// destructor
	
	
	// Phenotype imputation (Gaussian model only -- Student-t crashes due to a correlation of q and x_imp), improper prior
	void update(const Grp &mu, const SigmaI &SigIm, const gsl_rng *r);
	// imputation with a proper Gaussian prior
	void update(const Grp &mu, const SigmaI &SigIm, const SigmaI &SigIp, const gsl_rng *r);
	
	size_t nMissP() const{return _misPhenInd.size(); };
	const vector<size_t> getMisPhen() const{ return _misPhenInd; };
};

class MVnormBeta : public MVnorm { // by itself to be used strictly for regressions with a single predictor (or ignoring the effects of other predictors)
protected:
	gsl_vector_view _X; // vector view of predictor values (typically column of a predictor matrix that is in the encapsulating group class).  The matrix this points to can be modified through this, as when I scale the predictor
	double _scale;      // Scale -- typically XtX, but not always (special cases dealt with in derived classes)
	size_t _N;          // length of &_X.vector
	
	const size_t *_upLevel; // index into the vector of priors
	
public:
	//	constructors, similar to MVnormMu
	MVnormBeta();
	MVnormBeta(const size_t d);
	MVnormBeta(gsl_vector *b, const gsl_vector *sd, gsl_matrix *pred, const size_t &iCl, const gsl_rng *r);
	MVnormBeta(gsl_vector *b, const gsl_matrix *Sig, gsl_matrix *pred, const size_t &iCl, const gsl_rng *r);
	MVnormBeta(const gsl_matrix *resp, gsl_matrix *pred, const size_t &iCl, const gsl_matrix *Sig, const gsl_rng *r, gsl_matrix *bet, const size_t &iRw);
	MVnormBeta(const Grp &resp, gsl_matrix *pred, const size_t &iCl, const gsl_matrix *Sig, const gsl_rng *r, gsl_matrix *bet, const size_t &iRw);
	MVnormBeta(gsl_matrix *pred, const size_t &iCl, gsl_matrix *bet, const size_t &iRw); // deterministic constructor that just points to whatever is already in *bet
	
	MVnormBeta(gsl_vector *b, const gsl_vector *sd, gsl_matrix *pred, const size_t &iCl, const gsl_rng *r, const size_t &up);
	MVnormBeta(gsl_vector *b, const gsl_matrix *Sig, gsl_matrix *pred, const size_t &iCl, const gsl_rng *r, const size_t &up);
	MVnormBeta(const gsl_matrix *resp, gsl_matrix *pred, const size_t &iCl, const gsl_matrix *Sig, const gsl_rng *r, const size_t &up, gsl_matrix *bet, const size_t &iRw);
	MVnormBeta(const Grp &resp, gsl_matrix *pred, const size_t &iCl, const gsl_matrix *Sig, const gsl_rng *r, const size_t &up, gsl_matrix *bet, const size_t &iRw);
	MVnormBeta(gsl_matrix *pred, const size_t &iCl, const size_t &up, gsl_matrix *bet, const size_t &iRw);  // deterministic constructor that just points to whatever is already in *bet
	
	MVnormBeta(const MVnormBeta &); // copy constructor
	MVnormBeta & operator=(const MVnormBeta &);
	
	virtual ~MVnormBeta(); 	// destructor
	
	// legacy C/GSL methods
	void update(const gsl_matrix *resp, const SigmaI &SigIb, const gsl_rng *r);
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

class MVnormBetaCmp : public MVnormBeta { // for simple SNP regressions on complete SNP data
protected:
	// we multiply &_X.vector inherited from MVnormBeta by 1/_scale at initialization so that we can do the generic regression (XtX)^-1Xt faster.  This modifies the matrix _X is pointing to
	double _mhl(const SigmaI &SigI);
public:
	MVnormBetaCmp() : MVnormBeta() {};
	MVnormBetaCmp(gsl_matrix *pred, const size_t &iCl, gsl_matrix *bet, const size_t &iRw);
	MVnormBetaCmp(const gsl_matrix *resp, gsl_matrix *pred, const size_t &iCl, const gsl_matrix *Sig, const gsl_rng *r, gsl_matrix *bet, const size_t &iRw) : MVnormBeta(resp, pred, iCl, Sig, r, bet, iRw) {gsl_vector_scale(&_X.vector, 1.0/_scale); };
	MVnormBetaCmp(const Grp &resp, gsl_matrix *pred, const size_t &iCl, const gsl_matrix *Sig, const gsl_rng *r, gsl_matrix *bet, const size_t &iRw) : MVnormBeta(resp, pred, iCl, Sig, r, bet, iRw) {gsl_vector_scale(&_X.vector, 1.0/_scale); };
	
	MVnormBetaCmp(const MVnormBetaCmp &); // copy constructor
	MVnormBetaCmp & operator=(const MVnormBetaCmp &);
	
	~MVnormBetaCmp() {}; 	// destructor
	
	double mhl(const MVnorm *x, const SigmaI &SigI){return gsl_vector_get(&_vec.vector, _d - 1); };
	double mhl(const SigmaI &SigI){return gsl_vector_get(&_vec.vector, _d - 1); }; // distance from 0
	double mhl(const gsl_vector *x, const SigmaI &SigI){return gsl_vector_get(&_vec.vector, _d - 1); };
	double mhl(const MVnorm *x, const SigmaI &SigI) const{return gsl_vector_get(&_vec.vector, _d - 1); };
	double mhl(const SigmaI &SigI) const{return gsl_vector_get(&_vec.vector, _d - 1); }; // distance from 0
	double mhl(const gsl_vector *x, const SigmaI &SigI) const{return gsl_vector_get(&_vec.vector, _d - 1); };
	
	void update(const Grp &dat, const SigmaI &SigIb, const gsl_rng *r);
};

class MVnormBetaMiss : public MVnormBeta {
protected:
	vector<size_t> _presInd;
	
	double _mhl(const SigmaI &SigI);
public:
	MVnormBetaMiss() : MVnormBeta() {};
	MVnormBetaMiss(vector<double> &pred, const vector<size_t> &pres, gsl_matrix *bet, const size_t &iRw);
	MVnormBetaMiss(const gsl_matrix *resp, gsl_matrix *pred, const size_t &iCl, const gsl_matrix *Sig, const gsl_rng *r, const vector<size_t> &pres, gsl_matrix *bet, const size_t &iRw);
	MVnormBetaMiss(const Grp &resp, gsl_matrix *pred, const size_t &iCl, const gsl_matrix *Sig, const gsl_rng *r, const vector<size_t> &pres, gsl_matrix *bet, const size_t &iRw);
	
	MVnormBetaMiss(const MVnormBetaMiss &); // copy constructor
	MVnormBetaMiss & operator=(const MVnormBetaMiss &);
	
	~MVnormBetaMiss() {}; 	// destructor
	
	double mhl(const MVnorm *x, const SigmaI &SigI){return gsl_vector_get(&_vec.vector, _d - 1); };
	double mhl(const SigmaI &SigI){return gsl_vector_get(&_vec.vector, _d - 1); }; // distance from 0
	double mhl(const gsl_vector *x, const SigmaI &SigI){return gsl_vector_get(&_vec.vector, _d - 1); };
	double mhl(const MVnorm *x, const SigmaI &SigI) const{return gsl_vector_get(&_vec.vector, _d - 1); };
	double mhl(const SigmaI &SigI) const{return gsl_vector_get(&_vec.vector, _d - 1); }; // distance from 0
	double mhl(const gsl_vector *x, const SigmaI &SigI) const{return gsl_vector_get(&_vec.vector, _d - 1); };
	
	void update(const Grp &dat, const SigmaI &SigIb, const gsl_rng *r);
	
};

class MVnormBetaFt: public MVnormBeta {
protected:
	gsl_matrix_view _fitted;
	
public:
	MVnormBetaFt(): MVnormBeta() {};
	MVnormBetaFt(gsl_vector *b, const gsl_vector *sd, gsl_matrix *pred, const size_t &iCl, gsl_matrix *allFt, const size_t &begRw, const gsl_rng *r) : MVnormBeta(b, sd, pred, iCl, r) {_fitted = gsl_matrix_submatrix(allFt, begRw, 0, pred->size1, _d); };
	MVnormBetaFt(gsl_vector *b, const gsl_matrix *Sig, gsl_matrix *pred, const size_t &iCl, gsl_matrix *allFt, const size_t &begRw, const gsl_rng *r) : MVnormBeta(b, Sig, pred, iCl, r) {_fitted = gsl_matrix_submatrix(allFt, begRw, 0, pred->size1, _d); };
	MVnormBetaFt(const gsl_matrix *resp, gsl_matrix *pred, const size_t &iCl, gsl_matrix *allFt, const size_t &begRw, const gsl_matrix *Sig, const gsl_rng *r, gsl_matrix *bet, const size_t &iRw) : MVnormBeta(resp, pred, iCl, Sig, r, bet, iRw) {_fitted = gsl_matrix_submatrix(allFt, begRw, 0, pred->size1, _d); };
	MVnormBetaFt(const gsl_matrix *resp, gsl_matrix *pred, const size_t &iCl, vector<double> &eaFt, const gsl_matrix *Sig, const gsl_rng *r, gsl_matrix *bet, const size_t &iRw) : MVnormBeta(resp, pred, iCl, Sig, r, bet, iRw) {_fitted = gsl_matrix_view_array(eaFt.data(), pred->size1, _d); };
	MVnormBetaFt(const Grp &resp, gsl_matrix *pred, const size_t &iCl, gsl_matrix *allFt, const size_t &begRw, const gsl_matrix *Sig, const gsl_rng *r, gsl_matrix *bet, const size_t &iRw) : MVnormBeta(resp, pred, iCl, Sig, r, bet, iRw) {_fitted = gsl_matrix_submatrix(allFt, begRw, 0, pred->size1, _d); }; // deprecated
	MVnormBetaFt(const Grp &resp, gsl_matrix *pred, const size_t &iCl, vector<double> &eaFt, const gsl_matrix *Sig, const gsl_rng *r, gsl_matrix *bet, const size_t &iRw) : MVnormBeta(resp, pred, iCl, Sig, r, bet, iRw) {_fitted = gsl_matrix_view_array(eaFt.data(), pred->size1, _d); };
	MVnormBetaFt(gsl_matrix *pred, const size_t &iCl, vector<double> &eaFt, gsl_matrix *bet, const size_t &iRw) : MVnormBeta(pred, iCl, bet, iRw) {_fitted = gsl_matrix_view_array(eaFt.data(), pred->size1, _d); };
	
	MVnormBetaFt(gsl_vector *b, const gsl_vector *sd, gsl_matrix *pred, const size_t &iCl, gsl_matrix *allFt, const size_t &begRw, const gsl_rng *r, const size_t &up) : MVnormBeta(b, sd, pred, iCl, r, up) {_fitted = gsl_matrix_submatrix(allFt, begRw, 0, pred->size1, _d); };
	MVnormBetaFt(gsl_vector *b, const gsl_matrix *Sig, gsl_matrix *pred, const size_t &iCl, gsl_matrix *allFt, const size_t &begRw, const gsl_rng *r, const size_t &up) : MVnormBeta(b, Sig, pred, iCl, r, up) {_fitted = gsl_matrix_submatrix(allFt, begRw, 0, pred->size1, _d); };
	MVnormBetaFt(const gsl_matrix *resp, gsl_matrix *pred, const size_t &iCl, gsl_matrix *allFt, const size_t &begRw, const gsl_matrix *Sig, const gsl_rng *r, const size_t &up, gsl_matrix *bet, const size_t &iRw) : MVnormBeta(resp, pred, iCl, Sig, r, up, bet, iRw) {_fitted = gsl_matrix_submatrix(allFt, begRw, 0, pred->size1, _d); };
	MVnormBetaFt(const gsl_matrix *resp, gsl_matrix *pred, const size_t &iCl, vector<double> &eaFt, const gsl_matrix *Sig, const gsl_rng *r, const size_t &up, gsl_matrix *bet, const size_t &iRw) : MVnormBeta(resp, pred, iCl, Sig, r, up, bet, iRw) {_fitted = gsl_matrix_view_array(eaFt.data(), pred->size1, _d); };
	MVnormBetaFt(const Grp &resp, gsl_matrix *pred, const size_t &iCl, gsl_matrix *allFt, const size_t &begRw, const gsl_matrix *Sig, const gsl_rng *r, const size_t &up, gsl_matrix *bet, const size_t &iRw) : MVnormBeta(resp, pred, iCl, Sig, r, up, bet, iRw) {_fitted = gsl_matrix_submatrix(allFt, begRw, 0, pred->size1, _d); };
	MVnormBetaFt(const Grp &resp, gsl_matrix *pred, const size_t &iCl, vector<double> &eaFt, const gsl_matrix *Sig, const gsl_rng *r, const size_t &up, gsl_matrix *bet, const size_t &iRw) : MVnormBeta(resp, pred, iCl, Sig, r, up, bet, iRw) {_fitted = gsl_matrix_view_array(eaFt.data(), pred->size1, _d); };
	MVnormBetaFt(gsl_matrix *pred, const size_t &iCl, vector<double> &eaFt, const size_t &up, gsl_matrix *bet, const size_t &iRw) : MVnormBeta(pred, iCl, up, bet, iRw) {_fitted = gsl_matrix_view_array(eaFt.data(), pred->size1, _d); };

	MVnormBetaFt(const MVnormBetaFt &); // copy constructor
	MVnormBetaFt & operator=(const MVnormBetaFt &);
	
	virtual ~MVnormBetaFt() {}; 	// destructor
	
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

class MVnormBetaPC : public MVnormBetaFt {
protected:
	double _lam; // the eigenvalue
	
public:
	MVnormBetaPC() : MVnormBetaFt() {};
	MVnormBetaPC(gsl_vector *b, const gsl_vector *sd, gsl_matrix *pred, const size_t &iCl, gsl_matrix *allFt, const size_t &begRw, const gsl_rng *r, const double &lm) : MVnormBetaFt(b, sd, pred, iCl, allFt, begRw, r) { _lam = lm; };
	MVnormBetaPC(gsl_vector *b, const gsl_matrix *Sig, gsl_matrix *pred, const size_t &iCl, gsl_matrix *allFt, const size_t &begRw, const gsl_rng *r, const double &lm) : MVnormBetaFt(b, Sig, pred, iCl, allFt, begRw, r) { _lam = lm; };
	MVnormBetaPC(const gsl_matrix *resp, gsl_matrix *pred, const size_t &iCl, gsl_matrix *allFt, const size_t &begRw, const gsl_matrix *Sig, const gsl_rng *r, gsl_matrix *bet, const size_t &iRw, const double &lm) : MVnormBetaFt(resp, pred, iCl, allFt, begRw, Sig, r, bet, iRw) { _lam = lm; };
	MVnormBetaPC(const Grp &resp, gsl_matrix *pred, const size_t &iCl, gsl_matrix *allFt, const size_t &begRw, const gsl_matrix *Sig, const gsl_rng *r, gsl_matrix *bet, const size_t &iRw, const double &lm) : MVnormBetaFt(resp, pred, iCl, allFt, begRw, Sig, r, bet, iRw) { _lam = lm; };
	MVnormBetaPC(const Grp &resp, gsl_matrix *pred, const size_t &iCl, vector<double> &eaFt, const gsl_matrix *Sig, const gsl_rng *r, gsl_matrix *bet, const size_t &iRw, const double &lm) : MVnormBetaFt(resp, pred, iCl, eaFt, Sig, r, bet, iRw) { _lam = lm; };
	MVnormBetaPC(const gsl_matrix *resp, gsl_matrix *pred, const size_t &iCl, vector<double> &eaFt, const gsl_matrix *Sig, const gsl_rng *r, gsl_matrix *bet, const size_t &iRw, const double &lm) : MVnormBetaFt(resp, pred, iCl, eaFt, Sig, r, bet, iRw) { _lam = lm; };
	
	MVnormBetaPC(gsl_vector *b, const gsl_vector *sd, gsl_matrix *pred, const size_t &iCl, gsl_matrix *allFt, const size_t &begRw, const gsl_rng *r, const size_t &up, const double &lm) : MVnormBetaFt(b, sd, pred, iCl, allFt, begRw, r, up) { _lam = lm; };
	MVnormBetaPC(gsl_vector *b, const gsl_matrix *Sig, gsl_matrix *pred, const size_t &iCl, gsl_matrix *allFt, const size_t &begRw, const gsl_rng *r, const size_t &up, const double &lm) : MVnormBetaFt(b, Sig, pred, iCl, allFt, begRw, r, up) { _lam = lm; };
	MVnormBetaPC(const gsl_matrix *resp, gsl_matrix *pred, const size_t &iCl, gsl_matrix *allFt, const size_t &begRw, const gsl_matrix *Sig, const gsl_rng *r, const size_t &up, gsl_matrix *bet, const size_t &iRw, const double &lm) : MVnormBetaFt(resp, pred, iCl, allFt, begRw, Sig, r, up, bet, iRw) { _lam = lm; };
	MVnormBetaPC(const Grp &resp, gsl_matrix *pred, const size_t &iCl, gsl_matrix *allFt, const size_t &begRw, const gsl_matrix *Sig, const gsl_rng *r, const size_t &up, gsl_matrix *bet, const size_t &iRw, const double &lm) : MVnormBetaFt(resp, pred, iCl, allFt, begRw, Sig, r, up, bet, iRw) { _lam = lm; };
	MVnormBetaPC(const Grp &resp, gsl_matrix *pred, const size_t &iCl, vector<double> &eaFt, const gsl_matrix *Sig, const gsl_rng *r, const size_t &up, gsl_matrix *bet, const size_t &iRw, const double &lm) : MVnormBetaFt(resp, pred, iCl, eaFt, Sig, r, up, bet, iRw) { _lam = lm; };
	MVnormBetaPC(const gsl_matrix *resp, gsl_matrix *pred, const size_t &iCl, vector<double> &eaFt, const gsl_matrix *Sig, const gsl_rng *r, const size_t &up, gsl_matrix *bet, const size_t &iRw, const double &lm) : MVnormBetaFt(resp, pred, iCl, eaFt, Sig, r, up, bet, iRw) { _lam = lm; };
	
	MVnormBetaPC(const MVnormBetaPC &); // copy constructor
	MVnormBetaPC & operator=(const MVnormBetaPC &);
	
	~MVnormBetaPC() {}; 	// destructor
	
	// 0-mean prior
	// Gaussian model, Gaussian prior
	void update(const Grp &resp, const SigmaI &SigIb, const SigmaI &SigIp, const gsl_rng *r);
	// Gaussian model, Student-t prior
	void update(const Grp &resp, const SigmaI &SigIb, const double &qPr, const SigmaI &SigIp, const gsl_rng *r);
	// Student-t model, Gaussian prior
	void update(const Grp &resp, const Qgrp &q, const SigmaI &SigIb, const SigmaI &SigIp, const gsl_rng *r);
	// Student-t model, Student-t prior
	void update(const Grp &resp, const Qgrp &q, const SigmaI &SigIb, const double &qPr, const SigmaI &SigIp, const gsl_rng *r);
	
	// non-0-mean prior
	// Gaussian model, Gaussian prior
	void update(const Grp &resp, const SigmaI &SigIb, const Grp &muPr, const SigmaI &SigIp, const gsl_rng *r);
	// Gaussian model, Student-t prior
	void update(const Grp &resp, const SigmaI &SigIb, const Grp &muPr, const double &qPr, const SigmaI &SigIp, const gsl_rng *r);
	// Student-t model, Gaussian prior
	void update(const Grp &resp, const Qgrp &q, const SigmaI &SigIb, const Grp &muPr, const SigmaI &SigIp, const gsl_rng *r);
	// Student-t model, Student-t prior
	void update(const Grp &resp, const Qgrp &q, const SigmaI &SigIb, const Grp &muPr, const double &qPr, const SigmaI &SigIp, const gsl_rng *r);
	
};

class MVnormMuBlk : public MVnorm {
protected:
	const size_t *_upLevel;                  // pointer to the upper-level (prior) index -- the same for all blocks
	const vector< size_t > *_blkStart;       // starts of vector pieces; pointer to the vector that's in the corresponding Grp class
	vector< gsl_vector_view > _eachVec;      // views of vector pieces, each corresponding to a block of variables
	const vector< vector<size_t> > *_eachLL; // lowLevel (LL) for each block
	
public:
	MVnormMuBlk() : MVnorm() {};
	MVnormMuBlk(gsl_matrix *mn, const size_t &iRw, const vector<size_t> &blkStart, const size_t &up);
	MVnormMuBlk(gsl_matrix *mn, const size_t &iRw, const gsl_vector *sd, const gsl_rng *r, const vector<size_t> &blkStart, const vector< vector<size_t> > &eachLL, const size_t &up);
	
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
/*
	random index class that allows for mixture models.  Can be used as a deterministic index if there is deterministic initiation and no updating
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
	Index class that keeps track of the variable selection process.
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
	Parameter-expansion variable for 0-mean MVnormMu
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
	RanIndex *_lowLevel;   // points to the level below.  For Betas will correspond to the _fittedValues; if unset will be 0
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
	void update(const Grp &dat, const SigmaI &SigIm, const Grp &muPr, const SigmaI &SigIp);
};

class MuGrpMiss : public MuGrp {
protected:
	vector<size_t> _misInd;
	
public:
	MuGrpMiss(const string &datFlNam, const string &misMatFlNam, const string &misVecFlNam, RanIndex &up, const size_t &d);
	
	~MuGrpMiss() {};
	
	MuGrpMiss(const MuGrpMiss &mG); // copy constructor
	MuGrpMiss & operator=(const MuGrpMiss &mG);
	
	// only Gaussian updates work for imputation (in general)
	void update(const Grp &mu, const SigmaI &SigIm);
	void update(const Grp &mu, const SigmaI &SigIm, const SigmaI &SigIp);
	
	size_t nMis() const{return _misInd.size();};
	size_t nMis(){return _misInd.size();};
};

class BetaGrpFt : public Grp {
protected:
	vector<vector<double> > _fittedEach; // each member of the outer vector stores the element-specific fitted matrix as vector in the row-major format, to be accessed as a matrix_view of an array
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
	
	const gsl_matrix *fMat() const{return _fittedAll; };
	void save(const SigmaI &SigI);
	void dump();
	
	double lnOddsRat(const Grp &y, const SigmaI &SigI, const size_t i) const;
	
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

class BetaGrpPC : public BetaGrpFt {
protected:
	
public:
	BetaGrpPC() : BetaGrpFt() {};
	BetaGrpPC(const Grp &rsp, const string &predFlNam, const string &evFlNam, const size_t &Npred, const int &nThr);
	BetaGrpPC(const Grp &rsp, const string &predFlNam, const string &evFlNam, const size_t &Npred, RanIndex &up, const int &nThr);
	BetaGrpPC(const Grp &rsp, const string &predFlNam, const string &evFlNam, const size_t &Npred, RanIndex &low, RanIndex &up, const int &nThr);
	BetaGrpPC(const Grp &rsp, const string &predFlNam, const string &evFlNam, const size_t &Npred, const string &outFlNam, const int &nThr);
	BetaGrpPC(const Grp &rsp, const string &predFlNam, const string &evFlNam, const size_t &Npred, RanIndex &up, const string &outFlNam, const int &nThr);
	BetaGrpPC(const Grp &rsp, const string &predFlNam, const string &evFlNam, const size_t &Npred, RanIndex &low, RanIndex &up, const string &outFlNam, const int &nThr);
	
	~BetaGrpPC() {};
	
	BetaGrpPC(const BetaGrpPC &mG); // copy constructor
	BetaGrpPC & operator=(const BetaGrpPC &mG);
	
	//update methods taken care of by BetaGrpFt because the constructor makes pointers to individual MVnormBetaPC
};

class BetaGrpSnp : public MuGrp {
protected:
	gsl_matrix *_Xmat;
	gsl_matrix *_fakeFmat;
	
	int _nThr;
	double _numSaves;
public:
	BetaGrpSnp();
	BetaGrpSnp(const string &predFlNam, const string &outFlNam, const size_t &Ndat, const size_t &Npred, const size_t &d, const int &Nthr);
	BetaGrpSnp(const string &predFlNam, const string &outFlNam, RanIndex &low, const size_t &Npred, const size_t &d, const int &Nthr);
	
	~BetaGrpSnp();
	
	void dump();
	
	BetaGrpSnp(const BetaGrpSnp &mG); // copy constructor
	BetaGrpSnp & operator=(const BetaGrpSnp &mG);
	
	const gsl_matrix *fMat() const{return _fakeFmat; };

	void update(const Grp &dat, const SigmaI &SigIm);
};

class BetaGrpSnpMiss : public MuGrp {
protected:
	vector<vector<double> > _Xmat; // a rugged array that only stores present genotypes
	gsl_matrix *_fakeFmat;
	
	int _nThr;
	double _numSaves;
public:
	BetaGrpSnpMiss() : _numSaves(0.0), _nThr(1), MuGrp() {};
	BetaGrpSnpMiss(const string &predFlNam, const string &outFlNam, const size_t &Ndat, const size_t &Npred, const size_t &d, const int &Nthr, const double &absLab);
	BetaGrpSnpMiss(const string &predFlNam, const string &outFlNam, RanIndex &low, const size_t &Npred, const size_t &d, const int &Nthr, const double &absLab);
	
	~BetaGrpSnpMiss();
	
	BetaGrpSnpMiss(const BetaGrpSnpMiss &mG); // copy constructor
	BetaGrpSnpMiss & operator=(const BetaGrpSnpMiss &mG);
	
	void dump();
	
	const gsl_matrix *fMat() const{return _fakeFmat; };
	
	void update(const Grp &dat, const SigmaI &SigIm);
	
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

class MuBlk : public MuGrp {
protected:
	vector< size_t > _blkStart;                  // vector of indexes into beginnings of each block
	vector< vector< vector<size_t> > > _blkLow;  // lowLevel indexes, separate for each block
	gsl_matrix *_expandedVM;                     // the _valuMat exapnded according to _blkLow
	
	void _updateExp();
public:
	MuBlk() : MuGrp() {};
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
	void save(const char *how = "a");
	// save function, taking file name, appending by default
	void save(const string &fileNam, const char *how = "a");
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

	~SigmaIblk(){};
	
	// Gaussian
	void update(const Grp &dat);
	void update(const Grp &dat, const Grp &mu);
	
};

class StTq {	// the weigting parameter for Student-t sampling
private:
	double _q;
	double _nu; // Student-t degrees of freedom
	size_t _locInd;
public:
	// constructors
	StTq() : _nu(3.0), _q(1.0), _locInd(0) {};
	StTq(const double &nu) : _nu(nu), _q(1.0), _locInd(0) {};
	StTq(const double &nu, const int &d, gsl_rng *r);
	StTq(const double &q, const double &nu, const int &d, gsl_rng *r);
	StTq(const double &nu, const size_t &ind) : _nu(nu), _q(1.0), _locInd(ind) {};
	StTq(const double &q, const double &nu, const size_t &ind, const int &d, gsl_rng *r);
	
	~StTq() {};
	
	double getVal() const {return _q; };
	double getNu() const {return _nu; };
	void   setNu(const double nu){_nu = nu; };
	
	void update(const Grp &dat, const Grp &mu, const SigmaI &SigI, const gsl_rng *r);
	void update(const Grp &dat, const SigmaI &SigI, const gsl_rng *r);
	
	// save function, taking file name, appending by default
	void save(const string &fileNam, const char *how = "a");
	// save function, taking file stream name
	void save(FILE *fileStr);
};

class Qgrp {
private:
	vector<StTq> _qVec;
	vector<size_t> _presInd; // vector of indexes of samples with no missing phenotypes.  These are the ones we can update q for
	gsl_rng *_r;
	
public:
	Qgrp() {};
	Qgrp(const size_t &N);
	Qgrp(const size_t &N, const double &nu);
	Qgrp(const size_t &N, const double &nu, const string &misVecFlNam);
	
	~Qgrp() {};
	
	double operator[](const size_t i) const{return _qVec[i].getVal(); };
	double operator[](const size_t i) {return _qVec[i].getVal(); };
	
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
