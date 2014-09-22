//
//  myEMMAX.cpp
//  MuGenLib
//
//  Created by ajg67 on 9/18/14.
//  Copyright (c) 2014 SEELE. All rights reserved.
//

/*
 *	to compile on the Mezey servers:
 *  g++ myEMMAX.cpp -o myEMMAX -lz -Wl,--start-group /opt/intel/composer_xe_2013.3.163/mkl/lib/intel64/libmkl_intel_lp64.a /opt/intel/composer_xe_2013.3.163/mkl/lib/intel64/libmkl_core.a /opt/intel/composer_xe_2013.3.163/mkl/lib/intel64/libmkl_sequential.a /home/local/CORNELL/ajg67/gslLibs/lib/libgsl.a /opt/intel/composer_xe_2013.3.163/mkl/lib/intel64/libmkl_scalapack_ilp64.a -Wl,--end-group -lpthread -lm -ldl -fopenmp -DHAVE_INLINE -DGSL_RANGE_CHECK_OFF -O3 -march=native
 */
#include <omp.h>
#include <R_ext/Lapack.h> // have to make sure that there is a link to R_ext and Rconfig.h in the appropriate include dir
#include <cstdio>
#include <cstddef>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_statistics.h>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_min.h>

using std::cout;
using std::cerr;
using std::endl;
using std::flush;
using std::ofstream;
using std::stringstream;
using std::string;
using std::vector;
using std::advance;

void colCenter(gsl_matrix *inplace);
double fR(double delta, void *params);

struct twoVec {
	gsl_vector *xi;
	gsl_vector *etaSq;
};

void colCenter(gsl_matrix *inplace){
	gsl_vector *cl  = gsl_vector_alloc(inplace->size1);
	double clMn;
	for (size_t jCl = 0; jCl < inplace->size2; jCl++) {
		gsl_vector_view vCl = gsl_matrix_column(inplace, jCl);
		gsl_matrix_get_col(cl, inplace, jCl);
		clMn = gsl_stats_mean(cl->data, 1, cl->size);
		gsl_vector_add_constant(&vCl.vector, -clMn);
	}
	gsl_vector_free(cl);
}

double fR(double delta, void *params){
	struct twoVec *p = (struct twoVec *) params;
	double frac  = 0.0;
	double xiDel = 0.0;
	for (size_t iVec = 0; iVec < (p->xi)->size; iVec++) {
		xiDel += log(gsl_vector_get(p->xi, iVec) + delta);
		frac  += gsl_vector_get(p->etaSq, iVec)/(gsl_vector_get(p->xi, iVec) + delta);
	}
	return static_cast<double>((p->xi)->size)*log(frac) + xiDel;
}

int main(int argc, char *argv[]){
	const size_t N    = 4500;
	const size_t Nln  = 750;
	const size_t d    = 10;
	const size_t Nsnp = 669958;
	
	string KflNam("MESA_K.gbin");
	string ZflNam("simZ.gbin");
	string YflNam("Yout");
	string XflNam("MESAsnp.gbin");
	string pValFlNam("pVal");
	
	string sim;
	string her;
	int nThr = 4;
	
	// command line processing
	bool sOn = false;
	bool hOn = false;
	bool cOn = false;
	
	for (int iArg = 1; iArg < argc; iArg++) {
		char *pchar = argv[iArg];
		switch (pchar[0]) {
			case '-':{
				if (!pchar[1]) {
					cerr << "ERROR: forgot character after flag" << endl;
					exit(-1);
				}
				switch (pchar[1]) {
					case 's':
						sOn = true;
						break;
					case 'h':
						hOn = true;
						break;
					case 'c':
						cOn = true;
						break;
					default:
						cerr << "ERROR: unrecognized flag " << pchar[1] << endl;
						exit(-1);
						break;
				}
			}
			break;
				
			default:{
				if (sOn) {
					sOn = false;
					sim = pchar;
				}
				else if (hOn){
					hOn = false;
					her = pchar;
				}
				else {
					cOn = false;
					nThr = atoi(pchar);
				}
			}
			break;
			}
	}
	YflNam = YflNam + her + sim + ".gbin";
	pValFlNam = pValFlNam + her + sim + ".gbin";
	
	gsl_matrix *K    = gsl_matrix_alloc(Nln, Nln);
	gsl_matrix *Z    = gsl_matrix_alloc(N, Nln);
	gsl_matrix *Y    = gsl_matrix_alloc(N, d);
	gsl_matrix *Hinv = gsl_matrix_alloc(N, N);
	gsl_matrix *Ut   = gsl_matrix_calloc(N, N);
	gsl_vector *xi   = gsl_vector_calloc(N);  // xi of Kang et al
	gsl_vector *eta  = gsl_vector_alloc(N);
	
	
	FILE *Kin = fopen(KflNam.c_str(), "r");
	gsl_matrix_fread(Kin, K);
	fclose(Kin);

	FILE *Zin = fopen(ZflNam.c_str(), "r");
	gsl_matrix_fread(Zin, Z);
	fclose(Zin);

	FILE *Yin = fopen(YflNam.c_str(), "r");
	gsl_matrix_fread(Yin, Y);
	fclose(Yin);
	colCenter(Y); // so as not to worry about the intercept
	
	cout << "Read in data; phenotype file: " << YflNam << endl;
	gsl_vector_view Kdiag = gsl_matrix_diagonal(K);
	gsl_vector_add_constant(&Kdiag.vector, 1e-6);
	gsl_linalg_cholesky_decomp(K);
	gsl_blas_dtrmm(CblasRight, CblasLower, CblasNoTrans, CblasNonUnit, 1.0, K, Z);
	gsl_matrix_free(K);
	
	// everything has to be initiated to 0.0, otherwise SVD fails
	const int Nvt = 1;
	vector<double>vt(1, 0.0);
	int resSVD = 0;
	int m = static_cast<int>(N);
	int n = static_cast<int>(Nln);
	int Nw = 707250; // got this from a pre-run as recommended in the LAPACK manual
	vector<double>workArr(Nw, 0.0);
	
	gsl_matrix *A = gsl_matrix_alloc(Nln, N);
	gsl_matrix_transpose_memcpy(A, Z);
	dgesvd_("A", "N", &m, &n, A->data, &m, xi->data, Ut->data, &m, vt.data(), &Nvt, workArr.data(), &Nw, &resSVD);
	
	gsl_matrix_free(A);
	gsl_matrix_free(Z);

	for (size_t iXi = 0; iXi < Nln; iXi++) {
		double xiSq = gsl_pow_2(gsl_vector_get(xi, iXi));
		gsl_vector_set(xi, iXi, xiSq);
	}
	cout << "Pre-calculated the singular value matrix and vector" << endl;
	
	gsl_matrix *X = gsl_matrix_alloc(Nln, Nsnp);
	
	FILE *Xin = fopen(XflNam.c_str(), "r");
	gsl_matrix_fread(Xin, X);
	fclose(Xin);
	colCenter(X); // so as not to worry about the intercept
	cout << "Read the SNP matrix" << endl;
	
	gsl_matrix *pVal = gsl_matrix_calloc(Nsnp, d);
	for (size_t Yj = 0; Yj < d; Yj++) {
		gsl_vector_view y = gsl_matrix_column(Y, Yj);
		gsl_blas_dgemv(CblasNoTrans, 1.0, Ut, &y.vector, 0.0, eta);
		
		for (size_t iEta = 0; iEta < eta->size; iEta++) {
			double etaSq = gsl_pow_2(gsl_vector_get(eta, iEta));
			gsl_vector_set(eta, iEta, etaSq);
		}
		
		// minimizer set-up to find delta
		int status  = 0;
		int iter    = 0;
		int maxIter = 1e4;
		double deltaOpt = 1.0;
		double lowBound = 1e-9;
		double upBound  = 1e9;
		
		const gsl_min_fminimizer_type *T = gsl_min_fminimizer_brent;
		gsl_min_fminimizer *s = gsl_min_fminimizer_alloc(T);
		struct twoVec etaXi = {xi, eta};
		
		gsl_function reml;
		reml.function = &fR;
		reml.params   = &etaXi;
		gsl_min_fminimizer_set(s, &reml, deltaOpt, lowBound, upBound);
		while (iter < maxIter) {
			iter++;
			
			status   = gsl_min_fminimizer_iterate(s);
			deltaOpt = gsl_min_fminimizer_x_minimum(s);
			lowBound = gsl_min_fminimizer_x_lower(s);
			upBound  = gsl_min_fminimizer_x_upper(s);
			
			status = gsl_min_test_interval(lowBound, upBound, 1e-6, 0.0);
			if (status == GSL_SUCCESS) {
				break;
			}
			
		}
		gsl_min_fminimizer_free(s);
		cout << "Trait " << Yj + 1 << ": solved the random effect" << endl;
		
		gsl_matrix *locUt = gsl_matrix_alloc(N, N);
		gsl_matrix_memcpy(locUt, Ut);
		for (size_t iRw = 0; iRw < Ut->size1; iRw++) {
			gsl_vector_view row = gsl_matrix_row(locUt, iRw);
			gsl_vector_scale(&row.vector, 1.0/(sqrt(gsl_vector_get(xi, iRw) + deltaOpt)));
		}
		
		gsl_blas_dsyrk(CblasLower, CblasTrans, 1.0, locUt, 0.0, Hinv);
		gsl_blas_dsymv(CblasLower, 1.0, Hinv, &y.vector, 0.0, eta); // this will be re-used for each SNP
		gsl_matrix_free(locUt);
		
		cout << "Trait " << Yj + 1 << ": doing GWA" << endl;

#pragma omp parallel num_threads(nThr)
		{
			gsl_vector *XcolLong = gsl_vector_alloc(N);
			gsl_vector *rsd      = gsl_vector_alloc(N);
			gsl_vector *tmp      = gsl_vector_alloc(N);
			double w;
			double bet;
			double var;
			const double v2 = static_cast<double>(N - 1);
#pragma omp for
			for (size_t Xj = 0; Xj < Nsnp; Xj++) {
				gsl_vector_view Xcol = gsl_matrix_column(X, Xj);
				for (size_t iRep = 0; iRep < N/Nln; iRep++) {
					gsl_vector_view subLong = gsl_vector_subvector(XcolLong, iRep*Nln, Nln);
					gsl_vector_memcpy(&subLong.vector, &Xcol.vector);
				}
				gsl_blas_dsymv(CblasLower, 1.0, Hinv, XcolLong, 0.0, tmp);
				gsl_blas_ddot(XcolLong, tmp, &w);
				
				gsl_blas_ddot(XcolLong, eta, &bet);
				bet = bet/w;
				gsl_vector_memcpy(rsd, &y.vector);
				gsl_blas_daxpy(-bet, XcolLong, rsd);
				
				gsl_blas_dsymv(CblasLower, 1.0, Hinv, rsd, 0.0, tmp);
				gsl_blas_ddot(rsd, tmp, &var);
				var = var/(w*v2);
				
				double Fstat = v2/(v2 + gsl_pow_2(bet)/var);
				double pV = gsl_cdf_beta_P(Fstat, v2/2.0, 0.5);
				gsl_matrix_set(pVal, Xj, Yj, -log10(pV));
			}
			
			gsl_vector_free(tmp);
			gsl_vector_free(rsd);
			gsl_vector_free(XcolLong);
		}
		cout << "Trait " << Yj + 1 << ": GWA finished" << endl;
	}
	FILE *PVout = fopen(pValFlNam.c_str(), "w");
	gsl_matrix_fwrite(PVout, pVal);
	fclose(PVout);
	
	gsl_matrix_free(pVal);
	gsl_matrix_free(Ut);
	gsl_matrix_free(Y);
	gsl_vector_free(xi);
	gsl_vector_free(eta);
	gsl_matrix_free(X);
	gsl_matrix_free(Hinv);
}
