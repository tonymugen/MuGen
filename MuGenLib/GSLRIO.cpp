/*
 *  GSLRIO.cpp
 *  libMuGen
 *
 *  Copyright (c) 2015 Anthony J. Greenberg
 *
 *  This file is part of the MuGen library.
 *
 *  MuGen is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  MuGen is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with MuGen.  If not, see <http://www.gnu.org/licenses/>.
 *
 *  A set of functions to read and save GSL matrices in binary from R
 *  to compile: R CMD SHLIB GSLRIO.cpp -lgsl -O3 -DHAVE_INLINE -DGSL_RANGE_CHECK_OFF
 */

#include <R.h>
#include <string>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>

using namespace std;

/*
	note that to pass a double-quoted chracter to .C(), the **fileNam syntax is essential
*/

extern "C"{
	void GSLmatSave(const char **fileNam, const double *vec, const int *nRows, const int *nCols){
		gsl_matrix *mat = gsl_matrix_alloc(*nRows, *nCols);
		
		for (int iRow = 0; iRow < *nRows; iRow++){  // because R matrices are column-major
			for (int jCl = 0; jCl < *nCols; jCl++){
				gsl_matrix_set(mat, iRow, jCl, vec[iRow + (*nRows)*jCl]);
			}
		}
		
		FILE *fl = fopen(*fileNam, "w");
		gsl_matrix_fwrite(fl, mat);
		fclose(fl);
		
		gsl_matrix_free(mat);
	}
	
}
/*
	for saving vectors used as indeces
*/

extern "C"{
	void GSLvecSaveInt(const char **fileNam, const int *vec, const int *len){
		gsl_vector_int *svVec = gsl_vector_int_alloc(*len);
		
		for (int i = 0; i < *len; i++){
			gsl_vector_int_set(svVec, i, vec[i]);
		}
		
		FILE *fl = fopen(*fileNam, "w");
		gsl_vector_int_fwrite(fl, svVec);
		fclose(fl);
		
		gsl_vector_int_free(svVec);
	}
	
}
 
/*
	for saving double vectors
*/
extern "C"{
	void GSLvecSave(const char **fileNam, const double *vec, const int *len){
		gsl_vector *svVec = gsl_vector_alloc(*len);
		
		for (int i = 0; i < *len; i++){
			gsl_vector_set(svVec, i, vec[i]);
		}
		
		FILE *fl = fopen(*fileNam, "w");
		gsl_vector_fwrite(fl, svVec);
		fclose(fl);
		
		gsl_vector_free(svVec);
	}
	
}

/*
	for saving double vectors by appending to an existing file
*/
extern "C"{
	void GSLvecAppend(const char **fileNam, const double *vec, const int *len){
		gsl_vector *svVec = gsl_vector_alloc(*len);
		
		for (int i = 0; i < *len; i++){
			gsl_vector_set(svVec, i, vec[i]);
		}
		
		FILE *fl = fopen(*fileNam, "a");
		gsl_vector_fwrite(fl, svVec);
		fclose(fl);
		
		gsl_vector_free(svVec);
	}
	
}
/*
	loading int vectors
 */
extern "C" {
	void GSLiVecLoad(const char **fileNam, const int *len, int *arr){
		
		gsl_vector_int *vec = gsl_vector_int_alloc(*len);
		
		FILE *fl = fopen(*fileNam, "r");
		gsl_vector_int_fread(fl, vec);
		fclose(fl);
		
		for (int el = 0; el < *len; el++){
			arr[el] = gsl_vector_int_get(vec, el);
		}
		
		gsl_vector_int_free(vec);
	}
	
}

/*
	the matrix will be passed to R as a vector anyway, so I am reading the matrix in as a vector
	the matrix is stored by row, so make sure to use matrix(..., byrow = T) in the R code
*/

extern "C" {
	void GSLmatLoad(const char **fileNam, const int *nRows, const int *nCols, double *arr){
		int len = (*nRows)*(*nCols);
		
		gsl_vector *vec = gsl_vector_alloc(len);
		
		FILE *fl = fopen(*fileNam, "r");
		gsl_vector_fread(fl, vec);
		fclose(fl);
		
		for (int el = 0; el < len; el++){
			arr[el] = gsl_vector_get(vec, el);
		}
		
		gsl_vector_free(vec);
	}
	
}

