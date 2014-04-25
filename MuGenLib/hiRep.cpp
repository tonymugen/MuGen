//
//  hiRep.cpp
//  MuGenLib
//
//  Created by ajg67 on 1/7/14.
//  Copyright (c) 2014 SEELE. All rights reserved.
//

/*
	analysis of Diane's data
 */

#include <omp.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <vector>
#include <string>
#include <sstream>
#include <iostream>
#include <ctime>
#include <cmath>
#include <cctype>
#include <cstdio>

#include "MuGenLib.h"

using std::vector;
using std::string;
using std::cout;
using std::cerr;
using std::endl;
using std::flush;
using std::stringstream;
using std::remove;

int main(int argc, char *argv[]){
	
	bool sOn = false;
	bool bOn = false;
	bool tOn = false;
	bool nOn = false;
	bool NOn = false;
	bool cOn = false;
	bool fOn = false;
	bool iOn = false;
	
	const size_t d   = 9;
	const size_t Nln = 36;
	
	double nuDat = 3.0;
	int Nbnin    = 10;
	int Nsamp    = 10;
	int Nthin    = 1;
	size_t Ntot  = 694;
	string cNum;
	string fNum;
	string iNum;
	
	string datFlNam("hiRepData_");
	string naMatFlNam("hiRepMatNAind_");
	string naVecFlNam("hiRepTotNAind_");
	string lnIndFlNam("hiRepLnInd_");
	string outFlNam("hiRepMn_");
	
	for (int iArg = 1; iArg < argc; iArg++) {
		char *pchar = argv[iArg];
		switch (pchar[0]) {
			case '-':{
				if (!pchar[1]) {
					cerr << "ERROR: forgot character after flag" << endl;
					exit(0);
				}
				else
				if (pchar[2]) {
					cerr << "ERROR: unrecognized flag of more than one character" << endl;
					exit(0);
				}
				
				switch (pchar[1]) {
					case 's':
						sOn = true;
						break;
						
					case 'b':
						bOn = true;
						break;
						
					case 't':
						tOn = true;
						break;
						
					case 'n':
						nOn = true;
						break;
						
					case 'N':
						NOn = true;
						break;
						
					case 'f':
						fOn = true;
						break;
						
					case 'i':
						iOn = true;
						break;
						
					case 'c':
						cOn = true;
						break;
						
					default:
						cerr << "ERROR: unrecognized flag " << pchar[1] << endl;
						exit(0);
				}
			}
				break;
				
			default:{
				if (sOn) {
					sOn   = false;
					Nsamp = atoi(pchar);
				}
				else
				if (bOn) {
					bOn   = false;
					Nbnin = atoi(pchar);
				}
				else
				if (tOn) {
					tOn   = false;
					Nthin = atoi(pchar);
				}
				else
				if (NOn) {
					NOn   = false;
					Ntot = atoi(pchar);
				}
				else
				if (cOn) {
					cOn = false;
					cNum = pchar;
				}
				else
				if (fOn) {
					fOn = false;
					fNum = pchar;
				}
				else
				if (iOn) {
					iOn = false;
					iNum = pchar;
				}
				else
				{
					nOn = false;
					nuDat = atof(pchar);
				}
			}
				break;
		}
	}
	
	stringstream nuOSn;
	nuOSn << nuDat;
	string nString = nuOSn.str();
	size_t perPs = nString.find('.');
	if (perPs != string::npos) {
		nString.erase(perPs, 1);
	}
	
	datFlNam   = datFlNam   + fNum + "_" + iNum + ".gbin";
	naMatFlNam = naMatFlNam + fNum + "_" + iNum + ".gbin";
	naVecFlNam = naVecFlNam + fNum + "_" + iNum + ".gbin";
	lnIndFlNam = lnIndFlNam + fNum + "_" + iNum + ".gbin";
	outFlNam   = outFlNam   + fNum + "_" + iNum + "_" + nString + "_" + cNum + ".gbin";
		
	//cout << "Index file: " << lnIndFlNam << "; Ntot = " << Ntot << endl;
	RanIndex dat2ln(Ntot, Nln, lnIndFlNam);
	RanIndex ln2mu(Nln);
	RanIndex mu2pr;
	
	//cout << "Data file: " << datFlNam << "; naMat: " << naMatFlNam << "; naVec: " << naVecFlNam << endl;
	MuGrpMiss dataI(datFlNam, naMatFlNam, naVecFlNam, dat2ln, d);
	Grp &data = dataI;
	
	MuGrp muLnI(data, dat2ln, ln2mu);
	Grp &muLn = muLnI;
	
	MuGrp muI(muLn, ln2mu, mu2pr);
	Grp &mu = muI;
	
	MuGrp datDevI = data - muLn;
	Grp &datDev   = datDevI;
	
	MuGrp lnDevI = muLn - mu;
	Grp &lnDev   = lnDevI;

	SigmaI SigIe(datDev, 1.0, 2.0);
	SigmaI SigIln(lnDev, 1.0, 2.0);
	SigmaI SigIpr(d, 1e-6);
	
	Qgrp qDat(Ntot, nuDat, naVecFlNam);
	
	for (int iBn = 0; iBn < Nbnin; iBn++) {
		data.update(muLn, SigIe);
		muLn.update(data, qDat, SigIe, mu, SigIln);
		mu.update(muLn, SigIln, SigIpr);
		
		datDevI = muLn - data;
		lnDevI  = muLn - mu;
		
		SigIe.update(datDev, qDat);
		qDat.update(datDev, SigIe);
		SigIln.update(lnDev);
		
		cout << "+" << flush;
		
	}
	cout << endl;
	
	gsl_matrix *res = gsl_matrix_calloc(lnDev.dMat()->size1, lnDev.dMat()->size2);
	
	for (int iSm = 0; iSm < Nsamp; iSm++) {
		data.update(muLn, SigIe);
		muLn.update(data, qDat, SigIe, mu, SigIln);
		mu.update(muLn, SigIln, SigIpr);
		
		datDevI = muLn - data;
		lnDevI  = muLn - mu;
		
		SigIe.update(datDev, qDat);
		qDat.update(datDev, SigIe);
		SigIln.update(lnDev);
		
		if ((iSm + 1) % Nthin) {
			cout << "." << flush;
		}
		else {
			gsl_matrix_add(res, lnDev.dMat());
			cout << "|" << flush;
		}
		
	}
	cout << endl;
	gsl_matrix_scale(res, static_cast<double>(Nthin)/static_cast<double>(Nsamp));
	
	FILE *resOut = fopen(outFlNam.c_str(), "w");
	gsl_matrix_fwrite(resOut, res);
	fclose(resOut);
	
	gsl_matrix_free(res);
}

