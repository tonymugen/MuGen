//
//  tstFTR.cpp
//  MuGenLib
//
//  Created by ajg67 on 3/18/14.
//  Copyright (c) 2014 SEELE. All rights reserved.
//

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

#include <MuGen.h>

using std::vector;
using std::string;
using std::cout;
using std::cerr;
using std::endl;
using std::flush;
using std::stringstream;
using std::remove;

int main(int argc, char *argv[]){
	const size_t N  = 750;
	const size_t d  = 10;
	const int Nb    = atoi(argv[1]);
	const int Ns    = atoi(argv[2]);
	const int Nt    = atoi(argv[3]);
	const size_t Nx = 50 + atoi(argv[4]);
	const double nu = atof(argv[5]);
	
	string inXfnam("nullTrX/nullTrX");
	inXfnam = inXfnam + argv[4] + ".gbin";
	string inYfnam("nullTrY/Ymn");
	inYfnam = inYfnam + argv[6] + ".gbin";
	string outBnam("nullTrB/betRes");
	outBnam = outBnam + argv[6] + "_" + argv[4] + "_" + argv[5] + "_" + argv[7] + ".gbin";
	//string outSnam("nullTrB/SigBtRes");
	//outSnam = outSnam + argv[6] + "_" + argv[4] + "_" + argv[5] + "_" + argv[7] + ".gbin";
	
	RanIndex ln2mu = RanIndex(N);
	RanIndex b2pr  = RanIndex(Nx);
	RanIndex mu2pr = RanIndex();
	
	MuGrp yI = MuGrp(inYfnam, ln2mu, d);
	Grp &y   = yI;
	
	MuGrp muI = MuGrp(y, ln2mu, mu2pr);
	Grp &mu   = muI;
	
	MuGrp bPredI = y - mu;
	Grp &bPred   = bPredI;
	
	BetaGrpFt betI = BetaGrpFt(bPred, inXfnam, Nx, b2pr, outBnam, 4);
	Grp &bet       = betI;
	
	MuGrp mPrdI = y - bet;
	Grp &mPrd   = mPrdI;
	
	MuGrp rsdI = y - mu - bet;
	Grp &rsd   = rsdI;
	
	SigmaI SigIe(rsd, 1.0, 2.0);
	SigmaI SigIb(d, 1e-2, 750.0);
	//SigmaI SigIb(d, 1e-2, 2.0, outSnam);
	SigmaI SigIpr(d, 0.000001);
	
	Qgrp qB(Nx, nu);
	
	for (int iBn = 0; iBn < Nb; iBn++) {
		mPrdI = y - bet;
		mu.update(mPrd, SigIe, SigIpr);
		bPredI = y - mu;
		bet.update(bPred, SigIe, qB, SigIb);
		
		rsdI = y - mu - bet;
		SigIe.update(rsd);
		SigIb.update(bet, qB);
		qB.update(bet, SigIb);
		cout << "*" << flush;
	}
	cout << endl;
	for (int iSm = 0; iSm < Ns; iSm++) {
		mPrdI = y - bet;
		mu.update(mPrd, SigIe, SigIpr);
		bPredI = y - mu;
		bet.update(bPred, SigIe, qB, SigIb);
		
		rsdI = y - mu - bet;
		SigIe.update(rsd);
		SigIb.update(bet, qB);
		qB.update(bet, SigIb);
		
		if ((iSm + 1) % Nt) {
			cout << "." << flush;
		}
		else {
			bet.save();
			cout << "|" << flush;
		}
	}
	cout << endl;
	
}



