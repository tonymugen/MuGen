/*
*  MVstTHBgwa.cpp
*  MuGen
*
*  Created by ajg67 on 11/8/12.
*   Copyright (c) 2012 SEELE. All rights reserved.
*
*
*
*/
#include <vector>
#include <string>
#include <sstream>
#include <iostream>
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

void finishFlNam(string &, const vector<string> &);
string db2str(const double &);

void finishFlNam(string &start, const vector<string> &add){
	for (vector<string>::const_iterator stIt = add.begin(); stIt != add.end(); ++stIt) {
		start += "_";
		start += *stIt;
	}
	start += ".gbin";
}

string db2str(const double &val){
	stringstream nuOS;
	nuOS << val;
	string vString = nuOS.str();
	size_t perPs = vString.find('.');
	if (perPs != string::npos) {
		vString.replace(perPs, 1, "-");
	}
	return vString;
}

int main(int argc, char *argv[]){
	/*
	 *	Parsing command line flags
	 */
	bool sOn = false;
	bool bOn = false;
	bool tOn = false;
	bool nOn = false;
	bool gOn = false;
	bool pOn = false;
	bool lOn = false;
	bool fOn = false;
	bool cOn = false;
	bool COn = false;
	bool SOn = false;
	bool mOn = false;
	bool MOn = false;
	bool rOn = false;
	bool dOn = false;
	bool aOn = false;
	
	// intitialize the flag variables with default values
	double nuE   = 3.0;
	double nuG   = 3.0;
	double nuB   = 3.0;
	double ldCt  = 0.75;
	double prABF = 0.0;
	int Nbnin   = 10;
	int Nsamp   = 10;
	int Nthin   = 1;
	size_t Nsnp = 10000;
	double Nmul = 5.0;
	int nThr    = 4;
	size_t d    = 10;
	string cNum;
	bool miss  = false;
	string model("SM"); // default is single marker; other possibilities: VS == variable selection; RS == rank selection; flag -r controls this
	
	string dfNam("Y.gbin");
	string mMatFnam("simMatNAind.gbin");
	string mVecFnam("simTotNAind.gbin");
	string BVout("BV.gbin");
	string SgEout("Sig_e.gbin");
	string SgRout("Sig_rep.gbin");
	string SgAout("Sig_a.gbin");
	string betOut("betSNP.gbin");
	string pepOut("PEP.gbin");
	string snpIn("MESAsnp.gbin");
	string missF("");
	
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
						
					case 'M':
						MOn = true;
						break;
						
					case 'n':
						nOn = true;
						break;
						
					case 'g':
						gOn = true;
						break;
						
					case 'p':
						pOn = true;
						break;
						
					case 'l':
						lOn = true;
						break;
						
					case 'f':
						fOn = true;
						break;
						
					case 'c':
						cOn = true;
						break;
						
					case 'C':
						COn = true;
						break;
						
					case 'S':
						SOn = true;
						break;
					
					case 'm':
						mOn  = true;
						miss = true;
						break;
						
					case 'd':
						dOn = true;
						break;
						
					case 'a':
						aOn = true;
						break;
						
					case 'r':
						rOn  = true;
						break;
						
					case 'h':
						cerr << "usage:\n"
						<< "MVstTHBgwaNS [-s MCMC_sample_number] [-b burnin_length] [-t thinning_ratio] [-M multiple_of_sample#_to_keep] [-c chain_ID] [-C thread_number] [-m missing_indicator_file_tag] [-d phenotype_number] [-n Student_t_nu_for_reps] [-g Student_t_nu_for_PCs] [-S number_of_snps] [-f response_data_file_ID]" << endl;
						exit(0);
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
					if (MOn) {
						MOn  = false;
						Nmul = atof(pchar);
					}
				else
					if (cOn) {
						cOn  = false;
						cNum = pchar;
					}
				else
					if (COn) {
						COn  = false;
						nThr = atoi(pchar);
					}
				else
					if (SOn) {
						SOn  = false;
						Nsnp = atoi(pchar);
					}
				else
					if (dOn) {
						dOn = false;
						d   = atoi(pchar);
					}
				else
					if (aOn) {
						aOn = false;
						prABF = atof(pchar);
					}
				else
					if (mOn) {
						mOn = false;
						missF = pchar;
						
						mMatFnam = mMatFnam.substr(0,11);
						mMatFnam += pchar;
						mMatFnam += ".gbin";
						
						mVecFnam = mVecFnam.substr(0,11);
						mVecFnam += pchar;
						mVecFnam += ".gbin";
					}
				else
					if (fOn) {
						fOn = false;
						dfNam = dfNam[0];
						dfNam += pchar;
						dfNam += ".gbin";
						
						BVout = BVout.substr(0, 2);
						BVout += pchar;
						
						SgEout = SgEout.substr(0, 5);
						SgEout += pchar;
						
						SgRout = SgRout.substr(0, 7);
						SgRout += pchar;
						
						SgAout = SgAout.substr(0, 5);
						SgAout += pchar;
						
						betOut = betOut.substr(0, 6);
						betOut += pchar;
						
						pepOut = pepOut.substr(0, 3);
						pepOut += pchar;
						
					}
				else
					if (rOn){
						rOn = false;
						model = pchar;
					}
				else
					if (lOn){
						lOn = false;
						ldCt = atof(pchar);
					}
				else
					if (nOn){
						nOn = false;
						nuE = atof(pchar);
				}
				else
					if (gOn){
						gOn = false;
						nuG = atof(pchar);
					}
				else {
					pOn = false;
					nuB = atof(pchar);
				}
			}
				break;
		}
	}
	
	if ((model != "SM") && (model != "SMP") && (model != "VS") && (model != "RS")) {
		cerr << "ERROR: Specified SNP model (" << model << ") not supported." << endl;
		exit(-1);
	}
	
	if (d > 10) {
		cerr << "ERROR: d = " << d << " is greater than 10." << endl;
		exit(-1);
	}
	
	vector<string> addOn;
	addOn.push_back(db2str(nuE));
	addOn.push_back(db2str(nuG));
	addOn.push_back(db2str(nuB));
	
	if (miss) {
		addOn.push_back(missF);
	}
	addOn.push_back(model);
	if ((model == "SM") || (model == "SMP")) {
		if (prABF) {
			addOn.push_back("BF");
		}
		else {
			addOn.push_back("P");
		}
	}
	addOn.push_back(cNum);
	
	finishFlNam(BVout, addOn);
	finishFlNam(SgRout, addOn);
	finishFlNam(SgEout, addOn);
	finishFlNam(SgAout, addOn);
	finishFlNam(betOut, addOn);
	finishFlNam(pepOut, addOn);
	
	const size_t Ntt = 4500;
	const size_t Nrp = 2250;
	const size_t Nln = 750;
	const size_t Npc = 749;
	
	
	/*
	 *	Initialization using constructors
	 */
	
	RanIndex e2rp  = RanIndex(Ntt, Nrp, string("RPind.gbin"));
	RanIndex rp2ln = RanIndex(Nrp, Nln, string("LNind.gbin"));
	RanIndex rp2mu = RanIndex(Nrp);
	RanIndex ln2mu = RanIndex(Nln);
	RanIndex gm2pr = RanIndex(Npc);
	RanIndex mu2pr = RanIndex();
	
	
	Grp *datI;
	if (miss) {
		datI = new MuGrpMiss(dfNam, mMatFnam, mVecFnam, e2rp, d);
	}
	else {
		datI = new MuGrp(dfNam, e2rp, d);
	}
	
	Grp &dat = *datI;
	
	MuGrp muRepI(dat, e2rp, rp2ln);
	Grp &muRep = muRepI;
	
	MuGrp datDevI = dat - muRep;
	Grp &datDev   = datDevI;
	
	Grp *gammaI;
	if (nuG <= 10.0) {
		gammaI = new BetaGrpPC(muRep, string("MESAevec.gbin"), string("MESAeval.gbin"), Npc, rp2ln, gm2pr, nThr);
	}
	else {
		gammaI = new BetaGrpPCpex(muRep, string("MESAevec.gbin"), string("MESAeval.gbin"), Npc, 1e-6, rp2ln, gm2pr, nThr);
	}
	Grp &gamma = *gammaI;
	
	MuGrp muRspI = muRep - gamma;
	Grp &muRsp   = muRspI;
	
	MuGrp muI(muRsp, rp2mu, mu2pr);
	Grp &mu = muI;
	
	MuGrp bvI = mu + gamma;
	Grp &bv   = bvI;
	
	MuGrp repDevI = muRep - mu;
	Grp &repDev   = repDevI;
	
	MuGrp gmRspI = muRep - mu;
	Grp &gmRsp   = gmRspI;
	
	MuGrp gsI = gamma;
	Grp &gs   = gsI;
	
	SigmaI SigIe(datDev, SgEout, 1.0, 2.0);
	SigmaI SigIrep(repDev, SgRout, 1.0, 2.0);
	SigmaI SigIa(gamma, SgAout, 1.0, 2.0);
	SigmaI SigIpr(d, 0.000001);
	
	Qgrp qE;
	
	if (miss) {
		qE = Qgrp(Ntt, nuE, mVecFnam);
	}
	else {
		qE = Qgrp(Ntt, nuE);
	}
	Qgrp qG(Npc, nuG);
	cout << "Initial set-up done.  Starting SNP initialization..." << endl;
	
	for (int iPre = 0; iPre < floor(0.6*Nbnin); iPre++) {
		if (miss) {
			dat.update(muRep, SigIe);
		}
		
		muRep.update(dat, qE, SigIe, bv, SigIrep);
		datDevI = dat - muRep;
		
		SigIe.update(datDev, qE);
		qE.update(datDev, SigIe);
		
		bvI     = mu + gamma;
		repDevI = muRep - bv;
		
		SigIrep.update(repDev);
		
		gmRspI = muRep - mu;
		gamma.update(gmRsp, SigIrep, qG, SigIa);
		
		muRspI = muRep - gamma;
		mu.update(muRsp, SigIrep, SigIpr);
		
		SigIa.update(gamma, qG);
		qG.update(gamma, SigIa);
		
		cout << "*" << flush;

	}
	cout << endl;
	
	vector<double> initP(2, 0.5);
	vector<double> aPr(2);
	aPr[0] = 1.0;
	aPr[1] = Nmul - 1;
	MixP kappa(initP, aPr);
	
	RanIndexVS snp2prI(Nsnp, pepOut, kappa);
	RanIndex &snp2pr = snp2prI;
	
	cout << "Initializing SNP regression for model " << model << endl;
	Grp *snpBetI;
	if (model == "VS") {
		snpBetI = new BetaGrpBVSR(repDev, SigIrep, snpIn, Nmul, ldCt, rp2ln, snp2pr, betOut, nThr);
	}
	else
		if (model == "RS"){
			snpBetI = new BetaGrpFt(repDev, SigIrep, snpIn, Nsnp, Nmul, ldCt, rp2ln, snp2pr, betOut, nThr);
		}
		else if (model == "SMP"){
				if (d == 1) {
					cerr << "WARNING: trying to do partial regression with only one trait." << endl;
				}
				snpBetI = new BetaGrpPSR(snpIn, betOut, rp2ln, Nsnp, d, nThr, prABF); // if prABF is 0.0, it will output p-values
			}
		else {
			snpBetI = new BetaGrpSnp(snpIn, betOut, rp2ln, Nsnp, d, nThr, prABF); // if prABF is 0.0, it will output p-values
		}
	
	Grp &snpBet = *snpBetI;
	
	SigmaI SigIbet(d, 1.0, 20.0);
	Qgrp qB(Nln*Nmul, nuB);
	
	MuGrp snpDevI = repDev;
	Grp &snpDev   = snpDevI;
	
	if ( (model != "SM") && (model != "SMP") ) {
		cout << "Pre-burnin..." << endl;
		for (int iSt = 0; iSt < 100; iSt++) {
			if (miss) {
				dat.update(muRep, SigIe);
			}
			
			muRep.update(dat, qE, SigIe, bv, SigIrep);
			datDevI = dat - muRep;
			
			SigIe.update(datDev, qE);
			qE.update(datDev, SigIe);
			
			gsI     = gamma + snpBet;
			bvI     = mu + gs;
			repDevI = muRep - bv;
			
			SigIrep.update(repDev);
			
			gmRspI = muRep - mu - snpBet;
			gamma.update(gmRsp, SigIrep, qG, SigIa);
			
			snpDevI = muRep - mu - gamma;
			snpBet.update(snpDev, SigIrep, SigIbet);
			kappa.update(snp2pr);
			
			muRspI = muRep - gs;
			mu.update(muRsp, SigIrep, SigIpr);
			
			SigIa.update(gamma, qG);
			qG.update(gamma, SigIa);
			
			cout << "#" << flush;
			
		}
		cout<< endl;

	}
	
	if ( (model == "SM") || (model == "SMP") ) {
		cout << "Single marker model." << endl;
		cout << "Starting burnin..." << endl;
		for (int iBnin = 0; iBnin < floor(0.4*Nbnin); iBnin++) {
			if (miss) {
				dat.update(muRep, SigIe);
			}
			
			muRep.update(dat, qE, SigIe, bv, SigIrep);
			datDevI = dat - muRep;
			
			SigIe.update(datDev, qE);
			qE.update(datDev, SigIe);
			
			bvI     = mu + gamma;
			repDevI = muRep - bv;
			
			SigIrep.update(repDev);
			
			gmRspI = muRep - mu;
			gamma.update(gmRsp, SigIrep, qG, SigIa);
						
			muRspI = muRep - gamma;
			mu.update(muRsp, SigIrep, SigIpr);
			
			SigIa.update(gamma, qG);
			qG.update(gamma, SigIa);
			
			cout << "+" << flush;
		}
		cout << endl;
		cout << "Burn-in finished, starting the sampling..." << endl;
		
		for (int iSam = 0; iSam < Nsamp; iSam++) {
			if (miss) {
				dat.update(muRep, SigIe);
			}
			
			muRep.update(dat, qE, SigIe, bv, SigIrep);
			datDevI = dat - muRep;
			
			SigIe.update(datDev, qE);
			qE.update(datDev, SigIe);
			
			bvI     = mu + gamma;
			repDevI = muRep - bv;
			
			SigIrep.update(repDev);
			
			gmRspI = muRep - mu;
			gamma.update(gmRsp, SigIrep, qG, SigIa);
			
			muRspI = muRep - gamma;
			mu.update(muRsp, SigIrep, SigIpr);
			
			SigIa.update(gamma, qG);
			qG.update(gamma, SigIa);
			
			if ((iSam + 1) % Nthin) {
				cout << "." << flush;
			}
			else {
				bv.save(BVout);
				
				snpBet.update(repDev, SigIrep);
				
				SigIe.save();
				SigIrep.save();
				SigIa.save();
				
				cout << "|" << flush;
			}
		}
		cout << endl;
		snpBet.dump();
	}
	else {
		cout << "Multi-marker model." << endl;
		cout << "Starting burnin..." << endl;
		for (int iBnin = 0; iBnin < floor(0.4*Nbnin); iBnin++) {
			if (miss) {
				dat.update(muRep, SigIe);
			}
			
			muRep.update(dat, qE, SigIe, bv, SigIrep);
			datDevI = dat - muRep;
			
			SigIe.update(datDev, qE);
			qE.update(datDev, SigIe);
			
			gsI     = gamma + snpBet;
			bvI     = mu + gs;
			repDevI = muRep - bv;
			
			SigIrep.update(repDev);
			
			gmRspI = muRep - mu - snpBet;
			gamma.update(gmRsp, SigIrep, qG, SigIa);
			
			snpDevI = muRep - mu - gamma;
			if (model == "VS") {
				snpBet.update(snpDev, SigIrep, SigIbet);
				SigIbet.update(snpBet);
				kappa.update(snp2pr);
			}
			else{
				snpBet.update(snpDev, SigIrep, qB, SigIbet);
				SigIbet.update(snpBet, qB);
				qB.update(snpBet, SigIbet);
			}
			
			muRspI = muRep - gs;
			mu.update(muRsp, SigIrep, SigIpr);
			
			SigIa.update(gamma, qG);
			qG.update(gamma, SigIa);
			
			cout << "+" << flush;
		}
		cout << endl;
		cout << "Burn-in finished, starting the sampling..." << endl;
		
		for (int iSam = 0; iSam < Nsamp; iSam++) {
			if (miss) {
				dat.update(muRep, SigIe);
			}
			
			muRep.update(dat, qE, SigIe, bv, SigIrep);
			datDevI = dat - muRep;
			
			SigIe.update(datDev, qE);
			qE.update(datDev, SigIe);
			
			gsI     = gamma + snpBet;
			bvI     = mu + gs;
			repDevI = muRep - bv;
			
			SigIrep.update(repDev);
			
			gmRspI = muRep - mu - snpBet;
			gamma.update(gmRsp, SigIrep, qG, SigIa);
			
			snpDevI = muRep - mu - gamma;
			if (model == "VS") {
				snpBet.update(snpDev, SigIrep, SigIbet);
				SigIbet.update(snpBet);
				kappa.update(snp2pr);
			}
			else{
				snpBet.update(snpDev, SigIrep, qB, SigIbet);
				SigIbet.update(snpBet, qB);
				qB.update(snpBet, SigIbet);
			}
			
			muRspI = muRep - gs;
			mu.update(muRsp, SigIrep, SigIpr);
			
			SigIa.update(gamma, qG);
			qG.update(gamma, SigIa);
			
			if ((iSam + 1) % Nthin) {
				cout << "." << flush;
			}
			else {
				if (model == "VS") {
					snpBet.save(snpDev, SigIrep);
				}
				else {
						snpBet.save(SigIrep);
						snp2pr.save(snpDev, snpBet, SigIrep);
					}
				
				bv.save(BVout);
				
				SigIe.save();
				SigIrep.save();
				SigIa.save();
				
				cout << "|" << flush;
			}
		}
		cout << endl;
		snpBet.dump();
		snp2pr.dump();
	}
	
	delete snpBetI;
	
}
