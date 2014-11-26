//
//  ionoMV.cpp
//  MuGenLib
//
//  Created by ajg67 on 11/17/14.
//  Copyright (c) 2014 SEELE. All rights reserved.
//

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
		//vString.erase(perPs, 1);
		vString.replace(perPs, 1, "-");
	}
	return vString;
}

int main(int argc, char *argv[]){
	
	// parsing command-line flags
	bool sOn = false; // samplng
	bool bOn = false; // burnin-in
	bool tOn = false; // thinning
	bool nOn = false; // nuRep
	bool gOn = false; // nuG
	bool pOn = false; // pop data set
	bool dOn = false; // trait data set
	bool mOn = false; // SNP model
	bool cOn = false; // chain ID
	bool TOn = false; // # of threads
	
	double nuE   = 3.0;
	double nuG   = 3.0;
	double prABF = 1e3;
	int Nbnin    = 10;
	int Nsamp    = 10;
	int Nthin    = 1;
	size_t Nsnp  = 678742;
	int nThr     = 4;
	
	string snpModel("SM");
	string trtSet;
	char trtTok;
	string popID;
	char popTok;
	string cNam;
	
	string dfNam("Y");
	string sdfNam("SD");
	string eIndFnam("errInd");
	string cvFlNam("PCecov");
	string mMatFnam("MatNAind");
	string mVecFnam("TotNAind");
	string pcEvecFlNam("ionHDRA_EVC");
	string pcEvalFlNam("ionHDRA_EVL");
	string lnFacVec("LnInd");
	string batchFacVec("BatchInd");
	string LNout("LN");
	string BVout("BV");
	string SgEout("Sig_e");
	string betOut("SNP");
	string snpIn("snpHDRA");
	
	for (int iArg = 1; iArg < argc; iArg++) {
		char *pchar = argv[iArg];
		switch (pchar[0]) {
			case '-': {
				if (!pchar[1]) {
					cerr << "ERROR: forgot character after flag" << endl;
					exit(-1);
				}
				else if (pchar[2]) {
						cerr << "ERROR: unrecognized flag of more than one character" << endl;
						exit(-1);
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
						
					case 'g':
						gOn = true;
						break;
						
					case 'p':
						pOn = true;
						break;
						
					case 'd':
						dOn = true;
						break;
						
					case 'm':
						mOn = true;
						break;
						
					case 'c':
						cOn = true;
						break;
						
					case 'T':
						TOn = true;
						break;
						
					default:
						cerr << "ERROR: unrecognized flag " << pchar[1] << endl;
						exit(-1);
						break;
				}
			}
				break;
			
			default: {
				if (sOn) {
					sOn   = false;
					Nsamp = atoi(pchar);
				}
				else if (bOn) {
					bOn   = false;
					Nbnin = atoi(pchar);
				}
				else if (tOn) {
					tOn   = false;
					Nthin = atoi(pchar);
				}
				else if (nOn) {
					nOn = false;
					nuE = atof(pchar);
				}
				else if (gOn) {
					gOn = false;
					nuG = atof(pchar);
				}
				else if (pOn){
					pOn = false;
					switch (pchar[0]) {
						case 'A':
							popID  = "all";
							popTok = 'A';
							break;
						case 'J':
							popID  = "jap";
							popTok = 'J';
							break;
						case 'I':
							popID  = "ind";
							popTok = 'I';
							break;
						default:
							cerr << "ERROR: unknown population token " << pchar[0] << endl;
							exit(-1);
							break;
					}
				}
				else if (dOn){
					dOn = false;
					switch (pchar[0]) {
						case 'R':
							trtSet   = "root";
							trtTok   = 'R';
							cvFlNam += "Root";
							eIndFnam = "phenoData/" + eIndFnam + "R.gbin";
							break;
						case 'S':
							trtSet   = "shoot";
							trtTok   = 'S';
							cvFlNam += "Shoot";
							eIndFnam = "phenoData/" + eIndFnam + "S.gbin";
							break;
						case 'B':
							trtSet   = "all"; // 'B' for 'both' or 'A' for 'all' works
							trtTok   = 'A';
							cvFlNam += "All";
							eIndFnam = "phenoData/" + eIndFnam + "A.gbin";
							break;
						case 'A':
							trtSet   = "all";
							trtTok   = 'A';
							cvFlNam += "All";
							eIndFnam = "phenoData/" + eIndFnam + "A.gbin";
							break;
						default:
							cerr << "ERROR: unkown trait set token: " << pchar[0] << endl;
							exit(-1);
							break;
					}
				}
				else if (mOn){
					mOn = false;
					switch (pchar[0]) {
						case 'S':
							// leave it as "SM"
							break;
						case 'C':
							snpModel += "C";
							break;
						case 'P':
							snpModel += "P";
							break;
						default:
							cerr << "ERROR: unknown SNP model " << pchar[0] << endl;
							exit(-1);
							break;
					}
				}
				else if (cOn){
					cOn = false;
					cNam = pchar;
				}
				else if (TOn){
					TOn = false;
					nThr = atoi(pchar);
				}
				break;
			}
		}
	}
	if (popID == "all") {
		dfNam       = "phenoData/" + dfNam + trtSet + ".gbin";
		sdfNam      = "phenoData/" + sdfNam + trtSet + ".gbin";
		cvFlNam     = "phenoData/" + cvFlNam + ".gbin";
		mMatFnam    = "phenoData/" + trtSet + mMatFnam + ".gbin";
		mVecFnam    = "phenoData/" + trtSet + mVecFnam + ".gbin";
		lnFacVec    = "phenoData/" + trtSet + lnFacVec + ".gbin";
		batchFacVec = "phenoData/" + trtSet + batchFacVec + ".gbin";
	}
	else{
		dfNam       = "phenoData/" + dfNam + trtSet + "_" + popID + ".gbin";
		cvFlNam     = "phenoData/" + cvFlNam + "_" + popID + ".gbin";
		mMatFnam    = "phenoData/" + trtSet + mMatFnam + "_" + popID + ".gbin";
		mVecFnam    = "phenoData/" + trtSet + mVecFnam + "_" + popID + ".gbin";
		lnFacVec    = "phenoData/" + trtSet + lnFacVec + "_" + popID + ".gbin";
		batchFacVec = "phenoData/" + trtSet + batchFacVec + "_" + popID + ".gbin";
	}
	pcEvecFlNam = pcEvecFlNam + popID + ".gbin";
	pcEvalFlNam = pcEvalFlNam + popID + ".gbin";
	snpIn       = snpIn + popID + ".gbin";
	
	LNout  = trtSet + LNout;
	BVout  = trtSet + BVout;
	SgEout = trtSet + SgEout;
	betOut = trtSet + betOut;
	
	vector<string> addOn;
	addOn.push_back(snpModel);
	addOn.push_back(db2str(nuE));
	addOn.push_back(db2str(nuG));
	addOn.push_back(cNam);
	
	finishFlNam(LNout, addOn);
	finishFlNam(BVout, addOn);
	finishFlNam(SgEout, addOn);
	finishFlNam(betOut, addOn);
	
	size_t N;
	size_t Nln;
	size_t Nbatch;
	size_t Ncv = 25;
	size_t d;
	trtSet == "all" ? d = 50 : d = 25;
	
	switch (popTok) {
		case 'I':
			Nln  = 131;
			Nsnp = 586773;
			switch (trtTok) {
				case 'R':
					N      = 495;
					Nbatch = 21;
					break;
					
				case 'S':
					N      = 494;
					Nbatch = 19;
					break;
				
				default: // is all
					N      = 496;
					Nbatch = 21;
					break;
			}
			break;
			
		case 'J':
			Nln  = 232;
			Nsnp = 578377;
			Ncv  = 24;
			switch (trtTok) {
				case 'R':
					N      = 846;
					Nbatch = 22;
					break;
					
				case 'S':
					N      = 863;
					Nbatch = 23;
					break;
					
				default: // is all
					N      = 865;
					Nbatch = 23;
					break;
			}
			break;
		
		default:
			Nln = 373;
			switch (trtTok) {
				case 'R':
					N      = 1377;
					Nbatch = 22;
					break;
					
				case 'S':
					N      = 1393;
					Nbatch = 23;
					break;
					
				default: // is all
					N      = 1397;
					Nbatch = 23;
					break;
			}
			
			break;
	}
	
	const size_t Npc = Nln - 1;
	
	RanIndex e2ln(N, Nln, lnFacVec);
	RanIndex e2e(N, N);
	RanIndex b2pr(Nbatch);
	RanIndex ln2mu(Nln);
	RanIndex ln2bv(Nln, Nln);
	RanIndex gm2pr(Npc);
	RanIndex cv2pr(Ncv);
	RanIndex mu2pr;
	
	RanIndex *e2bch;
	
	//MuGrpMiss dataI(dfNam, mMatFnam, mVecFnam, e2e, d);
	MuGrpEEmiss dataI(dfNam, sdfNam, eIndFnam, mMatFnam, mVecFnam, e2e, d);
	Grp &data = dataI;
	
	MuGrp muLnI(data, e2ln, ln2bv, LNout);
	Grp &muLn = muLnI;
	
	MuGrp muI(muLn, ln2mu, mu2pr);
	Grp &mu = muI;
	
	MuGrp muDevI = muLn;
	Grp &muDev   = muDevI;
	
	MuGrp cvPredI = data - muLn;
	Grp &cvPred   = cvPredI;
	
	MuGrp batchPredI = data;
	Grp &batchPred   = batchPredI;
	
	Grp *batchI;
	if (trtSet == "all") {
		batchI = new MuBlk(batchPred, batchFacVec, Nbatch, b2pr, "phenoData/allBlkInd.gbin");
	}
	else {
		e2bch  = new RanIndex(N, Nbatch, batchFacVec);
		batchI = new MuGrp(batchPred, *e2bch, b2pr);
	}
	Grp &batch = *batchI;
	
	MuGrp corrMuLnI = muLn + cvPred + batch;
	Grp &corrMuLn   = corrMuLnI;
	
	MuGrp eDevI = data - muLn;
	Grp &eDev   = eDevI;
	
	MuGrp muLnPredI = data - cvPred;
	Grp &muLnPred   = muLnPredI;
	
	BetaGrpFt betaCvI(cvPred, cvFlNam, Ncv, cv2pr, nThr);
	Grp &betaCv = betaCvI;
	
	MuGrp gmPredI = muLn - mu;
	Grp &gmPred   = gmPredI;
	
	BetaGrpPC gammaI(gmPred, pcEvecFlNam, pcEvalFlNam, Npc, gm2pr, nThr);
	Grp &gamma = gammaI;
	
	MuGrp bvI = mu + gamma;
	Grp &bv   = bvI;
	
	MuGrp scaDevI = muLn - bv;
	Grp &scaDev   = scaDevI;
	
	SigmaI SigIe(eDev, SgEout, 1.0, 2.0);
	SigmaI SigIb(batch, 1.0, 25.0);
	SigmaI SigIs(scaDev, 1.0, 2.0);
	SigmaI SigIa(gamma, 1.0, 2.0);
	SigmaI SigIpr(d, 1e-6);
	
	Qgrp qG(Npc, nuG);
	Qgrp qE(N, nuE, mVecFnam);
	
	cout << "Burn-in..." << endl;
	for (int iBnin = 0; iBnin < Nbnin; iBnin++) {
		data.update(corrMuLn, SigIe);
		cvPredI = data - muLn - batch;
		
		betaCv.update(cvPred, qE, SigIe, SigIpr);
		
		batchPredI = data - muLn - betaCv;
		batch.update(batchPred, qE, SigIe, SigIb);
		SigIb.update(batch);
		
		muLnPredI = data - betaCv - batch;
		muLn.update(muLnPred, qE, SigIe, bv, SigIs);
		
		corrMuLnI = muLn + betaCv + batch;
		eDevI     = data - corrMuLn;
		
		SigIe.update(eDev, qE);
		qE.update(eDev, SigIe);
		
		gmPredI = muLn - mu;
		gamma.update(gmPred, SigIe, qG, SigIa);
		SigIa.update(gamma, qG);
		qG.update(gamma, SigIa);
		
		bvI     = mu + gamma;
		scaDevI = muLn - bv;
		
		SigIs.update(scaDev);
		
		muDevI = muLn - gamma;
		mu.update(muDev, SigIe, SigIpr);
		
		
		cout << "+" << flush;
	}
	cout << endl;
	
	cout << "Initializing SNP model..." << endl;
	
	Grp *snpBetI;
	if (snpModel == "SMC"){
		snpBetI = new BetaGrpSnpMissCV(snpIn, betOut, Nln, Nsnp, d, nThr, prABF, -9);
	}
	else if (snpModel == "SMP"){
		snpBetI = new BetaGrpPSRmiss(snpIn, betOut, Nln, Nsnp, d, nThr, prABF, -9);
	}
	else {
		snpBetI = new BetaGrpSnpMiss(snpIn, betOut, Nln, Nsnp, d, nThr, prABF, -9);
	}
	Grp &snpBet = *snpBetI;
	
	for (int iSam = 0; iSam < Nsamp; iSam++) {
		data.update(corrMuLn, SigIe);
		cvPredI = data - muLn - batch;
		betaCv.update(cvPred, qE, SigIe, SigIpr);
		
		batchPredI = data - muLn - betaCv;
		batch.update(batchPred, qE, SigIe, SigIb);
		SigIb.update(batch);
		
		muLnPredI = data - betaCv - batch;
		muLn.update(muLnPred, qE, SigIe, bv, SigIs);
		
		corrMuLnI = muLn + betaCv + batch;
		eDevI     = data - corrMuLn;
		
		SigIe.update(eDev, qE);
		qE.update(eDev, SigIe);
		
		gmPredI = muLn - mu;
		gamma.update(gmPred, SigIe, qG, SigIa);
		SigIa.update(gamma, qG);
		qG.update(gamma, SigIa);
		
		bvI     = mu + gamma;
		scaDevI = muLn - bv;
		
		SigIs.update(scaDev);
		
		muDevI = muLn - gamma;
		mu.update(muDev, SigIe, SigIpr);
		
		if ((iSam + 1) % Nthin) {
			cout << "." << flush;
		}
		else {
			cout << "|" << flush;
			muLn.save();
			bv.save(BVout);
			snpBet.update(scaDev, SigIs);
			SigIe.save();
			
		}
	}
	cout << endl;
	cout << "Saving GWA results..." << endl;
	
	snpBet.dump();
	
	delete snpBetI;
}

