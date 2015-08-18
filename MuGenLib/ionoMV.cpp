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
	bool sOn = false; // sampling
	bool bOn = false; // burnin-in
	bool tOn = false; // thinning
	bool nOn = false; // nuRep
	bool gOn = false; // nuG
	bool pOn = false; // pop data set
	bool dOn = false; // trait data set
	bool mOn = false; // SNP model
	bool SOn = false; // kind of phenotypic data
	bool TOn = false; // # of threads
	bool aOn = false; // ABF prior
	bool cOn = false; // chain ID
	
	double nuE   = 3.0;
	double nuG   = 3.0;
	double prABF = 0.0;
	int Nbnin    = 10;
	int Nsamp    = 10;
	int Nthin    = 1;
	size_t Nsnp  = 678742;
	int nThr     = 4;
	
	string snpModel("SM");
	string trtSet;
	char trtTok;
	char pssTok = 'A';
	string popID;
	char popTok;
	string cNam;
	
	string dataFolder("phenoData/");
	string dfNam("Y");
	string sdfNam("VAR");
	string eIndFnam("errInd");
	string cvFlNam("PCecov");
	string mMatFnam("MatNAind");
	string mVecFnam("TotNAind");
	string pcEvecFlNam("ionHDRA_EVC");
	string pcEvalFlNam("ionHDRA_EVL");
	string expFacVec("ExpInd");
	string lnFacVec("LnLInd");
	string batchFacVec("BatchInd");
	string bPrFacVec("BprInd");
	string tubFacVec("TubInd");
	string tubPrFacVec("TubPrInd");
	string LNout("LN");
	string SDout("SD");
	string BVout("BV");
	string TBout("TB");
	string BTout("BT");
	string EXout("EX");
	string SgEout("Sig_e");
	string SgEXout("Sig_exp");
	string SgSout("Sig_s");
	string SgBout("Sig_bt");
	string SgAout("Sig_a");
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
						
					case 'T':
						TOn = true;
						break;
						
					case 'a':
						aOn = true;
						break;
						
					case 'c':
						cOn = true;
						break;
						
					case 'S':
						SOn = true;
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
				else if (SOn){
					SOn = false;
					switch (pchar[0]) {
						case 'n':
							dataFolder = "phenoData/";
							break;
						case 's':
							dataFolder = "phenoDataSC/";
							break;
						case 'r':
							pssTok = 'R'; // "reduced" number of traits
							dataFolder = "phenoDataSM/";
							break;
						default:
							cerr << "ERROR: unknown phenotype scaling token " << pchar[0] << endl;
							exit(-1);
							break;
					}

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
							break;
						case 'S':
							trtSet   = "shoot";
							trtTok   = 'S';
							cvFlNam += "Shoot";
							break;
						case 'B':
							trtSet   = "all"; // 'B' for 'both' or 'A' for 'all' works
							trtTok   = 'A';
							cvFlNam += "All";
							
							break;
						case 'A':
							trtSet   = "all";
							trtTok   = 'A';
							cvFlNam += "All";
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
				else if (TOn){
					TOn = false;
					nThr = atoi(pchar);
				}
				else if (aOn){
					aOn = false;
					prABF = atof(pchar);
				}
				else if (cOn){
					cOn = false;
					cNam = pchar;
				}
				break;
			}
		}
	}
	if (popID == "all") {
		dfNam       = dataFolder + dfNam + trtSet + ".gbin";
		sdfNam      = dataFolder + sdfNam + trtSet + ".gbin";
		cvFlNam     = "phenoData/" + cvFlNam + ".gbin";
		mMatFnam    = dataFolder + trtSet + mMatFnam + ".gbin";
		mVecFnam    = dataFolder + trtSet + mVecFnam + ".gbin";
		expFacVec   = "phenoData/" + trtSet + expFacVec + ".gbin";
		lnFacVec    = "phenoData/" + trtSet + lnFacVec + ".gbin";
		batchFacVec = "phenoData/" + trtSet + batchFacVec + ".gbin";
		bPrFacVec   = "phenoData/" + trtSet + bPrFacVec + ".gbin";
		tubFacVec   = "phenoData/" + trtSet + tubFacVec + ".gbin";
		tubPrFacVec = "phenoData/" + trtSet + tubPrFacVec + ".gbin";
	}
	else{
		dfNam       = dataFolder + dfNam + trtSet + "_" + popID + ".gbin";
		sdfNam      = dataFolder + sdfNam + trtSet + "_" + popID + ".gbin";
		cvFlNam     = "phenoData/" + cvFlNam + "_" + popID + ".gbin";
		mMatFnam    = dataFolder + trtSet + mMatFnam + "_" + popID + ".gbin";
		mVecFnam    = dataFolder + trtSet + mVecFnam + "_" + popID + ".gbin";
		expFacVec   = "phenoData/" + trtSet + expFacVec + "_" + popID + ".gbin";
		lnFacVec    = "phenoData/" + trtSet + lnFacVec + "_" + popID + ".gbin";
		batchFacVec = "phenoData/" + trtSet + batchFacVec + "_" + popID + ".gbin";
		bPrFacVec   = "phenoData/" + trtSet + bPrFacVec + "_" + popID + ".gbin";
		tubFacVec   = "phenoData/" + trtSet + tubFacVec + "_" + popID + ".gbin";
		tubPrFacVec = "phenoData/" + trtSet + tubPrFacVec + "_" + popID + ".gbin";
	}
	pcEvecFlNam = pcEvecFlNam + popID + ".gbin";
	pcEvalFlNam = pcEvalFlNam + popID + ".gbin";
	snpIn       = snpIn + popID + ".gbin";
	
	LNout   = trtSet + LNout;
	SDout   = trtSet + SDout;
	BVout   = trtSet + BVout;
	BTout   = trtSet + BTout;
	TBout   = trtSet + TBout;
	EXout   = trtSet + EXout;
	SgEout  = trtSet + SgEout;
	SgEXout = trtSet + SgEXout;
	SgSout  = trtSet + SgSout;
	SgBout  = trtSet + SgBout;
	SgAout  = trtSet + SgAout;
	betOut  = trtSet + betOut;
	
	vector<string> addOn;
	addOn.push_back(popID);
	addOn.push_back(snpModel);
	addOn.push_back(db2str(nuE));
	addOn.push_back(db2str(nuG));
	addOn.push_back(cNam);
	
	finishFlNam(LNout, addOn);
	finishFlNam(SDout, addOn);
	finishFlNam(BVout, addOn);
	finishFlNam(BTout, addOn);
	finishFlNam(TBout, addOn);
	finishFlNam(EXout, addOn);
	finishFlNam(SgEout, addOn);
	finishFlNam(SgEXout, addOn);
	finishFlNam(SgSout, addOn);
	finishFlNam(SgBout, addOn);
	finishFlNam(SgAout, addOn);
	finishFlNam(betOut, addOn);
	
	const size_t Nel  = 2;
	const size_t Ntub = 8;
	
	size_t N;
	size_t Nxp;
	size_t Nln;
	size_t Nbatch;
	size_t Ncv;
	size_t d;
	if (pssTok == 'R') {
		trtSet == "all" ? d = 46 : d = 23;
	}
	else {
		trtSet == "all" ? d = 50 : d = 25;
	}
	
	switch (popTok) {
		case 'I':
			Nln  = 131;
			Nsnp = 586773;
			switch (trtTok) {
				case 'R':
					N      = 475;
					Nxp    = 250;
					Nbatch = 20;
					Ncv    = 17;
					
					eIndFnam = dataFolder + eIndFnam + "R.txt";
					break;
					
				case 'S':
					N      = 481;
					Nxp    = 250;
					Nbatch = 23;
					Ncv    = 17;
					
					eIndFnam = dataFolder + eIndFnam + "S.txt";
					break;
				
				default: // is all
					N      = 481;
					Nxp    = 250;
					Nbatch = 23;
					Ncv    = 17;
					
					eIndFnam = dataFolder + eIndFnam + "A.txt";
					break;
			}
			break;
			
		case 'J':
			Nln  = 232;
			Nsnp = 578377;
			switch (trtTok) {
				case 'R':
					N      = 862;
					Nxp    = 454;
					Nbatch = 22;
					Ncv    = 18;
					
					eIndFnam = dataFolder + eIndFnam + "R.txt";
					break;
					
				case 'S':
					N      = 873;
					Nxp    = 455;
					Nbatch = 23;
					Ncv    = 18;
					
					eIndFnam = dataFolder + eIndFnam + "S.txt";
					break;
					
				default: // is all
					N      = 876;
					Nxp    = 455;
					Nbatch = 23;
					Ncv    = 18;
					
					eIndFnam = dataFolder + eIndFnam + "A.txt";
					break;
			}
			break;
		
		default:
			Nln = 373;
			switch (trtTok) {
				case 'R':
					N      = 1377;
					Nxp    = 724;
					Nbatch = 22;
					Ncv    = 18;
					
					eIndFnam = dataFolder + eIndFnam + "R.txt";
					break;
					
				case 'S':
					N      = 1393;
					Nxp    = 725;
					Nbatch = 23;
					Ncv    = 19;
					
					eIndFnam = dataFolder + eIndFnam + "S.txt";
					break;
					
				default: // is all
					N      = 1397;
					Nxp    = 725;
					Nbatch = 23;
					Ncv    = 18;
					
					eIndFnam = dataFolder + eIndFnam + "A.txt";
					break;
			}
			
			break;
	}
	
	const size_t Npc = Nln - 1;
	
	RanIndex e2exp(N, Nxp, expFacVec);
	RanIndex e2ln(N, Nln, lnFacVec);
	RanIndex e2e(N, N);
	RanIndex e2tub(N, Ntub, tubFacVec);
	RanIndex tub2pr(Ntub);
	RanIndex el2pr(Nel);
	RanIndex ln2mu(Nln);
	RanIndex ln2ln(Nln, Nln);
	RanIndex exp2mu(Nxp);
	RanIndex gm2pr(Npc);
	RanIndex cv2pr(Ncv);
	RanIndex mu2pr;
	
	RanIndex *e2bch;
	RanIndex bch2pr(Nbatch);
	
	MuGrpEEmiss dataI(dfNam, sdfNam, eIndFnam, mMatFnam, mVecFnam, e2e, d);
	//MuGrpMiss dataI(dfNam, mMatFnam, mVecFnam, e2e, d);
	Grp &data = dataI;
	
	MuGrpPEX sExpI(data, e2exp, exp2mu, EXout, SgEXout, 1e-6, nThr);
	Grp &sExp = sExpI;
	
	MuGrp cvPredI = data - sExp;
	Grp &cvPred   = cvPredI;
	
	BetaGrpFt betaCvI(cvPred, cvFlNam, Ncv, cv2pr, BTout, nThr);
	Grp &betaCv = betaCvI;
	
	MuGrp tubI(cvPred, e2tub, tub2pr);
	Grp &tub = tubI;
	
	Grp *batchI;
	if (trtSet == "all") {
		batchI = new MuBlk(cvPred, batchFacVec, Nbatch, bch2pr, "phenoData/allBlkInd.gbin");
	}
	else {
		e2bch = new RanIndex(N, Nbatch, batchFacVec);
		batchI = new MuGrp(cvPred, *e2bch, bch2pr);
	}
	Grp &batch = *batchI;
	
	MuGrp tBchI = tub + batch;
	Grp &tBch   = tBchI;
	
	MuGrp beTbchI = betaCv + tBch;
	Grp &beTbch   = beTbchI;
	
	MuGrp bbI = betaCv + batch;
	Grp &bb   = bbI;
	
	MuGrp beTI = betaCv + tub;
	Grp &beT   = beTI;
	
	MuGrp muLnPredI = data - beTbch;
	Grp &muLnPred   = muLnPredI;
	
	MuGrp muLnI(muLnPred, e2ln, ln2ln, LNout);
	Grp &muLn = muLnI;
	
	MuGrp dLnI = data - muLn;
	Grp &dLn   = dLnI;
	
	MuGrp dLnExpI = dLn - sExp;
	Grp &dLnExp   = dLnExpI;
	
	MuGrp lnExpI = muLn + sExp;
	Grp &lnExp   = lnExpI;
	
	MuGrp sExpPredI = dLn;
	Grp &sExpPred   = sExpPredI;
	
	MuGrp dataPrI = lnExp + betaCv;
	Grp &dataPr   = dataPrI;
	
	MuGrp tubPredI = dLn;
	Grp &tubPred   = tubPredI;
	
	MuGrp bchPredI = dLn;
	Grp &bchPred   = bchPredI;
	
	MuGrp muI(muLn, ln2mu, mu2pr);
	Grp &mu = muI;
	
	MuGrp gamPredI = muLn - mu;
	Grp &gamPred   = gamPredI;
	
	BetaGrpPCpex gammaI(gamPred, pcEvecFlNam, pcEvalFlNam, Npc, 1e-6, gm2pr, nThr);
	Grp &gamma = gammaI;
	
	MuGrp bvI = mu + gamma;
	Grp &bv   = bvI;
	
	MuGrp muPredI = muLn - gamma;
	Grp &muPred   = muPredI;
	
	MuGrp eDevI = data - dataPr;
	Grp &eDev   = eDevI;
	
	MuGrp scaDevI = muLn - bv;
	Grp &scaDev   = scaDevI;
	
	SigmaI SigIe(eDev, SgEout, 1.0, 2.0);
	SigmaI SigItub(tub, 1.0, static_cast<double>(d));
	SigmaI SigIbch(batch, SgBout, 1.0, static_cast<double>(d));
	SigmaI SigIexp(sExp, 1.0, 2.0);
	SigmaI SigIs(scaDev, SgSout, 1.0, 2.0);
	SigmaI SigIa(gamma, SgAout, 1.0, 2.0);
	SigmaI SigIpr(d, 1e-6);
	
	Qgrp qE(N, nuE, mVecFnam);
	Qgrp qG(Npc, nuG);
	
	cout << "Burn-in..." << endl;
	for (int iBnin = 0; iBnin < Nbnin; iBnin++) {
		data.update(dataPr, qE, SigIe);
		//data.update(dataPr, SigIe);
		
		cvPredI = dLnExp - tBch;
		betaCv.update(cvPred, qE, SigIe, SigIpr);
		
		tubPredI = dLnExp - bb;
		tub.update(tubPred, qE, SigIe, SigItub);
		SigItub.update(tub);
		
		beTI     = betaCv + tub;
		bchPredI = dLnExp - beT;
		batch.update(bchPred, qE, SigIe, SigIbch);
		SigIbch.update(batch);
		
		tBchI   = tub + batch;
		beTbchI = betaCv + tBch;
		bbI     = betaCv + batch;
		
		sExpPredI = dLn - beTbch;
		sExp.update(sExpPred, qE, SigIe, SigIexp);
		SigIexp.update(sExp);
		
		muLnPredI = data - sExp - beTbch;
		muLn.update(muLnPred, qE, SigIe, bv, SigIs);
		
		dLnI    = data - muLn;
		dLnExpI = dLn - sExp;
		lnExpI  = muLn + sExp;
		dataPrI = lnExp + beTbch;
		eDevI   = data - dataPr;
		SigIe.update(eDev, qE);
		qE.update(eDev, SigIe);
		
		gamPredI = muLn - mu;
		gamma.update(gamPred, SigIs, qG, SigIa);
		SigIa.update(gamma, qG);
		qG.update(gamma, SigIa);
		
		muPredI = muLn - gamma;
		mu.update(muPred, SigIs, SigIpr);
		
		bvI     = mu + gamma;
		scaDevI = muLn - bv;
		
		SigIs.update(scaDev);
		
		cout << "+" << flush;
	}
	cout << endl;
	
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
		
	cout << "Sampling..." << endl;
	for (int iSam = 0; iSam < Nsamp; iSam++) {
		data.update(dataPr, qE, SigIe);
		//data.update(dataPr, SigIe);
		
		cvPredI = dLnExp - tBch;
		betaCv.update(cvPred, qE, SigIe, SigIpr);
		
		tubPredI = dLnExp - bb;
		tub.update(tubPred, qE, SigIe, SigItub);
		SigItub.update(tub);
		
		beTI     = betaCv + tub;
		bchPredI = dLnExp - beT;
		batch.update(bchPred, qE, SigIe, SigIbch);
		SigIbch.update(batch);
		
		tBchI   = tub + batch;
		beTbchI = betaCv + tBch;
		bbI     = betaCv + batch;
		
		sExpPredI = dLn - beTbch;
		sExp.update(sExpPred, qE, SigIe, SigIexp);
		SigIexp.update(sExp);
		
		muLnPredI = data - sExp - beTbch;
		muLn.update(muLnPred, qE, SigIe, bv, SigIs);
		
		dLnI    = data - muLn;
		dLnExpI = dLn - sExp;
		lnExpI  = muLn + sExp;
		dataPrI = lnExp + beTbch;
		eDevI   = data - dataPr;
		SigIe.update(eDev, qE);
		qE.update(eDev, SigIe);
		
		gamPredI = muLn - mu;
		gamma.update(gamPred, SigIs, qG, SigIa);
		SigIa.update(gamma, qG);
		qG.update(gamma, SigIa);
		
		muPredI = muLn - gamma;
		mu.update(muPred, SigIs, SigIpr);
		
		bvI     = mu + gamma;
		scaDevI = muLn - bv;
		
		SigIs.update(scaDev);
		
		if ((iSam + 1) % Nthin) {
			cout << "." << flush;
		}
		else {
			cout << "|" << flush;
			//tubSI = tub - muT;
			tub.save(TBout);
			betaCv.save();
			SigIbch.save();
			scaDev.save(SDout);
			bv.save(BVout);
			muLn.save(LNout);
			sExp.save(SigIexp);
			SigIe.save();
			SigIs.save();
			SigIa.save();
			snpBet.update(scaDev, SigIs);
			
		}
	}
	cout << endl;
	cout << "Saving GWA results..." << endl;
	
	snpBet.dump();
	
	delete snpBetI;
	delete batchI;
	delete e2bch;
}

