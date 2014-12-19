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
	string expFacVec("ExpInd");
	string lnFacVec("LnInd");
	string batchFacVec("BatchInd");
	string bPrFacVec("BprInd");
	string tubFacVec("TubInd");
	string tubPrFacVec("TubPrInd");
	string LNout("LN");
	string BVout("BV");
	string SgEout("Sig_e");
	string SgEXout("Sig_exp");
	string SgSout("Sig_s");
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
							eIndFnam = "phenoData/" + eIndFnam + "R.txt";
							break;
						case 'S':
							trtSet   = "shoot";
							trtTok   = 'S';
							cvFlNam += "Shoot";
							eIndFnam = "phenoData/" + eIndFnam + "S.txt";
							break;
						case 'B':
							trtSet   = "all"; // 'B' for 'both' or 'A' for 'all' works
							trtTok   = 'A';
							cvFlNam += "All";
							eIndFnam = "phenoData/" + eIndFnam + "A.txt";
							break;
						case 'A':
							trtSet   = "all";
							trtTok   = 'A';
							cvFlNam += "All";
							eIndFnam = "phenoData/" + eIndFnam + "A.txt";
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
		dfNam       = "phenoData/" + dfNam + trtSet + ".gbin";
		sdfNam      = "phenoData/" + sdfNam + trtSet + ".gbin";
		cvFlNam     = "phenoData/" + cvFlNam + ".gbin";
		mMatFnam    = "phenoData/" + trtSet + mMatFnam + ".gbin";
		mVecFnam    = "phenoData/" + trtSet + mVecFnam + ".gbin";
		expFacVec   = "phenoData/" + trtSet + expFacVec + ".gbin";
		lnFacVec    = "phenoData/" + trtSet + lnFacVec + ".gbin";
		batchFacVec = "phenoData/" + trtSet + batchFacVec + ".gbin";
		bPrFacVec   = "phenoData/" + trtSet + bPrFacVec + ".gbin";
		tubFacVec   = "phenoData/" + trtSet + tubFacVec + ".gbin";
		tubPrFacVec = "phenoData/" + trtSet + tubPrFacVec + ".gbin";
	}
	else{
		dfNam       = "phenoData/" + dfNam + trtSet + "_" + popID + ".gbin";
		sdfNam      = "phenoData/" + sdfNam + trtSet + "_" + popID + ".gbin";
		cvFlNam     = "phenoData/" + cvFlNam + "_" + popID + ".gbin";
		mMatFnam    = "phenoData/" + trtSet + mMatFnam + "_" + popID + ".gbin";
		mVecFnam    = "phenoData/" + trtSet + mVecFnam + "_" + popID + ".gbin";
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
	BVout   = trtSet + BVout;
	SgEout  = trtSet + SgEout;
	SgEXout = trtSet + SgEXout;
	SgSout  = trtSet + SgSout;
	betOut  = trtSet + betOut;
	
	vector<string> addOn;
	addOn.push_back(snpModel);
	addOn.push_back(db2str(nuE));
	addOn.push_back(db2str(nuG));
	addOn.push_back(cNam);
	
	finishFlNam(LNout, addOn);
	finishFlNam(BVout, addOn);
	finishFlNam(SgEout, addOn);
	finishFlNam(SgEXout, addOn);
	finishFlNam(SgSout, addOn);
	finishFlNam(betOut, addOn);
	
	const size_t Nel  = 2;
	const size_t Ntub = 8;
	
	size_t N;
	size_t Nxp;
	size_t Nln;
	size_t Nbatch;
	size_t Ncv;
	size_t d;
	trtSet == "all" ? d = 50 : d = 25;
	
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
					break;
					
				case 'S':
					N      = 481;
					Nxp    = 250;
					Nbatch = 23;
					Ncv    = 17;
					break;
				
				default: // is all
					N      = 481;
					Nxp    = 250;
					Nbatch = 23;
					Ncv    = 17;
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
					break;
					
				case 'S':
					N      = 873;
					Nxp    = 455;
					Nbatch = 23;
					Ncv    = 18;
					break;
					
				default: // is all
					N      = 876;
					Nxp    = 455;
					Nbatch = 23;
					Ncv    = 18;
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
					break;
					
				case 'S':
					N      = 1393;
					Nxp    = 725;
					Nbatch = 23;
					Ncv    = 19;
					break;
					
				default: // is all
					N      = 1397;
					Nxp    = 725;
					Nbatch = 23;
					Ncv    = 18;
					break;
			}
			
			break;
	}
	
	const size_t Npc = Nln - 1;
	
	RanIndex e2exp(N, Nxp, expFacVec);
	RanIndex exp2ln(Nxp, Nln, lnFacVec);
	RanIndex e2e(N, N);
	RanIndex e2tub(N, Ntub, tubFacVec);
	//RanIndex tub2pr(Ntub, Nel, tubPrFacVec);
	RanIndex tub2pr(Ntub);
	RanIndex el2pr(Nel);
	RanIndex ln2mu(Nln);
	RanIndex exp2mu(Nxp);
	RanIndex gm2pr(Npc);
	RanIndex cv2pr(Ncv);
	RanIndex mu2pr;
	
	RanIndex *e2bch;
	//RanIndex *bch2pr;
	RanIndex bch2pr(Nbatch);
	//RanIndex bch2pr(Nbatch, Nel, bPrFacVec);
	
	MuGrpEEmiss dataI(dfNam, sdfNam, eIndFnam, mMatFnam, mVecFnam, e2e, d);
	Grp &data = dataI;
	
	MuGrp muExpI(data, e2exp, exp2ln);
	Grp &muExp = muExpI;
	
	MuGrp cvPredI = data - muExp;
	Grp &cvPred   = cvPredI;
	
	BetaGrpFt betaCvI(cvPred, cvFlNam, Ncv, cv2pr, nThr);
	Grp &betaCv = betaCvI;
	
	MuGrp tubI(data, e2tub, tub2pr);
	//MuGrpPEX tubI(data, e2tub, tub2pr, 1e-6, nThr);
	Grp &tub = tubI;
	
	MuGrp muTI(tub, tub2pr, mu2pr);
	Grp &muT = muTI;
	
	/*
	MuGrp tubPrI(tub, tub2pr, el2pr);
	Grp &tubPr = tubPrI;
	
	MuGrp tubHPI(tubPr, el2pr, mu2pr);
	Grp &tubHP = tubHPI;
	*/
	Grp *batchI;
	if (trtSet == "all") {
		batchI = new MuBlk(data, batchFacVec, Nbatch, bch2pr, "phenoData/allBlkInd.gbin");
	}
	else {
		e2bch = new RanIndex(N, Nbatch, batchFacVec);
		batchI = new MuGrp(data, *e2bch, bch2pr);
		//batchI = new MuGrp(data, *e2bch, bch2pr);
	}
	Grp &batch = *batchI;
	
	MuGrp muBI(batch, bch2pr, mu2pr);
	Grp &muB = muBI;
	/*
	MuGrp bchPrI(batch, bch2pr, el2pr);
	Grp &bchPr = bchPrI;
	
	MuGrp bchHPI(bchPr, el2pr, mu2pr);
	Grp &bchHP = bchHPI;
	*/
	MuGrp dExpI = data - muExp;
	Grp &dExp   = dExpI;
	
	MuGrp tBchI = tub + batch;
	Grp &tBch   = tBchI;
	
	MuGrp muExpPredI = data - betaCv - tBch;
	Grp &muExpPred   = muExpPredI;
	
	MuGrp beTbchI = betaCv + tBch;
	Grp &beTbch   = beTbchI;
	
	MuGrp dataPrI = muExp + beTbch;
	Grp &dataPr   = dataPrI;
	
	MuGrp tubPredI = data - muExp - betaCv;
	Grp &tubPred   = tubPredI;
	
	MuGrp bchPredI = dExp - betaCv - tub;
	Grp &bchPred   = bchPredI;
	
	MuGrp muI(muExp, exp2mu, mu2pr);
	Grp &mu = muI;
	
	MuGrp scaPredI = muExp - mu;
	Grp &scaPred   = scaPredI;
	
	//MuGrpPEX scaI(scaPred, exp2ln, ln2mu, 1e-6, nThr);
	MuGrp scaI(scaPred, exp2ln, ln2mu);
	Grp &sca = scaI;
	
	//MuGrp muSI(sca, ln2mu, mu2pr);
	//Grp &muS = muSI;
	
	MuGrp gamPredI = muExp - mu;
	Grp &gamPred   = gamPredI;
	
	//BetaGrpPC gammaI(gamPred, pcEvecFlNam, pcEvalFlNam, Npc, exp2ln, gm2pr, nThr);
	BetaGrpPCpex gammaI(gamPred, pcEvecFlNam, pcEvalFlNam, Npc, 1e-6, exp2ln, gm2pr, nThr);
	Grp &gamma = gammaI;
	
	MuGrp bvI = mu + gamma;
	Grp &bv   = bvI;
	
	MuGrp gsI = gamma + sca;
	Grp &gs   = gsI;
	
	MuGrp msI = mu + sca;
	Grp &ms   = msI;
	
	MuGrp muLnI = bv + sca;
	Grp &muLn = muLnI;
	
	MuGrp muPredI = muExp - sca;
	Grp &muPred   = muPredI;
	
	MuGrp eDevI = data - dataPr;
	Grp &eDev   = eDevI;
	
	printMat(muExp.fMat());
	cout << "#####" << endl;
	
	printMat(muLn.fMat());
	cout << "#####" << endl;
	
	MuGrp expDevI = muExp - muLn;
	Grp &expDev   = expDevI;
	
	printMat(expDev.fMat());
	cout << "#####" << endl;
	
	SigmaIpex SigIeI(eDev, SgEout, 1.0, 2.0);
	SigmaI &SigIe = SigIeI;
	//SigmaI SigIe(eDev, SgEout, 1.0, 2.0);
	//SigmaI SigIe(eDev, SgEout, 1.0, 2000.0);
	//SigmaI SigItub(tub, 1.0, 2.0);
	SigmaI SigItub(tub, 1.0, static_cast<double>(d));
	//SigmaI SigItpr(tubPr, 1.0, static_cast<double>(d));
	//SigmaI SigIbch(batch, 1.0, 2.0);
	SigmaI SigIbch(batch, 1.0, static_cast<double>(d));
	//SigmaI SigIbpr(bchPr, 1.0, static_cast<double>(d));
	SigmaI SigIexp(expDev, SgEXout, 1.0, 2.0);
	//SigmaI SigIexp(expDev, SgEXout, 1.0, 2000.0);
	SigmaI SigIs(sca, SgSout, 1.0, 2.0);
	SigmaI SigIa(gamma, 1.0, 2.0);
	SigmaI SigIpr(d, 1e-6);
	
	QgrpPEX qEI(N, nuE, mVecFnam);
	Qgrp &qE = qEI;
	//Qgrp qE(N, nuE, mVecFnam);
	Qgrp qG(Npc, nuG);
	
	cout << "Burn-in..." << endl;
	for (int iBnin = 0; iBnin < Nbnin; iBnin++) {
		data.update(dataPr, qE, SigIe);
		dExpI = data - muExp;
		
		cvPredI = dExp - tBch;
		betaCv.update(cvPred, qE, SigIe, SigIpr);
		
		tubPredI = dExp - betaCv - batch;
		//tub.update(tubPred, qE, SigIe, tubPr, SigItub);
		//tub.update(tubPred, qE, SigIe, SigItub);
		tub.update(tubPred, qE, SigIe, muT, SigItub);
		muT.update(tub, SigItub, SigIpr);
		//tubPr.update(tub, SigItub, tubHP, SigItpr);
		//SigItub.update(tub, tubPr);
		//SigItub.update(tub);
		SigItub.update(tub, muT);
		//tubHP.update(tubPr, SigItpr, SigIpr);
		//SigItpr.update(tubPr, tubHP);
		
		bchPredI = dExp - betaCv - tub;
		//batch.update(bchPred, qE, SigIe, bchPr, SigIbch);
		//bchPr.update(batch, SigIbch, bchHP, SigIbpr);
		//batch.update(bchPred, qE, SigIe, SigIbch);
		batch.update(bchPred, qE, SigIe, muB, SigIbch);
		muB.update(batch, SigIbch, SigIpr);
		SigIbch.update(batch, muB);
		
		tBchI   = tub + batch;
		beTbchI = betaCv + tBch;
		
		muExpPredI = data - beTbch;
		muExp.update(muExpPred, qE, SigIe, muLn, SigIexp);
		
		dataPrI = muExp + beTbch;
		eDevI   = data - dataPr;
		SigIe.update(eDev, qE);
		qE.update(eDev, SigIe);
		
		scaPredI = muExp - bv;
		//sca.update(scaPred, SigIexp, muS, SigIs);
		sca.update(scaPred, SigIexp, SigIs);
		//muS.update(sca, SigIs, SigIpr);
		//SigIs.update(sca, muS);
		SigIs.update(sca);
		
		msI      = mu + sca;
		gamPredI = muExp - ms;
		gamma.update(gamPred, SigIexp, qG, SigIa);
		SigIa.update(gamma, qG);
		qG.update(gamma, SigIa);
		
		bvI     = mu + gamma;
		gsI     = gamma + sca;
		muPredI = muExp - gs;
		mu.update(muPred, SigIexp, SigIpr);
		
		muLnI   = bv + sca;
		expDevI = muExp - muLn;
		SigIexp.update(expDev);
		
		cout << "+" << flush;
	}
	cout << endl;
	
	cout << "Sampling..." << endl;
	for (int iSam = 0; iSam < Nsamp; iSam++) {
		data.update(dataPr, qE, SigIe);
		dExpI = data - muExp;
		
		cvPredI = dExp - tBch;
		betaCv.update(cvPred, qE, SigIe, SigIpr);
		
		tubPredI = dExp - betaCv - batch;
		//tub.update(tubPred, qE, SigIe, tubPr, SigItub);
		//tub.update(tubPred, qE, SigIe, SigItub);
		tub.update(tubPred, qE, SigIe, muT, SigItub);
		muT.update(tub, SigItub, SigIpr);
		//tubPr.update(tub, SigItub, tubHP, SigItpr);
		//SigItub.update(tub, tubPr);
		//SigItub.update(tub);
		SigItub.update(tub, muT);
		//tubHP.update(tubPr, SigItpr, SigIpr);
		//SigItpr.update(tubPr, tubHP);
		
		bchPredI = dExp - betaCv - tub;
		//batch.update(bchPred, qE, SigIe, bchPr, SigIbch);
		//bchPr.update(batch, SigIbch, bchHP, SigIbpr);
		//batch.update(bchPred, qE, SigIe, SigIbch);
		batch.update(bchPred, qE, SigIe, muB, SigIbch);
		muB.update(batch, SigIbch, SigIpr);
		SigIbch.update(batch, muB);
		//SigIbch.update(batch);
		//SigIbch.update(batch, bchPr);
		//bchPr.update(batch, SigIbch, bchHP, SigIbpr);
		//bchHP.update(bchPr, SigIbpr, SigIpr);
		//SigIbpr.update(bchPr, bchHP);
		
		tBchI   = tub + batch;
		beTbchI = betaCv + tBch;
		
		muExpPredI = data - beTbch;
		muExp.update(muExpPred, qE, SigIe, muLn, SigIexp);
		
		dataPrI = muExp + beTbch;
		eDevI   = data - dataPr;
		SigIe.update(eDev, qE);
		qE.update(eDev, SigIe);
		
		scaPredI = muExp - bv;
		//sca.update(scaPred, SigIexp, muS, SigIs);
		sca.update(scaPred, SigIexp, SigIs);
		//muS.update(sca, SigIs, SigIpr);
		//SigIs.update(sca, muS);
		SigIs.update(sca);
		
		msI      = mu + sca;
		gamPredI = muExp - ms;
		gamma.update(gamPred, SigIexp, qG, SigIa);
		SigIa.update(gamma, qG);
		qG.update(gamma, SigIa);
		
		bvI     = mu + gamma;
		gsI     = gamma + sca;
		muPredI = muExp - gs;
		mu.update(muPred, SigIexp, SigIpr);
		
		muLnI   = bv + sca;
		expDevI = muExp - muLn;
		SigIexp.update(expDev);
		
		if ((iSam + 1) % Nthin) {
			cout << "." << flush;
		}
		else {
			cout << "|" << flush;
			muLn.save(LNout);
			bv.save(BVout);
			SigIe.save();
			SigIexp.save();
			SigIs.save();
			
		}
	}
	cout << endl;
	//cout << "Saving GWA results..." << endl;
	
	/*
	 Grp *snpBetI;
	 if (snpModel == "SMC"){
		snpBetI = new BetaGrpSnpMissCV(snpIn, betOut, e2ln, Nsnp, d, nThr, prABF, -9);
	 }
	 else if (snpModel == "SMP"){
		snpBetI = new BetaGrpPSRmiss(snpIn, betOut, e2ln, Nsnp, d, nThr, prABF, -9);
	 }
	 else {
		snpBetI = new BetaGrpSnpMiss(snpIn, betOut, e2ln, Nsnp, d, nThr, prABF, -9);
	 }
	 Grp &snpBet = *snpBetI;
	 
	 MuGrp snpPredI = data - bv - bb;
	 Grp &snpPred   = snpPredI;

	 */
	//snpBet.dump();
	
	//delete snpBetI;
}

