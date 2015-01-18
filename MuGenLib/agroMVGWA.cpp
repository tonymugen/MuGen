/*
*  agroMVGWA.cpp
*  HBqgen
*
*  Created by ajg67 on 1/9/13.
*   Copyright (c) 2013 SEELE. All rights reserved.
*
*	to compile on CAC:
*	icc -I/home/fs01/ajg67/gslLibs/include -L/home/fs01/ajg67/gslLibs/lib agroMVGWA.cpp HBqgen.cpp -o agroMVGWA -march=native -O3 -lgsl -lz /opt/intel/composer_xe_2011_sp1.6.233/mkl/lib/intel64/libmkl_solver_lp64.a -Wl,--start-group /opt/intel/composer_xe_2011_sp1.6.233/mkl/lib/intel64/libmkl_intel_lp64.a /opt/intel/composer_xe_2011_sp1.6.233/mkl/lib/intel64/libmkl_intel_thread.a /opt/intel/composer_xe_2011_sp1.6.233/mkl/lib/intel64/libmkl_core.a -Wl,--end-group -DHAVE_INLINE -DGSL_RANGE_CHECK_OFF -fopenmp -lpthread -lifcore -static
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

#include "libMuGen.h"

using std::vector;
using std::string;
using std::cout;
using std::cerr;
using std::endl;
using std::flush;
using std::stringstream;
using std::remove;


int main(int argc, char *argv[]){
	/*
		Parsing command line flags
	*/
	bool sOn = false;
	bool bOn = false;
	bool tOn = false;
	bool nOn = false;
	bool rOn = false;
	bool gOn = false;
	bool pOn = false;
	bool fOn = false;
	bool cOn = false;
	bool COn = false;
	bool MOn = false;
	bool mOn = false;
	bool lOn = false;
	
	// intitialize the flag variables with default values
	double nuDat = 3.0;
	double nuR   = 3.0;
	double nuG   = 3.0;
	double nuB   = 3.0;
	double ldCt  = 0.75;
	int Nbnin    = 10;
	int Nsamp    = 10;
	int Nthin    = 1;
	int Nmul     = 5;
	int nThr     = 1;
	
	string model("SM"); // default is single marker; other possibilities: VS == variable selection; RS == rank selection; flag -r controls this
	string cNum;
	
	string snpFnam("SNP.gbin");
	string pcUfNam("agroSNP_U.gbin");
	string evRfNam("agroSNP_EV.gbin");
	string SgDout("Sig_dt.gbin");
	string SgRout("Sig_rp.gbin");
	string SgYout("Sig_yr.gbin");
	string SgSout("Sig_sca.gbin");
	string SgAout("Sig_a.gbin");
	string LNout("LN.gbin");
	string betOut("betSNP.gbin");
	string pepOut("pepSNP.gbin");
	
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
						
					case 'm':
						mOn = true;
						break;
						
					case 'l':
						lOn = true;
						break;
						
					case 'n':
						nOn = true;
						break;
						
					case 'r':
						rOn = true;
						break;
						
					case 'g':
						gOn = true;
						break;
						
					case 'p':
						pOn = true;
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
						
					case 'h':
						cerr << "usage:\n"
						<< "agroMVGWA [-s MCMC_sample_number] [-b burnin_length] [-t thinning_ratio] [-c chain_ID] [-C number_of_threads] [-M multiple_of_sample_size_to_include_SNPs] [-m model_ID] [-n repStudent_t_nu] [-g gStudent_t_nu] [-f SNP_file_tag]" << endl;
						exit(-1);
					default:
						cerr << "ERROR: unrecognized flag " << pchar[1] << endl;
						exit(-1);
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
								Nmul = atoi(pchar);
							}
							else
								if (mOn) {
									mOn   = false;
									model = pchar;
								}
							else
								if (cOn) {
									cOn = false;
									cNum = pchar;
								}
								else
									if (COn) {
										COn = false;
										nThr = atoi(pchar);
									}
								else
									if (fOn) {
										fOn = false;
										size_t chPos;
										snpFnam = snpFnam.substr(0,3);
										snpFnam += pchar;
										snpFnam += ".gbin";
										
										pcUfNam = pcUfNam.substr(0,7);
										pcUfNam += pchar;
										chPos = pcUfNam.find("tr");
										if (chPos == string::npos) {
											chPos = pcUfNam.find("p");
											pcUfNam.erase(pcUfNam.begin() + chPos, pcUfNam.end());
										}
										else {
											pcUfNam.erase(pcUfNam.end() - 2, pcUfNam.end());
										}
										pcUfNam += "_U.gbin";
										
										evRfNam = evRfNam.substr(0,7);
										evRfNam += pchar;
										chPos = evRfNam.find("tr");
										if (chPos == string::npos) {
											chPos = evRfNam.find("p");
											evRfNam.erase(evRfNam.begin() + chPos, evRfNam.end());
										}
										else {
											evRfNam.erase(evRfNam.end() - 2, evRfNam.end());
										}
										evRfNam += "_EV.gbin";
										
										SgDout = SgDout.substr(0, 6);
										SgDout += "_";
										SgDout += pchar;
										SgDout += "_";
										
										SgRout = SgRout.substr(0, 6);
										SgRout += "_";
										SgRout += pchar;
										SgRout += "_";
										
										SgYout = SgYout.substr(0, 6);
										SgYout += "_";
										SgYout += pchar;
										SgYout += "_";
										
										SgSout = SgSout.substr(0, 7);
										SgSout += "_";
										SgSout += pchar;
										SgSout += "_";
										
										SgAout = SgAout.substr(0, 5);
										SgAout += "_";
										SgAout += pchar;
										SgAout += "_";
										
										LNout = LNout.substr(0, 2);
										LNout += "_";
										LNout += pchar;
										LNout += "_";
										
										betOut = betOut.substr(0, 6);
										betOut += "_";
										betOut += pchar;
										betOut += "_";
										
										pepOut = pepOut.substr(0, 6);
										pepOut += "_";
										pepOut += pchar;
										pepOut += "_";
									}
									else
										if (gOn) {
											gOn = false;
											nuG = atof(pchar);
										}
									else
										if (pOn) {
											pOn = false;
											nuB = atof(pchar);
										}
									else
										if (lOn) {
											lOn = false;
											ldCt = atof(pchar);
										}
									else
										if (rOn) {
											rOn = false;
											nuR = atof(pchar);
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
	
	stringstream nuOSd;
	nuOSd << nuDat;
	string dString = nuOSd.str();
	size_t perPs = dString.find('.');
	if (perPs != string::npos) {
		dString.erase(perPs, 1);
	}
	
	stringstream nuOSr;
	nuOSr << nuR;
	string rString = nuOSr.str();
	perPs = rString.find('.');
	if (perPs != string::npos) {
		rString.erase(perPs, 1);
	}
	
	stringstream nuOSg;
	nuOSg << nuG;
	string gString = nuOSg.str();
	perPs = gString.find('.');
	if (perPs != string::npos) {
		gString.erase(perPs, 1);
	}
		
	stringstream nuOSb;
	nuOSb << nuB;
	string bString = nuOSb.str();
	perPs = bString.find('.');
	if (perPs != string::npos) {
		bString.erase(perPs, 1);
	}
	
	SgDout += model;
	SgDout += "_";
	SgDout += dString;
	SgDout += "_";
	SgDout += rString;
	SgDout += "_";
	SgDout += gString;
	SgDout += "_";
	SgDout += bString;
	SgDout += "_";
	SgDout += cNum;
	SgDout += ".gbin";
	
	SgRout += model;
	SgRout += "_";
	SgRout += dString;
	SgRout += "_";
	SgRout += rString;
	SgRout += "_";
	SgRout += gString;
	SgRout += "_";
	SgRout += bString;
	SgRout += "_";
	SgRout += cNum;
	SgRout += ".gbin";
	
	SgYout += model;
	SgYout += "_";
	SgYout += dString;
	SgYout += "_";
	SgYout += rString;
	SgYout += "_";
	SgYout += gString;
	SgYout += "_";
	SgYout += bString;
	SgYout += "_";
	SgYout += cNum;
	SgYout += ".gbin";
	
	SgSout += model;
	SgSout += "_";
	SgSout += dString;
	SgSout += "_";
	SgSout += rString;
	SgSout += "_";
	SgSout += gString;
	SgSout += "_";
	SgSout += bString;
	SgSout += "_";
	SgSout += cNum;
	SgSout += ".gbin";
	
	SgAout += model;
	SgAout += "_";
	SgAout += dString;
	SgAout += "_";
	SgAout += rString;
	SgAout += "_";
	SgAout += gString;
	SgAout += "_";
	SgAout += bString;
	SgAout += "_";
	SgAout += cNum;
	SgAout += ".gbin";
	
	LNout += model;
	LNout += "_";
	LNout += dString;
	LNout += "_";
	LNout += rString;
	LNout += "_";
	LNout += gString;
	LNout += "_";
	LNout += bString;
	LNout += "_";
	LNout += cNum;
	LNout += ".gbin";
	
	betOut += model;
	betOut += "_";
	betOut += dString;
	betOut += "_";
	betOut += rString;
	betOut += "_";
	betOut += gString;
	betOut += "_";
	betOut += bString;
	betOut += "_";
	betOut += cNum;
	betOut += ".gbin";
	
	pepOut += model;
	pepOut += "_";
	pepOut += dString;
	pepOut += "_";
	pepOut += rString;
	pepOut += "_";
	pepOut += gString;
	pepOut += "_";
	pepOut += bString;
	pepOut += "_";
	pepOut += cNum;
	pepOut += ".gbin";
	
	/*
		defining some more constants
	*/
	gsl_vector_int *consts = gsl_vector_int_alloc(7);
	
	FILE *CstsIF = fopen("agroConsts.gbin", "r");
	gsl_vector_int_fread(CstsIF, consts);
	fclose(CstsIF);
	
	
	const size_t d    = gsl_vector_int_get(consts, 0);
	const size_t Ntot = gsl_vector_int_get(consts, 1);
	const size_t Nrep = gsl_vector_int_get(consts, 2);
	const size_t Nyr  = gsl_vector_int_get(consts, 3);
	const size_t Nln  = gsl_vector_int_get(consts, 4);
	const size_t Npc  = gsl_vector_int_get(consts, 5);
	const size_t Nsnp = gsl_vector_int_get(consts, 6);
	
	gsl_vector_int_free(consts);
	
	
	/*
		Initializing variables
	 */
	
	RanIndex dat2rep(Ntot, Nrep, string("agroLnYrRepInd.gbin"));
	RanIndex rep2yr(Nrep, Nyr, string("agroLnYrInd.gbin"));
	RanIndex yr2ln(Nyr, Nln, string("agroLnInd.gbin"));
	RanIndex ln2mu(Nln);
	RanIndex yr2mu(Nyr);
	RanIndex gm2pr(Npc);
	RanIndex mu2pr;
	
	MuGrpMiss dataI(string("agroData.gbin"), string("agroMatNAind.gbin"), string("agroTotNAind.gbin"), dat2rep, d);
	Grp &data = dataI;
	
	MuGrp muRepI(data, dat2rep, rep2yr);
	Grp &muRep = muRepI;
	
	MuGrp muYrI(muRep, rep2yr, yr2ln);
	Grp &muYr = muYrI;
	
	MuGrp muRspI = muYr;
	Grp &muRsp   = muRspI;
	
	MuGrp muI(muRsp, yr2mu, mu2pr);
	Grp &mu = muI;
	
	MuGrp scaRspI = muYr - mu;
	Grp &scaRsp   = scaRspI;
	
	MuGrpPEX scaI(scaRsp, yr2ln, ln2mu, 1.0, nThr);
	Grp &sca = scaI;
	
	MuGrp muScaI(sca, ln2mu, mu2pr);
	Grp &muSca = muScaI;
	
	MuGrp datDevI = data - muRep;
	Grp &datDev   = datDevI;
	
	MuGrp repDevI = muRep - muYr;
	Grp &repDev   = repDevI;
	
	MuGrp yrDevI = muYr - mu;
	Grp &yrDev   = yrDevI;
	
	MuGrp muLnI = mu + sca;
	Grp &muLn = muLnI;
	
	SigmaI SigIdat(datDev, 1.0, 2.0);
	SigmaI SigIrep(repDev, 1.0, 2.0);
	SigmaI SigIyr(yrDev, 1.0, 2.0);
	SigmaI SigIsc(sca, 1.0, 2.0);
	SigmaI SigIpr(d, 1e-6);
	
	Qgrp qDat(Ntot, nuDat, string("agroTotNAind.gbin"));
	Qgrp qR(Nrep, nuR);
	
	/*
		Pre-run
	 */
	
	for (int iPr = 0; iPr < 200; iPr++) {
		data.update(muRep, SigIdat);
		muRep.update(data, qDat, SigIdat, muYr, qR, SigIrep);
		
		datDevI = muRep - data;
		SigIdat.update(datDev, qDat);
		qDat.update(datDev, SigIdat);
		
		muYr.update(muRep, qR, SigIrep, muLn, SigIyr);
		
		repDevI = muYr - muRep;
		SigIrep.update(repDev, qR);
		qR.update(repDev, SigIrep);
		
		scaRspI = muYr - mu;
		sca.update(scaRsp, SigIyr, muSca, SigIsc);
		muSca.update(sca, SigIsc, SigIpr);
		
		yrDevI = muYr - muLn;
		SigIyr.update(yrDev);
		SigIsc.update(sca, muSca);
		
		muRspI = muYr - sca;
		mu.update(muRsp, SigIyr, SigIpr);
		
		cout << "#" << flush;
	}
	cout << endl;
	
	MuGrp gmRspI = muYr - mu - sca;
	Grp &gmRsp   = gmRspI;
	
	BetaGrpPC gammaI(gmRsp, pcUfNam, evRfNam, Npc, yr2ln, gm2pr, nThr);
	Grp &gamma = gammaI;
	
	scaI.setApr(1e-6);
	
	MuGrp gsI = gamma + sca;
	Grp &gs   = gsI;
	
	SigmaI SigIa(gamma, 1.0, 4.0);
	Qgrp qG(Npc, nuG);
	
	cout << "...initialization done" << endl;
	cout << "Starting pre-SNP burnin..." << endl;
	
	/*
		Burnin-in
	*/
	for (int iBn = 0; iBn < floor(0.6*Nbnin); iBn++) {
		data.update(muRep, SigIdat);
		muRep.update(data, qDat, SigIdat, muYr, qR, SigIrep);
		
		datDevI = muRep - data;
		SigIdat.update(datDev, qDat);
		qDat.update(datDev, SigIdat);
		
		muYr.update(muRep, qR, SigIrep, muLn, SigIyr);
		
		repDevI = muYr - muRep;
		SigIrep.update(repDev, qR);
		qR.update(repDev, SigIrep);
		
		scaRspI = muYr - mu - gamma;
		sca.update(scaRsp, SigIyr, muSca, SigIsc);
		muSca.update(sca, SigIsc, SigIpr);
		
		gmRspI = muYr - mu - sca;
		gamma.update(gmRsp, SigIyr, qG, SigIa);
		
		SigIa.update(gamma, qG);
		qG.update(gamma, SigIa);
		
		gsI    = gamma + sca;
		muLnI  = mu + gs;
		yrDevI = muYr - muLn;
		SigIyr.update(yrDev);
		SigIsc.update(sca, muSca);
		
		muRspI = muYr - gs;
		mu.update(muRsp, SigIyr, SigIpr);
		
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
		snpBetI = new BetaGrpBVSR(yrDev, SigIyr, snpFnam, Nmul, ldCt, -9, yr2ln, snp2pr, betOut, nThr);
	}
	else
		if (model == "RS"){
			snpBetI = new BetaGrpFt(yrDev, SigIyr, snpFnam, Nsnp, Nmul, ldCt, -9, yr2ln, snp2pr, betOut, nThr);
		}
		else {
			snpBetI = new BetaGrpSnpMiss(snpFnam, betOut, yr2ln, Nsnp, d, nThr, -9);
		}
	
	Grp &snpBet = *snpBetI;
	cout << "                         ...done" << endl;
	
	
	SigmaI SigIbet(d, 1.0, 20.0);
	//SigmaI SigIbet(d, 1e-6, 20.0);
	Qgrp qB(Nln*Nmul, nuB);
	
	MuGrp snpDevI = yrDev;
	Grp &snpDev   = snpDevI;
	
	if (model == "SM") {
		cout << "Single marker model" << endl;
		cout << "Continuing burn-in as before..." << endl;
		for (int iBn = 0; iBn < ceil(0.4*Nbnin); iBn++) {
			data.update(muRep, SigIdat);
			muRep.update(data, qDat, SigIdat, muYr, qR, SigIrep);
			
			datDevI = muRep - data;
			SigIdat.update(datDev, qDat);
			qDat.update(datDev, SigIdat);
			
			muYr.update(muRep, qR, SigIrep, muLn, SigIyr);
			
			repDevI = muYr - muRep;
			SigIrep.update(repDev, qR);
			qR.update(repDev, SigIrep);
			
			scaRspI = muYr - mu - gamma;
			sca.update(scaRsp, SigIyr, muSca, SigIsc);
			muSca.update(sca, SigIsc, SigIpr);
			
			gmRspI = muYr - mu - sca;
			gamma.update(gmRsp, SigIyr, qG, SigIa);
			
			SigIa.update(gamma, qG);
			qG.update(gamma, SigIa);
			
			gsI    = gamma + sca;
			muLnI  = mu + gs;
			yrDevI = muYr - muLn;
			SigIyr.update(yrDev);
			SigIsc.update(sca, muSca);
			
			muRspI = muYr - gs;
			mu.update(muRsp, SigIyr, SigIpr);
			
			cout << "*" << flush;
			
		}
		cout << endl;
		
		cout << "Sampling..." << endl;
		for (int iSm = 0; iSm < Nsamp; iSm++) {
			data.update(muRep, SigIdat);
			muRep.update(data, qDat, SigIdat, muYr, qR, SigIrep);
			
			datDevI = muRep - data;
			SigIdat.update(datDev, qDat);
			qDat.update(datDev, SigIdat);
			
			muYr.update(muRep, qR, SigIrep, muLn, SigIyr);
			
			repDevI = muYr - muRep;
			SigIrep.update(repDev, qR);
			qR.update(repDev, SigIrep);
			
			scaRspI = muYr - mu - gamma;
			sca.update(scaRsp, SigIyr, muSca, SigIsc);
			muSca.update(sca, SigIsc, SigIpr);
			
			gmRspI = muYr - mu - sca;
			gamma.update(gmRsp, SigIyr, qG, SigIa);
			
			SigIa.update(gamma, qG);
			qG.update(gamma, SigIa);
			
			gsI    = gamma + sca;
			muLnI  = mu + gs;
			yrDevI = muYr - muLn;
			SigIyr.update(yrDev);
			SigIsc.update(sca, muSca);
			
			muRspI = muYr - gs;
			mu.update(muRsp, SigIyr, SigIpr);
			
			if ((iSm + 1) % Nthin) {
				cout << "." << flush;
			}
			else {
				snpBet.update(yrDev, SigIyr);
				muLn.save(LNout);
				SigIdat.save(SgDout);
				SigIrep.save(SgRout);
				SigIyr.save(SgYout);
				SigIsc.save(SgSout, scaI.getA());
				SigIa.save(SgAout);
				
				cout << "|" << flush;
			}
			
		}
		cout << endl;
		snpBet.dump();
	}
	else {
		cout << "Multi-marker model " << model << endl;
		cout << "Burnin with SNPs..." << endl;
		for (int iPbn = 0; iPbn < ceil(0.4*Nbnin); iPbn++) {
			data.update(muRep, SigIdat);
			muRep.update(data, qDat, SigIdat, muYr, qR, SigIrep);
			
			datDevI = muRep - data;
			SigIdat.update(datDev, qDat);
			qDat.update(datDev, SigIdat);
			
			muYr.update(muRep, qR, SigIrep, muLn, SigIyr);
			
			repDevI = muYr - muRep;
			SigIrep.update(repDev, qR);
			qR.update(repDev, SigIrep);
			
			//scaRspI = muYr - mu - gamma - snpBet;
			scaRspI = muYr - mu - gamma;
			sca.update(scaRsp, SigIyr, muSca, SigIsc);
			muSca.update(sca, SigIsc, SigIpr);
			
			//gmRspI = muYr - mu - sca - snpBet;
			gmRspI = muYr - mu - sca;
			gamma.update(gmRsp, SigIyr, qG, SigIa);
			
			SigIa.update(gamma, qG);
			qG.update(gamma, SigIa);
			
			//gsI    = gamma + sca + snpBet;
			gsI    = gamma + sca;
			muLnI  = mu + gs;
			yrDevI = muYr - muLn;
			SigIyr.update(yrDev);
			SigIsc.update(sca, muSca);
			
			snpDevI = muYr - mu - gamma - sca;
			if (model == "VS") {
				snpBet.update(snpDev, SigIyr, SigIbet);
				SigIbet.update(snpBet);
				kappa.update(snp2pr);
			}
			else{
				snpBet.update(yrDev, SigIyr, SigIbet);
				//snpBet.update(snpDev, SigIyr, SigIbet);
				//SigIbet.update(snpBet);
				//snpBet.update(snpDev, SigIyr, qB, SigIbet);
				//SigIbet.update(snpBet, qB);
				//qB.update(snpBet, SigIbet);
			}
			
			muRspI = muYr - gs;
			mu.update(muRsp, SigIyr, SigIpr);
			
			cout << "+" << flush;
		}
		cout << endl;
		
		cout << "Sampling..." << endl;
		for (int iSm = 0; iSm < Nsamp; iSm++) {
			data.update(muRep, SigIdat);
			muRep.update(data, qDat, SigIdat, muYr, qR, SigIrep);
			
			datDevI = muRep - data;
			SigIdat.update(datDev, qDat);
			qDat.update(datDev, SigIdat);
			
			muYr.update(muRep, qR, SigIrep, muLn, SigIyr);
			
			repDevI = muYr - muRep;
			SigIrep.update(repDev, qR);
			qR.update(repDev, SigIrep);
			
			//scaRspI = muYr - mu - gamma - snpBet;
			scaRspI = muYr - mu - gamma;
			sca.update(scaRsp, SigIyr, muSca, SigIsc);
			muSca.update(sca, SigIsc, SigIpr);
			
			//gmRspI = muYr - mu - sca - snpBet;
			gmRspI = muYr - mu - sca;
			gamma.update(gmRsp, SigIyr, qG, SigIa);
			
			SigIa.update(gamma, qG);
			qG.update(gamma, SigIa);
			
			//gsI    = gamma + sca + snpBet;
			gsI    = gamma + sca;
			muLnI  = mu + gs;
			yrDevI = muYr - muLn;
			SigIyr.update(yrDev);
			SigIsc.update(sca, muSca);
			
			snpDevI = muYr - mu - gamma - sca;
			if (model == "VS") {
				snpBet.update(snpDev, SigIyr, SigIbet);
				SigIbet.update(snpBet);
				kappa.update(snp2pr);
			}
			else{
				snpBet.update(yrDev, SigIyr, SigIbet);
				cout << "snpBet: " << gsl_matrix_get(snpBet.dMat(), 0, 0) << " " << gsl_matrix_get(snpBet.dMat(), 0, 1) << endl;
				//snpBet.update(snpDev, SigIyr, SigIbet);
				//SigIbet.update(snpBet);
				//snpBet.update(snpDev, SigIyr, qB, SigIbet);
				SigIbet.update(snpBet, qB);
				qB.update(snpBet, SigIbet);
			}
			
			muRspI = muYr - gs;
			mu.update(muRsp, SigIyr, SigIpr);
			
			if ((iSm + 1) % Nthin) {
				//cout << "." << flush;
			}
			else {
				if (model == "VS") {
					snpBet.save(snpDev, SigIyr);
				}
				else {
					snpBet.save(SigIyr);
					snp2pr.save(snpDev, snpBet, SigIyr);
				}
				
				muLn.save(LNout);
				SigIdat.save(SgDout);
				SigIrep.save(SgRout);
				SigIyr.save(SgYout);
				SigIsc.save(SgSout, scaI.getA());
				SigIa.save(SgAout);
				
				cout << "|" << flush;
			}
			
		}
		cout << endl;
		snpBet.dump();
		snp2pr.dump();
	}
	
	delete snpBetI;
}
