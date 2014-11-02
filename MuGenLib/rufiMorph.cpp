//
//  rufiMorph.cpp
//  MuGenLib
//
//  Created by ajg67 on 10/26/14.
//  Copyright (c) 2014 SEELE. All rights reserved.
//

#include <vector>
#include <string>
#include <sstream>
#include <iostream>
#include <cmath>
#include <cctype>
#include <cstdio>
#include <cstddef>

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
	
	const size_t d   = 29;
	const size_t N   = 437;
	const size_t Nln = 220;
	const size_t Npc = 219;
	
	const int Nbnin  = atoi(argv[1]);
	const int Nsamp  = atoi(argv[2]);
	const int Nthin  = atoi(argv[3]);
	const int Nthr   = atoi(argv[4]);
	const double nuG = atof(argv[5]);
	
	const char *chnID = argv[6];
	
	string dfNam("rufiMorphTraits.gbin");
	string mMatFnam("rufiMatNAind.gbin");
	string mVecFnam("rufiTotNAind.gbin");
	string LNout("LNrufi");
	LNout = LNout + chnID + ".gbin";
	string BVout("BVrufi");
	BVout = BVout + chnID + ".gbin";
	string SgEout("SigErufi");
	SgEout = SgEout + chnID + ".gbin";
	
	RanIndex e2ln(N, Nln, string("rufiLnInd.gbin"));
	RanIndex e2mu(N);
	RanIndex ln2mu(Nln);
	RanIndex gm2pr(Npc);
	RanIndex mu2pr = RanIndex();
	
	MuGrpMiss datI(dfNam, mMatFnam, mVecFnam, e2ln, d);
	Grp &dat = datI;
	
	MuGrp muRspI = dat;
	Grp &muRsp   = muRspI;
	
	MuGrp muI(muRsp, e2mu, mu2pr);
	Grp &mu = muI;
	
	MuGrp gmRspI = dat - mu;
	Grp &gmRsp   = gmRspI;
	
	//BetaGrpPC gammaI(gmRsp, string("rufiGBS_EVC.gbin"), string("rufiGBS_EVL.gbin"), Npc, gm2pr, Nthr);
	BetaGrpPCpex gammaI(gmRsp, string("rufiGBS_EVC.gbin"), string("rufiGBS_EVL.gbin"), Npc, 1e-6, e2ln, gm2pr, Nthr);
	Grp &gamma = gammaI;
	
	//MuGrp muGmI(gamma, gm2pr, mu2pr);
	//Grp &muGm = muGmI;
	
	MuGrp bvI = mu + gamma;
	Grp &bv   = bvI;
	
	MuGrp scaDevI = dat - mu - bv;
	Grp &scaDev   = scaDevI;
	
	MuGrpPEX scaI(scaDev, e2ln, ln2mu, 1e-6, Nthr);
	Grp &sca = scaI;
	
	MuGrp muLnI = bv + sca;
	Grp &muLn   = muLnI;
	
	MuGrp datDevI = dat - muLn;
	Grp &datDev   = datDevI;
	
	SigmaI SigIe(datDev, SgEout, 1.0, 2.0);
	SigmaI SigIs(sca, 1.0, 2.0);
	SigmaI SigIa(gamma, 1.0, 2.0);
	SigmaI SigIpr(d, 0.000001);
	
	Qgrp qG(Npc, nuG);
	
	for (int iBn = 0; iBn < Nbnin; iBn++) {
		dat.update(muLn, SigIe);
		
		bvI     = mu + gamma;
		scaDevI = dat - bv;
		sca.update(scaDev, SigIe, SigIs);
		SigIs.update(sca);
		
		muLnI   = bv + sca;
		datDevI = dat - muLn;
		SigIe.update(datDev);
		
		gmRspI = dat - mu - sca;
		gamma.update(gmRsp, SigIe, qG, SigIa);
		//muGm.update(gamma, qG, SigIa, SigIpr);
		SigIa.update(gamma, qG);
		qG.update(gamma, SigIa);
		
		muRspI = dat - gamma - sca;
		mu.update(muRsp, SigIe, SigIpr);
		
		cout << "+" << flush;
		
	}
	cout << endl;
	
	for (int iSm = 0; iSm < Nsamp; iSm++) {
		dat.update(muLn, SigIe);
		
		bvI     = mu + gamma;
		scaDevI = dat - bv;
		sca.update(scaDev, SigIe, SigIs);
		SigIs.update(sca);
		
		muLnI   = bv + sca;
		datDevI = dat - muLn;
		SigIe.update(datDev);
		
		gmRspI = dat - mu - sca;
		gamma.update(gmRsp, SigIe, qG, SigIa);
		//muGm.update(gamma, qG, SigIa, SigIpr);
		SigIa.update(gamma, qG);
		qG.update(gamma, SigIa);
		
		muRspI = dat - gamma - sca;
		mu.update(muRsp, SigIe, SigIpr);
		
		if ((iSm + 1) % Nthin) {
			cout << "." << flush;
		}
		else {
			bv.save(BVout);
			muLn.save(LNout);
			cout << "|" << flush;
		}
		
	}
	cout << endl;
	
}



