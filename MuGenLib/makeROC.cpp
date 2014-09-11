/*
*  makeROC.cpp
*  MuGenLib
*
*  Created by ajg67 on 2/6/13.
*   Copyright (c) 2013 SEELE. All rights reserved.

	Making ROC curves from simulated data analyses
*/

#include <omp.h>
#include <cstdio>
#include <cstddef>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <list>
#include <vector>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_sort_vector.h>
#include <gsl/gsl_sort.h>
#include <gsl/gsl_statistics.h>

using std::cout;
using std::cerr;
using std::endl;
using std::flush;
using std::ofstream;
using std::stringstream;
using std::string;
using std::list;
using std::vector;
using std::advance;

// test if a current miss hits another trait; if yes returns the index into the SNP ID (base 1), otherwise 0
size_t offHitID(const gsl_vector *curX, const int &curCHR, const int &curPOS, const size_t &curTrt, const gsl_matrix *xMat, const vector< vector<size_t> > &trueIND, const vector< vector<int> > &trueCHR, const vector< vector<int> > &truePOS, const double &ldCt, const int &distCt, const size_t &d, const size_t &mul);

size_t offHitID(const gsl_vector *curX, const int &curCHR, const int &curPOS, const size_t &curTrt, const gsl_matrix *xMat, const vector< vector<size_t> > &trueIND, const vector< vector<int> > &trueCHR, const vector< vector<int> > &truePOS, const double &ldCt, const int &distCt, const size_t &d, const size_t &mul){
	
	vector<size_t> htInd;
	vector<double> rSq;
	
	for (size_t iTrt = 0; iTrt < d - mul; iTrt++) { // stop comparison before reaching the multi-trait "trait"
		if (iTrt == curTrt) { // not interested in current trait
			continue;
		}
		
		for (size_t jX = 0; jX < trueIND[iTrt].size(); jX++) {
			if ( (trueCHR[iTrt][jX] == curCHR) && (abs(truePOS[iTrt][jX] - curPOS) <= distCt) ) {
				gsl_vector *tmpX = gsl_vector_alloc(curX->size);
				gsl_matrix_get_col(tmpX, xMat, trueIND[iTrt][jX]);
				double curRsq = gsl_pow_2(gsl_stats_correlation(curX->data, 1, tmpX->data, 1, curX->size));
				if (curRsq >= ldCt) {
					rSq.push_back(curRsq);
					htInd.push_back(trueIND[iTrt][jX]);
				}
				gsl_vector_free(tmpX);
			}
		}
		
	}
	if (rSq.size() > 1) {
		vector<size_t> rSqSortInd(rSq.size());
		gsl_sort_index(rSqSortInd.data(), rSq.data(), 1, rSq.size()); // will report the index with the greatest rSq, in case there are more than one
		
		return htInd[rSqSortInd.back()] + 1; // base-1 for R
	}
	else if (rSq.size() == 1){
		return htInd[0] + 1;
	}
	else {
		return 0;
	}
}

int main(int argc, char *argv[]){
	bool sOn = false;
	bool dOn = false;
	bool DOn = false;
	bool mOn = false;
	bool cOn = false;
	bool tOn = false;
	bool nOn = false;
	bool TOn = false; // transpose flag
	
	const size_t Nln  = 750;
	const size_t Nsnp = 669958;
	const int Dist    = 1e7;    // max distance between SNPs that still allows them to be "linked"
	const size_t Nts  = 100;    // total number of simulations
	const double LDct = 0.04;   // for testing LD to discard SNPs; the other cut-offs are for testing for positives
	
	string snpFlNam("MESAsnp.gbin");
	string chromIDflNam("MESAchromID.gbin");
	string chrPosflNam("MESAposition.gbin");
	string numTrSnpFlNam("numTrueSNP.gbin");
	string numEaSnpFlNam("numTrueSNPeach.gbin");
	string betFlNam("betSNP");
	string betDirNam;
	string hitsDirNam;
	
	string simNum;
	int simID;
	
	size_t d       = 11;
	int Ntst       = 2000;
	int Nthr       = 4;
	double LDctOff = 0.05;
	double LDall   = 0.20;
	size_t multi   = 0;      // is the last trait a multi-trait statistic? (yes == 1)
	bool trans     = false;  // transpose = false means Nsnp x d beta file; otherwise d x Nsnp, as in simple regression output
	
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
					case 't':
						tOn = true;
						break;
						
					case 'n':
						nOn = true;
						break;
						
					case 's':
						sOn = true;
						break;
						
					case 'd':
						dOn = true;
						break;
						
					case 'D':
						DOn = true;
						break;
						
					case 'm':
						mOn = true;
						multi = 1;
						break;
						
					case 'T':
						TOn = true;
						break;
						
					case 'c':
						cOn = true;
						break;
						
					case 'h':
						cerr << "usage:\n"
						<< "makeROC [-n hit_or_miss_SNP_number] [-s simulation_number] [-D directory_with_SNP_scores] [-d trait_number] [-m multi-trait test LD cut-off] [-c LD_cutoff] [-t n_threads] [-T transpose?]" << endl;
						exit(-1);
					default:
						cerr << "ERROR: inrecognized flag " << pchar[1] << endl;
						exit(-1);
				}
			}
				break;
				
			default:{
					if (nOn) {
						nOn   = false;
						Ntst = atoi(pchar);
					}
					else
						if (tOn) {
							tOn = false;
							Nthr = atoi(pchar);
						}
					else
						if (TOn) {
							TOn = false;
							trans = true;
						}
					else
						if (cOn) {
							cOn = false;
							LDctOff = atof(pchar);
						}
					else
						if (mOn) {
							mOn = false;
							LDall = atof(pchar);
						}
					else
						if (sOn) {
							sOn = false;
							simNum = pchar;
							simID  = atoi(pchar);
						}
					else
						if (dOn) {
							dOn = false;
							d = static_cast<size_t>(atoi(pchar));
						}
					else {
						DOn = false;
						betDirNam = pchar;
						cout << "dir name: " << betDirNam << endl;
					}
			}
				break;
		}
	}
	vector<double> LDctVec(d, LDctOff);
	if (multi) {
		LDctVec[d-1] = LDall;
	}
	vector<string> truSnpFlNam(d);
	vector<string> hitsFlNam(d);
	vector<string> offHitsFlNam(d); // will save the misses that hit other traits
	vector<string> truesFlNam(d);   // will dump the ID of the true SNP found
	hitsDirNam = "hits" + betDirNam.substr(3);  // adding everything from betDirNam but the first three letters (i.e., "bet")
	if (d < 1) {
		cerr << "ERROR: number of traits " << d << " is invalid" << endl;
		exit(-1);
	}
	if ((d == 1) && multi) {
		cerr << "ERROR: cannot have the only trait be multi-trait" << endl;
		exit(-1);
	}
	
	for (size_t iPhn = 0; iPhn < d - multi; iPhn++) {
		stringstream phnStrm;
		phnStrm << iPhn + 1;
		truSnpFlNam[iPhn] = "SNPtrueID/SNPtrue" + simNum + "_" + phnStrm.str() + ".gbin";
		hitsFlNam[iPhn]   = hitsDirNam + "/hits" + simNum + "_" + phnStrm.str() + ".tsv";
		remove(hitsFlNam[iPhn].c_str());
		offHitsFlNam[iPhn] = hitsDirNam + "/offHits" + simNum + "_" + phnStrm.str() + ".tsv";
		remove(offHitsFlNam[iPhn].c_str());
		truesFlNam[iPhn] = hitsDirNam + "/trueID" + simNum + "_" + phnStrm.str() + ".tsv";
		remove(truesFlNam[iPhn].c_str());
		phnStrm << flush;
	}
	if (multi) {  // no off-trait hits possible
		truSnpFlNam[d - 1] = "SNPtrueID/SNPtrue" + simNum + ".gbin";
		hitsFlNam[d - 1]   = hitsDirNam + "/hits" + simNum + ".tsv";
		remove(hitsFlNam[d - 1].c_str());
		truesFlNam[d - 1]   = hitsDirNam + "/trueID" + simNum + ".tsv";
		remove(truesFlNam[d - 1].c_str());
	}
	
	gsl_matrix *snpScore;
	gsl_matrix *curChnVl;
	if (trans) {
		cout << "Transposed score matrix..." << endl;
		snpScore = gsl_matrix_calloc(d, Nsnp);
		curChnVl = gsl_matrix_alloc(d, Nsnp);
	}
	else {
		snpScore = gsl_matrix_calloc(Nsnp, d);
		curChnVl = gsl_matrix_alloc(Nsnp, d);
	}
	betFlNam = betDirNam + "/" + betFlNam + simNum + "_";
	bool bFlExists = true;
	double chnNum = 0.0;
	while (bFlExists) {  // this set-up means that once we hit a non-exitent file, we don't look any more
		stringstream chnStrm;
		chnStrm << chnNum + 1.0;
		string curBetFlNam = betFlNam + chnStrm.str() + ".gbin";
		if (FILE *curBsnpFl = fopen(curBetFlNam.c_str(), "r")) {
			gsl_matrix_fread(curBsnpFl, curChnVl);
			fclose(curBsnpFl);
			chnNum += 1.0;
			gsl_matrix_add(snpScore, curChnVl);
			
		}
		else {
			bFlExists = false;
		}
	}
	if (chnNum == 0) {
		cerr << "No valid chains for SNP scores found!" << endl;
		exit(-1);
	}
	cout << "# of chains found: " << chnNum << endl;
	gsl_matrix_free(curChnVl);
	if (chnNum > 1.0) {
		gsl_matrix_scale(snpScore, 1.0/chnNum);
	}
	
	cout << "Opening SNP file..." << endl;
	gsl_matrix *snp = gsl_matrix_alloc(Nln, Nsnp);
	FILE *snpIn = fopen(snpFlNam.c_str(), "r");
	gsl_matrix_fread(snpIn, snp);
	fclose(snpIn);
	
	gsl_vector_int *chrID  = gsl_vector_int_alloc(Nsnp);
	gsl_vector_int *chrPos = gsl_vector_int_alloc(Nsnp);
	
	FILE *chrIDin = fopen(chromIDflNam.c_str(), "r");
	gsl_vector_int_fread(chrIDin, chrID);
	fclose(chrIDin);
	FILE *chrPSin = fopen(chrPosflNam.c_str(), "r");
	gsl_vector_int_fread(chrPSin, chrPos);
	fclose(chrPSin);
	
	gsl_vector_int *nTr = gsl_vector_int_alloc(Nts);
	FILE *trNMin = fopen(numTrSnpFlNam.c_str(), "r");
	gsl_vector_int_fread(trNMin, nTr);
	fclose(trNMin);
	const size_t NtruAll = gsl_vector_int_get(nTr, simID - 1); // that's the total number of true SNPs for all traits
	gsl_vector_int_free(nTr);
	
	gsl_matrix_int *nTrEa = gsl_matrix_int_alloc(100, 10);  // master file, so always 100 sims by 10 traits
	FILE *eaNMin = fopen(numEaSnpFlNam.c_str(), "r");
	gsl_matrix_int_fread(eaNMin, nTrEa);
	fclose(eaNMin);
	
	
	cout << "Populating true position list..." << endl;
	vector< vector<size_t> > trueX(d);
	vector< vector<int> > trueCHR(d);
	vector< vector<int> > truePOS(d);
	for (size_t iPhn = 0; iPhn < d - multi; iPhn++) {
		gsl_vector_int *curTruX = gsl_vector_int_alloc(gsl_matrix_int_get(nTrEa, simID - 1, iPhn)); // not all true vectors are of the same length
		FILE *trIn = fopen(truSnpFlNam[iPhn].c_str(), "r");
		gsl_vector_int_fread(trIn, curTruX);
		fclose(trIn);
		
		for (size_t iTr = 0; iTr < curTruX->size; iTr++) {
			trueX[iPhn].push_back(gsl_vector_int_get(curTruX, iTr));
			trueCHR[iPhn].push_back(gsl_vector_int_get(chrID, gsl_vector_int_get(curTruX, iTr)));
			truePOS[iPhn].push_back(gsl_vector_int_get(chrPos, gsl_vector_int_get(curTruX, iTr)));
		}
		
		gsl_vector_int_free(curTruX);
	}
	gsl_matrix_int_free(nTrEa);
	
	// Now the "all-trait" true positions
	if (multi) {
		gsl_vector_int *allTruX = gsl_vector_int_alloc(NtruAll);
		FILE *allTrIn = fopen(truSnpFlNam[d-1].c_str(), "r");
		gsl_vector_int_fread(allTrIn, allTruX);
		fclose(allTrIn);
		for (size_t iTr = 0; iTr < allTruX->size; iTr++) {
			trueX[d-1].push_back(gsl_vector_int_get(allTruX, iTr));
		}
		gsl_vector_int_free(allTruX);
	}
	
	vector<int> tossLD(d); // store here the number of SNPs eliminated b/c of LD
	
	cout << "Starting the processing by trait..." << endl;
#pragma omp parallel for num_threads(Nthr)
	for (size_t phI = 0; phI < d; phI++) {
		
		vector<size_t> hits;
		vector<size_t> offHits;
		vector<size_t> trueIDs;
		gsl_vector *betRow = gsl_vector_alloc(Nsnp); // it's a row in the SNP table, but could be a column in the score file
		if (trans) {
			gsl_matrix_get_row(betRow, snpScore, phI);
		}
		else {
			gsl_matrix_get_col(betRow, snpScore, phI);
		}
		
		gsl_permutation *snpRank = gsl_permutation_alloc(betRow->size);
		gsl_sort_vector_index(snpRank, betRow);
		gsl_permutation_reverse(snpRank); // sort is in the order of increase
		
		gsl_vector *curSNP = gsl_vector_alloc(Nln);
		size_t edge        = Ntst; // index of the permutation that is one beyond the current "edge" of the candidate picks
		list<size_t> candidates;
		
		for (size_t iEl = 0; iEl < Ntst; iEl++) {
			candidates.push_back(gsl_permutation_get(snpRank, iEl));
		}
		
		vector<gsl_vector *> fpSNP;  // store here the false positive SNPs; new false positives to be tested for LD with these
		vector<int> fpCHR;           // chromosome IDs of the FPs
		vector<int> fpPOS;           // chromosome positions of the FPs
		
		
		vector<size_t> trSNPid;
		vector<int> trCHR;
		vector<int> trPOS;
		for (size_t trxI = 0; trxI < trueX[phI].size(); trxI++) {
			trSNPid.push_back(trueX[phI][trxI]);
			trCHR.push_back(gsl_vector_int_get(chrID, trueX[phI][trxI]));
			trPOS.push_back(gsl_vector_int_get(chrPos, trueX[phI][trxI]));
		}
		
		while ((hits.size() < Ntst) && (edge < betRow->size) && (trSNPid.size())) { // either we get the number of data points we need, run out of candidate SNPs, or run out of true SNPs
			list<size_t>::iterator curCand = candidates.begin();
			gsl_matrix_get_col(curSNP, snp, *curCand);
			
			// figure out the LD of the current top SNP with the true ones, then sort that
			vector<double> crRsq;
			gsl_vector *tmpX = gsl_vector_alloc(Nln);
			for (size_t trI = 0; trI < trSNPid.size(); trI++) {
				// only the SNPs meeting closeness and same-chromosome criteria are assigned a non-0 rSq
				if ( (gsl_vector_int_get(chrID, *curCand) != trCHR[trI]) && abs(gsl_vector_int_get(chrPos, *curCand) - trPOS[trI]) > Dist ) {
					crRsq.push_back(0.0);
				}
				else {
					gsl_matrix_get_col(tmpX, snp, trSNPid[trI]);
					crRsq.push_back(gsl_pow_2(gsl_stats_correlation(curSNP->data, 1, tmpX->data, 1, Nln)));
				}
			}
			gsl_vector_free(tmpX);
			
			vector<size_t> crRsqSortInd(crRsq.size());
			gsl_sort_index(crRsqSortInd.data(), crRsq.data(), 1, crRsq.size());
			
			if (crRsq[crRsqSortInd.back()] < LDctVec[phI]) { // max LD does not meet cut-off: false positive; sort is in INCREASING order, so going backwards
				// see if it's in LD with some that are already on the FP list
				if (fpSNP.size()) {
					int nDrp = 0;
					size_t fpID = 0;
					for (size_t iFP = 0; iFP < fpSNP.size(); iFP++) {
						if ( (gsl_vector_int_get(chrID, *curCand) == fpCHR[iFP]) && (abs(gsl_vector_int_get(chrPos, *curCand) - fpPOS[iFP]) <= Dist) ) {  // only need to test for LD if the choromosome and distance are right
							double rSq = gsl_pow_2(gsl_stats_correlation(curSNP->data, 1, (fpSNP[iFP])->data, 1, Nln));
							if (rSq >= LDct) {
								curCand = candidates.erase(curCand);
								nDrp++;
								fpID = iFP;
								break;
							}
						}
						
					}
					if (nDrp) {
						while (nDrp) {
							size_t ind = gsl_permutation_get(snpRank, edge);
							gsl_vector *tstX = gsl_vector_alloc(Nln);
							gsl_matrix_get_col(tstX, snp, ind);
							if (gsl_pow_2(gsl_stats_correlation(curSNP->data, 1, tstX->data, 1, Nln)) >= LDct) { // make sure the new one is not in LD with the current positive
								edge++;
							}
							else if (gsl_pow_2(gsl_stats_correlation((fpSNP[fpID])->data, 1, tstX->data, 1, Nln)) >= LDct){ // make sure the new one is not in LD with the current false positive
								edge++;
							}
							else {
								candidates.push_back(ind);
								edge++;
								nDrp--;
							}
							gsl_vector_free(tstX);
						}
						
					}
					else {
						hits.push_back(0);
						offHits.push_back(offHitID(curSNP, gsl_vector_int_get(chrID, *curCand), gsl_vector_int_get(chrPos, *curCand), phI, snp, trueX, trueCHR, truePOS, LDctVec[phI], Dist, d, multi));
						
						size_t crSz = fpSNP.size() + 1;
						fpSNP.resize(crSz);
						fpSNP[crSz - 1] = gsl_vector_alloc(Nln);
						gsl_vector_memcpy(fpSNP[crSz - 1], curSNP);
						fpCHR.push_back(gsl_vector_int_get(chrID, *curCand));
						fpPOS.push_back(gsl_vector_int_get(chrPos, *curCand));
						
						curCand = candidates.erase(curCand);
						
						continue;
					}
				}
				else {
					hits.push_back(0);
					offHits.push_back(offHitID(curSNP, gsl_vector_int_get(chrID, *curCand), gsl_vector_int_get(chrPos, *curCand), phI, snp, trueX, trueCHR, truePOS, LDctVec[phI], Dist, d, multi));
					
					size_t crSz = fpSNP.size() + 1;
					fpSNP.resize(crSz);
					fpSNP[crSz - 1] = gsl_vector_alloc(Nln);
					gsl_vector_memcpy(fpSNP[crSz - 1], curSNP);
					fpCHR.push_back(gsl_vector_int_get(chrID, *curCand));
					fpPOS.push_back(gsl_vector_int_get(chrPos, *curCand));
					
					curCand = candidates.erase(curCand);
					
					continue;
				}
				
			}
			else { // in LD with at least one true SNP
				hits.push_back((*curCand) + 1); // make base-1 for R
				trueIDs.push_back(trSNPid[crRsqSortInd.back()] + 1);
				
				gsl_vector *crTrSNP = gsl_vector_alloc(Nln);  // the true SNP to be tested for LD with other candidates
				gsl_matrix_get_col(crTrSNP, snp, crRsqSortInd.back());
				int crTrCHRid = trCHR[crRsqSortInd.back()];
				int crTrPOSid = trPOS[crRsqSortInd.back()];
				int curCHRid  = gsl_vector_int_get(chrID, *curCand);
				int curPOSid  = gsl_vector_int_get(chrPos, *curCand);
				
				int nDrp = 0;
				gsl_vector *tstX = gsl_vector_alloc(Nln);
				
				curCand = candidates.erase(curCand); // erase the current, test the rest
				while (curCand != candidates.end()) {
					if ( (curCHRid == gsl_vector_int_get(chrID, *curCand)) && (abs(curPOSid - gsl_vector_int_get(chrPos, *curCand)) <= Dist) ) { // test for LD of the current positive with other candidates
						gsl_matrix_get_col(tstX, snp, *curCand);
						double rSq = gsl_pow_2(gsl_stats_correlation(curSNP->data, 1, tstX->data, 1, Nln));
						if (rSq >= LDct) {
							tossLD[phI]++;
							curCand = candidates.erase(curCand);
							nDrp++;
						}
						else {
							curCand++;
						}
					}
					else if ( (crTrCHRid == gsl_vector_int_get(chrID, *curCand)) && (abs(crTrPOSid - gsl_vector_int_get(chrPos, *curCand)) <= Dist) ){ // test for LD of the current true with other candidates
						gsl_matrix_get_col(tstX, snp, *curCand);
						double rSq = gsl_pow_2(gsl_stats_correlation(crTrSNP->data, 1, tstX->data, 1, Nln));
						if (rSq >= LDct) {
							tossLD[phI]++;
							curCand = candidates.erase(curCand);
							nDrp++;
						}
						else {
							curCand++;
						}
					}
					else {
						curCand++;
					}
				}
				while (nDrp) {
					size_t ind = gsl_permutation_get(snpRank, edge);
					gsl_matrix_get_col(tstX, snp, ind);
					if (gsl_pow_2(gsl_stats_correlation(curSNP->data, 1, tstX->data, 1, Nln)) >= LDct) { // make sure the new one is not in LD with the current positive
						edge++;
					}
					else if (gsl_pow_2(gsl_stats_correlation(crTrSNP->data, 1, tstX->data, 1, Nln)) >= LDct){ // make sure the new one is not in LD with the current true
						edge++;
					}
					else {
						candidates.push_back(ind);
						edge++;
						nDrp--;
					}
					
				}
				gsl_vector_free(crTrSNP);
				gsl_vector_free(tstX);
				trSNPid.erase(trSNPid.begin() + crRsqSortInd.back()); // erase the hit from the vector of trues so that it can't be hit more than once; base-0 counts in crRsqSortInd, so no need to subtract 1
				trCHR.erase(trCHR.begin() + crRsqSortInd.back());
				trPOS.erase(trPOS.begin() + crRsqSortInd.back());
			}
		}
		
		ofstream outHits(hitsFlNam[phI].c_str());
		for (vector<size_t>::iterator hI = hits.begin(); hI != hits.end(); ++hI) {
			outHits << *hI << " " << flush;
		}
		
		outHits << endl;
		outHits.close();
		
		ofstream outMtrues(truesFlNam[phI].c_str());
		for (vector<size_t>::iterator hI = trueIDs.begin(); hI != trueIDs.end(); ++hI) {
			outMtrues << *hI << " " << flush;
		}
		
		outMtrues << endl;
		outMtrues.close();
		
		ofstream outOffhts(offHitsFlNam[phI].c_str());
		for (vector<size_t>::iterator hI = offHits.begin(); hI != offHits.end(); ++hI) {
			outOffhts << *hI << " " << flush;
		}
		
		outOffhts << endl;
		outOffhts.close();
		
		gsl_permutation_free(snpRank);
		gsl_vector_free(betRow);
		gsl_vector_free(curSNP);
		
		for (vector<gsl_vector *>::iterator fpsIt = fpSNP.begin(); fpsIt != fpSNP.end(); ++fpsIt) {
			
			gsl_vector_free(*fpsIt);
		}
	}
	
	cout << "number of SNPs tossed for LD: " << flush;
	for (size_t iP = 0; iP < d; iP++) {
		cout << tossLD[iP] << " " << flush;
	}
	cout << endl;
	
	gsl_matrix_free(snpScore);
	gsl_matrix_free(snp);
	
	gsl_vector_int_free(chrID);
	gsl_vector_int_free(chrPos);
}


