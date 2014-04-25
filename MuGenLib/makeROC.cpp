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

void ldTest(const gsl_vector *mySNP, const int &myCHR, const int &myPOS, const double &ldCt, const int &dist, const vector<gsl_vector *> &tstSNP, const vector<int> &tstCHR, const vector<int> &tstPOS, size_t &disc);

void ldTest(const gsl_vector *mySNP, const int &myCHR, const int &myPOS, const double &ldCt, const int &dist, const vector<gsl_vector *> &tstSNP, const vector<int> &tstCHR, const vector<int> &tstPOS, size_t &disc){
	for (size_t iFP = 0; iFP < tstSNP.size(); iFP++) {
		if ( (myCHR == tstCHR[iFP]) && (abs(myPOS - tstPOS[iFP]) > dist) ) {  // only need to test for LD if the choromosome and distance are right
			double rSq = gsl_pow_2(gsl_stats_correlation(mySNP->data, 1, (tstSNP[iFP])->data, 1, mySNP->size));
			if (rSq >= ldCt) {
				disc++;
				break;
			}
		}
		
	}

}

int main(int argc, char *argv[]){
	bool sOn = false;
	bool dOn = false;
	bool cOn = false;
	bool tOn = false;
	bool nOn = false;
	bool TOn = false; // transpose flag
	
	const size_t d    = 11;
	const size_t Nln  = 750;
	const size_t Ntru = 500;
	const size_t Nsnp = 669958;
	const int Dist    = 1e6;    // max distance to between SNPs that still allows them to be "linked"
	const size_t Nts  = 100;    // total number of simulations
	
	string snpFlNam("MESAsnp.gbin");
	vector<string>truSnpFlNam(d);
	string chromIDflNam("MESAchromID.gbin");
	string chrPosflNam("MESAposition.gbin");
	string numTrSnpFlNam("numTrueSNP.gbin");
	string betFlNam("betSNP");
	string betDirNam;
	string hitsDirNam;
	vector<string> hitsFlNam(d);
	
	string simNum;
	int simID;
	
	int Ntst       = 2000;
	int Nthr       = 4;
	double LDctOff = 0.05;
	bool trans     = false; // transpose = false means Nsnp x d beta file; otherwise d x Nsnp, as in simple regression output
	
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
						
					case 'T':
						dOn = true;
						break;
						
					case 'c':
						cOn = true;
						break;
						
					case 'h':
						cerr << "usage:\n"
						<< "makeROC [-n hit_or_miss_SNP_number] [-s simulation_number] [-d directory_with_SNP_scores] [-c LD_cutoff] [-t n_threads] [-T transpose?]" << endl;
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
						if (sOn) {
							sOn = false;
							simNum = pchar;
							simID  = atoi(pchar);
						}
					else {
						dOn = false;
						betDirNam = pchar;
					}
			}
				break;
		}
	}
	hitsDirNam = "hits" + betDirNam.substr(3);
	
	for (size_t iPhn = 0; iPhn < d - 1; iPhn++) {
		stringstream phnStrm;
		phnStrm << iPhn + 1;
		truSnpFlNam[iPhn] = "SNPtrueID/SNPtrue" + simNum + "_" + phnStrm.str() + ".gbin";
		hitsFlNam[iPhn]   = hitsDirNam + "/hits" + simNum + "_" + phnStrm.str() + ".tsv";
		remove(hitsFlNam[iPhn].c_str());
		phnStrm << flush;
	}
	truSnpFlNam[d -1] = "SNPtrueID/SNPtrue" + simNum + ".gbin";
	hitsFlNam[d - 1]   = hitsDirNam + "/hits" + simNum + ".tsv";
	remove(hitsFlNam[d-1].c_str());
	
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
	gsl_matrix_scale(snpScore, 1.0/chnNum);
	
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
	const size_t NtruAll = gsl_vector_int_get(nTr, simID - 1); // that's the total number of true SNPs for all phenotypes
	gsl_vector_int_free(nTr);
	
	cout << "Populating true position list..." << endl;
	vector< vector<size_t> > trueX(d);
	gsl_vector_int *curTruX = gsl_vector_int_alloc(Ntru);
	for (size_t iPhn = 0; iPhn < d - 1; iPhn++) {
		FILE *trIn = fopen(truSnpFlNam[iPhn].c_str(), "r");
		gsl_vector_int_fread(trIn, curTruX);
		fclose(trIn);
		
		for (size_t iTr = 0; iTr < curTruX->size; iTr++) {
			trueX[iPhn].push_back(gsl_vector_int_get(curTruX, iTr));
		}
	}
	gsl_vector_int_free(curTruX);
	
	// Now the "all-phenotype" true positions
	gsl_vector_int *allTruX = gsl_vector_int_alloc(NtruAll);
	FILE *allTrIn = fopen(truSnpFlNam[d-1].c_str(), "r");
	gsl_vector_int_fread(allTrIn, allTruX);
	fclose(allTrIn);
	for (size_t iTr = 0; iTr < allTruX->size; iTr++) {
		trueX[d-1].push_back(gsl_vector_int_get(allTruX, iTr));
	}
	gsl_vector_int_free(allTruX);
	
	vector<int> tossLD(d); // store here the number of SNPs eliminated b/c of LD
	
	cout << "Starting the processing by phenotype..." << endl;
#pragma omp parallel for num_threads(Nthr)
	for (size_t phI = 0; phI < d; phI++) {
		
		vector<size_t> hits;
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
			size_t nDiscard = 0;
			
			// figure out the LD of the current top SNP with the true ones, then sort that
			vector<double> crRsq;
			gsl_vector *tmpX = gsl_vector_alloc(Nln);
			for (vector<size_t>::iterator truIt = trSNPid.begin(); truIt != trSNPid.end(); ++truIt) {
				gsl_matrix_get_col(tmpX, snp, *truIt);
				crRsq.push_back(gsl_pow_2(gsl_stats_correlation(curSNP->data, 1, tmpX->data, 1, Nln)));
			}
			gsl_vector_free(tmpX);
			
			vector<size_t> crRsqSortInd(crRsq.size());
			gsl_sort_index(crRsqSortInd.data(), crRsq.data(), 1, crRsq.size());
			
			if (crRsq[crRsqSortInd.back()] < LDctOff) { // max LD does not meet cut-off: false positive; sort is in INCREASING order, so going backwards
				// see if it's in LD with some that are already on the FP list
				if (fpSNP.size()) {
					ldTest(curSNP, gsl_vector_int_get(chrID, *curCand), gsl_vector_int_get(chrPos, *curCand), LDctOff, Dist, fpSNP, fpCHR, fpPOS, nDiscard);
				}
				
				if (nDiscard) {  // tossing it if it's in LD with an already IDed FP (don't want to count those multiple times)
					curCand = candidates.erase(curCand);
					candidates.push_back(gsl_permutation_get(snpRank, edge)); // just appending the extra candidate to the end.  No need to check if it's in LD with any exsting FP, since it still might be a hit, and if not that LD will be checked
					edge++;
					tossLD[phI]++;
					continue;
				}
				else {
					hits.push_back(0);
					
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
				if ( (gsl_vector_int_get(chrID, *curCand) == trCHR[crRsqSortInd.back()]) && (abs(gsl_vector_int_get(chrPos, *curCand) - trPOS[crRsqSortInd.back()]) <= Dist) ) { // both on the same chromosome and close enough: true positive (hit)
					hits.push_back(*curCand);
					
					gsl_matrix_get_col(curSNP, snp, crRsqSortInd.back()); // switch over the curSNP to the true SNP to be tested for LD with other candidates
					int curCHRid = trCHR[crRsqSortInd.back()];
					int curPOSid = trPOS[crRsqSortInd.back()];
					
					int nDrp = 0;
					gsl_vector *tstX = gsl_vector_alloc(Nln);
					
					curCand = candidates.erase(curCand); // erase the current, test the rest
					while (curCand != candidates.end()) {
						if ( (curCHRid == gsl_vector_int_get(chrID, *curCand)) && (abs(curPOSid - gsl_vector_int_get(chrPos, *curCand)) <= Dist) ) {
							gsl_matrix_get_col(tstX, snp, *curCand);
							double rSq = gsl_pow_2(gsl_stats_correlation(curSNP->data, 1, tstX->data, 1, Nln));
							if (rSq >= LDctOff) {
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
						if (gsl_pow_2(gsl_stats_correlation(curSNP->data, 1, tstX->data, 1, Nln)) >= LDctOff) { // make sure the new one is not in LD with the current
							edge++;
						}
						else {
							candidates.push_back(ind);
							edge++;
							nDrp--;
						}
						
					}
					gsl_vector_free(tstX);
					trSNPid.erase(trSNPid.begin() + crRsqSortInd.back()); // erase the hit from the vector of trues so that it can't be hit more than once; base-0 counts in crRsqSortInd, so no need to subtract 1
					trCHR.erase(trCHR.begin() + crRsqSortInd.back());
					trPOS.erase(trPOS.begin() + crRsqSortInd.back());
					continue;
				}
				else {  // if high enough LD, but location wrong, look to see if lesser LD guys meet the criteria
					bool isFP = true;
					vector<size_t>::reverse_iterator ldIt = crRsqSortInd.rbegin(); // stepping backwards from the end
					ldIt++;
					for (; ldIt != crRsqSortInd.rend(); ++ldIt) {
						if (crRsq[*ldIt] < LDctOff) {
							break;
						}
						else {
							if ( (gsl_vector_int_get(chrID, *curCand) == trCHR[*ldIt]) && (abs(gsl_vector_int_get(chrPos, *curCand) - trPOS[*ldIt]) <= Dist) ) { // found a hit
								isFP = false;
								hits.push_back(*curCand);
								
								gsl_matrix_get_col(curSNP, snp, *ldIt); // switch over the curSNP to the true SNP to be tested for LD with other candidates
								int curCHRid = trCHR[*ldIt];
								int curPOSid = trPOS[*ldIt];
								
								int nDrp = 0;
								gsl_vector *tstX = gsl_vector_alloc(Nln);
								
								curCand = candidates.erase(curCand);  // erase the current, test the rest
								while (curCand != candidates.end()) {
									if ( (curCHRid == gsl_vector_int_get(chrID, *curCand)) && (abs(curPOSid - gsl_vector_int_get(chrPos, *curCand)) <= Dist) ) {
										gsl_matrix_get_col(tstX, snp, *curCand);
										double rSq = gsl_pow_2(gsl_stats_correlation(curSNP->data, 1, tstX->data, 1, Nln));
										if (rSq >= LDctOff) {
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
									if (gsl_pow_2(gsl_stats_correlation(curSNP->data, 1, tstX->data, 1, Nln)) >= LDctOff) { // make sure the new one is not in LD with the current
										edge++;
									}
									else {
										candidates.push_back(ind);
										edge++;
										nDrp--;
									}
									
								}
								gsl_vector_free(tstX);
								trSNPid.erase(trSNPid.begin() + (*ldIt)); // erase the hit from the vector of trues so that it can't be hit more than once; *ldIt has base-0 positions, so no need to subtract 1
								trCHR.erase(trCHR.begin() + (*ldIt));
								trPOS.erase(trPOS.begin() + (*ldIt));
								break;
							}
							else {
								continue;
							}
						}
					}
					if (isFP) {
						if (fpSNP.size()) {
							ldTest(curSNP, gsl_vector_int_get(chrID, *curCand), gsl_vector_int_get(chrPos, *curCand), LDctOff, Dist, fpSNP, fpCHR, fpPOS, nDiscard);
						}
						
						if (nDiscard) {  // tossing it if it's in LD with an already IDed FP (don't want to count those multiple times)
							curCand = candidates.erase(curCand);
							candidates.push_back(gsl_permutation_get(snpRank, edge)); // just appending the extra candidate to the end.  No need to check if it's in LD with any exsting FP, since it still might be a hit, and if not that LD will be checked
							edge++;
							tossLD[phI]++;
							continue;
						}
						else {
							hits.push_back(0);
							
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
				}
			}
		}
		
		ofstream outHits(hitsFlNam[phI].c_str());
		for (vector<size_t>::iterator hI = hits.begin(); hI != hits.end(); ++hI) {
			outHits << *hI << " " << flush;
		}
		
		outHits << endl;
		outHits.close();
		
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


