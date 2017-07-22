/*
 * defs.h
 *
 *  Created on: Mar 27, 2017
 *      Author: mahmoud
 */

#ifndef DEFS_H_
#define DEFS_H_

#include "mltaln.h"

//int TLS commonAlloc1 = 0; // what is TLS ?
//int TLS commonAlloc2 = 0;
//int TLS **commonIP = NULL;
//int TLS **commonJP = NULL;
//int nthread = 1;
//int randomseed = 0;
//int parallelizationstrategy = BAATARI1;
extern int TLS commonAlloc1; // what is TLS ?
extern int TLS commonAlloc2;
extern int TLS **commonIP;
extern int TLS **commonJP;
extern int nthread;
extern int randomseed;
extern int parallelizationstrategy;

char modelname[500];
int njob, nlenmax;
int amino_n[0x100]; //256 chars
char amino_grp[0x100]; //256 chars
//int amino_dis[0x100][0x100];
extern int **amino_dis;
extern double **n_disLN;
//double amino_dis_consweight_multi[0x100][0x100];
extern double **amino_dis_consweight_multi;
extern int **n_dis;
extern int **n_disFFT;
extern double **n_dis_consweight_multi;
unsigned char amino[0x100]; //256 chars
double polarity[0x100]; //256 chars
double volume[0x100]; //256 chars
int ribosumdis[37][37];

int ppid;
double thrinter;
double fastathreshold;
int pslocal, ppslocal;
int constraint;
int divpairscore;
int fmodel; // 1-> fmodel 0->default -1->raw
int nblosum; // 45, 50, 62, 80
int kobetsubunkatsu; //kobetsubunkatsu = individual division
int bunkatsu;
extern int dorp; // arguments de shitei suruto, tbfast -> pairlocalalign no yobidashi de futsugou
int niter;
int contin;
int calledByXced;
int devide;
int scmtd;
int weight;
int utree;
int tbutree;
int refine;
int check;
double cut;
int cooling;
extern int trywarp;
int penalty, ppenalty, penaltyLN;
int penalty_dist, ppenalty_dist;
int RNApenalty, RNAppenalty;
int RNApenalty_ex, RNAppenalty_ex;
int penalty_ex, ppenalty_ex, penalty_exLN;
int penalty_EX, ppenalty_EX;
int penalty_OP, ppenalty_OP;
int penalty_shift, ppenalty_shift;
extern double penalty_shift_factor;
int RNAthr, RNApthr;
int offset, poffset, offsetLN, offsetFFT;
int scoremtx;
int TMorJTT;
char use_fft;
char force_fft;
int nevermemsave;
int fftscore;
int fftWinSize; //may be means fft window size
int fftThreshold;
int fftRepeatStop;
int fftNoAnchStop;
int divWinSize;
int divThreshold;
int disp;
extern int outgap;
char alg;
int cnst;
int mix;
int tbitr;
int tbweight;
int tbrweight;
int disopt;
int pamN;
int checkC;
double geta2;
int treemethod;
int kimuraR;
char *swopt;
int fftkeika;
int score_check;
int makedistmtx;
char *inputfile;
char *addfile;
extern int addprofile;
int rnakozo;
char rnaprediction;
extern int scoreout;
extern int spscoreout;
extern int outnumber;
extern int legacygapcost;
extern double minimumweight;
extern int nwildcard;

char *signalSM;
FILE *prep_g;
FILE *trap_g;
char **seq_g;
char **res_g;

extern double consweight_multi;
extern double consweight_rna;
extern char RNAscoremtx;

extern char TLS *newgapstr;

extern int nalphabets;
extern int nscoredalphabets;

extern double specificityconsideration;
extern int ndistclass;
extern int maxdistclass;

extern int gmsg;

extern double sueff_global;

double lenfaca, lenfacb, lenfacc, lenfacd;
int maxl, tsize;

void initglobalvariables();


#endif /* DEFS_H_ */
