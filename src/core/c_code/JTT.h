/*
 * JTT.h
 *
 *  Created on: Mar 27, 2017
 *      Author: mahmoud
 */

#ifndef JTT_H_
#define JTT_H_

#define DEFAULTGOP_J -1530
#define DEFAULTGEP_J   -00
#define DEFAULTOFS_J  -123  /* +10 -- -50  teido ka ? */
#define DEFAULTPAMN  200

void JTTmtx( double **rsr, double *freq, unsigned char locamino[26], char locgrp[26], int isTM );

#endif /* JTT_H_ */
