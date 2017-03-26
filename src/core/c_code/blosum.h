/*
int locpenaltyd = -2400;
int locoffsetd = -60;
char locaminod[26] = "GASTPLIMVDNEQFYWKRHCXXX.-U";
char locaminod[] = "ARNDCQEGHILKMFPSTWYVBZX.-U";
char locgrpd[] = 
{
	0, 3, 2, 2, 5, 2, 2, 0, 3, 1, 1, 3, 1, 4, 0, 0, 0, 4, 4, 1, 2, 2,
	6, 6, 6, 6,
};
*/

#ifndef BLOSUM_H_
#define BLOSUM_H_

#define DEFAULTGOP_B -1530
#define DEFAULTGEP_B   -00
#define DEFAULTOFS_B  -123   /* +10 -- -50  teido ka ? */

void BLOSUMmtx( int n, double **matrix, double *freq, unsigned char *amino, char *amino_grp );

void extendedmtx( double **matrix, double *freq, unsigned char *amino, char *amino_grp );

#endif /* BLOSUM_H_ */
