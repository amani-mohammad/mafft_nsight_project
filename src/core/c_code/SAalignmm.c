//SA alignment
#include "mltaln.h"
#include "dp.h"

#define DEBUG 0

//fill match array
static void match_calc( double *match, double **cpmx1, double **cpmx2, int i1, int lgth2, double **doublework, int **intwork, int initialize )
{
	int j, k, l;
//	double scarr[26];
	double **cpmxpd = doublework;
	int **cpmxpdn = intwork;
	int count = 0;
	double *scarr;
	scarr = calloc( nalphabets, sizeof( double ) );

	if( initialize )
	{
		for( j=0; j<lgth2; j++ )
		{
			count = 0;
			for( l=0; l<nalphabets; l++ )
			{
				if( cpmx2[l][j] )
				{
					cpmxpd[count][j] = cpmx2[l][j];
					cpmxpdn[count][j] = l;
					count++;
				}
			}
			cpmxpdn[count][j] = -1;
		}
	}

	for( l=0; l<nalphabets; l++ )
	{
		scarr[l] = 0.0;
		for( k=0; k<nalphabets; k++ )
			scarr[l] += n_dis[k][l] * cpmx1[k][i1]; //n_dis defined in defs.c
	}
	for( j=0; j<lgth2; j++ )
	{
		match[j] = 0;
		for( k=0; cpmxpdn[k][j] > -1;  k++ )
			match[j] += scarr[cpmxpdn[k][j]] * cpmxpd[k][j];
	} 
	free( scarr );
}

//updates mseq1, mseq2 and ijp values based on some calculations
static double Atracking( double *lasthorizontalw, double *lastverticalw, 
						char **seq1, char **seq2, 
                        char **mseq1, char **mseq2, 
                        double **cpmx1, double **cpmx2, 
                        int **ijp, int icyc, int jcyc )
{
	int i, j, k, l, iin, jin, ifi, jfi, lgth1, lgth2;
//	char gap[] = "-";
	char *gap;
	double wm;
	gap = newgapstr;
	lgth1 = strlen( seq1[0] );
	lgth2 = strlen( seq2[0] );

#if DEBUG
	for( i=0; i<lgth1; i++ ) 
	{
		fprintf( stderr, "lastverticalw[%d] = %f\n", i, lastverticalw[i] );
	}
#endif
 
	if( outgap == 1 )
		;
	else
	{
		wm = lastverticalw[0];
		for( i=0; i<lgth1; i++ )
		{
			if( lastverticalw[i] >= wm )
			{
				wm = lastverticalw[i];
				iin = i; jin = lgth2-1;
				ijp[lgth1][lgth2] = +( lgth1 - i );
			}
		}
		for( j=0; j<lgth2; j++ )
		{
			if( lasthorizontalw[j] >= wm )
			{
				wm = lasthorizontalw[j];
				iin = lgth1-1; jin = j;
				ijp[lgth1][lgth2] = -( lgth2 - j );
			}
		}
	}

    for( i=0; i<lgth1+1; i++ ) 
    {
        ijp[i][0] = i + 1;
    }
    for( j=0; j<lgth2+1; j++ ) 
    {
        ijp[0][j] = -( j + 1 );
    }

	for( i=0; i<icyc; i++ )
	{
		mseq1[i] += lgth1+lgth2;
		*mseq1[i] = 0;
	}
	for( j=0; j<jcyc; j++ )
	{
		mseq2[j] += lgth1+lgth2;
		*mseq2[j] = 0;
	}
	iin = lgth1; jin = lgth2;
	for( k=0; k<=lgth1+lgth2; k++ ) 
	{
		if( ijp[iin][jin] < 0 ) 
		{
			ifi = iin-1; jfi = jin+ijp[iin][jin];
		}
		else if( ijp[iin][jin] > 0 )
		{
			ifi = iin-ijp[iin][jin]; jfi = jin-1;
		}
		else
		{
			ifi = iin-1; jfi = jin-1;
		}
		l = iin - ifi;
		while( --l ) 
		{
			for( i=0; i<icyc; i++ )
				*--mseq1[i] = seq1[i][ifi+l];
			for( j=0; j<jcyc; j++ ) 
				*--mseq2[j] = *gap;
			k++;
		}
		l= jin - jfi;
		while( --l )
		{
			for( i=0; i<icyc; i++ ) 
				*--mseq1[i] = *gap;
			for( j=0; j<jcyc; j++ ) 
				*--mseq2[j] = seq2[j][jfi+l];
			k++;
		}
		if( iin <= 0 || jin <= 0 ) break;
		for( i=0; i<icyc; i++ ) 
			*--mseq1[i] = seq1[i][ifi];
		for( j=0; j<jcyc; j++ ) 
			*--mseq2[j] = seq2[j][jfi];
		k++;
		iin = ifi; jin = jfi;
	}
	return( 0.0 );
}

//apply specific algorithm to align seq1 and seq2
double Aalign( char **seq1, char **seq2, double *eff1, double *eff2, int icyc, int jcyc, int alloclen )
/* score no keisan no sai motokaraaru gap no atukai ni mondai ga aru */
{
	register int i, j; //register suggests to the compiler to allocate these variable in registers, if possible
	int lasti;                      /* outgap == 0 -> lgth1, outgap == 1 -> lgth1+1 */
	int lgth1, lgth2;
	int resultlen;
	double wm = 0.0;   /* int ?????? */
	double g;
	double x;
	static TLS double mi, *m;
	static TLS int **ijp;
	static TLS int mpi, *mp;
	static TLS double *currentw;
	static TLS double *previousw;
	static TLS double *match;
	static TLS double *initverticalw;    /* kufuu sureba iranai */
	static TLS double *lastverticalw;    /* kufuu sureba iranai */
	static TLS char **mseq1;
	static TLS char **mseq2;
	static TLS char **mseq;
	static TLS double **cpmx1;
	static TLS double **cpmx2;
	static TLS int **intwork;
	static TLS double **doublework;
	static TLS int orlgth1 = 0, orlgth2 = 0;

#if DEBUG
	fprintf( stderr, "eff in SA+++align\n" );
	for( i=0; i<icyc; i++ ) fprintf( stderr, "eff1[%d] = %f\n", i, eff1[i] );
#endif
	if( orlgth1 == 0 )
	{
		mseq1 = AllocateCharMtx( njob, 1 ); 
		mseq2 = AllocateCharMtx( njob, 1 ); /* by J. Thompson */
	}

	lgth1 = strlen( seq1[0] );
	lgth2 = strlen( seq2[0] );

	if( lgth1 > orlgth1 || lgth2 > orlgth2 )
	{
		int ll1, ll2;

		if( orlgth1 > 0 && orlgth2 > 0 )
		{
			FreeFloatVec( currentw );
			FreeFloatVec( previousw );
			FreeFloatVec( match );
			FreeFloatVec( initverticalw );
			FreeFloatVec( lastverticalw );

			FreeFloatVec( m );
			FreeIntVec( mp );

			FreeCharMtx( mseq );

			FreeFloatMtx( cpmx1 );
			FreeFloatMtx( cpmx2 );

			FreeFloatMtx( doublework );
			FreeIntMtx( intwork );
		}

		ll1 = MAX( (int)(1.1*lgth1), orlgth1 ) + 100;
		ll2 = MAX( (int)(1.1*lgth2), orlgth2 ) + 100;

		fprintf( stderr, "\ntrying to allocate (%d+%d)xn matrices ... ", ll1, ll2 );

		currentw = AllocateFloatVec( ll2+2 );
		previousw = AllocateFloatVec( ll2+2 );
		match = AllocateFloatVec( ll2+2 );

		initverticalw = AllocateFloatVec( ll1+2 );
		lastverticalw = AllocateFloatVec( ll1+2 );

		m = AllocateFloatVec( ll2+2 );
		mp = AllocateIntVec( ll2+2 );

		mseq = AllocateCharMtx( njob, ll1+ll2 );

		cpmx1 = AllocateFloatMtx( nalphabets, ll1+2 );
		cpmx2 = AllocateFloatMtx( nalphabets, ll2+2 );

		doublework = AllocateFloatMtx( nalphabets, MAX( ll1, ll2 )+2 ); 
		intwork = AllocateIntMtx( nalphabets, MAX( ll1, ll2 )+2 ); 

		fprintf( stderr, "succeeded\n" );

		orlgth1 = ll1;
		orlgth2 = ll2;
	}

	for( i=0; i<icyc; i++ ) mseq1[i] = mseq[i];
	for( j=0; j<jcyc; j++ ) mseq2[j] = mseq[icyc+j];


	if( orlgth1 > commonAlloc1 || orlgth2 > commonAlloc2 ) //commonAlloc1 and commonAlloc2 defined in defs.c and = 0
	{
		int ll1, ll2;

		if( commonAlloc1 && commonAlloc2 ) //if commonAlloc1 != 0 & commonAlloc2 != 0
		{
			FreeIntMtx( commonIP ); //commonIP is defined in defs.c
		}

		ll1 = MAX( orlgth1, commonAlloc1 );
		ll2 = MAX( orlgth2, commonAlloc2 );

		fprintf( stderr, "\n\ntrying to allocate %dx%d matrices ... ", ll1+1, ll2+1 );

		commonIP = AllocateIntMtx( ll1+10, ll2+10 );

		fprintf( stderr, "succeeded\n\n" );

		commonAlloc1 = ll1;
		commonAlloc2 = ll2;
	}
	ijp = commonIP;

	cpmx_calc( seq1, cpmx1, eff1, strlen( seq1[0] ), icyc ); //defined in tddis.c. fill in cpmx1 matrix based on amino_n and eff1 values.
	cpmx_calc( seq2, cpmx2, eff2, strlen( seq2[0] ), jcyc ); //fill cpmx2

	match_calc( initverticalw, cpmx2, cpmx1, 0, lgth1, doublework, intwork, 1 ); //defined here. fills match array based on other params values
	match_calc( currentw, cpmx1, cpmx2, 0, lgth2, doublework, intwork, 1 ); //fill currentw array

	if( outgap == 1 ) //defined in defs.c and = 1
	{
		for( i=1; i<lgth1+1; i++ )
		{
			initverticalw[i] += penalty * 0.5; //penalty defined in defs.h and set in constants.c
		}
		for( j=1; j<lgth2+1; j++ )
		{
			currentw[j] += penalty * 0.5;
		}
	}

	for( j=0; j<lgth2+1; ++j ) 
	{
		m[j] = currentw[j-1] + penalty * 0.5; mp[j] = 0;
	}

	lastverticalw[0] = currentw[lgth2-1];

	if( outgap ) lasti = lgth1+1; else lasti = lgth1;

	for( i=1; i<lasti; i++ )
	{

		doublencpy( previousw, currentw, lgth2+1 ); //defined in mltaln9.c. copy (lgth2+1) numbers from currentw to previousw
		previousw[0] = initverticalw[i-1];

		match_calc( currentw, cpmx1, cpmx2, i, lgth2, doublework, intwork, 0 ); //defined here. fill currentw array without initialization step
		currentw[0] = initverticalw[i];

		mi = previousw[0] + penalty * 0.5; mpi = 0;
		for( j=1; j<lgth2+1; j++ )
		{
			wm = previousw[j-1];
			ijp[i][j] = 0;

			g = penalty * 0.5;
			x = mi + g;
			if( x > wm )
			{
				wm = x;
				ijp[i][j] = -( j - mpi );
			}
			g = penalty * 0.5;
			x = previousw[j-1] + g;
			if( mi <= x )
			{
				mi = x;
				mpi = j-1;
			}

			g = penalty * 0.5;
			x = m[j] + g;
			if( x > wm )
			{
				wm = x;
				ijp[i][j] = +( i - mp[j] );
			}
			g = penalty * 0.5;
			x = previousw[j-1] + g;
			if( m[j] <= x )
			{
				m[j] = x;
				mp[j] = i-1;
			}
			currentw[j] += wm;
		}
		lastverticalw[i] = currentw[lgth2-1];
	} //I think this loop fills lastverticalw based on currentw values and some calculations
	/*
	fprintf( stderr, "\n" );
	for( i=0; i<icyc; i++ ) fprintf( stderr,"%s\n", seq1[i] );
	fprintf( stderr, "#####\n" );
	for( j=0; j<jcyc; j++ ) fprintf( stderr,"%s\n", seq2[j] );
	fprintf( stderr, "====>" );
	for( i=0; i<icyc; i++ ) strcpy( mseq1[i], seq1[i] );
	for( j=0; j<jcyc; j++ ) strcpy( mseq2[j], seq2[j] );
	*/
	Atracking( currentw, lastverticalw, seq1, seq2, mseq1, mseq2, cpmx1, cpmx2, ijp, icyc, jcyc ); //defined here.
	//updates mseq1, mseq2 and ijp values based on some calculations

	resultlen = strlen( mseq1[0] );
	if( alloclen < resultlen || resultlen > N ) //if resultlen > allocated length or max length of sequence
	{
		fprintf( stderr, "alloclen=%d, resultlen=%d, N=%d\n", alloclen, resultlen, N );
		ErrorExit( "LENGTH OVER!\n" );
	}

	//copy new sequences to original ones
	for( i=0; i<icyc; i++ ) strcpy( seq1[i], mseq1[i] );
	for( j=0; j<jcyc; j++ ) strcpy( seq2[j], mseq2[j] );
	/*
	fprintf( stderr, "\n" );
	for( i=0; i<icyc; i++ ) fprintf( stderr, "%s\n", mseq1[i] );
	fprintf( stderr, "#####\n" );
	for( j=0; j<jcyc; j++ ) fprintf( stderr, "%s\n", mseq2[j] );
	*/
	return( wm );
}
