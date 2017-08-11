//This file contains the code of pairwise alignment. It is so important and needs more study
#include "mltaln.h"

#define DEBUG 0
#define IODEBUG 0
#define SCOREOUT 0
#define SHISHAGONYU 0 // for debug

#define NODIST -9999

static char *whereispairalign;
static char *laraparams;
static char foldalignopt[1000];
static int stdout_align;
static int stdout_dist;
static int store_localhom;
static int store_dist;
static int nadd;
static int laste;
static int lastm;
static int lastsubopt;
static int lastonce;
static int usenaivescoreinsteadofalignmentscore;
static int specifictarget;

typedef struct _lastres
{
	int score;
	int start1;
	int start2;
	char *aln1;
	char *aln2;
} Lastres;

typedef struct _reg
{
	int start;
	int end;
} Reg;

typedef struct _aln
{
	int nreg;
	Reg *reg1;
	Reg *reg2;
} Aln;

typedef struct _lastresx
{
	int score;
	int naln;
	Aln *aln;
} Lastresx;

#ifdef enablemultithread //i need to know where it is defined?
typedef struct _jobtable
{
	int i;
	int j;
} Jobtable;

typedef struct _thread_arg
{
	int thread_no;
	int njob;
	Jobtable *jobpospt;
	char **name;
	char **seq;
	char **dseq;
	int *thereisxineachseq;
	LocalHom **localhomtable;
	double **distancemtx;
	double *selfscore;
	char ***bpp;
	Lastresx **lastresx;
	int alloclen;
	int *targetmap;
	pthread_mutex_t *mutex_counter;
	pthread_mutex_t *mutex_stdout;
} thread_arg_t;
#endif

typedef struct _lastcallthread_arg
{
	int nq, nd;
	char **dseq;
	char **qseq;
	Lastresx **lastresx;
#ifdef enablemultithread
	int thread_no;
	int *kshare;
	pthread_mutex_t *mutex;
#endif
} lastcallthread_arg_t;

static void t2u( char *seq ) //convert capital letters to small and t to u
{
	while( *seq )
	{
		if     ( *seq == 'A' ) *seq = 'a';
		else if( *seq == 'a' ) *seq = 'a';
		else if( *seq == 'T' ) *seq = 'u';
		else if( *seq == 't' ) *seq = 'u';
		else if( *seq == 'U' ) *seq = 'u';
		else if( *seq == 'u' ) *seq = 'u';
		else if( *seq == 'G' ) *seq = 'g';
		else if( *seq == 'g' ) *seq = 'g';
		else if( *seq == 'C' ) *seq = 'c';
		else if( *seq == 'c' ) *seq = 'c';
		else *seq = 'n';
		seq++;
	}
}

static int removex( char *d, char *m ) //create new sequence 'd' without x and return x count in the input sequence
{
	int val = 0;
	while( *m != 0 )
	{
		if( *m == 'X' || *m == 'x' ) 
		{
			m++;
			val++;
		}
		else 
		{
			*d++ = *m++;
		}
	}
	*d = 0;
	return( val );
}

//This method calculates some scores based on s1, s2 and amino_n values and stores them in localhompt
//It also updates some references in localhompt based on other lastresx values.
static void putlocalhom_last( char *s1, char *s2, LocalHom *localhompt, Lastresx *lastresx, char korh )
{
	char *pt1, *pt2;
	int naln, nreg;
	int iscore;
	int isumscore;
	int sumoverlap;
	LocalHom *tmppt = localhompt; //LocalHom is a structure defined in mltaln.h.
	LocalHom *tmppt2;
	LocalHom *localhompt0;
	Reg *rpt1, *rpt2; //Reg is a structure defined here.
	Aln *apt; //Aln is a structure defined here.
	int nlocalhom = 0;
	int len;

//	fprintf( stderr, "s1=%s\n", s1 );
//	fprintf( stderr, "s2=%s\n", s2 );


	naln = lastresx->naln;
	apt = lastresx->aln;

	if( naln == 0 ) return;
	while( naln-- )
	{
		rpt1 = apt->reg1;
		rpt2 = apt->reg2;
		nreg = apt->nreg;
		isumscore = 0;
		sumoverlap = 0;
		while( nreg-- )
		{
			if( nlocalhom++ > 0 )
			{
//				fprintf( stderr, "reallocating ...\n" );
				tmppt->next = (LocalHom *)calloc( 1, sizeof( LocalHom ) );
//				fprintf( stderr, "done\n" );
				tmppt = tmppt->next;
				tmppt->next = NULL;
			}
			tmppt->start1 = rpt1->start;
			tmppt->start2 = rpt2->start;
			tmppt->end1   = rpt1->end;
			tmppt->end2   = rpt2->end;
			tmppt->korh   = 'h';
			if( rpt1 == apt->reg1 ) localhompt0 = tmppt; // ?
	
//			fprintf( stderr, "in putlocalhom, reg1: %d-%d (nreg=%d)\n", rpt1->start, rpt1->end, lastresx->nreg );
//			fprintf( stderr, "in putlocalhom, reg2: %d-%d (nreg=%d)\n", rpt2->start, rpt2->end, lastresx->nreg );
	
			len = tmppt->end1 - tmppt->start1 + 1;
	
//			fprintf( stderr, "tmppt->start1=%d\n", tmppt->start1 );
//			fprintf( stderr, "tmppt->start2=%d\n", tmppt->start2 );

//			fprintf( stderr, "s1+tmppt->start1=%*.*s\n", len, len, s1+tmppt->start1 );
//			fprintf( stderr, "s2+tmppt->start2=%*.*s\n", len, len, s2+tmppt->start2 );
	
			pt1 = s1 + tmppt->start1;
			pt2 = s2 + tmppt->start2;
			iscore = 0;
			while( len-- )
			{
				iscore += n_dis[(int)amino_n[(unsigned char)*pt1++]][(int)amino_n[(unsigned char)*pt2++]]; // - offset はいらないかも
//				fprintf( stderr, "len=%d, %c-%c, iscore(0) = %d\n", len, *(pt1-1), *(pt2-1), iscore );
			}
	
			if( divpairscore ) //defined in defs.h.
			{
				tmppt->overlapaa   = tmppt->end2-tmppt->start2+1;
				tmppt->opt = (double)iscore / tmppt->overlapaa * 5.8 / 600;
			}
			else
			{
				isumscore += iscore;
				sumoverlap += tmppt->end2-tmppt->start2+1;
			}
			rpt1++;
			rpt2++;
		}
#if 0
		fprintf( stderr, "iscore (1)= %d\n", iscore );
		fprintf( stderr, "al1: %d - %d\n", start1, end1 );
		fprintf( stderr, "al2: %d - %d\n", start2, end2 );
#endif

		if( !divpairscore )
		{
			for( tmppt2=localhompt0; tmppt2; tmppt2=tmppt2->next )
			{
				tmppt2->overlapaa = sumoverlap;
				tmppt2->opt = (double)isumscore * 5.8 / ( 600 * sumoverlap );
//				fprintf( stderr, "tmpptr->opt = %f\n", tmppt->opt );
			}
		}
		apt++;
	}
}

static int countcomma( char *s ) //get count of comma in sequence s
{
	int v = 0;
	while( *s ) if( *s++ == ',' ) v++;
	return( v );
}

//read fold align file -> align mseq1 and mseq2 based on it -> calculate alignment score -> copy new aligned sequences to mseq1 and mseq2 -> return alignment score
static double recallpairfoldalign( char **mseq1, char **mseq2, int m1, int m2, int *of1pt, int *of2pt, int alloclen )
{
	static FILE *fp = NULL;
	double value;
	char *aln1;
	char *aln2;
	int of1tmp, of2tmp;

	if( fp == NULL )
	{
		fp = fopen( "_foldalignout", "r" );
		if( fp == NULL )
		{
			fprintf( stderr, "Cannot open _foldalignout\n" );
			exit( 1 );
		}
	}

	aln1 = calloc( alloclen, sizeof( char ) );
	aln2 = calloc( alloclen, sizeof( char ) );

	//fill aln1 and aln2 based on alignment values read from fold align file and processed using other args and some calcs.
	readpairfoldalign( fp, *mseq1, *mseq2, aln1, aln2, m1, m2, &of1tmp, &of2tmp, alloclen ); //defined in io.c.

	if( strstr( foldalignopt, "-global") ) //i think this choose between local or global alignment - yeaaaa right :D, this argument is set to '-global' from command
	{
		fprintf( stderr, "Calling G__align11\n" );
		//Calculates distance between mseq1 and mseq2 based on specific algo.
		value = G__align11( n_dis_consweight_multi, mseq1, mseq2, alloclen, outgap, outgap );  //defined in Galign11.c.
		*of1pt = 0;
		*of2pt = 0;
	}
	else
	{
		fprintf( stderr, "Calling L__align11\n" );
		//This methods finds the matchings between mseq1 and mseq2 and returns the score of matching them.
		//It is similar to G__align11 with some small differences.
		value = L__align11( n_dis_consweight_multi, 0.0, mseq1, mseq2, alloclen, of1pt, of2pt ); //defined in Lalign11.c.
	}

//	value = (double)naivepairscore11( *mseq1, *mseq2, penalty ); // nennnotame

	if( aln1[0] == 0 )
	{
		fprintf( stderr, "FOLDALIGN returned no alignment between %d and %d.  Sequence alignment is used instead.\n", m1+1, m2+1 );
	}
	else
	{
		strcpy( *mseq1, aln1 ); //copy alignment result to mseq
		strcpy( *mseq2, aln2 );
		*of1pt = of1tmp;
		*of2pt = of2tmp;
	}

//	value = naivepairscore11( *mseq1, *mseq2, penalty ); // v6.511 ha kore wo tsukau, global nomi dakara.

//	fclose( fp ); // saigo dake yatta houga yoi.

//	fprintf( stderr, "*mseq1 = %s\n", *mseq1 );
//	fprintf( stderr, "*mseq2 = %s\n", *mseq2 );


	free( aln1 );
	free( aln2 );

	return( value );
}

//update reg1 and reg2 data based on block info
static void block2reg( char *block, Reg *reg1, Reg *reg2, int start1, int start2 )
{
	Reg *rpt1, *rpt2;
	char *tpt, *npt;
	int pos1, pos2;
	int len, glen1, glen2;
	pos1 = start1;
	pos2 = start2;
	rpt1 = reg1;
	rpt2 = reg2;
	while( block )
	{
		block++;
//		fprintf( stderr, "block = %s\n", block );
		tpt = strchr( block, ':' );
		npt = strchr( block, ',' );
		if( !tpt || tpt > npt )
		{
			len = atoi( block );
			reg1->start = pos1;
			reg2->start = pos2;
			pos1 += len - 1;
			pos2 += len - 1;
			reg1->end = pos1;
			reg2->end = pos2;
//			fprintf( stderr, "in loop reg1: %d-%d\n", reg1->start, reg1->end );
//			fprintf( stderr, "in loop reg2: %d-%d\n", reg2->start, reg2->end );
			reg1++;
			reg2++;
		}
		else
		{
			sscanf( block, "%d:%d", &glen1, &glen2 );
			pos1 += glen1 + 1;
			pos2 += glen2 + 1;
		}
		block = npt;

	}
	reg1->start = reg1->end = reg2->start = reg2->end = -1;
	
	while( rpt1->start != -1 )
	{
//		fprintf( stderr, "reg1: %d-%d\n", rpt1->start, rpt1->end );
//		fprintf( stderr, "reg2: %d-%d\n", rpt2->start, rpt2->end );
		rpt1++;
		rpt2++;
	}
//	*apt1 = *apt2 = 0;
//	fprintf( stderr, "aln1 = %s\n", aln1 );
//	fprintf( stderr, "aln2 = %s\n", aln2 );
}


static void readlastresx_singleq( FILE *fp, int n1, int nameq, Lastresx **lastresx )
{
	char *gett;
	Aln *tmpaln;
	int prevnaln, naln, nreg;
#if 0
	int i, pstart, pend, end1, end2;
#endif
	int score, name1, start1, alnSize1, seqSize1;
	int        name2, start2, alnSize2, seqSize2;
	char strand1, strand2;
	int includeintoscore;
	gett = calloc( 10000, sizeof( char ) );

//	fprintf( stderr, "seq2[0] = %s\n", seq2[0] );
//	fprintf( stderr, "seq1[0] = %s\n", seq1[0] );

	while( 1 )
	{
		fgets( gett, 9999, fp );
		if( feof( fp ) ) break;
		if( gett[0] == '#' ) continue;
//		fprintf( stdout, "gett = %s\n", gett );
		if( gett[strlen(gett)-1] != '\n' )
		{
			fprintf( stderr, "Too long line?\n" );
			exit( 1 );
		}

		sscanf( gett, "%d %d %d %d %c %d %d %d %d %c %d", 
					&score, &name1, &start1, &alnSize1, &strand1, &seqSize1,
					        &name2, &start2, &alnSize2, &strand2, &seqSize2 );

		if( alg == 'R' && name2 <= name1 ) continue;
		if( name2 != nameq )
		{
			fprintf( stderr, "BUG!!!\n" );
			exit( 1 );
		}

//		if( lastresx[name1][name2].score ) continue; // dame!!!!


		prevnaln = lastresx[name1][name2].naln;
#if 0
		for( i=0; i<prevnaln; i++ )
		{
			nreg = lastresx[name1][name2].aln[i].nreg;

			pstart = lastresx[name1][name2].aln[i].reg1[0].start + 0;
			pend   = lastresx[name1][name2].aln[i].reg1[nreg-1].end - 0;
			end1 = start1 + alnSize1;
//			fprintf( stderr, "pstart = %d, pend = %d\n", pstart, pend );
			if( pstart <= start1 && start1 <= pend && pend - start1 > 1 ) break;
			if( pstart <= end1   && end1   <= pend && end1 - pstart > 1 ) break;

			pstart = lastresx[name1][name2].aln[i].reg2[0].start + 0;
			pend   = lastresx[name1][name2].aln[i].reg2[nreg-1].end - 0;
			end2 = start2 + alnSize2;
//			fprintf( stderr, "pstart = %d, pend = %d\n", pstart, pend );
			if( pstart <= start2 && start2 <= pend && pend - start2 > 1 ) break;
			if( pstart <= end2   && end2   <= pend && end2 - pstart > 1 ) break;
		}
		includeintoscore = ( i == prevnaln );
#else
		if( prevnaln ) includeintoscore = 0;
		else includeintoscore = 1;
#endif
		if( !includeintoscore && !lastsubopt )
			continue;

		naln = prevnaln + 1;
		lastresx[name1][name2].naln = naln;
//		fprintf( stderr, "OK! add this alignment to hat3, %d-%d, naln = %d->%d\n", name1, name2, prevnaln, naln );

		if( ( tmpaln = (Aln *)realloc( lastresx[name1][name2].aln, (naln) * sizeof( Aln ) ) ) == NULL ) // yoyu nashi
		{
			fprintf( stderr, "Cannot reallocate lastresx[][].aln\n" );
			exit( 1 );
		}
		else
			lastresx[name1][name2].aln = tmpaln;

		nreg = countcomma( gett )/2 + 1;
		lastresx[name1][name2].aln[prevnaln].nreg = nreg;
//		lastresx[name1][name2].aln[naln].nreg = -1;
//		lastresx[name1][name2].aln[naln].reg1 = NULL;
//		lastresx[name1][name2].aln[naln].reg2 = NULL;
//		fprintf( stderr, "name1=%d, name2=%d, nreg=%d, prevnaln=%d\n", name1, name2, nreg, prevnaln );

		if( ( lastresx[name1][name2].aln[prevnaln].reg1 = (Reg *)calloc( nreg+1, sizeof( Reg ) ) ) == NULL ) // yoyu nashi
		{
			fprintf( stderr, "Cannot reallocate lastresx[][].reg2\n" );
			exit( 1 );
		}

		if( ( lastresx[name1][name2].aln[prevnaln].reg2 = (Reg *)calloc( nreg+1, sizeof( Reg ) ) ) == NULL ) // yoyu nashi
		{
			fprintf( stderr, "Cannot reallocate lastresx[][].reg2\n" );
			exit( 1 );
		}

//		lastresx[name1][name2].aln[prevnaln].reg1[0].start = -1; // iranai?
//		lastresx[name1][name2].aln[prevnaln].reg2[0].start = -1; // iranai?
		block2reg( strrchr( gett, '\t' ), lastresx[name1][name2].aln[prevnaln].reg1, lastresx[name1][name2].aln[prevnaln].reg2, start1, start2 );

		if( includeintoscore )
		{
			if( lastresx[name1][name2].score ) score += penalty;
			lastresx[name1][name2].score += score;
		}

//		fprintf( stderr, "score(%d,%d) = %d\n", name1, name2, lastresx[name1][name2].score );
	}
	free( gett );
}

#ifdef enablemultithread
#if 0
static void readlastresx_group( FILE *fp, Lastresx **lastresx )
{
	char *gett;
	Aln *tmpaln;
	int prevnaln, naln, nreg;
#if 0
	int i, pstart, pend, end1, end2;
#endif
	int score, name1, start1, alnSize1, seqSize1;
	int        name2, start2, alnSize2, seqSize2;
	char strand1, strand2;
	int includeintoscore;
	gett = calloc( 10000, sizeof( char ) );

//	fprintf( stderr, "seq2[0] = %s\n", seq2[0] );
//	fprintf( stderr, "seq1[0] = %s\n", seq1[0] );

	while( 1 )
	{
		fgets( gett, 9999, fp );
		if( feof( fp ) ) break;
		if( gett[0] == '#' ) continue;
//		fprintf( stdout, "gett = %s\n", gett );
		if( gett[strlen(gett)-1] != '\n' )
		{
			fprintf( stderr, "Too long line?\n" );
			exit( 1 );
		}

		sscanf( gett, "%d %d %d %d %c %d %d %d %d %c %d", 
					&score, &name1, &start1, &alnSize1, &strand1, &seqSize1,
					        &name2, &start2, &alnSize2, &strand2, &seqSize2 );

		if( alg == 'R' && name2 <= name1 ) continue;

//		if( lastresx[name1][name2].score ) continue; // dame!!!!

		prevnaln = lastresx[name1][name2].naln;
#if 0
		for( i=0; i<prevnaln; i++ )
		{
			nreg = lastresx[name1][name2].aln[i].nreg;

			pstart = lastresx[name1][name2].aln[i].reg1[0].start + 0;
			pend   = lastresx[name1][name2].aln[i].reg1[nreg-1].end - 0;
			end1 = start1 + alnSize1;
//			fprintf( stderr, "pstart = %d, pend = %d\n", pstart, pend );
			if( pstart <= start1 && start1 <= pend && pend - start1 > 3 ) break;
			if( pstart <= end1   && end1   <= pend && end1 - pstart > 3 ) break;

			pstart = lastresx[name1][name2].aln[i].reg2[0].start + 0;
			pend   = lastresx[name1][name2].aln[i].reg2[nreg-1].end - 0;
			end2 = start2 + alnSize2;
//			fprintf( stderr, "pstart = %d, pend = %d\n", pstart, pend );
			if( pstart <= start2 && start2 <= pend && pend - start2 > 3 ) break;
			if( pstart <= end2   && end2   <= pend && end2 - pstart > 3 ) break;
		}
		includeintoscore = ( i == prevnaln );
#else
		if( prevnaln ) includeintoscore = 0;
		else includeintoscore = 1;
#endif
		if( !includeintoscore && !lastsubopt )
			continue;

		naln = prevnaln + 1;
		lastresx[name1][name2].naln = naln;
//		fprintf( stderr, "OK! add this alignment to hat3, %d-%d, naln = %d->%d\n", name1, name2, prevnaln, naln );

		if( ( tmpaln = (Aln *)realloc( lastresx[name1][name2].aln, (naln) * sizeof( Aln ) ) ) == NULL ) // yoyu nashi
		{
			fprintf( stderr, "Cannot reallocate lastresx[][].aln\n" );
			exit( 1 );
		}
		else
			lastresx[name1][name2].aln = tmpaln;



		nreg = countcomma( gett )/2 + 1;
		lastresx[name1][name2].aln[prevnaln].nreg = nreg;
//		lastresx[name1][name2].aln[naln].nreg = -1;
//		lastresx[name1][name2].aln[naln].reg1 = NULL;
//		lastresx[name1][name2].aln[naln].reg2 = NULL;
//		fprintf( stderr, "name1=%d, name2=%d, nreg=%d, prevnaln=%d\n", name1, name2, nreg, prevnaln );

		if( ( lastresx[name1][name2].aln[prevnaln].reg1 = (Reg *)calloc( nreg+1, sizeof( Reg ) ) ) == NULL ) // yoyu nashi
		{
			fprintf( stderr, "Cannot reallocate lastresx[][].reg2\n" );
			exit( 1 );
		}

		if( ( lastresx[name1][name2].aln[prevnaln].reg2 = (Reg *)calloc( nreg+1, sizeof( Reg ) ) ) == NULL ) // yoyu nashi
		{
			fprintf( stderr, "Cannot reallocate lastresx[][].reg2\n" );
			exit( 1 );
		}

//		lastresx[name1][name2].aln[prevnaln].reg1[0].start = -1; // iranai?
//		lastresx[name1][name2].aln[prevnaln].reg2[0].start = -1; // iranai?
		block2reg( strrchr( gett, '\t' ), lastresx[name1][name2].aln[prevnaln].reg1, lastresx[name1][name2].aln[prevnaln].reg2, start1, start2 );

		if( includeintoscore )
		{
			if( lastresx[name1][name2].score ) score += penalty;
			lastresx[name1][name2].score += score;
		}

//		fprintf( stderr, "score(%d,%d) = %d\n", name1, name2, lastresx[name1][name2].score );
	}
	free( gett );
}
#endif
#endif

static void readlastresx( FILE *fp, int n1, int n2, Lastresx **lastresx, char **seq1, char **seq2 )
{
	char *gett;
	Aln *tmpaln;
	int prevnaln, naln, nreg;
#if 0
	int i, pstart, pend, end1, end2;
#endif
	int score, name1, start1, alnSize1, seqSize1;
	int        name2, start2, alnSize2, seqSize2;
	char strand1, strand2;
	int includeintoscore;
	gett = calloc( 10000, sizeof( char ) );

//	fprintf( stderr, "seq2[0] = %s\n", seq2[0] );
//	fprintf( stderr, "seq1[0] = %s\n", seq1[0] );

	while( 1 )
	{
		fgets( gett, 9999, fp ); //read line from fp file into gett
		if( feof( fp ) ) break;
		if( gett[0] == '#' ) continue; //if comment, go to next iteration
//		fprintf( stdout, "gett = %s\n", gett );
		if( gett[strlen(gett)-1] != '\n' ) //if line is longer than max chars allowed, exit
		{
			fprintf( stderr, "Too long line?\n" );
			exit( 1 );
		}

		sscanf( gett, "%d %d %d %d %c %d %d %d %d %c %d", 
					&score, &name1, &start1, &alnSize1, &strand1, &seqSize1,
					        &name2, &start2, &alnSize2, &strand2, &seqSize2 ); //read all these values from line read

		if( alg == 'R' && name2 <= name1 ) continue;

//		if( lastresx[name1][name2].score ) continue; // dame!!!!

		prevnaln = lastresx[name1][name2].naln;
#if 0
		for( i=0; i<prevnaln; i++ )
		{
			nreg = lastresx[name1][name2].aln[i].nreg;

			pstart = lastresx[name1][name2].aln[i].reg1[0].start + 0;
			pend   = lastresx[name1][name2].aln[i].reg1[nreg-1].end - 0;
			end1 = start1 + alnSize1;
//			fprintf( stderr, "pstart = %d, pend = %d\n", pstart, pend );
			if( pstart <= start1 && start1 <= pend && pend - start1 > 3 ) break;
			if( pstart <= end1   && end1   <= pend && end1 - pstart > 3 ) break;

			pstart = lastresx[name1][name2].aln[i].reg2[0].start + 0;
			pend   = lastresx[name1][name2].aln[i].reg2[nreg-1].end - 0;
			end2 = start2 + alnSize2;
//			fprintf( stderr, "pstart = %d, pend = %d\n", pstart, pend );
			if( pstart <= start2 && start2 <= pend && pend - start2 > 3 ) break;
			if( pstart <= end2   && end2   <= pend && end2 - pstart > 3 ) break;
		}
		includeintoscore = ( i == prevnaln );
#else
		if( prevnaln ) includeintoscore = 0;
		else includeintoscore = 1;
#endif
		if( !includeintoscore && !lastsubopt )
			continue;

		naln = prevnaln + 1;
		lastresx[name1][name2].naln = naln;
//		fprintf( stderr, "OK! add this alignment to hat3, %d-%d, naln = %d->%d\n", name1, name2, prevnaln, naln );

		if( ( tmpaln = (Aln *)realloc( lastresx[name1][name2].aln, (naln) * sizeof( Aln ) ) ) == NULL ) // yoyu nashi
		{
			fprintf( stderr, "Cannot reallocate lastresx[][].aln\n" );
			exit( 1 );
		}
		else
			lastresx[name1][name2].aln = tmpaln;



		nreg = countcomma( gett )/2 + 1;
		lastresx[name1][name2].aln[prevnaln].nreg = nreg;
//		lastresx[name1][name2].aln[naln].nreg = -1;
//		lastresx[name1][name2].aln[naln].reg1 = NULL;
//		lastresx[name1][name2].aln[naln].reg2 = NULL;
//		fprintf( stderr, "name1=%d, name2=%d, nreg=%d, prevnaln=%d\n", name1, name2, nreg, prevnaln );

		if( ( lastresx[name1][name2].aln[prevnaln].reg1 = (Reg *)calloc( nreg+1, sizeof( Reg ) ) ) == NULL ) // yoyu nashi
		{
			fprintf( stderr, "Cannot reallocate lastresx[][].reg2\n" );
			exit( 1 );
		}

		if( ( lastresx[name1][name2].aln[prevnaln].reg2 = (Reg *)calloc( nreg+1, sizeof( Reg ) ) ) == NULL ) // yoyu nashi
		{
			fprintf( stderr, "Cannot reallocate lastresx[][].reg2\n" );
			exit( 1 );
		}

//		lastresx[name1][name2].aln[prevnaln].reg1[0].start = -1; // iranai?
//		lastresx[name1][name2].aln[prevnaln].reg2[0].start = -1; // iranai?
		//defined here. update reg1 and reg2 data based on first parameter - block -.
		block2reg( strrchr( gett, '\t' ), lastresx[name1][name2].aln[prevnaln].reg1, lastresx[name1][name2].aln[prevnaln].reg2, start1, start2 );

		if( includeintoscore )
		{
			if( lastresx[name1][name2].score ) score += penalty;
			lastresx[name1][name2].score += score;
		}

//		fprintf( stderr, "score(%d,%d) = %d\n", name1, name2, lastresx[name1][name2].score );
	}
	free( gett );
}

#ifdef enablemultithread
#if 0
static void *lastcallthread_group( void *arg )
{
	lastcallthread_arg_t *targ = (lastcallthread_arg_t *)arg;
	int k, i;
	int nq = targ->nq;
	int nd = targ->nd;
#ifdef enablemultithread
	int thread_no = targ->thread_no;
	int *kshare = targ->kshare; 
#endif
	Lastresx **lastresx = targ->lastresx;
	char **dseq = targ->dseq;
	char **qseq = targ->qseq;
	char command[5000];
	FILE *lfp;
	int msize;
	int klim;
	int qstart, qend, shou, amari;
	char kd[1000];

	if( nthread )
	{
		shou = nq / nthread;
		amari = nq - shou * nthread;
		fprintf( stderr, "shou: %d, amari: %d\n", shou, amari );

		qstart = thread_no * shou;
		if( thread_no - 1 < amari ) qstart += thread_no;
		else qstart += amari;

		qend = qstart + shou - 1;
		if( thread_no < amari ) qend += 1;
		fprintf( stderr, "%d: %d-%d\n", thread_no, qstart, qend );
	}
	k = -1;
	while( 1 )
	{
		if( nthread )
		{
			if( qstart > qend ) break;
			if( k == thread_no ) break;
			fprintf( stderr, "\n%d-%d / %d (thread %d)                    \n", qstart, qend, nq, thread_no );
			k = thread_no;
		}
		else
		{
			k++;
			if( k == nq ) break;
			fprintf( stderr, "\r%d / %d                    \r", k, nq );
		}

		if( alg == 'R' ) // if 'r' -> calllast_fast
		{
			fprintf( stderr, "Not supported\n" );
			exit( 1 );
		}
		else // 'r'
		{
			kd[0] = 0;
		}
		
		sprintf( command, "_q%d", k );
		lfp = fopen( command, "w" );
		if( !lfp )
		{
			fprintf( stderr, "Cannot open %s", command );
			exit( 1 );
		}
		for( i=qstart; i<=qend; i++ )
			fprintf( lfp, ">%d\n%s\n", i, qseq[i] );
		fclose( lfp );
	
//		if( alg == 'R' ) msize = MAX(10,k+nq);
//			else msize = MAX(10,nd+nq);
		if( alg == 'R' ) msize = MAX(10,k*lastm);
			else msize = MAX(10,nd*lastm);

//		fprintf( stderr, "Calling lastal from lastcallthread, msize = %d, k=%d\n", msize, k );
//		sprintf( command, "grep '>' _db%sd", kd );
//		system( command );
		sprintf( command, "%s/lastal -m %d -e %d -f 0 -s 1 -p _scoringmatrixforlast -a %d -b %d _db%sd _q%d > _lastres%d", whereispairalign, msize, laste, -penalty, -penalty_ex, kd, k, k );
		if( system( command ) ) exit( 1 );
	
		sprintf( command, "_lastres%d", k );
		lfp = fopen( command, "r" );
		if( !lfp )
		{
			fprintf( stderr, "Cannot read _lastres%d", k );
			exit( 1 );
		}
//		readlastres( lfp, nd, nq, lastres, dseq, qseq );
//		fprintf( stderr, "Reading lastres\n" );
		readlastresx_group( lfp, lastresx );
		fclose( lfp );
	}
	return( NULL );
}
#endif
#endif

static void *lastcallthread( void *arg )
{
	lastcallthread_arg_t *targ = (lastcallthread_arg_t *)arg;
	int k, i;
	int nq = targ->nq;
	int nd = targ->nd;
#ifdef enablemultithread
	int thread_no = targ->thread_no;
	int *kshare = targ->kshare; 
#endif
	Lastresx **lastresx = targ->lastresx;
	char **dseq = targ->dseq;
	char **qseq = targ->qseq;
	char command[5000];
	FILE *lfp;
	int msize;
	int klim;
	char kd[1000];

	k = -1;
	while( 1 )
	{

#ifdef enablemultithread
		if( nthread )
		{
			pthread_mutex_lock( targ->mutex );
			k = *kshare;
			if( k == nq )
			{
				pthread_mutex_unlock( targ->mutex );
				break;
			}
			fprintf( stderr, "\r%d / %d (thread %d)                    \r", k, nq, thread_no );
			++(*kshare);
			pthread_mutex_unlock( targ->mutex );
		}
		else
#endif
		{
			k++;
			if( k == nq ) break;
			fprintf( stderr, "\r%d / %d                    \r", k, nq );
		}

		if( alg == 'R' ) // if 'r' -> calllast_fast
		{
			klim = MIN( k, njob-nadd );
//			klim = k; // dochira demo yoi
			if( klim == k ) 
			{
				sprintf( command, "_db%dd", k );
				lfp = fopen( command, "w" );
				if( !lfp )
				{
					fprintf( stderr, "Cannot open _db." );
					exit( 1 );
				}
				for( i=0; i<klim; i++ ) fprintf( lfp, ">%d\n%s\n", i, dseq[i] );
				fclose( lfp );

//				sprintf( command, "md5sum _db%dd > /dev/tty", k );
//				system( command );

				if( dorp == 'd' ) 
					sprintf( command, "%s/lastdb _db%dd _db%dd", whereispairalign, k, k );
				else
					sprintf( command, "%s/lastdb -p _db%dd _db%dd", whereispairalign, k, k );
				system( command );
				sprintf( kd, "%d", k );
			}
			else // calllast_fast de tsukutta nowo riyou
			{
				kd[0] = 0;
//				fprintf( stderr, "klim=%d, njob=%d, nadd=%d, skip!\n", klim, njob, nadd );
			}
		}
		else // 'r'
		{
			kd[0] = 0;
		}
		
		sprintf( command, "_q%d", k );
		lfp = fopen( command, "w" );
		if( !lfp )
		{
			fprintf( stderr, "Cannot open %s", command );
			exit( 1 );
		}
		fprintf( lfp, ">%d\n%s\n", k, qseq[k] );
		fclose( lfp );
	
//		if( alg == 'R' ) msize = MAX(10,k+nq);
//			else msize = MAX(10,nd+nq);
		if( alg == 'R' ) msize = MAX(10,k*lastm);
			else msize = MAX(10,nd*lastm);

//		fprintf( stderr, "Calling lastal from lastcallthread, msize = %d, k=%d\n", msize, k );
//		sprintf( command, "grep '>' _db%sd", kd );
//		system( command );
		sprintf( command, "%s/lastal -m %d -e %d -f 0 -s 1 -p _scoringmatrixforlast -a %d -b %d _db%sd _q%d > _lastres%d", whereispairalign, msize, laste, -penalty, -penalty_ex, kd, k, k );
		if( system( command ) ) exit( 1 );
	
		sprintf( command, "_lastres%d", k );
		lfp = fopen( command, "r" );
		if( !lfp )
		{
			fprintf( stderr, "Cannot read _lastres%d", k );
			exit( 1 );
		}
//		readlastres( lfp, nd, nq, lastres, dseq, qseq );
//		fprintf( stderr, "Reading lastres\n" );
		readlastresx_singleq( lfp, nd, k, lastresx ); //defined here. read data from _lastres%d file and fill in lastresx matrix
		fclose( lfp );
	}
	return( NULL );
}


static void calllast_fast( int nd, char **dseq, int nq, char **qseq, Lastresx **lastresx )
{
	int i, j;
	FILE *lfp;
	char command[1000];

	lfp = fopen( "_scoringmatrixforlast", "w" ); //open scoring matrix for last file
	if( !lfp )
	{
		fprintf( stderr, "Cannot open _scoringmatrixforlast" );
		exit( 1 );
	}
	if( dorp == 'd' )
	{
		fprintf( lfp, "      " );
		for( j=0; j<4; j++ ) fprintf( lfp, " %c ", amino[j] );
		fprintf( lfp, "\n" );
		for( i=0; i<4; i++ )
		{
			fprintf( lfp, "%c ", amino[i] );
			for( j=0; j<4; j++ ) fprintf( lfp, " %d ", n_dis[i][j] );
			fprintf( lfp, "\n" );
		} //previous two loops print DNA nucleotides distances to scoringmatrixforlast file
	}
	else
	{
		fprintf( lfp, "      " );
		for( j=0; j<20; j++ ) fprintf( lfp, " %c ", amino[j] );
		fprintf( lfp, "\n" );
		for( i=0; i<20; i++ )
		{
			fprintf( lfp, "%c ", amino[i] );
			for( j=0; j<20; j++ ) fprintf( lfp, " %d ", n_dis[i][j] );
			fprintf( lfp, "\n" );
		} //previous two loops print proteins distances to scoringmatrixforlast file
	}
	fclose( lfp );

//	if( alg == 'r' ) // if 'R' -> lastcallthread, kokonoha nadd>0 no toki nomi shiyou
	{
		sprintf( command, "_dbd" );
		lfp = fopen( command, "w" );
		if( !lfp )
		{
			fprintf( stderr, "Cannot open _dbd" );
			exit( 1 );
		}
		if( alg == 'R' )
			j = njob-nadd;
		else
			j = nd;
		for( i=0; i<j; i++ ) fprintf( lfp, ">%d\n%s\n", i, dseq[i] ); //print sequences to _dbd file

		fclose( lfp );
		if( dorp == 'd' ) 
			sprintf( command, "%s/lastdb _dbd _dbd", whereispairalign );
		else
			sprintf( command, "%s/lastdb -p _dbd _dbd", whereispairalign );
		system( command ); //execute command, i.e. run lastdb command
	}

#ifdef enablemultithread
	if( nthread )
	{
		pthread_t *handle;
		pthread_mutex_t mutex;
		lastcallthread_arg_t *targ;
		int *ksharept;
		targ = (lastcallthread_arg_t *)calloc( nthread, sizeof( lastcallthread_arg_t ) );
		handle = calloc( nthread, sizeof( pthread_t ) );
		ksharept = calloc( 1, sizeof(int) );
		*ksharept = 0;
		pthread_mutex_init( &mutex, NULL );
		for( i=0; i<nthread; i++ )
		{
			targ[i].thread_no = i;
			targ[i].kshare = ksharept;
			targ[i].nq = nq;
			targ[i].nd = nd;
			targ[i].dseq = dseq;
			targ[i].qseq = qseq;
			targ[i].lastresx = lastresx;
			targ[i].mutex = &mutex;
			pthread_create( handle+i, NULL, lastcallthread, (void *)(targ+i) );
		}

		for( i=0; i<nthread; i++ )
		{
			pthread_join( handle[i], NULL );
		}
		pthread_mutex_destroy( &mutex );
		free( handle );
		free( targ );
		free( ksharept );
	}
	else
#endif
	{
		lastcallthread_arg_t *targ; //lastcallthread_arg_t is a structure defined here
		targ = (lastcallthread_arg_t *)calloc( 1, sizeof( lastcallthread_arg_t ) );
		targ[0].nq = nq;
		targ[0].nd = nd;
		targ[0].dseq = dseq;
		targ[0].qseq = qseq;
		targ[0].lastresx = lastresx;
		lastcallthread( targ ); //defined here. I think this calls last service also but on threads/separated steps, not like calllast_once method
		free( targ );
	}

}

//this method calls lastal - which finds similar regions between sequences - and save its result in lastresx
static void calllast_once( int nd, char **dseq, int nq, char **qseq, Lastresx **lastresx )
{
	int i, j;
	char command[5000];
	FILE *lfp;
	int msize;
	int res;

	fprintf( stderr, "nq=%d\n", nq );

	lfp = fopen( "_db", "w" );
	if( !lfp )
	{
		fprintf( stderr, "Cannot open _db" );
		exit( 1 );
	}
	for( i=0; i<nd; i++ ) fprintf( lfp, ">%d\n%s\n", i, dseq[i] ); //print sequences in _db file
	fclose( lfp );

	if( dorp == 'd' ) //DNA
	{
		sprintf( command, "%s/lastdb _db _db", whereispairalign ); //whereispairalign is set from input argument.
		system( command ); //execute previous command. I think this command sends sequences to last service to be aligned in next command
		lfp = fopen( "_scoringmatrixforlast", "w" ); //open scoring matrix for last file
		if( !lfp )
		{
			fprintf( stderr, "Cannot open _scoringmatrixforlast" );
			exit( 1 );
		}
		fprintf( lfp, "      " );
		for( j=0; j<4; j++ ) fprintf( lfp, " %c ", amino[j] );
		fprintf( lfp, "\n" );
		for( i=0; i<4; i++ )
		{
			fprintf( lfp, "%c ", amino[i] );
			for( j=0; j<4; j++ ) fprintf( lfp, " %d ", n_dis[i][j] );
			fprintf( lfp, "\n" );
		}
		fclose( lfp ); //previous two loops print DNA nucleotides distances to scoringmatrixforlast file
#if 0
		sprintf( command, "lastex -s 2 -a %d -b %d -p _scoringmatrixforlast -E 10000 _db.prj _db.prj > _lastex", -penalty, -penalty_ex );
		system( command );
		lfp = fopen( "_lastex", "r" );
		fgets( command, 4999, lfp );
		fgets( command, 4999, lfp );
		fgets( command, 4999, lfp );
		fgets( command, 4999, lfp );
		laste = atoi( command );
		fclose( lfp );
		fprintf( stderr, "laste = %d\n", laste );
		sleep( 10 );
#else
//		laste = 5000;
#endif
	}
	else
	{
		sprintf( command, "%s/lastdb -p _db _db", whereispairalign ); //whereispairalign is set from input argument
		system( command ); //execute previous command
		lfp = fopen( "_scoringmatrixforlast", "w" ); //open scoring matrix for last file
		if( !lfp )
		{
			fprintf( stderr, "Cannot open _scoringmatrixforlast" );
			exit( 1 );
		}
		fprintf( lfp, "      " );
		for( j=0; j<20; j++ ) fprintf( lfp, " %c ", amino[j] );
		fprintf( lfp, "\n" );
		for( i=0; i<20; i++ )
		{
			fprintf( lfp, "%c ", amino[i] );
			for( j=0; j<20; j++ ) fprintf( lfp, " %d ", n_dis[i][j] );
			fprintf( lfp, "\n" );
		}
		fclose( lfp ); //previous two loops print proteins distances to scoringmatrixforlast file
//		fprintf( stderr, "Not written yet\n" );
	}

	lfp = fopen( "_q", "w" );
	if( !lfp )
	{
		fprintf( stderr, "Cannot open _q" );
		exit( 1 );
	}
	for( i=0; i<nq; i++ )
	{
		fprintf( lfp, ">%d\n%s\n", i, qseq[i] ); //print each sequence number followed by the sequence to _q file
	}
	fclose( lfp );

	msize = MAX(10,nd*lastm); //MAX is defined in fft.h. lastm default = 3, otherwise set from arguments

//	fprintf( stderr, "Calling lastal from calllast_once, msize=%d\n", msize );
	sprintf( command, "%s/lastal -v -m %d -e %d -f 0 -s 1 -p _scoringmatrixforlast -a %d -b %d _db _q > _lastres", whereispairalign, msize, laste, -penalty, -penalty_ex );
//	sprintf( command, "lastal -v -m %d -e %d -f 0 -s 1 -p _scoringmatrixforlast -a %d -b %d _db _q > _lastres", 1, laste, -penalty, -penalty_ex );
//	sprintf( command, "lastal -v -e 40 -f 0 -s 1 -p _scoringmatrixforlast -a %d -b %d _db _q > _lastres", -penalty, -penalty_ex );
	res = system( command ); //execute lastal command
	if( res )
	{
		fprintf( stderr, "LAST aborted\n" );
		exit( 1 );
	}

	lfp = fopen( "_lastres", "r" );
	if( !lfp )
	{
		fprintf( stderr, "Cannot read _lastres" );
		exit( 1 );
	}
//	readlastres( lfp, nd, nq, lastres, dseq, qseq );
	fprintf( stderr, "Reading lastres\n" );
	readlastresx( lfp, nd, nq, lastresx, dseq, qseq ); //defined here. read data from _lastres file and fill in lastres matrix
	fclose( lfp );
}

//this method calls foldalign command with sequences as input and output is saved in _foldalignout file
static void callfoldalign( int nseq, char **mseq )
{
	FILE *fp;
	int i;
	int res;
	static char com[10000];

	for( i=0; i<nseq; i++ )
		t2u( mseq[i] ); //defined here. convert capital letters to small and t to u

	fp = fopen( "_foldalignin", "w" );
	if( !fp )
	{
		fprintf( stderr, "Cannot open _foldalignin\n" );
		exit( 1 );
	}
	for( i=0; i<nseq; i++ ) //print sequences and their numbers to _foldalignin file
	{
		fprintf( fp, ">%d\n", i+1 );
		fprintf( fp, "%s\n", mseq[i] );
	}
	fclose( fp );

	sprintf( com, "env PATH=%s  foldalign210 %s _foldalignin > _foldalignout ", whereispairalign, foldalignopt );
	res = system( com ); //execute foldalign command and save output in _foldalignout
	if( res )
	{
		fprintf( stderr, "Error in foldalign\n" );
		exit( 1 );
	}

}

//this method executes mafft_lara command with sequences as input and output saved to _laraout
static void calllara( int nseq, char **mseq, char *laraarg )
{
	FILE *fp;
	int i;
	int res;
	static char com[10000];

//	for( i=0; i<nseq; i++ )

	fp = fopen( "_larain", "w" );
	if( !fp )
	{
		fprintf( stderr, "Cannot open _larain\n" );
		exit( 1 );
	}
	for( i=0; i<nseq; i++ ) //print sequences and their numbers to _larain file
	{
		fprintf( fp, ">%d\n", i+1 );
		fprintf( fp, "%s\n", mseq[i] );
	}
	fclose( fp );


//	fprintf( stderr, "calling LaRA\n" );
	sprintf( com, "env PATH=%s:/bin:/usr/bin mafft_lara -i _larain -w _laraout -o _lara.params %s", whereispairalign, laraarg );
	res = system( com ); //execute lara command with input in _larain and output to _laraout
	if( res )
	{
		fprintf( stderr, "Error in lara\n" );
		exit( 1 );
	}
}

//read sequences from lara file and if matched with given ones, calculate naive score and return it.
static double recalllara( char **mseq1, char **mseq2, int alloclen )
{
	static FILE *fp = NULL;
	static char *ungap1;
	static char *ungap2;
	static char *ori1;
	static char *ori2;
//	int res;
	static char com[10000];
	double value;


	if( fp == NULL )
	{
		fp = fopen( "_laraout", "r" );
		if( fp == NULL )
		{
			fprintf( stderr, "Cannot open _laraout\n" );
			exit( 1 );
		}
		ungap1 = AllocateCharVec( alloclen );
		ungap2 = AllocateCharVec( alloclen );
		ori1 = AllocateCharVec( alloclen );
		ori2 = AllocateCharVec( alloclen );
	}


	strcpy( ori1, *mseq1 );
	strcpy( ori2, *mseq2 );

	// I think these five lines from lara file read three lines from the file and fill seq1 and seq2 with last two lines.
	fgets( com, 999, fp );
	myfgets( com, 9999, fp ); //defined in io.c. read line from fp into com
	strcpy( *mseq1, com );
	myfgets( com, 9999, fp );
	strcpy( *mseq2, com );

	gappick0( ungap1, *mseq1 ); //defined in io.c. copy mseq1 chars to ungap1 without gaps chars
	gappick0( ungap2, *mseq2 ); //copy mseq2 chars to ungap2 without gaps chars
	t2u( ungap1 ); //defined here. converts capital letters to small and t to u
	t2u( ungap2 );
	t2u( ori1 );
	t2u( ori2 );

	if( strcmp( ungap1, ori1 ) || strcmp( ungap2, ori2 ) ) //compare original sequences and those read from lara file. if changed, exit, else continue.
	{
		fprintf( stderr, "SEQUENCE CHANGED!!\n" );
		fprintf( stderr, "*mseq1  = %s\n", *mseq1 );
		fprintf( stderr, "ungap1  = %s\n", ungap1 );
		fprintf( stderr, "ori1    = %s\n", ori1 );
		fprintf( stderr, "*mseq2  = %s\n", *mseq2 );
		fprintf( stderr, "ungap2  = %s\n", ungap2 );
		fprintf( stderr, "ori2    = %s\n", ori2 );
		exit( 1 );
	}

	value = (double)naivepairscore11( *mseq1, *mseq2, penalty ); //defined in mltaln9.c. calculates score between seq1 and seq2 based on penal and amino_dis values

//	fclose( fp ); // saigo dake yatta houga yoi.

	return( value );
}

//executes some commands then write mseq1 and mseq2 to file and align them, then apply naive score and return it.
static double calldafs_giving_bpp( char **mseq1, char **mseq2, char **bpp1, char **bpp2, int alloclen, int i, int j )
{
	FILE *fp;
	int res;
	char *com;
	double value;
	char *dirname;


	dirname = calloc( 100, sizeof( char ) );
	com = calloc( 1000, sizeof( char ) );
	sprintf( dirname, "_%d-%d", i, j );
	sprintf( com, "rm -rf %s", dirname );
	system( com ); //execute the given command
	sprintf( com, "mkdir %s", dirname );
	system( com );


	sprintf( com, "%s/_bpporg", dirname );
	fp = fopen( com, "w" ); //open file to write
	if( !fp )
	{
		fprintf( stderr, "Cannot write to %s/_bpporg\n", dirname );
		exit( 1 );
	}
	fprintf( fp, ">a\n" );
	while( *bpp1 ) //write bpp1 to the file
		fprintf( fp, "%s", *bpp1++ );

	fprintf( fp, ">b\n" );
	while( *bpp2 ) //write bpp2 to the file
		fprintf( fp, "%s", *bpp2++ );
	fclose( fp );

	sprintf( com, "tr -d '\\r' < %s/_bpporg > %s/_bpp", dirname, dirname );
	system( com ); // for cygwin, wakaran

	t2u( *mseq1 ); //defined here. convert capital letters to small and t to u.
	t2u( *mseq2 );

	sprintf( com, "%s/_dafsinorg", dirname );
	fp = fopen( com, "w" );
	if( !fp )
	{
		fprintf( stderr, "Cannot open %s/_dafsinorg\n", dirname );
		exit( 1 );
	}
	fprintf( fp, ">1\n" );
//	fprintf( fp, "%s\n", *mseq1 );
	write1seq( fp, *mseq1 ); //defined in io.c. write mseq1 to fp
	fprintf( fp, ">2\n" );
//	fprintf( fp, "%s\n", *mseq2 );
	write1seq( fp, *mseq2 ); //write mseq2 to fp
	fclose( fp );

	sprintf( com, "tr -d '\\r' < %s/_dafsinorg > %s/_dafsin", dirname, dirname );
	system( com ); // for cygwin, wakaran

	sprintf( com, "_dafssh%s", dirname );
	fp = fopen( com, "w" );
	fprintf( fp, "cd %s\n", dirname );
	fprintf( fp, "%s/dafs --mafft-in _bpp _dafsin > _dafsout 2>_dum\n", whereispairalign );
	fprintf( fp, "exit $tatus\n" );
	fclose( fp );

	sprintf( com, "tr -d '\\r' < _dafssh%s > _dafssh%s.unix", dirname, dirname );
	system( com ); // for cygwin, wakaran

	sprintf( com, "sh _dafssh%s.unix 2>_dum%s", dirname, dirname );
	res = system( com );
	if( res )
	{
		fprintf( stderr, "Error in dafs\n" );
		exit( 1 );
	}

	sprintf( com, "%s/_dafsout", dirname );

	fp = fopen( com, "r" );
	if( !fp )
	{
		fprintf( stderr, "Cannot open %s/_dafsout\n", dirname );
		exit( 1 );
	}

	myfgets( com, 999, fp ); // nagai kanousei ga arunode
	fgets( com, 999, fp );
	myfgets( com, 999, fp ); // nagai kanousei ga arunode
	fgets( com, 999, fp );
	load1SeqWithoutName_new( fp, *mseq1 ); //I think this reads sequence from fp into mseq1 - without name -.
	fgets( com, 999, fp );
	load1SeqWithoutName_new( fp, *mseq2 ); //I think this reads sequence from fp into mseq2 - without name -.

	fclose( fp );

//	fprintf( stderr, "*mseq1 = %s\n", *mseq1 );
//	fprintf( stderr, "*mseq2 = %s\n", *mseq2 );

	//calculates score between seq1 and seq2 based on penal and amino_dis values
	value = (double)naivepairscore11( *mseq1, *mseq2, penalty ); //defined in mltaln9.c.

#if 0
	sprintf( com, "rm -rf %s > /dev/null 2>&1", dirname );
	if( system( com ) )
	{
		fprintf( stderr, "retrying to rmdir\n" );
		usleep( 2000 );
		system( com );
	}
#endif

	free( dirname );
	free( com );


	return( value );
}

//executes some commands then write mseq1 and mseq2 to file and align them, then apply naive score and return it.
static double callmxscarna_giving_bpp( char **mseq1, char **mseq2, char **bpp1, char **bpp2, int alloclen, int i, int j )
{
	FILE *fp;
	int res;
	char *com;
	double value;
	char *dirname;


	dirname = calloc( 100, sizeof( char ) );
	com = calloc( 1000, sizeof( char ) );
	sprintf( dirname, "_%d-%d", i, j );
	sprintf( com, "rm -rf %s", dirname );
	system( com ); //execute the given command
	sprintf( com, "mkdir %s", dirname );
	system( com );


	sprintf( com, "%s/_bpporg", dirname );
	fp = fopen( com, "w" ); //open file to write
	if( !fp )
	{
		fprintf( stderr, "Cannot write to %s/_bpporg\n", dirname );
		exit( 1 );
	}
	fprintf( fp, ">a\n" );
	while( *bpp1 ) //write bpp1 to the file
		fprintf( fp, "%s", *bpp1++ );

	fprintf( fp, ">b\n" );
	while( *bpp2 ) //write bpp2 to the file
		fprintf( fp, "%s", *bpp2++ );
	fclose( fp );

	sprintf( com, "tr -d '\\r' < %s/_bpporg > %s/_bpp", dirname, dirname );
	system( com ); // for cygwin, wakaran

	t2u( *mseq1 ); //defined here. convert capital letters to small and t to u.
	t2u( *mseq2 );

	sprintf( com, "%s/_mxscarnainorg", dirname );
	fp = fopen( com, "w" );
	if( !fp )
	{
		fprintf( stderr, "Cannot open %s/_mxscarnainorg\n", dirname );
		exit( 1 );
	}
	fprintf( fp, ">1\n" );
//	fprintf( fp, "%s\n", *mseq1 );
	write1seq( fp, *mseq1 ); //defined in io.c. write mseq1 to fp
	fprintf( fp, ">2\n" );
//	fprintf( fp, "%s\n", *mseq2 );
	write1seq( fp, *mseq2 );  //write mseq2 to fp
	fclose( fp );

	sprintf( com, "tr -d '\\r' < %s/_mxscarnainorg > %s/_mxscarnain", dirname, dirname );
	system( com ); // for cygwin, wakaran

#if 0
	sprintf( com, "cd %s; %s/mxscarnamod -readbpp _mxscarnain > _mxscarnaout 2>_dum", dirname, whereispairalign );
#else
	sprintf( com, "_mxscarnash%s", dirname );
	fp = fopen( com, "w" );
	fprintf( fp, "cd %s\n", dirname );
	fprintf( fp, "%s/mxscarnamod -readbpp _mxscarnain > _mxscarnaout 2>_dum\n", whereispairalign );
	fprintf( fp, "exit $tatus\n" );
	fclose( fp );
//sleep( 10000 );

	sprintf( com, "tr -d '\\r' < _mxscarnash%s > _mxscarnash%s.unix", dirname, dirname );
	system( com ); // for cygwin, wakaran

	sprintf( com, "sh _mxscarnash%s.unix 2>_dum%s", dirname, dirname );
#endif
	res = system( com );
	if( res )
	{
		fprintf( stderr, "Error in mxscarna\n" );
		exit( 1 );
	}

	sprintf( com, "%s/_mxscarnaout", dirname );

	fp = fopen( com, "r" );
	if( !fp )
	{
		fprintf( stderr, "Cannot open %s/_mxscarnaout\n", dirname );
		exit( 1 );
	}

	fgets( com, 999, fp );
	load1SeqWithoutName_new( fp, *mseq1 ); //I think this reads sequence from fp into mseq1 - without name -.
	fgets( com, 999, fp );
	load1SeqWithoutName_new( fp, *mseq2 ); //I think this reads sequence from fp into mseq2 - without name -.

	fclose( fp );

//	fprintf( stderr, "*mseq1 = %s\n", *mseq1 );
//	fprintf( stderr, "*mseq2 = %s\n", *mseq2 );

	//calculates score between seq1 and seq2 based on penal and amino_dis values
	value = (double)naivepairscore11( *mseq1, *mseq2, penalty ); //defined in mltaln9.c.

#if 0
	sprintf( com, "rm -rf %s > /dev/null 2>&1", dirname );
	if( system( com ) )
	{
		fprintf( stderr, "retrying to rmdir\n" );
		usleep( 2000 );
		system( com );
	}
#endif

	free( dirname );
	free( com );


	return( value );
}

//read fp file into bpp matrix
static void readhat4( FILE *fp, char ***bpp )
{
	char oneline[1000];
	int bppsize;
	int onechar;
//	double prob;
//	int posi, posj;

	bppsize = 0;
//	fprintf( stderr, "reading hat4\n" );
	onechar = getc(fp);
//	fprintf( stderr, "onechar = %c\n", onechar );
	if( onechar != '>' )
	{
		fprintf( stderr, "Format error\n" );
		exit( 1 );
	}
	ungetc( onechar, fp );
	fgets( oneline, 999, fp ); //read line from hat4 file
	while( 1 )
	{
		onechar = getc(fp);
		ungetc( onechar, fp );
		if( onechar == '>' || onechar == EOF ) //I think this means end loop and return back to the caller
		{
//			fprintf( stderr, "Next\n" );
			*bpp = realloc( *bpp, (bppsize+2) * sizeof( char * ) );
			(*bpp)[bppsize] = NULL;
			break;
		}
		fgets( oneline, 999, fp );
//		fprintf( stderr, "oneline=%s\n", oneline );
//		sscanf( oneline, "%d %d %lf", &posi, &posj, &prob );
//		fprintf( stderr, "%d %d -> %f\n", posi, posj, prob );
		*bpp = realloc( *bpp, (bppsize+2) * sizeof( char * ) );
		(*bpp)[bppsize] = calloc( 100, sizeof( char ) );
		strcpy( (*bpp)[bppsize], oneline ); //copy line read from hat4 file to bpp
		bppsize++; //and increment bppsize by 1
	}
}

//read hat4 file into bpp matrix
static void preparebpp( int nseq, char ***bpp )
{
	FILE *fp;
	int i;

	fp = fopen( "hat4", "r" );
	if( !fp )
	{
		fprintf( stderr, "Cannot open hat4\n" );
		exit( 1 );
	}
	for( i=0; i<nseq; i++ )
		readhat4( fp, bpp+i ); //read hat4 file into bpp matrix
	fclose( fp );
}

static void pair_local_align_arguments( int argc, char *argv[] )
{
    int c;

	nthread = 1; //defined in defs.c
	laste = 5000; //defined here
	lastm = 3; //defined here
	nadd = 0; //defined here
	lastsubopt = 0; //defined here
	lastonce = 0; //defined here
	foldalignopt[0] = 0; //defined here
	laraparams = NULL; //defined here
	inputfile = NULL; //defined in defs.h
	fftkeika = 0; //defined in defs.h
	pslocal = -1000.0; //defined in defs.h
	constraint = 0; //defined in defs.h
	nblosum = 62; //defined in defs.h
	fmodel = 0; //defined in defs.h
	calledByXced = 0; //defined in defs.h
	devide = 0; //defined in defs.h
	use_fft = 0; //defined in defs.h
	fftscore = 1; //defined in defs.h
	fftRepeatStop = 0; //defined in defs.h
	fftNoAnchStop = 0; //defined in defs.h
    weight = 3; //defined in defs.h
    utree = 1; //defined in defs.h
	tbutree = 1; //defined in defs.h
    refine = 0; //defined in defs.h
    check = 1; //defined in defs.h
    cut = 0.0; //defined in defs.h
    disp = 0;  //defined in defs.h
    outgap = 1; //defined in defs.c
    alg = 'A'; //defined in defs.h
    mix = 0; //defined in defs.h
	tbitr = 0; //defined in defs.h
	scmtd = 5; //defined in defs.h
	tbweight = 0; //defined in defs.h
	tbrweight = 3; //defined in defs.h
	checkC = 0; //defined in defs.h
	treemethod = 'x'; //defined in defs.h
	contin = 0; //defined in defs.h
	scoremtx = 1; //defined in defs.h
	kobetsubunkatsu = 0; //defined in defs.h. It means Individual division
	divpairscore = 0; //defined in defs.h
	stdout_align = 0; //defined here
	stdout_dist = 0; //defined here
	store_dist = 1; //defined here
	store_localhom = 1; //defined here
//	dorp = NOTSPECIFIED;
	ppenalty = NOTSPECIFIED; //defined in defs.h
	ppenalty_OP = NOTSPECIFIED; //defined in defs.h
	ppenalty_ex = NOTSPECIFIED; //defined in defs.h
	ppenalty_EX = NOTSPECIFIED; //defined in defs.h
	penalty_shift_factor = 1000.0; //defined in defs.c
	poffset = NOTSPECIFIED; //defined in defs.h
	kimuraR = NOTSPECIFIED; //defined in defs.h
	pamN = NOTSPECIFIED; //defined in defs.h
	geta2 = GETA2; //defined in defs.h //GETA2 defined in mltaln.h
	fftWinSize = NOTSPECIFIED; //defined in defs.h
	fftThreshold = NOTSPECIFIED; //defined in defs.h
	RNAppenalty = NOTSPECIFIED; //defined in defs.h
	RNApthr = NOTSPECIFIED; //defined in defs.h
	specificityconsideration = 0.0; //defined in defs.c
	usenaivescoreinsteadofalignmentscore = 0; //defined here. use naive score instead of alignment score
	specifictarget = 0; //defined here
	nwildcard = 0; //defined in defs.c

//	reporterr( "argc=%d\n", argc );
//	reporterr( "*argv=%s\n", *argv );
//	reporterr( "(*argv)[0]=%c\n", (*argv)[0] );
    while( --argc > 0 && (*++argv)[0] == '-' )
	{
//		reporterr( "(*argv)[0] in while loop = %s\n", (*argv) );
        while ( ( c = *++argv[0] ) ) //parse arguments in the running command
		{
            switch( c )
            {
				case 'i':
					inputfile = *++argv;
//					fprintf( stderr, "inputfile = %s\n", inputfile );
					--argc;
					goto nextoption;
				case 'f':
					ppenalty = (int)( atof( *++argv ) * 1000 - 0.5 );
					--argc;
					goto nextoption;
				case 'g':
					ppenalty_ex = (int)( atof( *++argv ) * 1000 - 0.5 );
					--argc;
					goto nextoption;
				case 'O':
					ppenalty_OP = (int)( atof( *++argv ) * 1000 - 0.5 );
					--argc;
					goto nextoption;
				case 'E':
					ppenalty_EX = (int)( atof( *++argv ) * 1000 - 0.5 );
					--argc;
					goto nextoption;
				case 'Q':
					penalty_shift_factor = atof( *++argv );
					--argc;
					goto nextoption;
				case 'h':
					poffset = (int)( atof( *++argv ) * 1000 - 0.5 );
					--argc;
					goto nextoption;
				case 'k':
					kimuraR = myatoi( *++argv );
//					fprintf( stderr, "kimuraR = %d\n", kimuraR );
					--argc;
					goto nextoption;
				case 'b':
					nblosum = myatoi( *++argv );
					scoremtx = 1;
//					fprintf( stderr, "blosum %d\n", nblosum );
					--argc;
					goto nextoption;
				case 'j':
					pamN = myatoi( *++argv );
					scoremtx = 0;
					TMorJTT = JTT;
//					fprintf( stderr, "jtt %d\n", pamN );
					--argc;
					goto nextoption;
				case 'm':
					pamN = myatoi( *++argv );
					scoremtx = 0;
					TMorJTT = TM;
//					fprintf( stderr, "TM %d\n", pamN );
					--argc;
					goto nextoption;
#if 0
				case 'l':
					ppslocal = (int)( atof( *++argv ) * 1000 + 0.5 );
					pslocal = (int)( 600.0 / 1000.0 * ppslocal + 0.5);
//					fprintf( stderr, "ppslocal = %d\n", ppslocal );
//					fprintf( stderr, "pslocal = %d\n", pslocal );
					--argc;
					goto nextoption;
#else
				case 'l':
					if( atof( *++argv ) < 0.00001 ) store_localhom = 0;
					--argc;
					goto nextoption;
#endif
				case 'd':
					whereispairalign = *++argv;
					fprintf( stderr, "whereispairalign = %s\n", whereispairalign );
					--argc; 
					goto nextoption;
				case 'p':
					laraparams = *++argv;
					fprintf( stderr, "laraparams = %s\n", laraparams );
					--argc; 
					goto nextoption;
				case 'C':
					nthread = myatoi( *++argv );
//					fprintf( stderr, "nthread = %d\n", nthread );
					--argc; 
#ifndef enablemultithread
					nthread = 0;
#endif
					goto nextoption;
				case 'I':
					nadd = myatoi( *++argv );
//					fprintf( stderr, "nadd = %d\n", nadd );
					--argc;
					goto nextoption;
				case 'w':
					lastm = myatoi( *++argv );
					fprintf( stderr, "lastm = %d\n", lastm );
					--argc;
					goto nextoption;
				case 'e':
					laste = myatoi( *++argv );
					fprintf( stderr, "laste = %d\n", laste );
					--argc;
					goto nextoption;
				case 'u':
					specificityconsideration = (double)myatof( *++argv );
//					fprintf( stderr, "specificityconsideration = %f\n", specificityconsideration );
					--argc; 
					goto nextoption;
				case 'K': // Hontou ha iranai. disttbfast.c, tbfast.c to awaserutame.
					break;
				case 'c':
					stdout_dist = 1;
					break;
				case 'n':
					stdout_align = 1;
					break;
				case 'x':
					store_localhom = 0;
					store_dist = 0;
					break;
#if 1
				case 'a':
					fmodel = 1;
					break;
#endif
#if 0
				case 'r':
					fmodel = -1;
					break;
#endif
				case 'D':
					dorp = 'd';
					break;
				case 'P':
					dorp = 'p';
					break;
#if 0
				case 'e':
					fftscore = 0;
					break;
				case 'O':
					fftNoAnchStop = 1;
					break;
#endif
#if 0
				case 'Q':
					calledByXced = 1;
					break;
				case 'x':
					disp = 1;
					break;
				case 'a':
					alg = 'a';
					break;
				case 'S':
					alg = 'S';
					break;
#endif
				case 'U':
					lastonce = 1;
					break;
				case 'S':
					lastsubopt = 1;
					break;
				case 't':
					alg = 't';
					store_localhom = 0;
					break;
				case 'L':
					alg = 'L';
					break;
				case 'Y':
					alg = 'Y'; // nadd>0 no toki nomi. moto no hairetsu to atarashii hairetsuno alignmnt -> L;
					break;
				case 'Z':
					usenaivescoreinsteadofalignmentscore = 1;
					break;
				case 's':
					alg = 's';
					break;
				case 'G':
					alg = 'G';
					break;
				case 'B':
					alg = 'B';
					break;
				case 'T':
					alg = 'T';
					break;
				case 'H':
					alg = 'H';
					break;
				case 'M':
					alg = 'M';
					break;
				case 'R':
					alg = 'R';
					break;
				case 'r':
					alg = 'r'; // nadd>0 no toki nomi. moto no hairetsu to atarashii hairetsuno alignmnt -> R, last
					break;
				case 'N':
					alg = 'N';
					break;
				case 'A':
					alg = 'A';
					break;
				case 'V':
					alg = 'V';
					break;
				case 'F':
					use_fft = 1;
					break;
				case 'v':
					tbrweight = 3;
					break;
				case 'y':
					divpairscore = 1;
					break;
				case '=':
					specifictarget = 1;
					break;
				case ':':
					nwildcard = 1;
					break;
/* Modified 01/08/27, default: user tree */
				case 'J':
					tbutree = 0;
					break;
/* modification end. */
				case 'o':
//					foldalignopt = *++argv;
					strcat( foldalignopt, " " ); //append " " to foldalignopt
					strcat( foldalignopt, *++argv ); //then append argument value to it
					fprintf( stderr, "foldalignopt = %s\n", foldalignopt );
					--argc; 
					goto nextoption;
#if 0
				case 'z':
					fftThreshold = myatoi( *++argv );
					--argc; 
					goto nextoption;
				case 'w':
					fftWinSize = myatoi( *++argv );
					--argc;
					goto nextoption;
				case 'Z':
					checkC = 1;
					break;
#endif
                default:
                    fprintf( stderr, "illegal option %c\n", c );
                    argc = 0;
                    break;
            }
		}
		nextoption:
			;
	}
    if( argc == 1 )
    {
        cut = atof( (*argv) );
        argc--;
    }
    if( argc != 0 ) 
    {
        fprintf( stderr, "pairlocalalign options: Check source file !\n" );
        exit( 1 );
    }
	if( tbitr == 1 && outgap == 0 )
	{
		fprintf( stderr, "conflicting options : o, m or u\n" );
		exit( 1 );
	}
}

int pair_local_align_countamino( char *s, int end ) //count the number of amino acids in the sequence
{
	int val = 0;
	while( end-- )
		if( *s++ != '-' ) val++;
	return( val );
}

//returns value based on the three args values.
static double score2dist( double pscore, double selfscore1, double selfscore2)
{
	double val;
	double bunbo;
//	fprintf( stderr, "In score2dist\n" );

	if( (bunbo=MIN( selfscore1, selfscore2 )) == 0.0 )
		val = 2.0;
	else if( bunbo < pscore ) // mondai ari
		val = 0.0;
	else
		val = ( 1.0 - pscore / bunbo ) * 2.0;
	return( val );
}

#if enablemultithread
static void *athread( void *arg ) // alg='R', alg='r' -> tsukawarenai.
{
	thread_arg_t *targ = (thread_arg_t *)arg;
	int i, ilim, j, jst;
	int off1, off2, dum1, dum2, thereisx;
	int intdum;
	double pscore = 0.0; // by D.Mathog
	double *effarr1;
	double *effarr2;
	char **mseq1, **mseq2, **distseq1, **distseq2, **dumseq1, **dumseq2;
	char **aseq;
	double **dynamicmtx = NULL;
	double dist;
	double scoreoffset;

// thread_arg
	int thread_no = targ->thread_no;
	int njob = targ->njob;
	Jobtable *jobpospt = targ->jobpospt;
	char **name = targ->name;
	char **seq = targ->seq;
	char **dseq = targ->dseq;
	int *thereisxineachseq = targ->thereisxineachseq;
	LocalHom **localhomtable = targ->localhomtable;
	double **distancemtx = targ->distancemtx;
	double *selfscore = targ->selfscore;
	char ***bpp = targ->bpp;
	Lastresx **lastresx = targ->lastresx;
	int alloclen = targ->alloclen;
	int *targetmap = targ->targetmap;

//	fprintf( stderr, "thread %d start!\n", thread_no );

	effarr1 = AllocateDoubleVec( 1 );
	effarr2 = AllocateDoubleVec( 1 );
	mseq1 = AllocateCharMtx( njob, 0 );
	mseq2 = AllocateCharMtx( njob, 0 );
	if( alg == 'N' )
	{
		dumseq1 = AllocateCharMtx( 1, alloclen+10 );
		dumseq2 = AllocateCharMtx( 1, alloclen+10 );
	}
	distseq1 = AllocateCharMtx( 1, 0 );
	distseq2 = AllocateCharMtx( 1, 0 );
	aseq = AllocateCharMtx( 2, alloclen+10 );
	if( specificityconsideration > 0.0 ) dynamicmtx = AllocateDoubleMtx( nalphabets, nalphabets );

	if( alg == 'Y' || alg == 'r' ) ilim = njob - nadd;
	else ilim = njob - 1;


	while( 1 )
	{
		pthread_mutex_lock( targ->mutex_counter );
		j = jobpospt->j;
		i = jobpospt->i;
		j++;
		if( j == njob )
		{
			i++;

			if( alg == 'Y' || alg == 'r' ) jst = njob - nadd;
			else jst = i + 1;
			j = jst; 

			if( i == ilim )
			{
//				fprintf( stderr, "thread %d end!\n", thread_no );
				pthread_mutex_unlock( targ->mutex_counter );

				if( commonIP ) FreeIntMtx( commonIP );
				commonIP = NULL;
				if( commonJP ) FreeIntMtx( commonJP );
				commonJP = NULL;
				Falign( NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, 0, 0, 0, NULL, NULL, 0, NULL );
				G__align11( NULL, NULL, NULL, 0, 0, 0 ); // 20130603
				G__align11_noalign( NULL, 0, 0, NULL, NULL, 0 );
				L__align11( NULL, 0.0, NULL, NULL, 0, NULL, NULL );
				L__align11_noalign( NULL, NULL, NULL );
				genL__align11( NULL, NULL, NULL, 0, NULL, NULL );
				free( effarr1 );
				free( effarr2 );
				free( mseq1 );
				free( mseq2 );
				if( alg == 'N' )
				{
					FreeCharMtx( dumseq1 );
					FreeCharMtx( dumseq2 );
				}
				free( distseq1 );
				free( distseq2 );
				FreeCharMtx( aseq  );
				if( dynamicmtx ) FreeDoubleMtx( dynamicmtx  );
				return( NULL );
			}
		}
		jobpospt->j = j;
		jobpospt->i = i;
		pthread_mutex_unlock( targ->mutex_counter );


		if( j == i+1 || j % 100 == 0 ) 
		{
			fprintf( stderr, "% 5d / %d (by thread %3d) \r", i, njob-nadd, thread_no );
//			fprintf( stderr, "% 5d - %5d / %d (thread %d)\n", i, j, njob, thread_no );
		}


		if( strlen( seq[i] ) == 0 || strlen( seq[j] ) == 0 )
		{
			if( store_dist )
			{
				if( alg == 'Y' || alg == 'r' ) distancemtx[i][j-(njob-nadd)] = 3.0;
				else distancemtx[i][j-i] = 3.0;
			}
			if( stdout_dist) 
			{
				pthread_mutex_lock( targ->mutex_stdout );
				fprintf( stdout, "%d %d d=%.3f\n", i+1, j+1, 3.0 );
				pthread_mutex_unlock( targ->mutex_stdout );
			}
			continue;
		}

		strcpy( aseq[0], seq[i] );
		strcpy( aseq[1], seq[j] );
//		clus1 = conjuctionfortbfast( pair, i, aseq, mseq1, effarr1, effarr, indication1 );
//		clus2 = conjuctionfortbfast( pair, j, aseq, mseq2, effarr2, effarr, indication2 );
//		fprintf( stderr, "Skipping conjuction..\n" );

		effarr1[0] = 1.0;
		effarr2[0] = 1.0;
		mseq1[0] = aseq[0];
		mseq2[0] = aseq[1];

		thereisx = thereisxineachseq[i] + thereisxineachseq[j];
//		strcpy( distseq1[0], dseq[i] ); // nen no tame
//		strcpy( distseq2[0], dseq[j] ); // nen no tame
		distseq1[0] = dseq[i];
		distseq2[0] = dseq[j];

//		fprintf( stderr, "mseq1 = %s\n", mseq1[0] );
//		fprintf( stderr, "mseq2 = %s\n", mseq2[0] );
	
#if 0
		fprintf( stderr, "group1 = %.66s", indication1 );
		fprintf( stderr, "\n" );
		fprintf( stderr, "group2 = %.66s", indication2 );
		fprintf( stderr, "\n" );
#endif
//		for( l=0; l<clus1; l++ ) fprintf( stderr, "## STEP-eff for mseq1-%d %f\n", l, effarr1[l] );

		if( use_fft )
		{
			pscore = Falign( NULL, NULL, n_dis_consweight_multi, mseq1, mseq2, effarr1, effarr2, NULL, NULL, 1, 1, alloclen, &intdum, NULL, 0, NULL );
//			fprintf( stderr, "pscore (fft) = %f\n", pscore );
			off1 = off2 = 0;
		}
		else
		{
			switch( alg )
			{
				case( 'R' ):
					if( nadd && njob-nadd <= j && njob-nadd <= i ) // new sequence doushi ha mushi
						pscore = 0.0;
					else
						pscore = (double)lastresx[i][j].score; // all pair
					break;
				case( 'r' ):
					if( nadd == 0 || ( i < njob-nadd && njob-nadd <= j ) )
						pscore = (double)lastresx[i][j-(njob-nadd)].score;
					else
						pscore = 0.0;
					break;
				case( 'L' ):
					if( nadd && njob-nadd <= j && njob-nadd <= i ) // new sequence doushi ha mushi
						pscore = 0.0;
					else
					{
						if( usenaivescoreinsteadofalignmentscore )
						{
							L__align11( n_dis_consweight_multi, 0.0, mseq1, mseq2, alloclen, &off1, &off2 );
							pscore = (double)naivepairscore11( mseq1[0], mseq2[0], 0.0 ); // uwagaki
						}
						else
						{
//							if( store_localhom )
							if( store_localhom && ( targetmap[i] != -1 || targetmap[j] != -1 ) )
							{
								pscore = L__align11( n_dis_consweight_multi, 0.0, mseq1, mseq2, alloclen, &off1, &off2 );
								if( thereisx ) pscore = L__align11_noalign( n_dis_consweight_multi, distseq1, distseq2 ); // uwagaki
#if 1
								if( specificityconsideration > 0.0 )
								{
									dist = score2dist( pscore, selfscore[i], selfscore[j] );
									if( ( scoreoffset = dist2offset( dist ) ) < 0.0 )
									{
										makedynamicmtx( dynamicmtx, n_dis_consweight_multi, 0.5 * dist ); // upgma ni awaseru.
										strcpy( mseq1[0], seq[i] );
										strcpy( mseq2[0], seq[j] );
										L__align11( dynamicmtx, scoreoffset, mseq1, mseq2, alloclen, &off1, &off2 );
									}
								}
#endif
							}
							else
								pscore = L__align11_noalign( n_dis_consweight_multi, distseq1, distseq2 );
						}
					}
//					pscore = G__align11_noalign( n_dis_consweight_multi, penalty, penalty_ex, distseq1, distseq2, alloclen ); // CHUUI!!!!!!
					break;
				case( 'Y' ):
					if( nadd == 0 || ( i < njob-nadd && njob-nadd <= j ) ) // new sequence vs exiting sequence nomi keisan
					{
						if( usenaivescoreinsteadofalignmentscore )
						{
							L__align11( n_dis_consweight_multi, 0.0, mseq1, mseq2, alloclen, &off1, &off2 );
							pscore = (double)naivepairscore11( mseq1[0], mseq2[0], 0.0 ); // uwagaki
						}
						else
						{
							if( store_localhom )
							{
								pscore = L__align11( n_dis_consweight_multi, 0.0, mseq1, mseq2, alloclen, &off1, &off2 );
								if( thereisx ) pscore = L__align11_noalign( n_dis_consweight_multi, distseq1, distseq2 ); // uwagaki
							}
							else
								pscore = L__align11_noalign( n_dis_consweight_multi, distseq1, distseq2 );
						}
					}
					else
						pscore = 0.0;
					break;
				case( 'A' ):
					if( usenaivescoreinsteadofalignmentscore )
					{
						G__align11( n_dis_consweight_multi, mseq1, mseq2, alloclen, outgap, outgap );
						pscore = (double)naivepairscore11( mseq1[0], mseq2[0], 0.0 ); // uwagaki
					}
					else
					{
//						if( store_localhom )
						if( store_localhom && ( targetmap[i] != -1 || targetmap[j] != -1 ) )
						{
							pscore = G__align11( n_dis_consweight_multi, mseq1, mseq2, alloclen, outgap, outgap );
							if( thereisx ) pscore = G__align11_noalign( n_dis_consweight_multi, penalty, penalty_ex, distseq1, distseq2, alloclen ); // uwagaki
#if 1
							if( specificityconsideration > 0.0 )
							{
								dist = score2dist( pscore, selfscore[i], selfscore[j] );
//								dist = score2dist( L__align11_noalign( n_dis_consweight_multi, distseq1, distseq2 ), selfscore[i], selfscore[j] ); // 2014/Feb/20
								if( dist2offset( dist ) < 0.0 )
								{
									makedynamicmtx( dynamicmtx, n_dis_consweight_multi, 0.5 * dist ); // upgma ni awaseru.
									strcpy( mseq1[0], seq[i] );
									strcpy( mseq2[0], seq[j] );
									G__align11( dynamicmtx, mseq1, mseq2, alloclen, outgap, outgap );
					
								}
//								pscore = (double)naivepairscore11( *mseq1, *mseq2, 0.0 );
							}
#endif
						}
						else
							pscore = G__align11_noalign( n_dis_consweight_multi, penalty, penalty_ex, distseq1, distseq2, alloclen ); // uwagaki
					}
					off1 = off2 = 0;
					break;
				case( 'N' ):
					if( usenaivescoreinsteadofalignmentscore )
					{
						genL__align11( n_dis_consweight_multi, mseq1, mseq2, alloclen, &off1, &off2 );
						pscore = (double)naivepairscore11( mseq1[0], mseq2[0], 0.0 ); // uwagaki
					}
					else
					{
//						pscore = G__align11_noalign( n_dis_consweight_multi, penalty, penalty_ex, mseq1, mseq2, alloclen );
						pscore = genL__align11( n_dis_consweight_multi, mseq1, mseq2, alloclen, &off1, &off2 );
						if( thereisx )
						{
							strcpy( dumseq1[0], distseq1[0] );
							strcpy( dumseq2[0], distseq2[0] );
							pscore = genL__align11( n_dis_consweight_multi, dumseq1, dumseq2, alloclen, &dum1, &dum2 ); // uwagaki
						}
#if 1
						if( specificityconsideration > 0.0 )
						{
							dist = score2dist( pscore, selfscore[i], selfscore[j] );
							if( dist2offset( dist ) < 0.0 )
							{
								makedynamicmtx( dynamicmtx, n_dis_consweight_multi, 0.5 * dist ); // upgma ni awaseru.
								strcpy( mseq1[0], seq[i] );
								strcpy( mseq2[0], seq[j] );
								genL__align11( dynamicmtx, mseq1, mseq2, alloclen, &off1, &off2 );
							}
						}
#endif
					}
					break;
				case( 't' ):
					off1 = off2 = 0;
//					pscore = G__align11_noalign( n_dis_consweight_multi, penalty, penalty_ex, mseq1, mseq2, alloclen );
					pscore = G__align11_noalign( n_dis_consweight_multi, penalty, penalty_ex, distseq1, distseq2, alloclen ); // tsuneni distseq shiyou
					break;
				case( 's' ):
					pscore = callmxscarna_giving_bpp( mseq1, mseq2, bpp[i], bpp[j], alloclen, i, j );
					off1 = off2 = 0;
					break;
				case( 'G' ):
					pscore = calldafs_giving_bpp( mseq1, mseq2, bpp[i], bpp[j], alloclen, i, j );
					off1 = off2 = 0;
					break;
#if 0 
				case( 'a' ):
					pscore = Aalign( mseq1, mseq2, effarr1, effarr2, clus1, clus2, alloclen );
					off1 = off2 = 0;
					break;
				case( 'K' ):
					pscore = genG__align11( n_dis_consweight_multi, mseq1, mseq2, alloclen );
					off1 = off2 = 0;
					break;
				case( 'H' ):
					pscore = recallpairfoldalign( mseq1, mseq2, i, j, &off1, &off2, alloclen );
					break;
				case( 'B' ):
				case( 'T' ):
					pscore = recalllara( mseq1, mseq2, alloclen );
					off1 = off2 = 0;
					break;
				case( 'M' ):
					pscore = MSalign11( mseq1, mseq2, alloclen );
					break;
#endif
				default:
					ErrorExit( "\n\nERROR IN SOURCE FILE\n\n" );
			}
		}

		if( alg == 't' || ( mseq1[0][0] != 0 && mseq2[0][0] != 0  ) ) // 't' no jouken ha iranai to omou. if( ( mseq1[0][0] != 0 && mseq2[0][0] != 0  ) )
		{
#if SCOREOUT
			fprintf( stderr, "score = %10.2f (%d,%d)\n", pscore, i, j );
#endif
//			if( pscore > 0.0 && ( nadd == 0 || ( alg != 'Y' && alg != 'r' ) || ( i < njob-nadd && njob-nadd <= j ) ) ) x-ins-i de seido teika
			if( ( nadd == 0 || ( alg != 'Y' && alg != 'r' ) || ( i < njob-nadd && njob-nadd <= j ) ) )
			{
				if( !store_localhom )
					;
				else if( specifictarget && targetmap[i] == -1 && targetmap[j] == -1)
					;
				else if( alg == 'R' )
					putlocalhom_last( mseq1[0], mseq2[0], localhomtable[i]+j, lastresx[i]+j, 'h' );
				else if( alg == 'r' )
					putlocalhom_last( mseq1[0], mseq2[0], localhomtable[i]+j-(njob-nadd), lastresx[i]+j-(njob-nadd), 'h' );// ?????
				else if( alg == 'H' )
					putlocalhom_ext( mseq1[0], mseq2[0], localhomtable[i]+j, off1, off2, (int)pscore, strlen( mseq1[0] ), 'h' );
				else if( alg == 'Y' )
					putlocalhom2( mseq1[0], mseq2[0], localhomtable[i]+j-(njob-nadd), off1, off2, (int)pscore, strlen( mseq1[0] ), 'h' );
				else if( !specifictarget && alg != 'S' && alg != 'V' )
					putlocalhom2( mseq1[0], mseq2[0], localhomtable[i]+j-i, off1, off2, (int)pscore, strlen( mseq1[0] ), 'h' );
				else
//					putlocalhom2( mseq1[0], mseq2[0], localhomtable[i]+j, off1, off2, (int)pscore, strlen( mseq1[0] ) );
				{
					if( targetmap[i] != -1 && targetmap[j] != -1 )
					{
						putlocalhom2( mseq2[0], mseq1[0], localhomtable[targetmap[j]]+i, off2, off1, (int)pscore, strlen( mseq2[0] ), 'h' );
						putlocalhom2( mseq1[0], mseq2[0], localhomtable[targetmap[i]]+j, off1, off2, (int)pscore, strlen( mseq1[0] ), 'h' ); // sukoshi muda.
					}
					else if( targetmap[j] != -1 )
						putlocalhom2( mseq2[0], mseq1[0], localhomtable[targetmap[j]]+i, off2, off1, (int)pscore, strlen( mseq2[0] ), 'h' );
					else if( targetmap[i] != -1 )
						putlocalhom2( mseq1[0], mseq2[0], localhomtable[targetmap[i]]+j, off1, off2, (int)pscore, strlen( mseq1[0] ), 'h' );
#if 0
					if( targetmap[i] != -1 )
						putlocalhom2( mseq1[0], mseq2[0], localhomtable[targetmap[i]]+j, off1, off2, (int)pscore, strlen( mseq1[0] ), 'h' );
					
					else if( targetmap[j] != -1 )
						putlocalhom2( mseq2[0], mseq1[0], localhomtable[targetmap[j]]+i, off2, off1, (int)pscore, strlen( mseq2[0] ), 'h' );
#endif
					else
					{
						reporterr( "okashii\n" );
						exit( 1 );
					}
				}
			}
			pscore = score2dist( pscore, selfscore[i], selfscore[j] );

//			pscore = (double)naivepairscore11( mseq1[0], mseq2[0], 0.0 );
//			pscore = score2dist( pscore, selfscore[i], selfscore[j] );
//			reporterr( "->pscore = %f\n", pscore );

		}
		else
		{
			pscore = 2.0;
		}

#if 1 // mutex
		if( stdout_align )
		{
			pthread_mutex_lock( targ->mutex_stdout );
			if( alg != 't' )
			{
				fprintf( stdout, "sequence %d - sequence %d, pairwise distance = %10.5f\n", i+1, j+1, pscore );
				fprintf( stdout, ">%s\n", name[i] );
				write1seq( stdout, mseq1[0] );
				fprintf( stdout, ">%s\n", name[j] );
				write1seq( stdout, mseq2[0] );
				fprintf( stdout, "\n" );
			}
			pthread_mutex_unlock( targ->mutex_stdout );
		}
		if( stdout_dist )
		{
			pthread_mutex_lock( targ->mutex_stdout );
			if( j == i+1 ) fprintf( stdout, "%d %d d=%.3f\n", i+1, i+1, 0.0 );
			fprintf( stdout, "%d %d d=%.3f\n", i+1, j+1, pscore );
			pthread_mutex_unlock( targ->mutex_stdout );
		}
#endif // mutex
		if( store_dist )
		{
			if( alg == 'Y' || alg == 'r' ) distancemtx[i][j-(njob-nadd)] = pscore;
			else distancemtx[i][j-i] = pscore;
		}
	}
}
#endif

//Aligns seq with appropriate algorithm based on alg value.Here are the steps:
// 1. Initialize local variables.
// 2. Initialize targetmap and targetmapr. - save references to sequences to focus on in these arrays.
// 3. Initialize localhomtable
// 4. Initialize distance matrix.
// 5. If alg == R, call LAST web service to find similar regions between sequences and save result to lastresx
// 6. If alg == r,  call LAST web service to find similar regions between sequences with some difference than R and save the result to lastresx
// 7. If alg == H, call foldalign command with sequences as input and save output to _foldalignout file.
// 8. If alg == B, execute mafft_lara on input sequences and save the result to _laraout
// 9. If alg == T, execute lara with -s option on input sequences and save the result to _laraout
// 10. If alg == s, read hat4 file into bop 3d matrix then print Running MXSCARNA
// 11. If alg == G, read hat4 file into bop 3d matrix then print Running DAFS
// 12. Accumulate amino acids distances to pscore based on initial chars in seq input.
// 13. Initialize distance matrix value.
// 14. If use_fft == 1, call FAlign to align each two sequences with FFT algorithm. Else, switch on alg value.
// 15. If alg == t, call G__align11_noalign to align sequences.
// 16. If alg == A, call G__align11 / G__align11_noalign /  to align sequences and calculate score using naivepairscore11 / score2dist / makedynamicmtx scoring. All this based on the parameters used.
// 17. If alg == N, call genL__align11 and calculate score using naivepairscore11 / score2dist / makedynamicmtx to align sequences.
// 18. If alg == R, take score from lastresx - which was calculated before -.
// 19. If alg == r, take score from lastresx - which was calculated before -.
// 20. If alg == L, call L__align11 / L__align11_noalign to align sequences and use naivepairscore11 / score2dist / makedynamicmtx to calculate score.
// 21. If alg == Y, call L__align11 / L__align11_noalign to align sequences and use naivepairscore11to calculate score.
// 22. If alg == a, call Aalign to align sequences.
// 23. If alg == H, call recallpairfoldalign -> read fold align file that was filled in previous step -> align mseq1 and mseq2 based on it -> calculate alignment score -> copy new aligned sequences to mseq1 and mseq2 -> return alignment score
// 24. If alg == B or T, call recalllara which reads sequences from lara file and if matched with given ones then calculate naive score and return it.
// 25. If alg == s, call callmxscarna_giving_bpp to align sequences and use naive score to score them.
// 26. If alg == G, call calldafs_giving_bpp to align sequences then use naive score to score them.
// 27. If alg == M, call MSalign11 to align sequences
// 28. Update localhomtable values based on alg type
// 29. If store_dist, write name and distance matrix to hat2 file.
// 30. If store_localhom, write localhomtable to hat3 file.
// 31. Free all allocated matrices
static void pairalign( char **name, int *nlen, char **seq, char **aseq, char **dseq, int *thereisxineachseq, char **mseq1, char **mseq2, int alloclen, Lastresx **lastresx, double **distancemtx, LocalHom **localhomtable, int ngui )
{
	int i, j, ilim, jst, jj;
	int off1, off2, dum1, dum2, thereisx;
	double pscore = 0.0; // by D.Mathog
	FILE *hat2p, *hat3p;
//	double **distancemtx;
	double *selfscore;
	double *effarr1;
	double *effarr2;
	char *pt;
	char *hat2file = "hat2";
//	LocalHom **localhomtable = NULL, 
	LocalHom *tmpptr; //LocalHom is a structure defined in mltaln.h
	int intdum;
	char ***bpp = NULL; // mxscarna no toki dake
	char **distseq1, **distseq2;
	char **dumseq1, **dumseq2;
	double dist;
	double scoreoffset;
	int ntarget;
	int *targetmap, *targetmapr;

	//save references to sequences in focus in targetmap and targetmapr
	if( specifictarget ) //default = 0, and set to 1 if input argument inserted
	{
		targetmap = calloc( njob, sizeof( int ) );
		ntarget = 0;
		for( i=0; i<njob; i++ )
		{
			targetmap[i] = -1;
			if( !strncmp( name[i]+1, "_focus_", 7 ) ) //if name[i] contains _focus_
				targetmap[i] = ntarget++;
		}
		targetmapr = calloc( ntarget, sizeof( int ) );
		for( i=0; i<njob; i++ )
			if( targetmap[i] != -1 ) targetmapr[targetmap[i]] = i; //save indices of name array that contained target _focus_ in targetmapr

		if( ntarget == 0 )
		{
			reporterr( "\n\nAdd '>_focus_' to the title lines of the sequences to be focused on.\n\n" );
			exit( 1 );
		}
		else
		{
			reporterr( "nfocus = %d \n", ntarget );
		}
	}
	else
	{
		ntarget = njob; //target all sequences
		targetmap = calloc( njob, sizeof( int ) );
		targetmapr = calloc( njob, sizeof( int ) );
		for( i=0; i<njob; i++ )
			targetmap[i] = targetmapr[i] = i;
	}

#if 0
	for( i=0; i<njob; i++ )
		reporterr( "targetmap[%d] = %d\n", i, targetmap[i] );
	for( i=0; i<ntarget; i++ )
		reporterr( "targetmapr[%d] = %d\n", i, targetmapr[i] );
#endif

	if( store_localhom && localhomtable == NULL ) //store_localhom default = 1 and if set as from arguments = 0
	{
		if( alg == 'Y' || alg == 'r' )
		{
			ilim = njob - nadd;
			jst = nadd;
		}
		else
		{
			ilim = ntarget;
			jst = njob;
		}
		localhomtable = (LocalHom **)calloc( ilim, sizeof( LocalHom *) );
		for( i=0; i<ilim; i++)
		{
			localhomtable[i] = (LocalHom *)calloc( jst, sizeof( LocalHom ) );
			for( j=0; j<jst; j++)
			{
				localhomtable[i][j].start1 = -1;
				localhomtable[i][j].end1 = -1;
				localhomtable[i][j].start2 = -1; 
				localhomtable[i][j].end2 = -1; 
				localhomtable[i][j].opt = -1.0;
				localhomtable[i][j].next = NULL;
				localhomtable[i][j].nokori = 0;
				localhomtable[i][j].extended = -1;
				localhomtable[i][j].last = localhomtable[i]+j;
				localhomtable[i][j].korh = 'h';
			}
			if( !specifictarget && alg != 'Y' && alg != 'r' ) jst--;
		}
	}

	if( store_dist ) //default = 1, and if set from parameter = 0
	{
		if( ngui == 0 )
		{
			if( alg == 'Y' || alg == 'r' )
				distancemtx = AllocateDoubleMtx( njob, nadd );
			else
				distancemtx = AllocateDoubleHalfMtx( njob );
//				distancemtx = AllocateDoubleMtx( njob, njob );
		}
	}
	else distancemtx = NULL;

	if( alg == 'N' )
	{
		dumseq1 = AllocateCharMtx( 1, alloclen+10 );
		dumseq2 = AllocateCharMtx( 1, alloclen+10 );
	}
	distseq1 = AllocateCharMtx( 1, 0 ); // muda
	distseq2 = AllocateCharMtx( 1, 0 ); // muda

	selfscore = AllocateDoubleVec( njob );
	effarr1 = AllocateDoubleVec( njob );
	effarr2 = AllocateDoubleVec( njob );

#if 0
	fprintf( stderr, "##### fftwinsize = %d, fftthreshold = %d\n", fftWinSize, fftThreshold );
#endif

#if 0
	for( i=0; i<njob; i++ )
		fprintf( stderr, "TBFAST effarr[%d] = %f\n", i, effarr[i] );
#endif


//	writePre( njob, name, nlen, aseq, 0 );

	reporterr( "All-to-all alignment.\n" );
	if( alg == 'R' )
	{
		fprintf( stderr, "Calling last (http://last.cbrc.jp/)\n" ); //LAST is a web service that finds similar regions between sequences
		if( lastonce ) //default val = 0, and if argument set = 1
			calllast_once( njob, seq, njob, seq, lastresx ); //defined here. call lastal service and save result to lastresx
		else
			calllast_fast( njob, seq, njob, seq, lastresx ); //defined here. call lastal%d service and save result to lastresx
		fprintf( stderr, "done.\n" );
//		nthread = 0; // igo multithread nashi
	}
	if( alg == 'r' )
	{
		fprintf( stderr, "Calling last (http://last.cbrc.jp/)\n" ); //LAST is a web service that finds similar regions between sequences
		fprintf( stderr, "nadd=%d\n", nadd );
#if 1 // last_fast ha, lastdb ga muda
		if( lastonce ) //default val = 0, and if argument set = 1
			calllast_once( njob-nadd, seq, nadd, seq+njob-nadd, lastresx ); //defined here. call lastal service and save result to lastresx
		else
			calllast_fast( njob-nadd, seq, nadd, seq+njob-nadd, lastresx ); //defined here. call lastal%d service and save result to lastresx
#else
		calllast_once( njob-nadd, seq, nadd, seq+njob-nadd, lastresx );
#endif

		fprintf( stderr, "nadd=%d\n", nadd );
		fprintf( stderr, "done.\n" );
//		nthread = 0; // igo multithread nashi
	}

	if( alg == 'H' )
	{
		fprintf( stderr, "Calling FOLDALIGN with option '%s'\n", foldalignopt ); //foldalignopt is a constant char array, set from arguments
		callfoldalign( njob, seq ); //defined here. it calls foldalign command with sequences as input and output is saved in _foldalignout file
		fprintf( stderr, "done.\n" );
	}
	if( alg == 'B' )
	{
		fprintf( stderr, "Running LARA (Bauer et al. http://www.planet-lisa.net/)\n" ); //I couldn't find this site
		calllara( njob, seq, "" ); //defined here. calls mafft_lara and executes it on sequences then saves output to _laraout
		fprintf( stderr, "done.\n" );
	}
	if( alg == 'T' )
	{
		fprintf( stderr, "Running SLARA (Bauer et al. http://www.planet-lisa.net/)\n" );
		calllara( njob, seq, "-s" ); //defined here. calls mafft_lara and executes it on sequences with -s option then saves output to _laraout
		fprintf( stderr, "done.\n" );
	}
	if( alg == 's' )
	{
		fprintf( stderr, "Preparing bpp\n" );
//		bpp = AllocateCharCub( njob, nlenmax, 0 );
		bpp = calloc( njob, sizeof( char ** ) );
		preparebpp( njob, bpp ); //defined here. it reads hat4 file into bpp
		fprintf( stderr, "done.\n" );
		fprintf( stderr, "Running MXSCARNA (Tabei et al. http://www.ncrna.org/software/mxscarna)\n" );
		//mxscarna is a fast structural multiple alignment method for long RNA sequences
	}
	if( alg == 'G' )
	{
		fprintf( stderr, "Preparing bpp\n" );
//		bpp = AllocateCharCub( njob, nlenmax, 0 );
		bpp = calloc( njob, sizeof( char ** ) );
		preparebpp( njob, bpp ); //defined here. in reads hat4 file into bpp
		fprintf( stderr, "done.\n" );
		fprintf( stderr, "Running DAFS (Sato et al. http://www.ncrna.org/)\n" ); //ncrna Bioinformatics tools and databases for functional RNA analysis
	}

	for( i=0; i<njob; i++ )
	{
		pscore = 0.0;
		for( pt=seq[i]; *pt; pt++ )
			pscore += amino_dis[(unsigned char)*pt][(unsigned char)*pt]; //accumulate amino acids distances to pscore
		selfscore[i] = pscore;
//		fprintf( stderr, "selfscore[%d] = %f\n", i, selfscore[i] );
	}

#if enablemultithread
	if( nthread > 0 ) // alg=='r' || alg=='R' -> nthread:=0 (sukoshi ue)
	{
		Jobtable jobpos;
		pthread_t *handle;
		pthread_mutex_t mutex_counter;
		pthread_mutex_t mutex_stdout;
		thread_arg_t *targ;

		if( alg == 'Y' || alg == 'r' ) jobpos.j = njob - nadd - 1;
		else jobpos.j = 0;
		jobpos.i = 0;

		targ = calloc( nthread, sizeof( thread_arg_t ) );
		handle = calloc( nthread, sizeof( pthread_t ) );
		pthread_mutex_init( &mutex_counter, NULL );
		pthread_mutex_init( &mutex_stdout, NULL );

		for( i=0; i<nthread; i++ )
		{
			targ[i].thread_no = i;
			targ[i].njob = njob;
			targ[i].jobpospt = &jobpos;
			targ[i].name = name;
			targ[i].seq = seq;
			targ[i].dseq = dseq;
			targ[i].thereisxineachseq = thereisxineachseq;
			targ[i].localhomtable = localhomtable;
			targ[i].distancemtx = distancemtx;
			targ[i].selfscore = selfscore;
			targ[i].bpp = bpp; 
			targ[i].lastresx = lastresx;
			targ[i].alloclen = alloclen;
			targ[i].targetmap = targetmap;
			targ[i].mutex_counter = &mutex_counter;
			targ[i].mutex_stdout = &mutex_stdout;

//			athread( (void *)targ );
			pthread_create( handle+i, NULL, athread, (void *)(targ+i) );
//			pthread_create( handle+i, NULL, bthread, (void *)(targ+i) );
		}


		for( i=0; i<nthread; i++ )
		{
			pthread_join( handle[i], NULL );
		}
		pthread_mutex_destroy( &mutex_counter );
		pthread_mutex_destroy( &mutex_stdout );
		free( handle );
		free( targ );
	}
	else
#endif
	{
		double **dynamicmtx = NULL;
		if( specificityconsideration > 0.0 ) dynamicmtx = AllocateDoubleMtx( nalphabets, nalphabets );

		if( alg == 'Y' || alg == 'r' ) ilim = njob - nadd;
		else ilim = njob - 1;
		for( i=0; i<ilim; i++ ) 
		{
			if( stdout_dist) fprintf( stdout, "%d %d d=%.3f\n", i+1, i+1, 0.0 );
			fprintf( stderr, "% 5d / %d\r", i, njob-nadd );
			fflush( stderr );

			if( alg == 'Y' || alg == 'r' ) jst = njob - nadd;
			else jst = i + 1;
			for( j=jst; j<njob; j++ )
			{
	
				if( strlen( seq[i] ) == 0 || strlen( seq[j] ) == 0 )
				{
					if( store_dist ) 
					{
						if( alg == 'Y' || alg == 'r' ) distancemtx[i][j-(njob-nadd)] = 3.0;
						else distancemtx[i][j-i] = 3.0;
					}
					if( stdout_dist) fprintf( stdout, "%d %d d=%.3f\n", i+1, j+1, 3.0 );
					continue;
				}
	
				strcpy( aseq[0], seq[i] );
				strcpy( aseq[1], seq[j] );
//				clus1 = conjuctionfortbfast( pair, i, aseq, mseq1, effarr1, effarr, indication1 );
//				clus2 = conjuctionfortbfast( pair, j, aseq, mseq2, effarr2, effarr, indication2 );
//				fprintf( stderr, "Skipping conjuction..\n" );

				effarr1[0] = 1.0;
				effarr2[0] = 1.0;
				mseq1[0] = aseq[0];
				mseq2[0] = aseq[1];

				thereisx = thereisxineachseq[i] + thereisxineachseq[j];
//				strcpy( distseq1[0], dseq[i] ); // nen no tame
//				strcpy( distseq2[0], dseq[j] ); // nen no tame
				distseq1[0] = dseq[i];
				distseq2[0] = dseq[j];

	//			fprintf( stderr, "mseq1 = %s\n", mseq1[0] );
	//			fprintf( stderr, "mseq2 = %s\n", mseq2[0] );
		
#if 0
				fprintf( stderr, "group1 = %.66s", indication1 );
				fprintf( stderr, "\n" );
				fprintf( stderr, "group2 = %.66s", indication2 );
				fprintf( stderr, "\n" );
#endif
	//			for( l=0; l<clus1; l++ ) fprintf( stderr, "## STEP-eff for mseq1-%d %f\n", l, effarr1[l] );
	
				if( use_fft ) //default = 0, and if set from arguments = 1 //if using FFT, then call Falign
				{
					//defined in Falign.c.
					//Falign aligns mseq1 and mseq2 based on FFT and algorithm selected
					pscore = Falign( NULL, NULL, n_dis_consweight_multi, mseq1, mseq2, effarr1, effarr2, NULL, NULL, 1, 1, alloclen, &intdum, NULL, 0, NULL );
//					fprintf( stderr, "pscore (fft) = %f\n", pscore );
					off1 = off2 = 0;
				}
				else
				{
					switch( alg ) //defined in defs.h. I need to know what those algorithms stand for ?
					{
						case( 't' ):
//							pscore = G__align11_noalign( n_dis_consweight_multi, penalty, penalty_ex, mseq1, mseq2, alloclen );
							pscore = G__align11_noalign( n_dis_consweight_multi, penalty, penalty_ex, distseq1, distseq2, alloclen ); // tsuneni distseq shiyou
							//defined in Galign11.c. distseq1 and distseq2 are sequences without gaps and x char.
							//it is like G__align11 but without warp and some other details
							//aligns distseq1 and distseq2 with some specific algorithm
							off1 = off2 = 0;
							break;
						case( 'A' ):
							if( usenaivescoreinsteadofalignmentscore ) //use naive score instead of alignment score
							{
								//Calculates distance between mseq1 and mseq2 based on specific algo.
								G__align11( n_dis_consweight_multi, mseq1, mseq2, alloclen, outgap, outgap ); //defined in Galign11.c.
								pscore = (double)naivepairscore11( mseq1[0], mseq2[0], 0.0 ); // uwagaki //defined in mltaln9.c.
								//calculates score between mseq1[0] and mseq2[0] based on penal and amino_dis values
							}
							else
							{
//								if( store_localhom )
								if( store_localhom && ( targetmap[i] != -1 || targetmap[j] != -1 ) ) //store_localhom default = 1 and if set as from arguments = 0
								{
									//Calculates distance between mseq1 and mseq2 based on specific algo.
									pscore = G__align11( n_dis_consweight_multi, mseq1, mseq2, alloclen, outgap, outgap ); //defined in Galign11.c.
									if( thereisx ) pscore = G__align11_noalign( n_dis_consweight_multi, penalty, penalty_ex, distseq1, distseq2, alloclen ); // uwagaki
									//defined in Galign11.c. distseq1 and distseq2 are sequences without gaps and x char.
									//it is like G__align11 but without warp and some other details
									//aligns distseq1 and distseq2 with some specific algorithm
#if 1
									if( specificityconsideration > 0.0 ) //defined in defs.c
									{
										dist = score2dist( pscore, selfscore[i], selfscore[j] ); //defined here. returns value based on the three args values.
//										dist = score2dist( L__align11_noalign( n_dis_consweight_multi, distseq1, distseq2 ), selfscore[i], selfscore[j] ); // 2014/Feb/20
										if( dist2offset( dist ) < 0.0 ) //defined in mltaln9.c. returns value based on some calcs on dist value
										{
											//fill dynamicmtx matrix based on n_dis_consweight_multi matrix values and offset value
											makedynamicmtx( dynamicmtx, n_dis_consweight_multi, 0.5 * dist ); // upgma ni awaseru. //defined in mltaln9.c.
											strcpy( mseq1[0], seq[i] ); //copy seq to mseq
											strcpy( mseq2[0], seq[j] );
											//Calculates distance between mseq1 and mseq2 based on specific algo.
											G__align11( dynamicmtx, mseq1, mseq2, alloclen, outgap, outgap ); //defined in Galign11.c.
										}
//										pscore = (double)naivepairscore11( *mseq1, *mseq2, 0.0 );
									}
#endif
								}
								else
									pscore = G__align11_noalign( n_dis_consweight_multi, penalty, penalty_ex, distseq1, distseq2, alloclen ); // uwagaki
									//defined in Galign11.c. distseq1 and distseq2 are sequences without gaps and x char.
									//it is like G__align11 but without warp and some other details
									//aligns distseq1 and distseq2 with some specific algorithm
							}
							off1 = off2 = 0;
							break;
						case( 'N' ):
//							pscore = G__align11_noalign( n_dis_consweight_multi, penalty, penalty_ex, mseq1, mseq2, alloclen );
							if( usenaivescoreinsteadofalignmentscore ) //use naive score instead of alignment score
							{
								genL__align11( n_dis_consweight_multi, mseq1, mseq2, alloclen, &off1, &off2 ); //defined in genalign11.c.
								//Calculates distance between mseq1 and mseq2 based on specific algo.
								//It is like G__align11 with some small differences specially absence of warp
								pscore = (double)naivepairscore11( mseq1[0], mseq2[0], 0.0 ); // uwagaki  //defined in mltaln9.c.
								//calculates score between mseq1[0] and mseq2[0] based on penal and amino_dis values
							}
							else
							{
								//Calculates distance between mseq1 and mseq2 based on specific algo.
								//It is like G__align11 with some small differences specially absence of warp
								pscore = genL__align11( n_dis_consweight_multi, mseq1, mseq2, alloclen, &off1, &off2 ); //defined in genalign11.c.
								if( thereisx )
								{
									strcpy( dumseq1[0], distseq1[0] );
									strcpy( dumseq2[0], distseq2[0] );
									//Calculates distance between dumseq1 and dumseq2 based on specific algo.
									pscore = genL__align11( n_dis_consweight_multi, dumseq1, dumseq2, alloclen, &dum1, &dum2 ); // uwagaki
								}
#if 1
								if( specificityconsideration > 0.0 )
								{
//									fprintf( stderr, "dist = %f\n", score2dist( pscore, selfscore[i], selfscore[j] ) );
									dist = score2dist( pscore, selfscore[i], selfscore[j] ); //defined here. returns value based on the three args values.
									if( dist2offset( dist ) < 0.0 ) //defined in mltaln9.c. returns value based on some calcs on dist value
									{
										//fill dynamicmtx matrix based on n_dis_consweight_multi matrix values and offset value
										makedynamicmtx( dynamicmtx, n_dis_consweight_multi, 0.5 * dist ); // upgma ni awaseru.  //defined in mltaln9.c.
										strcpy( mseq1[0], seq[i] ); //copy seq to mseq
										strcpy( mseq2[0], seq[j] );
										//Calculates distance between dumseq1 and dumseq2 based on specific algo.
										genL__align11( dynamicmtx, mseq1, mseq2, alloclen, &off1, &off2 ); //defined in genalign11.c.
									}
								}
#endif
							}
							break;
						case( 'R' ):
							if( nadd && njob-nadd <= j && njob-nadd <= i ) // new sequence doushi ha mushi
								pscore = 0.0;
							else
								pscore = (double)lastresx[i][j].score; // all pair
							break;
						case( 'r' ):
							if( nadd == 0 || ( i < njob-nadd && njob-nadd <= j ) )
								pscore = (double)lastresx[i][j-(njob-nadd)].score;
							else
								pscore = 0.0;
							break;
						case( 'L' ):
							if( nadd && njob-nadd <= j && njob-nadd <= i ) // new sequence doushi ha mushi
								pscore = 0.0;
							else
							{
								if( usenaivescoreinsteadofalignmentscore ) //use naive score instead of alignment score
								{
									//This method finds the matchings between mseq1 and mseq2 and returns the score of matching them.
									//It is similar to G__align11 with some small differences.
									L__align11( n_dis_consweight_multi, 0.0, mseq1, mseq2, alloclen, &off1, &off2 ); //defined in lalign11.c.
									pscore = (double)naivepairscore11( mseq1[0], mseq2[0], 0.0 ); // uwagaki  //defined in mltaln9.c.
									//calculates score between mseq1[0] and mseq2[0] based on penal and amino_dis values
								}
								else
								{
//									if( store_localhom )
									if( store_localhom && ( targetmap[i] != -1 || targetmap[j] != -1 ) ) //store_localhom default = 1 and if set as from arguments = 0
									{
										//This method finds the matchings between mseq1 and mseq2 and returns the score of matching them.
										//It is similar to G__align11 with some small differences.
										pscore = L__align11( n_dis_consweight_multi, 0.0, mseq1, mseq2, alloclen, &off1, &off2 ); // all pair  //defined in lalign11.c.
										if( thereisx ) pscore = L__align11_noalign( n_dis_consweight_multi, distseq1, distseq2 ); // all pair  //defined in lalign11.c.
										//distseq1 and distseq2 are sequences without gaps and x char.
										//This method finds the matchings between distseq1 and distseq2 and returns the score of matching them.
										//It is similar to L__align11 with some small differences. It is simpler and doesn't contain warp
#if 1
										if( specificityconsideration > 0.0 )
										{
											dist = score2dist( pscore, selfscore[i], selfscore[j] ); //defined here. returns value based on the three args values.
											if( ( scoreoffset = dist2offset( dist ) ) < 0.0 ) //defined in mltaln9.c. returns value based on some calcs on dist value
											{
												//fill dynamicmtx matrix based on n_dis_consweight_multi matrix values and offset value
												makedynamicmtx( dynamicmtx, n_dis_consweight_multi, 0.5 * dist ); // upgma ni awaseru.  //defined in mltaln9.c.
												strcpy( mseq1[0], seq[i] );
												strcpy( mseq2[0], seq[j] );
												//This method finds the matchings between mseq1 and mseq2 and returns the score of matching them.
												//It is similar to G__align11 with some small differences.
												L__align11( dynamicmtx, scoreoffset, mseq1, mseq2, alloclen, &off1, &off2 ); //defined in lalign11.c.
											}
										}
#endif
									}
									else
										pscore = L__align11_noalign( n_dis_consweight_multi, distseq1, distseq2 ); // all pair  //defined in lalign11.c.
										//distseq1 and distseq2 are sequences without gaps and x char.
										//This method finds the matchings between distseq1 and distseq2 and returns the score of matching them.
										//It is similar to L__align11 with some small differences. It is simpler and doesn't contain warp.
								}
							}
//							pscore = G__align11_noalign( n_dis_consweight_multi, penalty, penalty_ex, distseq1, distseq2, alloclen ); // CHUUI!!!!!!
							break;
						case( 'Y' ):
							if( nadd == 0 || ( i < njob-nadd && njob-nadd <= j ) ) // new sequence vs exiting sequence nomi keisan
							{
								if( usenaivescoreinsteadofalignmentscore ) //use naive score instead of alignment score
								{
									//This method finds the matchings between mseq1 and mseq2 and returns the score of matching them.
									//It is similar to G__align11 with some small differences.
									L__align11( n_dis_consweight_multi, 0.0, mseq1, mseq2, alloclen, &off1, &off2 ); //defined in lalign11.c.
									pscore = (double)naivepairscore11( mseq1[0], mseq2[0], 0.0 ); // uwagaki  //defined in mltaln9.c.
									//calculates score between mseq1[0] and mseq2[0] based on penal and amino_dis values
								}
								else
								{
									if( store_localhom )
									{
										//This method finds the matchings between mseq1 and mseq2 and returns the score of matching them.
										//It is similar to G__align11 with some small differences.
										pscore = L__align11( n_dis_consweight_multi, 0.0, mseq1, mseq2, alloclen, &off1, &off2 ); //defined in lalign11.c.
										if( thereisx ) pscore = L__align11_noalign( n_dis_consweight_multi, distseq1, distseq2 ); // uwagaki   //defined in lalign11.c.
										//distseq1 and distseq2 are sequences without gaps and x char.
										//This method finds the matchings between distseq1 and distseq2 and returns the score of matching them.
										//It is similar to L__align11 with some small differences. It is simpler and doesn't contain warp.
									}
									else
										pscore = L__align11_noalign( n_dis_consweight_multi, distseq1, distseq2 );  //defined in lalign11.c.
										//distseq1 and distseq2 are sequences without gaps and x char.
										//This method finds the matchings between distseq1 and distseq2 and returns the score of matching them.
										//It is similar to L__align11 with some small differences. It is simpler and doesn't contain warp.
								}
							}
							else
								pscore = 0.0;
							break;
						case( 'a' ):
							//apply specific algorithm to align mseq1 and mseq2
							pscore = Aalign( mseq1, mseq2, effarr1, effarr2, 1, 1, alloclen ); //defined in SAalignmm.c.
							off1 = off2 = 0;
							break;
#if 0
						case( 'K' ):
							pscore = genG__align11( mseq1, mseq2, alloclen );
							off1 = off2 = 0;
							break;
#endif
						case( 'H' ):
							//read fold align file -> align mseq1 and mseq2 based on it -> calculate alignment score -> copy new aligned sequences to mseq1 and mseq2 -> return alignment score
							pscore = recallpairfoldalign( mseq1, mseq2, i, j, &off1, &off2, alloclen ); //defined here.
							break;
						case( 'B' ):
						case( 'T' ):
							//read sequences from lara file and if matched with given ones then calculate naive score and return it.
							pscore = recalllara( mseq1, mseq2, alloclen ); //defined here.
							off1 = off2 = 0;
							break;
						case( 's' ):
							//executes some commands then write mseq1 and mseq2 to file and align them, then apply naive score and return it.
							pscore = callmxscarna_giving_bpp( mseq1, mseq2, bpp[i], bpp[j], alloclen, i, j ); //defined here.
							off1 = off2 = 0;
							break;
						case( 'G' ):
							//executes some commands then write mseq1 and mseq2 to file and align them, then apply naive score and return it.
							pscore = calldafs_giving_bpp( mseq1, mseq2, bpp[i], bpp[j], alloclen, i, j ); //defined here.
							off1 = off2 = 0;
							break;
						case( 'M' ):
							//aligns seq1 and seq2
							pscore = MSalign11( mseq1, mseq2, alloclen );  //defined in MSalign11.c.
							break;
						default:
							ErrorExit( "ERROR IN SOURCE FILE" );
					}
				}
	

				if( alg == 't' || ( mseq1[0][0] != 0 && mseq2[0][0] != 0  ) ) // 't' no jouken ha iranai to omou. if( ( mseq1[0][0] != 0 && mseq2[0][0] != 0  ) )
				{
#if SCOREOUT
					fprintf( stderr, "score = %10.2f (%d,%d)\n", pscore, i, j );
#endif
//					if( pscore > 0.0 && ( nadd == 0 || ( alg != 'Y' && alg != 'r' ) || ( i < njob-nadd && njob-nadd <= j ) ) ) // x-ins-i de seido teika
					if( ( nadd == 0 || ( alg != 'Y' && alg != 'r' ) || ( i < njob-nadd && njob-nadd <= j ) ) )
					{
						if( !store_localhom )
							;
						else if( specifictarget && targetmap[i] == -1 && targetmap[j] == -1)
							;
						else if( alg == 'R' )
							//This method calculates some scores based on mseq1[0], mseq2[0] and amino_n values and stores them in localhomtable[i]+j
							//It also updates some references in localhomtable[i]+j based on other lastresx[i]+j values.
							putlocalhom_last( mseq1[0], mseq2[0], localhomtable[i]+j, lastresx[i]+j, 'h' ); //defined here.
						else if( alg == 'r' )
							//This method calculates some scores based on mseq1[0], mseq2[0] and amino_n values and stores them in localhomtable[i]+j-(njob-nadd)
							//It also updates some references in localhomtable[i]+j-(njob-nadd) based on other lastresx[i]+j-(njob-nadd) values.
							putlocalhom_last( mseq1[0], mseq2[0], localhomtable[i]+j-(njob-nadd), lastresx[i]+j-(njob-nadd), 'h' );// ?????
						else if( alg == 'H' )
							//This method is similar to putlocalhom2 except in localhomtable[i]+j->last assignment and localhomtable[i]+j->opt calculation
							putlocalhom_ext( mseq1[0], mseq2[0], localhomtable[i]+j, off1, off2, (int)pscore, strlen( mseq1[0] ), 'h' ); //defined in io.c.
						else if( alg == 'Y' )
							//I think this method scans the two sequences mseq1[0] and mseq2[1] and saves score and other alignment info in localhomtable[i]+j-(njob-nadd)
							//what i need to understand now is LocalHom mechanism and how it works exactly
							putlocalhom2( mseq1[0], mseq2[0], localhomtable[i]+j-(njob-nadd), off1, off2, (int)pscore, strlen( mseq1[0] ), 'h' ); //defined in io.c.
						else if( !specifictarget && alg != 'S' && alg != 'V' )
							//I think this method scans the two sequences mseq1[0] and mseq2[1] and saves score and other alignment info in localhomtable[i]+j-i
							//what i need to understand now is LocalHom mechanism and how it works exactly
							putlocalhom2( mseq1[0], mseq2[0], localhomtable[i]+j-i, off1, off2, (int)pscore, strlen( mseq1[0] ), 'h' ); //defined in io.c.
						else
						{
							if( targetmap[i] != -1 && targetmap[j] != -1 )
							{
								putlocalhom2( mseq1[0], mseq2[0], localhomtable[targetmap[i]]+j, off1, off2, (int)pscore, strlen( mseq1[0] ), 'h' );
								putlocalhom2( mseq2[0], mseq1[0], localhomtable[targetmap[j]]+i, off2, off1, (int)pscore, strlen( mseq2[0] ), 'h' ); // sukoshi muda.
							}
							else if( targetmap[j] != -1 )
								putlocalhom2( mseq2[0], mseq1[0], localhomtable[targetmap[j]]+i, off2, off1, (int)pscore, strlen( mseq2[0] ), 'h' );
							else if( targetmap[i] != -1 )
								putlocalhom2( mseq1[0], mseq2[0], localhomtable[targetmap[i]]+j, off1, off2, (int)pscore, strlen( mseq1[0] ), 'h' );
							else
							{
								reporterr( "okashii\n" );
								exit( 1 );
							}
						}
					}

					pscore = score2dist( pscore, selfscore[i], selfscore[j] ); //defined here. returns value based on the three args values.
				}
				else
				{
					pscore = 2.0;
				}
	
				if( stdout_align ) //defined here. set to 0 by default and to 1 if set from input arguments
				{
					if( alg != 't' )
					{
						fprintf( stdout, "sequence %d - sequence %d, pairwise distance = %10.5f\n", i+1, j+1, pscore );
						fprintf( stdout, ">%s\n", name[i] );
						write1seq( stdout, mseq1[0] ); //defined in io.c. writes mseq1[0] to stdout
						fprintf( stdout, ">%s\n", name[j] );
						write1seq( stdout, mseq2[0] ); //writes mseq2[0] to stdout
						fprintf( stdout, "\n" );
					}
				}
				if( stdout_dist ) fprintf( stdout, "%d %d d=%.3f\n", i+1, j+1, pscore ); //stdout_dist defined here. default = 0, from args = 1.
				if( store_dist) //defined here. default = 1, and if set from parameter = 0
				{
					if( alg == 'Y' || alg == 'r' ) distancemtx[i][j-(njob-nadd)] = pscore; //distancemtx is an argument to this method
					else distancemtx[i][j-i] = pscore;
				}
			}
		}
		if( dynamicmtx ) FreeDoubleMtx( dynamicmtx ); //free dynamicmtx if allocated
	}


	if( store_dist && ngui == 0 ) //store_dist defined here. default = 1, and if set from parameter = 0. ngui is an argument to this method
	{
		hat2p = fopen( hat2file, "w" ); //open hat2 file for write.
		if( !hat2p ) ErrorExit( "Cannot open hat2." );
		if( alg == 'Y' || alg == 'r' )
			WriteHat2_part_pointer( hat2p, njob, nadd, name, distancemtx ); //defined in io.c. write name and distancemtx content to hat2p file in specific format
		else
//			WriteHat2_pointer( hat2p, njob, name, distancemtx );
			WriteFloatHat2_pointer_halfmtx( hat2p, njob, name, distancemtx ); // jissiha double //defined in io.c. write name and distancemtx content to hat2p file in specific format
		fclose( hat2p );
	}

	hat3p = fopen( "hat3", "w" ); //open hat3 file for write.
	if( !hat3p ) ErrorExit( "Cannot open hat3." );
	if( store_localhom && ngui == 0 ) //store_localhom default = 1 and if set as from arguments = 0. ngui is an argument to this method
	{

		fprintf( stderr, "\n\n##### writing hat3\n" );
		if( alg == 'Y' || alg == 'r' )
			ilim = njob-nadd;	
		else if( specifictarget ) //defined here. default = 0 and if set from args = 1.
			ilim = ntarget;
		else
			ilim = njob-1;	
		for( i=0; i<ilim; i++ ) 
		{
			if( alg == 'Y' || alg == 'r' )
			{
				jst = njob-nadd;
				jj = 0;
			}
			else if( specifictarget )
			{
				jst = 0;
				jj = 0;
			}
			else
			{
				jst = i;
				jj = 0;
			}
			for( j=jst; j<njob; j++, jj++ )
			{
				for( tmpptr=localhomtable[i]+jj; tmpptr; tmpptr=tmpptr->next )
				{
//					fprintf( stderr, "j=%d, jj=%d\n", j, jj );
					if( tmpptr->opt == -1.0 ) continue; //go to next iteration
// tmptmptmptmptmp
//					if( alg == 'B' || alg == 'T' )
//						fprintf( hat3p, "%d %d %d %7.5f %d %d %d %d %p\n", i, j, tmpptr->overlapaa, 1.0, tmpptr->start1, tmpptr->end1, tmpptr->start2, tmpptr->end2, (void *)tmpptr->next ); 
//					else
					if( targetmap[j] == -1 || targetmap[i] < targetmap[j] ) //print values from targetmap and tmpptr to hat3 file
						fprintf( hat3p, "%d %d %d %7.5f %d %d %d %d h\n", targetmapr[i], j, tmpptr->overlapaa, tmpptr->opt, tmpptr->start1, tmpptr->end1, tmpptr->start2, tmpptr->end2 );
//						fprintf( hat3p, "%d %d %d %7.5f %d %d %d %d h\n", i, j, tmpptr->overlapaa, tmpptr->opt, tmpptr->start1, tmpptr->end1, tmpptr->start2+1, tmpptr->end2+1 ); // zettai dame!!!!
				}
			}
		}
//		if( ngui == 0 )
//		{
#if DEBUG
			fprintf( stderr, "calling FreeLocalHomTable\n" );
#endif
			if( alg == 'Y' || alg == 'r' )
				FreeLocalHomTable_part( localhomtable, (njob-nadd), nadd ); //free localhomtable
			else if( specifictarget )
				FreeLocalHomTable_part( localhomtable, ntarget, njob );
			else
				FreeLocalHomTable_half( localhomtable, njob );
#if DEBUG
			fprintf( stderr, "done. FreeLocalHomTable\n" );
#endif
//		}
	}
	fclose( hat3p );

	if( alg == 's' ) //free bpp matrix
	{
		char **ptpt;
		for( i=0; i<njob; i++ )
		{
			ptpt = bpp[i];
			while( 1 )
			{
				if( *ptpt ) free( *ptpt );
				else break;
				ptpt++;
			}
			free( bpp[i] );
		}
		free( bpp );
	}

	//free all other allocated memory
	free( selfscore );
	free( effarr1 );
	free( effarr2 );
	if( alg == 'N' )
	{
		FreeCharMtx( dumseq1 );
		FreeCharMtx( dumseq2 );
	}
	free( distseq1 );
	free( distseq2 );
	if( store_dist && ngui == 0 ) FreeDoubleHalfMtx( distancemtx, njob ); //defined in mtxutl.c. frees distancemtx memory

	free( targetmap );
	free( targetmapr ); //ok, i stopped here and now need to revise this whole method and write its steps then go back to its caller to continue.
}

// 1. Initialize variables
// 2. Read input arguments to this main method
// 3. If no sequences number is provided, get it from input file.
// 4. If alg == R or r, allocate lasersx matrix.
// 5. If sequences are not given, read them from input file.
// 6. Call constants.c to fill weighting and scoring matrices.
// 7. Inits prep_g and trap_g files for tracing
// 8. Save both sequences without gaps and sequences without gaps and x chars. Also save number of x chars in each sequence.
// 9. Then call the giant method pairalign to do all the alignment effort :D
// 10. Then finally deallocate all memory.
int pairlocalalign( int ngui, int lgui, char **namegui, char **seqgui, double **distancemtx, LocalHom **localhomtable, int argc, char **argv )
{
	int  *nlen, *thereisxineachseq; //there is x in each seq :D
	char **name, **seq;
	char **mseq1, **mseq2;
	char **aseq;
	char **bseq;
	char **dseq;
	int i, j, k;
	FILE *infp;
	char c;
	int alloclen;
	Lastresx **lastresx; //structure defined here.

//	reporterr( "argc=%d, argv[0]=%s\n", argc, argv[0] );

	pair_local_align_arguments( argc, argv ); //read input arguments


	if( !ngui ) //if no sequences number provided
	{
		if( inputfile )
		{
			infp = fopen( inputfile, "r" );
			if( !infp )
			{
				fprintf( stderr, "Cannot open %s\n", inputfile );
				exit( 1 );
			}
		}
		else
			infp = stdin;
	
		getnumlen( infp ); //defined in io.c. finds sequences count, max length and dna or protein from input file
		rewind( infp );
	
		if( njob < 2 )
		{
			fprintf( stderr, "At least 2 sequences should be input!\n"
							 "Only %d sequence found.\n", njob ); 
			exit( 1 );
		}
		if( njob > M ) //M is the maximum number of jobs(sequences) allowed. It is defined in mltaln.h and = 500,000
		{
			fprintf( stderr, "The number of sequences must be < %d\n", M );
			fprintf( stderr, "Please try the splittbfast program for such large data.\n" );
			exit( 1 );
		}
	}

	if( ( alg == 'r' || alg == 'R' ) && dorp == 'p' ) //R algorithm is not supported with protein yet
	{
		fprintf( stderr, "Not yet supported\n" );
		exit( 1 );
	}

	alloclen = nlenmax*2; //pairalign defined in defs.h.
	if( ngui ) 
	{
		seq = seqgui;
		name = namegui;
	}
	else
	{
		seq = AllocateCharMtx( njob, alloclen+10 );
		name = AllocateCharMtx( njob, B );
	}

	aseq = AllocateCharMtx( 2, alloclen+10 );
	bseq = AllocateCharMtx( njob, alloclen+10 );
	dseq = AllocateCharMtx( njob, alloclen+10 );
	mseq1 = AllocateCharMtx( njob, 0 );
	mseq2 = AllocateCharMtx( njob, 0 );
	nlen = AllocateIntVec( njob );
	thereisxineachseq = AllocateIntVec( njob ); //there is x in each sequence


	if( alg == 'R' )
	{
		lastresx = calloc( njob+1, sizeof( Lastresx * ) ); //Lastresx is a structure defined here
		for( i=0; i<njob; i++ ) 
		{
			lastresx[i] = calloc( njob+1, sizeof( Lastresx ) ); // muda   //muda = useless
			for( j=0; j<njob; j++ ) 
			{
				lastresx[i][j].score = 0;
				lastresx[i][j].naln = 0;
				lastresx[i][j].aln = NULL;
			}
			lastresx[i][njob].naln = -1;
		}
		lastresx[njob] = NULL;
	}
	else if( alg == 'r' )
	{
//		fprintf( stderr, "Allocating lastresx (%d), njob=%d, nadd=%d\n", njob-nadd+1, njob, nadd );
		lastresx = calloc( njob-nadd+1, sizeof( Lastresx * ) );
		for( i=0; i<njob-nadd; i++ )
		{
//			fprintf( stderr, "Allocating lastresx[%d]\n", i );
			lastresx[i] = calloc( nadd+1, sizeof( Lastresx ) );
			for( j=0; j<nadd; j++ ) 
			{
//				fprintf( stderr, "Initializing lastresx[%d][%d]\n", i, j );
				lastresx[i][j].score = 0;
				lastresx[i][j].naln = 0;
				lastresx[i][j].aln = NULL;
			}
			lastresx[i][nadd].naln = -1;
		}
		lastresx[njob-nadd] = NULL;
	}
	else
		lastresx = NULL;

#if 0
	Read( name, nlen, seq );
#else
	if( !ngui )  //if no sequences number provided
	{
		readData_pointer( infp, name, nlen, seq ); //defined in io.c. It reads sequences and their names in seq, name and nlen arrays.
		fclose( infp );
	}
#endif

	constants( njob, seq ); //defined in constants.c
	//after all this method, n_dis, ribosumdis, amino_dis, amino_dis_consweight_multi, n_dis_consweight_multi,
	//n_disLN, n_disFFT, polarity, volume arrays are initialized and some constants are set.

#if 0
	fprintf( stderr, "params = %d, %d, %d\n", penalty, penalty_ex, offset );
#endif

	initSignalSM(); //defined in io.c - inits signalSM value - which I think is not used at all - it exists only in not compiled code

	initFiles(); //defined in io.c - inits prep_g and trap_g files. I think these files are for tracing

//	WriteOptions( trap_g );

	c = seqcheck( seq ); //defined in mltaln9.c. It checks 'seq' characters for unusual characters
	if( c )
	{
		fprintf( stderr, "Illegal character %c\n", c );
		exit( 1 );
	}

//	writePre( njob, name, nlen, seq, 0 );


	for( i=0; i<njob; i++ ) 
	{
		gappick0( bseq[i], seq[i] ); //defined in mltaln9.c. copy seq[i] chars to bseq[i] without gaps
		thereisxineachseq[i] = removex( dseq[i], bseq[i] ); //defined here.
		//It fills dseq[i] with chars from bseq[i] without x char and assigns the count of 'x' to thereisxineachseq[i]
	}
	//so now 'seq' contains original sequences, 'bseq' contains sequences without gaps and 'dseq' contains sequences without gaps and x's

	//and this method aligns sequences based on algorithm type chosen
	pairalign( name, nlen, bseq, aseq, dseq, thereisxineachseq, mseq1, mseq2, alloclen, lastresx, distancemtx, localhomtable, ngui );

	fprintf( trap_g, "done.\n" ); //defined in defs.h.
#if DEBUG
	fprintf( stderr, "closing trap_g\n" );
#endif
	fclose( trap_g );
	fclose( prep_g ); //defined in defs.h.

//	writePre( njob, name, nlen, aseq, !contin );
#if 0
	writeData( stdout, njob, name, nlen, aseq );
#endif
#if IODEBUG
	fprintf( stderr, "OSHIMAI\n" );
#endif
	SHOWVERSION;

	if( stdout_dist && nthread > 1 )
	{
		fprintf( stderr, "\nThe order of distances is not identical to that in the input file, because of the parallel calculation.  Reorder them by yourself, using sort -n -k 2 | sort -n -k 1 -s\n" );
	}
	if( stdout_align && nthread > 1 )
	{
		fprintf( stderr, "\nThe order of pairwise alignments is not identical to that in the input file, because of the parallel calculation.  Reorder them by yourself.\n" );
	}

	//free all allocated memory
#if 1
	if( lastresx ) //free lastresx ** matrix
	{
		for( i=0; lastresx[i]; i++ ) 
		{
			for( j=0; lastresx[i][j].naln!=-1; j++ ) 
			{
				for( k=0; k<lastresx[i][j].naln; k++ )
				{
					free( lastresx[i][j].aln[k].reg1 );
					free( lastresx[i][j].aln[k].reg2 );
				}
				free( lastresx[i][j].aln );
			}
			free( lastresx[i] );
		}
		free( lastresx );
	}
#endif
	if( ngui == 0 ) 
	{
		FreeCharMtx( seq );
		FreeCharMtx( name );
	}
	FreeCharMtx( aseq );
	FreeCharMtx( bseq );
	FreeCharMtx( dseq );
	free( mseq1 );
	free( mseq2 );
	free( nlen );
	free( thereisxineachseq );
	freeconstants(); //free memory allocated in constants() call

	if( !ngui )
	{
		FreeCommonIP();
	}
	Falign( NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, 0, 0, 0, NULL, NULL, 0, NULL ); //defined in Falign.c. I see it doesn't do any thing, just returns 0
	G__align11( NULL, NULL, NULL, 0, 0, 0 ); // 20130603 //defined in Galign11.c. also only returns 0
	G__align11_noalign( NULL, 0, 0, NULL, NULL, 0 ); //defined in Galign11.c. also only returns 0
	L__align11( NULL, 0.0, NULL, NULL, 0, NULL, NULL ); //defined in Lalign11.c. also only returns 0
	L__align11_noalign( NULL, NULL, NULL ); //defined in Lalign11.c. also only returns 0
	genL__align11( NULL, NULL, NULL, 0, NULL, NULL ); //defined in genalign11.c. also only returns 0

#if SHISHAGONYU
	if( ngui )
	{
		char buf[100];
		for( i=0; i<njob-1; i++ ) for( j=i+1; j<njob; j++ )
		{
			sprintf( buf, "%5.3f", distancemtx[i][j-i] );
			distancemtx[i][j-i] = 0.0;
			sscanf( buf, "%lf", distancemtx[i]+j-i );
//			distancemtx[i][j-i] = 0.001 * (int)(distancemtx[i][j-i] * 1000 + 0.5);
		}

	}
#endif


	return( 0 );
}

