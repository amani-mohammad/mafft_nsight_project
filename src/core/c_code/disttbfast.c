//i think this files calculates distance between sequences
#include "mltaln.h"



#define DEBUG 0
#define IODEBUG 0
#define SCOREOUT 0
#define SKIP 1

#define END_OF_VEC -1

static int nadd;
static int treein;
static int topin;
static int treeout;
static int noalign;
static int distout;
static int tuplesize;
static int subalignment;
static int subalignmentoffset;
static int nguidetree;
static int sparsepickup;
static int keeplength;
static int ndeleted;
static int mapout;
static int smoothing;
static int compacttree = 0;
static double maxdistmtxsize;

#if 0
#define PLENFACA 0.0123
#define PLENFACB 10252
#define PLENFACC 10822
#define PLENFACD 0.5
#define DLENFACA 0.01
#define DLENFACB 2445
#define DLENFACC 2412
#define DLENFACD 0.1
#else
#define PLENFACA 0.01
#define PLENFACB 10000
#define PLENFACC 10000
#define PLENFACD 0.1
#define D6LENFACA 0.01
#define D6LENFACB 2500
#define D6LENFACC 2500
#define D6LENFACD 0.1
#define D10LENFACA 0.01
#define D10LENFACB 1000000
#define D10LENFACC 1000000
#define D10LENFACD 0.0
#endif

typedef struct _jobtable
{
    int i;  
    int j;  
} Jobtable;

typedef struct _msacompactdistmtxthread_arg
{
	int njob;
	int thread_no;
	int *selfscore;
	double **partmtx;
	char **seq;
	int **skiptable;
	double *mindist;
	int *mindistfrom;
 	int *jobpospt;
#ifdef enablemultithread
	pthread_mutex_t *mutex;
#endif
} msacompactdistmtxthread_arg_t;

typedef struct _compactdistmtxthread_arg
{
	int njob;
	int thread_no;
	int *nogaplen;
	int **pointt;
	int *selfscore;
	double **partmtx;
	int *jobpospt;
	double *mindist;
	int *mindistfrom;
#ifdef enablemultithread
	pthread_mutex_t *mutex;
#endif
} compactdistmtxthread_arg_t;

typedef struct _msadistmtxthread_arg
{
	int njob;
	int thread_no;
	int *selfscore;
	double **iscore;
	double **partmtx;
	char **seq;
	int **skiptable;
	Jobtable *jobpospt;
#ifdef enablemultithread
	pthread_mutex_t *mutex;
#endif
} msadistmtxthread_arg_t;

#ifdef enablemultithread
// ue futatsu ha singlethread demo tsukau
typedef struct _treebasethread_arg
{
	int thread_no;
	int njob;
	int *nrunpt;
	int *nlen;
	int *jobpospt;
	int ***topol;
	Treedep *dep;
	char **aseq;
	double *effarr;
	int *alloclenpt;
	int *fftlog;
	char *mergeoralign;
	double **newdistmtx;
	int *selfscore;
	pthread_mutex_t *mutex;
	pthread_cond_t *treecond;
} treebasethread_arg_t;

typedef struct _distancematrixthread_arg
{
	int thread_no;
	int njob;
	int *jobpospt;
	int **pointt;
	double **mtx;
	pthread_mutex_t *mutex;
} distancematrixthread_arg_t;
#endif


void dist_tb_fast_arguments( int argc, char *argv[] )
{
    int c;

	nthread = 1; //defined in defs.c.
	outnumber = 0; //defined in defs.h.
	topin = 0; //defined here.
	treein = 0; //defined here.
	treeout = 0; //defined here.
	distout = 0; //defined here.
	noalign = 0; //defined here.
	nevermemsave = 0; //defined in defs.h.
	inputfile = NULL; //defined in defs.h.
	nadd = 0; //defined here.
	addprofile = 1; //defined in defs.c.
	fftkeika = 0; //defined in defs.h.
	constraint = 0; //defined in defs.h.
	nblosum = 62; //defined in defs.h.
	fmodel = 0; //defined in defs.h.
	calledByXced = 0; //defined in defs.h.
	devide = 0; //defined in defs.h.
	use_fft = 0; //defined in defs.h.
	force_fft = 0; //defined in defs.h.
	fftscore = 1; //defined in defs.h.
	fftRepeatStop = 0; //defined in defs.h.
	fftNoAnchStop = 0; //defined in defs.h.
    weight = 3; //defined in defs.h.
    utree = 1; //defined in defs.h.
	tbutree = 1; //defined in defs.h.
    refine = 0; //defined in defs.h.
    check = 1; //defined in defs.h.
    cut = 0.0; //defined in defs.h.
    disp = 0; //defined in defs.h.
    outgap = 1; //defined in defs.c.
    alg = 'A'; //defined in defs.h.
    mix = 0; //defined in defs.h.
	tbitr = 0; //defined in defs.h.
	scmtd = 5; //defined in defs.h.
	tbweight = 0; //defined in defs.h.
	tbrweight = 3; //defined in defs.h.
	checkC = 0; //defined in defs.h.
	treemethod = 'X'; //defined in defs.h.
	sueff_global = 0.1; //defined in defs.c.
	contin = 0; //defined in defs.h.
	scoremtx = 1; //defined in defs.h.
	kobetsubunkatsu = 0; //defined in defs.h.
	dorp = NOTSPECIFIED; //defined in defs.c.
	ppenalty_dist = NOTSPECIFIED; //defined in defs.h.
	ppenalty = -1530; //defined in defs.h.
	ppenalty_ex = NOTSPECIFIED; //defined in defs.h.
	penalty_shift_factor = 1000.0; //defined in defs.c.
	poffset = -123; //defined in defs.h.
	kimuraR = NOTSPECIFIED; //defined in defs.h.
	pamN = NOTSPECIFIED; //defined in defs.h.
	geta2 = GETA2; //geta2 defined in defs.h and GETA2 defined in mltaln.h.
	fftWinSize = NOTSPECIFIED; //defined in defs.h.
	fftThreshold = NOTSPECIFIED; //defined in defs.h.
	TMorJTT = JTT; //TMorJTT defined in defs.h and JTT defined in mltaln.h.
	scoreout = 0; //defined in defs.c.
	spscoreout = 0; //defined in defs.h.
	tuplesize = 6; //defined here.
	subalignment = 0; //defined here.
	subalignmentoffset = 0; //defined here.
	legacygapcost = 0; //defined in defs.c.
	specificityconsideration = 0.0; //defined in defs.c.
	nguidetree = 1; //defined here.
	sparsepickup = 0; //defined here.
	keeplength = 0; //defined here.
	mapout = 0; //defined here.
	smoothing = 0; //defined here.
	nwildcard = 0; //defined in defs.c.

    while( --argc > 0 && (*++argv)[0] == '-' )
	{
        while ( (c = *++argv[0]) )
		{
            switch( c )
            {
				case 'i':
					inputfile = *++argv;
					reporterr(       "inputfile = %s\n", inputfile );
					--argc;
					goto nextoption;
				case 'I':
					nadd = myatoi( *++argv );
					reporterr(       "nadd = %d\n", nadd );
					--argc;
					goto nextoption;
				case 'V':
					ppenalty_dist = (int)( atof( *++argv ) * 1000 - 0.5 );
//					fprintf( stderr, "ppenalty = %d\n", ppenalty );
					--argc;
					goto nextoption;
				case 'f':
					ppenalty = (int)( atof( *++argv ) * 1000 - 0.5 );
//					reporterr(       "ppenalty = %d\n", ppenalty );
					--argc;
					goto nextoption;
				case 'Q':
					penalty_shift_factor = atof( *++argv );
					--argc;
					goto nextoption;
				case 'g':
					ppenalty_ex = (int)( atof( *++argv ) * 1000 - 0.5 );
					reporterr(       "ppenalty_ex = %d\n", ppenalty_ex );
					--argc;
					goto nextoption;
				case 'h':
					poffset = (int)( atof( *++argv ) * 1000 - 0.5 );
//					reporterr(       "poffset = %d\n", poffset );
					--argc;
					goto nextoption;
				case 'k':
					kimuraR = myatoi( *++argv );
					reporterr(       "kappa = %d\n", kimuraR );
					--argc;
					goto nextoption;
				case 'b':
					nblosum = myatoi( *++argv );
					scoremtx = 1;
//					reporterr(       "blosum %d / kimura 200 \n", nblosum );
					--argc;
					goto nextoption;
				case 'j':
					pamN = myatoi( *++argv );
					scoremtx = 0;
					TMorJTT = JTT;
					reporterr(       "jtt/kimura %d\n", pamN );
					--argc;
					goto nextoption;
				case 'm':
					pamN = myatoi( *++argv );
					scoremtx = 0;
					TMorJTT = TM;
					reporterr(       "tm %d\n", pamN );
					--argc;
					goto nextoption;
				case 'C':
					nthread = myatoi( *++argv );
					reporterr(       "nthread = %d\n", nthread );
					--argc; 
					goto nextoption;
				case 's':
					specificityconsideration = (double)myatof( *++argv );
//					reporterr(       "specificityconsideration = %f\n", specificityconsideration );
					--argc; 
					goto nextoption;
#if 1
				case 'a':
					fmodel = 1;
					break;
#endif
				case 'K':
					addprofile = 0;
					break;
				case 'y':
					distout = 1;
					break;
				case 't':
					treeout = 1;
					break;
				case 'T':
					noalign = 1;
					break;
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
				case 'L':
					legacygapcost = 1;
					break;
				case 'e':
					fftscore = 0;
					break;
				case 'H':
					subalignment = 1;
					subalignmentoffset = myatoi( *++argv );
					--argc;
					goto nextoption;
#if 0
				case 'R':
					fftRepeatStop = 1;
					break;
#endif
				case 'n' :
					outnumber = 1;
					break;
#if 0
				case 's':
					treemethod = 's';
					break;
				case 'q':
					treemethod = 'q'; // minimum
					break;
#endif
				case 'q':
					sparsepickup = myatoi( *++argv );
//					reporterr(       "sparsepickup = %d\n", sparsepickup );
					--argc; 
					goto nextoption;
				case 'X':
					treemethod = 'X';
					sueff_global = atof( *++argv );
//					fprintf( stderr, "sueff_global = %f\n", sueff_global );
					--argc;
					goto nextoption;
				case 'E':
					nguidetree = myatoi( *++argv );
//					reporterr(       "nguidetree = %d\n", nguidetree );
					--argc; 
					goto nextoption;
#if 0
				case 'a':
					alg = 'a';
					break;
				case 'H':
					alg = 'H';
					break;
				case 'R':
					alg = 'R';
					break;
#endif
				case 'A':
					alg = 'A';
					break;
				case '&':
					alg = 'a';
					break;
				case '@':
					alg = 'd';
					break;
				case 'N':
					nevermemsave = 1;
					break;
				case 'M':
					alg = 'M';
					break;
#if 0
				case 'S' :
					scoreout = 1; // for checking parallel calculation
					break;
#else
				case 'S' :
					spscoreout = 1; // 2014/Dec/30, sp score
					break;
#endif
				case 'B': // hitsuyou! memopt -M -B no tame
					break;
				case 'F':
					use_fft = 1;
					break;
				case 'G':
					use_fft = 1;
					force_fft = 1;
					break;
#if 0
				case 'V':
					topin = 1;
					break;
#endif
				case 'U':
					treein = 1;
					break;
				case 'u':
					weight = 0;
					tbrweight = 0;
					break;
				case 'v':
					tbrweight = 3;
					break;
#if 1
				case 'd':
					disp = 1;
					break;
#endif
#if 1
				case 'O':
					outgap = 0;
					break;
#else
				case 'O':
					fftNoAnchStop = 1;
					break;
#endif
				case 'J':
					tbutree = 0;
					break;
				case 'z':
					fftThreshold = myatoi( *++argv );
					--argc; 
					goto nextoption;
				case 'w':
					fftWinSize = myatoi( *++argv );
					--argc;
					goto nextoption;
				case 'W':
					tuplesize = myatoi( *++argv );
					--argc;
					goto nextoption;
#if 0
				case 'Z':
					checkC = 1;
					break;
#endif
				case 'Y':
					keeplength = 1;
					break;
				case 'Z':
					mapout = 1;
					break;
				case 'p':
					smoothing = 1;
					break;
				case ':':
					nwildcard = 1;
					break;
                default:
                    reporterr(       "illegal option %c\n", c );
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
        reporterr(       "options: Check source file !\n" );
        exit( 1 );
    }
	if( tbitr == 1 && outgap == 0 )
	{
		reporterr(       "conflicting options : o, m or u\n" );
		exit( 1 );
	}
}

//align sequences from seq using G__align11_noalign based on npick and other args
static int varpairscore( int nseq, int npick, int nlenmax, char **seq, int seed )
{
	int i, j, npair;
	int *slist;
	char **pickseq;
	double score;
	double scoreav;
	double scoreav2;
	double scorestd;
	double scorevar;
	slist = calloc( nseq, sizeof( int ) );
	pickseq = AllocateCharMtx( npick, nlenmax );
	reporterr( "nseq = %d, nlenmax=%d, seed=%d\n", nseq, nlenmax, seed );

	srand( seed );

	for( i=0; i<nseq; i++ ) slist[i] = i;
//	for( i=0; i<nseq; i++ ) reporterr( "slist[%d] = %d\n", i, slist[i] );

	stringshuffle( slist, nseq ); //defined in mltaln9.c. shuffle slist items
	for( i=0; i<npick; i++ ) gappick0( pickseq[i], seq[slist[i]] ); //defined in mltaln9.c. copy 'seq[slist[i]]' chars to 'pickseq[i]' without gaps chars

	scoreav = 0.0;
	scoreav2 = 0.0;
	npair = npick * (npick-1) / 2;
	for( i=1; i<npick; i++ ) 
	{
		reporterr( "%d / %d\r", i, npick );
		for( j=0; j<i; j++ ) 
		{
			//it is like G__align11 but without warp and some other details
			score = G__align11_noalign( n_dis_consweight_multi, -1200, -60, pickseq+i, pickseq+j, nlenmax ); //defined in Galign11.c.
			scoreav += score;
			scoreav2 += score * score;
			printf( "score = %d\n", (int)score );
		}
	}

	scoreav /= (double)npair;
	scoreav2 /= (double)npair;
	scorevar = ( scoreav2 - scoreav * scoreav )*npair/(npair-1);
	scorestd = sqrt( scorevar );
	printf( "av = %f\n", scoreav );
	printf( "stddev = %f\n", scorestd );
	printf( "cv = %f\n", scorestd/scoreav );

	FreeCharMtx( pickseq );

	if( scorestd/scoreav < 0.2 ) return( 's' );
	else return( 't' );
}

//pickup part of sequences based on sparsepickup value, and write selected sequences and their names to stdout, and write not used ones to notused file
static void pickup( int n, int *seqlen, int ***topol, char **name, char **seq ) // memsave ni mitaiou
{
	int i, j, k, m;
	int **longestseq;
	int **longestlen;
	int *select;
	char **nameout, **seqout;
	int *nlenout;
	char **namenotused, **seqnotused;
	int *nlennotused;
	FILE *notusedfp;

	longestseq = AllocateIntMtx( n-1, 2 );
	longestlen = AllocateIntMtx( n-1, 2 );
	select = AllocateIntVec( n );
	for( i=0; i<n; i++ ) select[i] = 0;
	nameout = AllocateCharMtx( n, 0 );
	seqout = AllocateCharMtx( n, 0 );
	nlenout = AllocateIntVec( n );
	namenotused = AllocateCharMtx( n, 0 );
	seqnotused = AllocateCharMtx( n, 0 );
	nlennotused = AllocateIntVec( n );

	for( i=0; i<n-1; i++ )
	{
//		reporterr( "STEP %d\n", i );
		longestlen[i][0] = -1;
		longestseq[i][0] = -1;
		for( j=0; (m=topol[i][0][j])!=-1; j++ ) // sukoshi muda
		{
			if( seqlen[m] > longestlen[i][0] )
			{
				longestlen[i][0] = seqlen[m];
				longestseq[i][0] = m;
			}
//			reporterr( "%d ", topol[i][0][j] );
		}
//		reporterr( "longest = %d (%d)\n", longestlen[i][0], longestseq[i][0] );


		longestlen[i][1] = -1;
		longestseq[i][1] = -1;
		for( j=0; (m=topol[i][1][j])!=-1; j++ ) // sukoshi muda
		{
			if( seqlen[m] > longestlen[i][1] )
			{
				longestlen[i][1] = seqlen[m];
				longestseq[i][1] = m;
			}
//			reporterr( "%d ", topol[i][1][j] );
		}
//		reporterr( "longest = %d (%d)\n", longestlen[i][1], longestseq[i][1] );
	}

	m = 1;
	for( i=n-2; i>-1; i-- )
	{
//		reporterr( "longest[%d][0] = %d (%d)\n", i, longestlen[i][0], longestseq[i][0] );
//		reporterr( "longest[%d][1] = %d (%d)\n", i, longestlen[i][1], longestseq[i][1] );
		select[longestseq[i][0]] = 1;
		select[longestseq[i][1]] = 1;
		m += 1;
		if( m >= sparsepickup ) break;
	}
	for( i=0, k=0, j=0; i<n; i++ ) 
	{
		if( select[i] )
		{
			nameout[k] = name[i];
			seqout[k] = seq[i];
			nlenout[k] = strlen( seqout[k] );
			k++;
		}
		else
		{
			namenotused[j] = name[i];
			seqnotused[j] = seq[i];
			nlennotused[j] = strlen( seqnotused[j] );
			j++;
		}
	}
	writeData_pointer( stdout, m, nameout, nlenout, seqout ); //defined in io.c. write selected sequences and their names to stdout.

	notusedfp = fopen( "notused", "w" ); //open 'notused' file for writing
	writeData_pointer( notusedfp, n-m, namenotused, nlennotused, seqnotused ); //write not used sequences and their names to 'notused' file.
	fclose( notusedfp );


	free( nameout );
	free( nlenout );
	free( seqout );
	free( namenotused );
	free( nlennotused );
	free( seqnotused );
	FreeIntMtx( longestseq );
	FreeIntMtx( longestlen );
	free( select );
}


static int nunknown = 0;

//fill grp with values from 'amino_grp' based on seq chars
void disttbfast_seq_grp_nuc( int *grp, char *seq )
{
	int tmp;
	int *grpbk = grp;
	while( *seq )
	{
		tmp = amino_grp[(int)*seq++];
		if( tmp < 4 )
			*grp++ = tmp;
		else
			nunknown++;
	}
	*grp = END_OF_VEC;
	if( grp - grpbk < tuplesize )
	{
//		reporterr(       "\n\nWARNING: Too short.\nPlease also consider use mafft-ginsi, mafft-linsi or mafft-ginsi.\n\n\n" );
//		exit( 1 );
		*grpbk = -1;
	}
}

//fill grp with values from 'amino_grp' based on seq chars
void disttbfast_seq_grp( int *grp, char *seq )
{
	int tmp;
	int *grpbk = grp;
	while( *seq )
	{
		tmp = amino_grp[(unsigned char)*seq++];
		if( tmp < 6 )
			*grp++ = tmp;
		else
			nunknown++;
	}
	*grp = END_OF_VEC;
	if( grp - grpbk < 6 )
	{
//		reporterr(       "\n\nWARNING: Too short.\nPlease also consider use mafft-ginsi, mafft-linsi or mafft-ginsi.\n\n\n" );
//		exit( 1 );
		*grpbk = -1;
	}
}

//update table values based on pointt values
void disttbfast_makecompositiontable_p( int *table, int *pointt )
{
	int point;

	while( ( point = *pointt++ ) != END_OF_VEC )
		table[point]++;
}

//fill pointt with values based on n and some calculations on it
void disttbfast_makepointtable_nuc_dectet( int *pointt, int *n )
{
	int point;
	register int *p;

	if( *n == -1 )
	{
		*pointt = -1;
		return;
	}

	p = n;
	point  = *n++ *262144;
	point += *n++ * 65536;
	point += *n++ * 16384;
	point += *n++ *  4096;
	point += *n++ *  1024;
	point += *n++ *   256;
	point += *n++ *    64;
	point += *n++ *    16;
	point += *n++ *     4;
	point += *n++;
	*pointt++ = point;

	while( *n != END_OF_VEC )
	{
		point -= *p++ *262144;
		point *= 4;
		point += *n++;
		*pointt++ = point;

	}
	*pointt = END_OF_VEC;
}

void disttbfast_makepointtable_nuc_octet( int *pointt, int *n )
{
	int point;
	register int *p;

	if( *n == -1 )
	{
		*pointt = -1;
		return;
	}

	p = n;
	point  = *n++ * 16384;
	point += *n++ *  4096;
	point += *n++ *  1024;
	point += *n++ *   256;
	point += *n++ *    64;
	point += *n++ *    16;
	point += *n++ *     4;
	point += *n++;
	*pointt++ = point;

	while( *n != END_OF_VEC )
	{
		point -= *p++ * 16384;
		point *= 4;
		point += *n++;
		*pointt++ = point;
	}
	*pointt = END_OF_VEC;
}

//fill pointt with values based on n and some calculations on it
void disttbfast_makepointtable_nuc( int *pointt, int *n )
{
	int point;
	register int *p;

	if( *n == -1 )
	{
		*pointt = -1;
		return;
	}

	p = n;
	point  = *n++ *  1024;
	point += *n++ *   256;
	point += *n++ *    64;
	point += *n++ *    16;
	point += *n++ *     4;
	point += *n++;
	*pointt++ = point;

	while( *n != END_OF_VEC )
	{
		point -= *p++ * 1024;
		point *= 4;
		point += *n++;
		*pointt++ = point;
	}
	*pointt = END_OF_VEC;
}

//fill pointt with values based on n and some calculations on it
void disttbfast_makepointtable( int *pointt, int *n )
{
	int point;
	register int *p;

	if( *n == -1 )
	{
		*pointt = -1;
		return;
	}

	p = n;
	point  = *n++ *  7776;
	point += *n++ *  1296;
	point += *n++ *   216;
	point += *n++ *    36;
	point += *n++ *     6;
	point += *n++;
	*pointt++ = point;

	while( *n != END_OF_VEC )
	{
		point -= *p++ * 7776;
		point *= 6;
		point += *n++;
		*pointt++ = point;
	}
	*pointt = END_OF_VEC;
}

//return value based on pos and ori subtraction
static double preferenceval( int ori, int pos, int max ) // for debug
{
	pos -= ori;
	if( pos < 0 ) pos += max;
	return( 0.00000000000001 * pos );
}

//update arg->partmtx, arg->jobpospt, arg->mindist, arg->mindistfrom values based on calculations on other args.
static void *compactdisthalfmtxthread( void *arg ) // enablemultithread == 0 demo tsukau
{
	compactdistmtxthread_arg_t *targ = (compactdistmtxthread_arg_t *)arg;
	int njob = targ->njob;
	int thread_no = targ->thread_no;
	int *selfscore = targ->selfscore;
	double **partmtx = targ->partmtx;
	int *nogaplen = targ->nogaplen;
	int **pointt = targ->pointt;
 	int *jobpospt = targ->jobpospt;
	double *mindist = targ->mindist;
	int *mindistfrom = targ->mindistfrom;
	int i, j;
	double tmpdist, preference, tmpdistx, tmpdisty;
	int *table1;

	while( 1 )
	{
#ifdef enablemultithread
		if( nthread ) 
		{
			pthread_mutex_lock( targ->mutex );
			i = *jobpospt;
			if( i == njob-1 )
			{
				pthread_mutex_unlock( targ->mutex );
				commonsextet_p( NULL, NULL );
				return( NULL );
			}
			*jobpospt = i+1;
			pthread_mutex_unlock( targ->mutex );
		}
		else
#endif 
		{
			i = *jobpospt;
			if( i == njob-1 )
			{
				commonsextet_p( NULL, NULL ); //defined in mltaln9.c. free allocated memory from previous call.
				return( NULL );
			}
			*jobpospt = i+1;
		}

		table1 = (int *)calloc( tsize, sizeof( int ) );
		if( !table1 ) ErrorExit( "Cannot allocate table1\n" );
		if( i % 100 == 0 )
		{
			if( nthread )
				reporterr(       "\r% 5d / %d (thread %4d)", i+1, njob, thread_no );
			else
				reporterr(       "\r% 5d / %d", i+1, njob );
		}
		disttbfast_makecompositiontable_p( table1, pointt[i] ); //defined here. update table1 values based on pointt[i] values

		for( j=i+1; j<njob; j++ ) 
		{

			tmpdist = distcompact( nogaplen[i], nogaplen[j], table1, pointt[j], selfscore[i], selfscore[j] ); //defined in mltaln9.c.
			//return value based on 'commonsextet_p' call for table1 and pointt[j] and calculations on this value and other args
			preference = preferenceval( i, j, njob ); //defined here. return value based on i and j subtraction
			tmpdistx = tmpdist + preference;
			if( tmpdistx < mindist[i] )
			{
				mindist[i] = tmpdistx;
				mindistfrom[i] = j;
			}

			preference = preferenceval( j, i, njob ); //defined here. return value based on j and i subtraction
			tmpdisty = tmpdist + preference;
			if( tmpdisty < mindist[j] )
			{
				mindist[j] = tmpdisty;
				mindistfrom[j] = i;
			}

			if( partmtx[i] ) partmtx[i][j] = tmpdist;
			if( partmtx[j] ) partmtx[j][i] = tmpdist;
		} 
		free( table1 );
	}
}


static void *msacompactdisthalfmtxthread( void *arg ) // enablemultithread == 0 demo tsukau
{
	msacompactdistmtxthread_arg_t *targ = (msacompactdistmtxthread_arg_t *)arg;
	int njob = targ->njob;
	int thread_no = targ->thread_no;
	int *selfscore = targ->selfscore;
	double **partmtx = targ->partmtx;
	char **seq = targ->seq;
	int **skiptable = targ->skiptable;
	double *mindist = targ->mindist;
	int *mindistfrom = targ->mindistfrom;
 	int *jobpospt = targ->jobpospt;
	double tmpdist, preference, tmpdistx, tmpdisty;
	int i, j;

	while( 1 )
	{
#ifdef enablemultithread
		if( nthread ) 
		{
			pthread_mutex_lock( targ->mutex );
			i = *jobpospt;
			if( i == njob-1 )
			{
				pthread_mutex_unlock( targ->mutex );
				return( NULL );
			}
			*jobpospt = i+1;
			pthread_mutex_unlock( targ->mutex );
		}
		else
#endif
		{
			i = *jobpospt;
			if( i == njob-1 )
			{
				return( NULL );
			}
			*jobpospt = i+1;
		}

		if( i % 100 == 0 ) 
		{
			if( nthread )
				fprintf( stderr, "\r% 5d / %d (thread %4d)", i, njob, thread_no );
			else
				fprintf( stderr, "\r% 5d / %d", i, njob );
		}

		for( j=i+1; j<njob; j++ ) 
		{
			tmpdist = distcompact_msa( seq[i], seq[j], skiptable[i], skiptable[j], selfscore[i], selfscore[j] ); // osoikedo,

			preference = preferenceval( i, j, njob );
			tmpdistx = tmpdist + preference;
			if( tmpdistx < mindist[i] )
			{
				mindist[i] = tmpdistx;
				mindistfrom[i] = j;
			}

			preference = preferenceval( j, i, njob );
			tmpdisty = tmpdist + preference;
			if( tmpdisty < mindist[j] )
			{
				mindist[j] = tmpdisty;
				mindistfrom[j] = i;
			}
			if( partmtx[i] ) partmtx[i][j] = tmpdist;
			if( partmtx[j] ) partmtx[j][i] = tmpdist;
		}
	}
}

#if 1
static void *msadistmtxthread( void *arg ) // enablemultithread == 0 demo tsukau
{
	msadistmtxthread_arg_t *targ = (msadistmtxthread_arg_t *)arg;
	int njob = targ->njob;
	int thread_no = targ->thread_no;
	int *selfscore = targ->selfscore;
	double **iscore = targ->iscore;
	char **seq = targ->seq;
	int **skiptable = targ->skiptable;
	Jobtable *jobpospt = targ->jobpospt;


	double ssi, ssj, bunbo, iscoretmp;
	int i, j;
	int nlim = njob-1;

	while( 1 )
	{
#ifdef enablemultithread
		if( nthread ) 
		{
			pthread_mutex_lock( targ->mutex );
			i = jobpospt->i; // (jobpospt-i)++ dato, shuuryou hantei no mae ni ++ surunode, tomaranakunaru.

			if( i == nlim )
			{
				pthread_mutex_unlock( targ->mutex );
				return( NULL );
			}
			jobpospt->i += 1;
			pthread_mutex_unlock( targ->mutex );
			if( i % 100 == 0 ) fprintf( stderr, "\r% 5d / %d (thread %4d)", i, njob, thread_no );
		}
		else
#endif
		{
			i = (jobpospt->i)++;
			if( i == nlim ) return( NULL );
			if( i % 100 == 0 ) fprintf( stderr, "\r% 5d / %d", i, njob );
		}

		ssi = selfscore[i];
		for( j=i+1; j<njob; j++ )
		{
			ssj = selfscore[j];
			bunbo = MIN( ssi, ssj );
//fprintf( stderr, "bunbo = %f\n", bunbo );
//fprintf( stderr, "naivepairscorefast() = %f\n", naivepairscorefast( seq[i], seq[j], skiptable[i], skiptable[j], penalty_dist ) );
			if( bunbo == 0.0 )
				iscoretmp = 2.0; // 2013/Oct/17
			else
			{
				iscoretmp = ( 1.0 - naivepairscorefast( seq[i], seq[j], skiptable[i], skiptable[j], penalty_dist ) / bunbo ) * 2.0; // 2014/Aug/15 fast 
				if( iscoretmp > 10 ) iscoretmp = 10.0;  // 2015/Mar/17
	
			}
			if( iscoretmp < 0.0 ) 
			{
				reporterr( "WARNING: negative distance, iscoretmp = %f\n", iscoretmp );
				iscoretmp = 0.0;
			}
			iscore[i][j-i] = iscoretmp;
//			printf( "i,j=%d,%d, iscoretmp=%f\n", i, j, iscoretmp );

		}
	}
}
#else
static void *msadistmtxthread( void *arg ) // enablemultithread == 0 demo tsukau
{
	msadistmtxthread_arg_t *targ = (msadistmtxthread_arg_t *)arg;
	int njob = targ->njob;
	int thread_no = targ->thread_no;
	int *selfscore = targ->selfscore;
	double **iscore = targ->iscore;
	char **seq = targ->seq;
	int **skiptable = targ->skiptable;
	Jobtable *jobpospt = targ->jobpospt;


	double ssi, ssj, bunbo, iscoretmp;
	int i, j;

	while( 1 )
	{
#ifdef enablemultithread
		if( nthread ) pthread_mutex_lock( targ->mutex );
#endif
		j = jobpospt->j;
		i = jobpospt->i;
		j++;
		if( j == njob )
		{
			i++;
			j = i + 1;
			if( i == njob-1 )
			{
#ifdef enablemultithread
				if( nthread ) pthread_mutex_unlock( targ->mutex );
#endif
				return( NULL );
			}
		}
		jobpospt->j = j;
		jobpospt->i = i;
#ifdef enablemultithread
		if( nthread ) pthread_mutex_unlock( targ->mutex );
#endif


		if( nthread )
		{
			if( j==i+1 && i % 10 == 0 ) fprintf( stderr, "\r% 5d / %d (thread %4d)", i, njob, thread_no );
		}
		else
		{
			if( j==i+1 && i % 10 == 0 ) fprintf( stderr, "\r% 5d / %d", i, njob );
		}
		ssi = selfscore[i];
		ssj = selfscore[j];
		bunbo = MIN( ssi, ssj );
//fprintf( stderr, "bunbo = %f\n", bunbo );
//fprintf( stderr, "naivepairscorefast() = %f\n", naivepairscorefast( seq[i], seq[j], skiptable[i], skiptable[j], penalty_dist ) );
		if( bunbo == 0.0 )
			iscoretmp = 2.0; // 2013/Oct/17
		else
		{
			iscoretmp = ( 1.0 - naivepairscorefast( seq[i], seq[j], skiptable[i], skiptable[j], penalty_dist ) / bunbo ) * 2.0; // 2014/Aug/15 fast 
			if( iscoretmp > 10 ) iscoretmp = 10.0;  // 2015/Mar/17

		}
		iscore[i][j-i] = iscoretmp;


	}
}
#endif

#ifdef enablemultithread
static void *distancematrixthread( void *arg )
{
	distancematrixthread_arg_t *targ = (distancematrixthread_arg_t *)arg;
	int thread_no = targ->thread_no;
	int njob = targ->njob;
	int *jobpospt = targ->jobpospt;
	int **pointt = targ->pointt;
	double **mtx = targ->mtx;

	int *table1;
	int i, j;

	while( 1 )
	{
		pthread_mutex_lock( targ->mutex );
		i = *jobpospt;
		if( i == njob )
		{
			pthread_mutex_unlock( targ->mutex );
			commonsextet_p( NULL, NULL );
			return( NULL );
		}
		*jobpospt = i+1;
		pthread_mutex_unlock( targ->mutex );

		table1 = (int *)calloc( tsize, sizeof( int ) );
		if( !table1 ) ErrorExit( "Cannot allocate table1\n" );
		if( i % 100 == 0 )
		{
			reporterr(       "\r% 5d / %d (thread %4d)", i+1, njob, thread_no );
		}
		disttbfast_makecompositiontable_p( table1, pointt[i] );

		for( j=i; j<njob; j++ ) 
		{
			mtx[i][j-i] = (double)commonsextet_p( table1, pointt[j] );
		} 
		free( table1 );
	}
}


static void *treebasethread( void *arg )
{
	treebasethread_arg_t *targ = (treebasethread_arg_t *)arg;
	int thread_no = targ->thread_no;
	int *nrunpt = targ->nrunpt;
	int njob = targ->njob;
	int *nlen = targ->nlen;
	int *jobpospt = targ->jobpospt;
	int ***topol = targ->topol;
	Treedep *dep = targ->dep;
	char **aseq = targ->aseq;
	double *effarr = targ->effarr;
	int *alloclen = targ->alloclenpt;
	int *fftlog = targ->fftlog;
	char *mergeoralign = targ->mergeoralign;
	double **newdistmtx = targ->newdistmtx;
	int *selfscore = targ->selfscore;

	char **mseq1, **mseq2;
	char **localcopy;
	int i, m, j, l;
	int immin, immax;
	int len1, len2;
	int clus1, clus2;
	double pscore, tscore;
	char *indication1, *indication2;
	double *effarr1 = NULL;
	double *effarr2 = NULL;
//	double dumfl = 0.0;
	double dumdb = 0.0;
	int ffttry;
	int m1, m2;
	double **dynamicmtx;
	int ssi, ssm, bunbo;
	int tm, ti;
	int **localmem = NULL;
	int posinmem;
#if SKIP
	int **skiptable1 = NULL, **skiptable2 = NULL;
#endif
#if 0
	int i, j;
#endif

	mseq1 = AllocateCharMtx( njob, 0 );
	mseq2 = AllocateCharMtx( njob, 0 );
	localcopy = calloc( njob, sizeof( char * ) );
	for( i=0; i<njob; i++ ) localcopy[i] = NULL;
	effarr1 = AllocateDoubleVec( njob );
	effarr2 = AllocateDoubleVec( njob );
	indication1 = AllocateCharVec( 150 );
	indication2 = AllocateCharVec( 150 );
	dynamicmtx = AllocateDoubleMtx( nalphabets, nalphabets );
	localmem = AllocateIntMtx( 2, njob+1 );


#if 0
	reporterr(       "##### fftwinsize = %d, fftthreshold = %d\n", fftWinSize, fftThreshold );
#endif

#if 0
	for( i=0; i<njob; i++ )
		reporterr(       "TBFAST effarr[%d] = %f\n", i, effarr[i] );
#endif

//	for( i=0; i<njob; i++ ) strcpy( aseq[i], seq[i] );


//	writePre( njob, name, nlen, aseq, 0 );

	while( 1 )
	{
		pthread_mutex_lock( targ->mutex );
		l = *jobpospt;
		if( l == njob-1 )
		{
			pthread_mutex_unlock( targ->mutex );
			if( commonIP ) FreeIntMtx( commonIP );
			commonIP = NULL;
			Falign( NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, 0, 0, 0, NULL, NULL, 0, NULL );
			Falign_udpari_long( NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, 0, 0, 0, NULL );
			A__align( NULL,  NULL, NULL, NULL, NULL, 0, 0, 0, NULL, NULL, NULL, NULL, NULL, NULL, NULL, 0, NULL, 0, 0, -1, -1 );
			D__align( NULL,  NULL, NULL, NULL, NULL, 0, 0, 0, NULL, NULL, NULL, NULL, NULL, NULL, NULL, 0, NULL, 0, 0 );
			G__align11( NULL, NULL, NULL, 0, 0, 0 ); // iru?
			free( mseq1 );
			free( mseq2 );
			free( localcopy );
			free( effarr1 );
			free( effarr2 );
			free( indication1 );
			free( indication2 );
			FreeDoubleMtx( dynamicmtx );
			FreeIntMtx( localmem );
			return( NULL );
		}
		*jobpospt = l+1;

		if( dep[l].child0 != -1 )
		{
			while( dep[dep[l].child0].done == 0 )
				pthread_cond_wait( targ->treecond, targ->mutex );
		}
		if( dep[l].child1 != -1 ) 
		{
			while( dep[dep[l].child1].done == 0 )
				pthread_cond_wait( targ->treecond, targ->mutex );
		}
		while( *nrunpt >= nthread )
			pthread_cond_wait( targ->treecond, targ->mutex );
		(*nrunpt)++;


		if( mergeoralign[l] == 'n' )
		{
//			reporterr(       "SKIP!\n" );
			dep[l].done = 1;
			(*nrunpt)--;
			pthread_cond_broadcast( targ->treecond );
//			free( topol[l][0] ); topol[l][0] = NULL;
//			free( topol[l][1] ); topol[l][1] = NULL;
//			free( topol[l] ); topol[l] = NULL;
			pthread_mutex_unlock( targ->mutex );
			continue;
		}



		m1 = topol[l][0][0];
		m2 = topol[l][1][0];

//		reporterr(       "\ndistfromtip = %f\n", dep[l].distfromtip );
//		makedynamicmtx( dynamicmtx, n_dis_consweight_multi, dep[l].distfromtip - 0.5 );
		makedynamicmtx( dynamicmtx, n_dis_consweight_multi, dep[l].distfromtip );

		len1 = strlen( aseq[m1] );
		len2 = strlen( aseq[m2] );
		if( *alloclen <= len1 + len2 )
		{
			reporterr(       "\nReallocating.." );
			*alloclen = ( len1 + len2 ) + 1000;
			ReallocateCharMtx( aseq, njob, *alloclen + 10  ); 
			reporterr(       "done. *alloclen = %d\n", *alloclen );
		}

		localmem[0][0] = -1;
		posinmem=0;
		topolorder( njob, localmem[0], &posinmem, topol, dep, l, 0 );
		localmem[1][0] = -1;
		posinmem=0;
		topolorder( njob, localmem[1], &posinmem, topol, dep, l, 1 );
		for( i=0; (j=localmem[0][i])!=-1; i++ )
		{
			localcopy[j] = calloc( *alloclen, sizeof( char ) );
			strcpy( localcopy[j], aseq[j] );
		}
		for( i=0; (j=localmem[1][i])!=-1; i++ )
		{
			localcopy[j] = calloc( *alloclen, sizeof( char ) );
			strcpy( localcopy[j], aseq[j] );
		}

		if( !nevermemsave && ( alg != 'M' && ( len1 > 30000 || len2 > 30000  ) ) )
		{
			reporterr(       "\nlen1=%d, len2=%d, Switching to the memsave mode\n", len1, len2 );
			alg = 'M';
		}

		if( alg == 'M' ) // hoka no thread ga M ni shitakamo shirenainode
		{
//			reporterr(       "Freeing commonIP (thread %d)\n", thread_no );
			if( commonIP ) FreeIntMtx( commonIP );
			commonIP = NULL;
			commonAlloc1 = 0;
			commonAlloc2 = 0;
		}

		pthread_mutex_unlock( targ->mutex );

#if 1 // CHUUI@@@@
		clus1 = fastconjuction_noname( localmem[0], localcopy, mseq1, effarr1, effarr, indication1, 0.0 );
		clus2 = fastconjuction_noname( localmem[1], localcopy, mseq2, effarr2, effarr, indication2, 0.0 );
#else
		clus1 = fastconjuction_noweight( topol[l][0], localcopy, mseq1, effarr1,  indication1 );
		clus2 = fastconjuction_noweight( topol[l][1], localcopy, mseq2, effarr2,  indication2 );
#endif


#if 0
    for( i=0; i<clus1; i++ ) 
    {
        if( strlen( mseq1[i] ) != len1 ) 
        {
            reporterr(       "i = %d / %d\n", i, clus1 ); 
            reporterr(       "hairetsu ga kowareta (in treebase, after conjuction) !\n" ); 
            exit( 1 );
        }
    }
    for( j=0; j<clus2; j++ )
    {
        if( strlen( mseq2[j] ) != len2 ) 
        {
            reporterr(       "j = %d / %d\n", j, clus2 ); 
            reporterr(       "hairetsu ga kowareta (in treebase, after conjuction) !\n" ); 
            exit( 1 );
        }
    }
#endif

#if 0
    for( i=0; i<clus1; i++ ) 
    {
        if( strlen( mseq1[i] ) != len1 ) 
        {
            reporterr(       "i = %d / %d\n", i, clus1 ); 
            reporterr(       "hairetsu ga kowareta (in treebase, after free topol) !\n" ); 
            exit( 1 );
        }
    }
    for( j=0; j<clus2; j++ )
    {
        if( strlen( mseq2[j] ) != len2 ) 
        {
            reporterr(       "j = %d / %d\n", j, clus2 ); 
            reporterr(       "hairetsu ga kowareta (in treebase, after free topol) !\n" ); 
            exit( 1 );
        }
    }
#endif


//		fprintf( trap_g, "\nSTEP-%d\n", l );
//		fprintf( trap_g, "group1 = %s\n", indication1 );
//		fprintf( trap_g, "group2 = %s\n", indication2 );

//		reporterr(       "\rSTEP % 5d / %d %d-%d", l+1, njob-1, clus1, clus2 );
		reporterr(       "\rSTEP % 5d / %d (thread %4d)", l+1, njob-1, thread_no );

#if 0
		reporterr(       "STEP %d /%d\n", l+1, njob-1 );
		reporterr(       "group1 = %.66s", indication1 );
		if( strlen( indication1 ) > 66 ) reporterr(       "..." );
		reporterr(       "\n" );
		reporterr(       "group2 = %.66s", indication2 );
		if( strlen( indication2 ) > 66 ) reporterr(       "..." );
		reporterr(       "\n" );
#endif

/*
		reporterr(       "before align all\n" );
		display( aseq, njob );
		reporterr(       "\n" );
		reporterr(       "before align 1 %s \n", indication1 );
		display( mseq1, clus1 );
		reporterr(       "\n" );
		reporterr(       "before align 2 %s \n", indication2 );
		display( mseq2, clus2 );
		reporterr(       "\n" );
*/



//		if( fftlog[m1] && fftlog[m2] ) ffttry = ( nlen[m1] > clus1 && nlen[m2] > clus2 );
		if( fftlog[m1] && fftlog[m2] ) ffttry = ( nlen[m1] > clus1 && nlen[m2] > clus2 && clus1 < 1000 && clus2 < 1000);
		else						   ffttry = 0;
//		ffttry = ( nlen[m1] > clus1 && nlen[m2] > clus2 && clus1 < 5000 && clus2 < 5000); // v6.708
//		reporterr(       "f=%d, len1/fftlog[m1]=%f, clus1=%d, len2/fftlog[m2]=%f, clus2=%d\n", ffttry, (double)len1/fftlog[m1], clus1, (double)len2/fftlog[m2], clus2 );

		if( force_fft || ( use_fft && ffttry ) )
		{
			reporterr(       "f" );
			if( alg == 'M' )
			{
				reporterr(       "m" );
				pscore = Falign_udpari_long( NULL, NULL, dynamicmtx, mseq1, mseq2, effarr1, effarr2, NULL, NULL, clus1, clus2, *alloclen, fftlog+m1 );
			}
			else
			{
				pscore = Falign( NULL, NULL, dynamicmtx, mseq1, mseq2, effarr1, effarr2, NULL, NULL, clus1, clus2, *alloclen, fftlog+m1, NULL, 0, NULL );
			}
		}
		else
		{
			reporterr(       "d" );
			fftlog[m1] = 0;
			switch( alg )
			{
				case( 'a' ):
					pscore = Aalign( mseq1, mseq2, effarr1, effarr2, clus1, clus2, *alloclen );
					break;
				case( 'M' ):
					reporterr(       "m" );
//					reporterr(       "%d-%d", clus1, clus2 );
					pscore = MSalignmm( dynamicmtx, mseq1, mseq2, effarr1, effarr2, clus1, clus2, *alloclen, NULL, NULL, NULL, NULL, NULL, 0, NULL, outgap, outgap );
					break;
				case( 'd' ):
					if( 1 && clus1 == 1 && clus2 == 1 )
					{
//						reporterr(       "%d-%d", clus1, clus2 );
						pscore = G__align11( dynamicmtx, mseq1, mseq2, *alloclen, outgap, outgap );
					}
					else
					{
//						reporterr(       "%d-%d", clus1, clus2 );
						pscore = D__align_ls( dynamicmtx, mseq1, mseq2, effarr1, effarr2, clus1, clus2, *alloclen, NULL, &dumdb, NULL, NULL, NULL, NULL, NULL, 0, NULL, outgap, outgap );
					}
					break;
				case( 'A' ):
					if( clus1 == 1 && clus2 == 1 )
					{
//						reporterr(       "%d-%d", clus1, clus2 );
						pscore = G__align11( dynamicmtx, mseq1, mseq2, *alloclen, outgap, outgap );
					}
					else
					{
//						reporterr(       "%d-%d", clus1, clus2 );
						pscore = A__align( dynamicmtx, mseq1, mseq2, effarr1, effarr2, clus1, clus2, *alloclen, NULL, &dumdb, NULL, NULL, NULL, NULL, NULL, 0, NULL, outgap, outgap, -1, -1 );
					}
					break;
				default:
					ErrorExit( "ERROR IN SOURCE FILE" );
			}
		}
#if SCOREOUT
		reporterr(       "score = %10.2f\n", pscore );
#endif
		tscore += pscore;
		nlen[m1] = 0.5 * ( nlen[m1] + nlen[m2] );

		if( disp ) display( localcopy, njob );

		if( newdistmtx ) // tsukawanai
		{
#if 0
			reporterr( "group1 = " );
			for( i=0; i<clus1; i++ ) reporterr( "%d ", topol[l][0][i] );
			reporterr( "\n" );
			reporterr( "group2 = " );
			for( m=0; m<clus2; m++ ) reporterr( "%d ", topol[l][1][m] );
			reporterr( "\n" );
#endif
#if SKIP
			skiptable1 = AllocateIntMtx( clus1, 0 );
			skiptable2 = AllocateIntMtx( clus2, 0 );
			makeskiptable( clus1, skiptable1, mseq1 ); // allocate suru.
			makeskiptable( clus2, skiptable2, mseq2 ); // allocate suru.
#endif
			for( i=0; i<clus1; i++ ) 
			{
				ti = localmem[0][i];
				ssi = selfscore[localmem[0][i]];
				for( m=0; m<clus2; m++ )
				{
					ssm = selfscore[localmem[1][m]];
					tm = localmem[1][m];
					if( ti<tm )
					{
						immin = ti;
						immax = tm;
					}
					else
					{
						immin = tm;
						immax = ti;
					}
					bunbo = MIN( ssi, ssm );
					if( bunbo == 0 )
						newdistmtx[immin][immax-immin] = 2.0; // 2013/Oct/17
					else
#if SKIP
						newdistmtx[immin][immax-immin] = ( 1.0 - naivepairscorefast( mseq1[i], mseq2[m], skiptable1[i], skiptable2[m], penalty_dist ) / bunbo ) * 2.0;
#else
						newdistmtx[immin][immax-immin] = ( 1.0 - naivepairscore11( mseq1[i], mseq2[m], penalty_dist ) / bunbo ) * 2.0;
#endif
				}
			}
#if SKIP
			FreeIntMtx( skiptable1 ); skiptable1 = NULL;
			FreeIntMtx( skiptable2 ); skiptable2 = NULL;
#endif
		}






		pthread_mutex_lock( targ->mutex );
		dep[l].done = 1;
		(*nrunpt)--;
		pthread_cond_broadcast( targ->treecond );

		for( i=0; (j=localmem[0][i])!=-1; i++ )
			strcpy( aseq[j], localcopy[j] );
		for( i=0; (j=localmem[1][i])!=-1; i++ )
			strcpy( aseq[j], localcopy[j] );

//		reporterr( "at step %d\n", l );
//		use_getrusage();

		pthread_mutex_unlock( targ->mutex );



		for( i=0; (j=localmem[0][i])!=-1; i++ )
		{
			if(localcopy[j] ) free( localcopy[j] );
			localcopy[j] = NULL;
		}
		for( i=0; (j=localmem[1][i])!=-1; i++ )
		{
			if( localcopy[j] ) free( localcopy[j] );
			localcopy[j] = NULL;
		}


//		if( topol[l][0] ) free( topol[l][0] );
//		topol[l][0] = NULL;
//		if( topol[l][1] ) free( topol[l][1] );
//		topol[l][1] = NULL;
//		if( topol[l] ) free( topol[l] );
//		topol[l] = NULL;


//		reporterr(       "\n" );
	}
#if SCOREOUT
	reporterr(       "totalscore = %10.2f\n\n", tscore );
#endif
}
#endif


//I think it is for progressive alignment step
static int disttbfast_treebase( int *nlen, char **aseq, int nadd, char *mergeoralign, char **mseq1, char **mseq2, int ***topol, Treedep *dep, double *effarr, double **newdistmtx, int *selfscore, int *alloclen, int (*callback)(int, int, char*) )
{
	int l, len1, len2, i, m, immin, immax;
	int len1nocommongap, len2nocommongap;
	int clus1, clus2;
	double pscore, tscore;
	char *indication1 = NULL, *indication2 = NULL;
	double *effarr1 = NULL;
	double *effarr2 = NULL;
	int *fftlog = NULL; // fixed at 2006/07/26
//	double dumfl = 0.0;
	double dumdb = 0.0;
	int ffttry;
	int m1, m2;
	int *gaplen = NULL;
	int *gapmap = NULL;
	int *alreadyaligned = NULL;
	double **dynamicmtx = NULL;
	double ssi, ssm, bunbo;
	int tm, ti;
	int gapmaplen;
	int **localmem = NULL;
	int posinmem;
#if SKIP
	int **skiptable1 = NULL, **skiptable2 = NULL;
#endif
#if 0
	int i, j;
#endif

//	reporterr( "treebase newdistmtx=%p\n", newdistmtx );

	if( effarr1 == NULL ) 
	{
		effarr1 = AllocateDoubleVec( njob );
		effarr2 = AllocateDoubleVec( njob );
		indication1 = AllocateCharVec( 150 );
		indication2 = AllocateCharVec( 150 );
		fftlog = AllocateIntVec( njob );
		gaplen = AllocateIntVec( *alloclen+10 );
		gapmap = AllocateIntVec( *alloclen+10 );
		alreadyaligned = AllocateIntVec( njob );
		dynamicmtx = AllocateDoubleMtx( nalphabets, nalphabets );
		localmem = AllocateIntMtx( 2, njob+1 );
	}
	for( i=0; i<njob-nadd; i++ ) alreadyaligned[i] = 1;
	for( i=njob-nadd; i<njob; i++ ) alreadyaligned[i] = 0;

	if( callback && callback( 0, 50, "Progressive alignment" ) ) goto chudan_tbfast;

	for( l=0; l<njob; l++ ) fftlog[l] = 1;
	localmem[0][0] = -1;
	localmem[1][0] = -1;
	clus1 = 1;// chain ni hitsuyou

#if 0
	reporterr(       "##### fftwinsize = %d, fftthreshold = %d\n", fftWinSize, fftThreshold );
#endif

#if 0
	for( i=0; i<njob; i++ )
		reporterr(       "TBFAST effarr[%d] = %f\n", i, effarr[i] );
#endif

//	for( i=0; i<njob; i++ ) strcpy( aseq[i], seq[i] );


//	writePre( njob, name, nlen, aseq, 0 );

	tscore = 0.0;
	for( l=0; l<njob-1; l++ )
	{
//		reporterr( " at the beginning of the loop, clus1,clus2=%d,%d\n", clus1, clus2 );

		if(  l > 0 && dep[l].child0 == l-1 && dep[l].child1 == -1 && dep[dep[l].child0].child1 == -1 )
		{
			localmem[0][clus1] = topol[l-1][1][0];
			localmem[0][clus1+1] = -1;

			localmem[1][0] = topol[l][1][0];
			localmem[1][1] = -1;
		}
		else
		{
			localmem[0][0] = -1;
			posinmem = 0;
			//recursion method. update localmem[0] and posinmem values based on other arguments
			topolorder( njob, localmem[0], &posinmem, topol, dep, l, 0 ); //defined in mltaln9.c.
			localmem[1][0] = -1;
			posinmem = 0;
			//recursion method. update localmem[1] and posinmem values based on other arguments
			topolorder( njob, localmem[1], &posinmem, topol, dep, l, 1 );
		}

		if( mergeoralign[l] == 'n' )
		{
//			reporterr(       "SKIP!\n" );
//			free( topol[l][0] ); topol[l][0] = NULL;
//			free( topol[l][1] ); topol[l][1] = NULL;
//			free( topol[l] ); topol[l] = NULL;
			continue;
		}

//		reporterr(       "\ndistfromtip = %f\n", dep[l].distfromtip );
		//fill out matrix based on n_dis_consweight_multi matrix values and dep[l].distfromtip value
		makedynamicmtx( dynamicmtx, n_dis_consweight_multi, dep[l].distfromtip ); //defined in mltaln9.c.
//		makedynamicmtx( dynamicmtx, n_dis_consweight_multi, ( dep[l].distfromtip - 0.2 ) * 3 );


		m1 = topol[l][0][0];
		m2 = topol[l][1][0];
		len1 = strlen( aseq[m1] );
		len2 = strlen( aseq[m2] );
		if( *alloclen < len1 + len2 )
		{
			reporterr(       "\nReallocating.." );
			*alloclen = ( len1 + len2 ) + 1000;
			ReallocateCharMtx( aseq, njob, *alloclen + 10  ); 
			gaplen = realloc( gaplen, ( *alloclen + 10 ) * sizeof( int ) );
			if( gaplen == NULL )
			{
				reporterr(       "Cannot realloc gaplen\n" );
				exit( 1 );
			}
			gapmap = realloc( gapmap, ( *alloclen + 10 ) * sizeof( int ) );
			if( gapmap == NULL )
			{
				reporterr(       "Cannot realloc gapmap\n" );
				exit( 1 );
			}
			reporterr(       "done. *alloclen = %d\n", *alloclen );
		}

#if 1 // CHUUI@@@@
		//update mseq1, effarr1, indication1 values based on other arguments values and calculations on them
		clus1 = fastconjuction_noname( localmem[0], aseq, mseq1, effarr1, effarr, indication1, 0.0 ); //defined in tddis.c.
		//update mseq2, effarr2, indication2 values based on other arguments values and calculations on them
		clus2 = fastconjuction_noname( localmem[1], aseq, mseq2, effarr2, effarr, indication2, 0.0 );
#else
		clus1 = fastconjuction_noname( topol[l][0], aseq, mseq1, effarr1, effarr, indication1, 0.0 );
		clus2 = fastconjuction_noname( topol[l][1], aseq, mseq2, effarr2, effarr, indication2, 0.0 );
//		clus1 = fastconjuction_noweight( topol[l][0], aseq, mseq1, effarr1,  indication1 );
//		clus2 = fastconjuction_noweight( topol[l][1], aseq, mseq2, effarr2,  indication2 );
#endif


		if( mergeoralign[l] == '1' || mergeoralign[l] == '2' )
		{
			newgapstr = "=";
		}
		else
			newgapstr = "-";

		len1nocommongap = len1;
		len2nocommongap = len2;
		if( mergeoralign[l] == '1' ) // nai
		{
			findcommongaps( clus2, mseq2, gapmap ); //defined in addfunctions.c. fill gapmap with values based on gaps positions in mseq2 and value of clus2
			commongappick( clus2, mseq2 ); //defined in mltaln9.c. update mseq2 values based on gaps positions in it and clus2 value
			len2nocommongap = strlen( mseq2[0] );
		}
		else if( mergeoralign[l] == '2' )
		{
			findcommongaps( clus1, mseq1, gapmap ); //fill gapmap with values based on gaps positions in mseq1 and value of clus1
			commongappick( clus1, mseq1 ); //update mseq1 values based on gaps positions in it and clus1 value
			len1nocommongap = strlen( mseq1[0] );
		}

#if 0
    for( i=0; i<clus1; i++ ) 
    {
        if( strlen( mseq1[i] ) != len1 ) 
        {
            reporterr(       "i = %d / %d\n", i, clus1 ); 
            reporterr(       "hairetsu ga kowareta (in treebase, after conjuction) !\n" ); 
            exit( 1 );
        }
    }
    for( j=0; j<clus2; j++ )
    {
        if( strlen( mseq2[j] ) != len2 ) 
        {
            reporterr(       "j = %d / %d\n", j, clus2 ); 
            reporterr(       "hairetsu ga kowareta (in treebase, after conjuction) !\n" ); 
            exit( 1 );
        }
    }
#endif


#if 0
    for( i=0; i<clus1; i++ ) 
    {
        if( strlen( mseq1[i] ) != len1 ) 
        {
            reporterr(       "i = %d / %d\n", i, clus1 ); 
            reporterr(       "hairetsu ga kowareta (in treebase, after free topol) !\n" ); 
            exit( 1 );
        }
    }
    for( j=0; j<clus2; j++ )
    {
        if( strlen( mseq2[j] ) != len2 ) 
        {
            reporterr(       "j = %d / %d\n", j, clus2 ); 
            reporterr(       "hairetsu ga kowareta (in treebase, after free topol) !\n" ); 
            exit( 1 );
        }
    }
#endif


//		fprintf( trap_g, "\nSTEP-%d\n", l );
//		fprintf( trap_g, "group1 = %s\n", indication1 );
//		fprintf( trap_g, "group2 = %s\n", indication2 );

//		reporterr(       "\rSTEP % 5d / %d %d-%d", l+1, njob-1, clus1, clus2 );
		reporterr(       "\rSTEP % 5d / %d ", l+1, njob-1 );
		if( callback && callback( 0, 50+50*l/(njob-1), "Progressive alignment" ) ) goto chudan_tbfast;

#if 0
		reporterr(       "STEP %d /%d\n", l+1, njob-1 );
		reporterr(       "group1 = %.66s", indication1 );
		if( strlen( indication1 ) > 66 ) reporterr(       "..." );
		reporterr(       "\n" );
		reporterr(       "group2 = %.66s", indication2 );
		if( strlen( indication2 ) > 66 ) reporterr(       "..." );
		reporterr(       "\n" );
#endif

/*
		reporterr(       "before align all\n" );
		display( aseq, njob );
		reporterr(       "\n" );
		reporterr(       "before align 1 %s \n", indication1 );
		display( mseq1, clus1 );
		reporterr(       "\n" );
		reporterr(       "before align 2 %s \n", indication2 );
		display( mseq2, clus2 );
		reporterr(       "\n" );
*/


		if( !nevermemsave && ( alg != 'M' && ( len1 > 30000 || len2 > 30000  ) ) )
		{
			reporterr(       "\nlen1=%d, len2=%d, Switching to the memsave mode\n", len1, len2 );
			alg = 'M';
			if( commonIP ) FreeIntMtx( commonIP );
			commonIP = NULL;
			commonAlloc1 = 0;
			commonAlloc2 = 0;
		}

//		if( fftlog[m1] && fftlog[m2] ) ffttry = ( nlen[m1] > clus1 && nlen[m2] > clus2 );
		if( fftlog[m1] && fftlog[m2] ) ffttry = ( nlen[m1] > clus1 && nlen[m2] > clus2 && clus1 < 1000 && clus2 < 1000);
		else						   ffttry = 0;
//		ffttry = ( nlen[m1] > clus1 && nlen[m2] > clus2 && clus1 < 5000 && clus2 < 5000); // v6.708
//		reporterr(       "f=%d, len1/fftlog[m1]=%f, clus1=%d, len2/fftlog[m2]=%f, clus2=%d\n", ffttry, (double)len1/fftlog[m1], clus1, (double)len2/fftlog[m2], clus2 );

		if( force_fft || ( use_fft && ffttry ) )
		{
			reporterr(       "f" );
			if( alg == 'M' )
			{
				reporterr(       "m" );
				//defined in Falign.c.
				//this method aligns mseq1 and mseq2 based on FFT and algorithm selected
				//it is nearly identical to Falign except some small details about the algorithm calculations, and also it supports only M algorithm
				//unlike Falign which supports  a, M, d and A
				pscore = Falign_udpari_long( NULL, NULL, dynamicmtx, mseq1, mseq2, effarr1, effarr2, NULL, NULL, clus1, clus2, *alloclen, fftlog+m1 );
			}
			else
			{
				//defined in Falign.c. this method aligns mseq1 and mseq2 based on FFT and algorithm selected
				pscore = Falign( NULL, NULL, dynamicmtx, mseq1, mseq2, effarr1, effarr2, NULL, NULL, clus1, clus2, *alloclen, fftlog+m1, NULL, 0, NULL );
//				reporterr(       "######### mseq1[0] = %s\n", mseq1[0] );
			}
		}
		else
		{
			reporterr(       "d" );
			fftlog[m1] = 0;
			switch( alg )
			{
				case( 'a' ):
					//apply specific algorithm to align mseq1 and mseq2
					pscore = Aalign( mseq1, mseq2, effarr1, effarr2, clus1, clus2, *alloclen ); //defined in SAalignmm.c.
					break;
				case( 'M' ):
					reporterr(       "m" );
//					reporterr(       "%d-%d", clus1, clus2 );
				//defined in MSalignmm.c. apply specific algorithm to align mseq1 and mseq2
					pscore = MSalignmm( dynamicmtx, mseq1, mseq2, effarr1, effarr2, clus1, clus2, *alloclen, NULL, NULL, NULL, NULL, NULL, 0, NULL, outgap, outgap );
					break;
				case( 'd' ):
					if( 1 && clus1 == 1 && clus2 == 1 )
					{
//						reporterr(       "%d-%d", clus1, clus2 );
						//apply specific algorithm to align mseq1 and mseq2
						pscore = G__align11( dynamicmtx, mseq1, mseq2, *alloclen, outgap, outgap ); //defined in Galign11.c.
					}
					else
					{
//						reporterr(       "%d-%d", clus1, clus2 );
						//defined in Dalignmm.c. apply specific algorithm to align mseq1 and mseq2
						pscore = D__align_ls( dynamicmtx, mseq1, mseq2, effarr1, effarr2, clus1, clus2, *alloclen, NULL, &dumdb, NULL, NULL, NULL, NULL, NULL, 0, NULL, outgap, outgap );
					}
					break;
				case( 'A' ):
					if( clus1 == 1 && clus2 == 1 )
					{
//						reporterr(       "%d-%d", clus1, clus2 );
						//apply specific algorithm to align mseq1 and mseq2
						pscore = G__align11( dynamicmtx, mseq1, mseq2, *alloclen, outgap, outgap ); //defined in Galign11.c.
					}
					else
					{
//						reporterr(       "\n\n %d - %d (%d-%d) : ", topol[l][0][0], topol[l][1][0], clus1, clus2 );
						//defined in Salignmm.c. apply specific algorithm to align mseq1 and mseq2
						pscore = A__align( dynamicmtx, mseq1, mseq2, effarr1, effarr2, clus1, clus2, *alloclen, NULL, &dumdb, NULL, NULL, NULL, NULL, NULL, 0, NULL, outgap, outgap, localmem[0][0], 1 );
					}
					break;
				default:
					ErrorExit( "ERROR IN SOURCE FILE" );
			}
		}
#if SCOREOUT
		reporterr(       "score = %10.2f\n", pscore );
#endif
		tscore += pscore;
		nlen[m1] = 0.5 * ( nlen[m1] + nlen[m2] );

//		writePre( njob, name, nlen, aseq, 0 );

		if( disp ) display( aseq, njob ); //defined in mltaln9.c. displays aseq chars to canvas
//		reporterr(       "\n" );

		if( mergeoralign[l] == '1' ) // jissainiha nai. atarashii hairetsu ha saigo dakara.
		{
			reporterr( "Check source!!!\n" );
			exit( 1 );
		}
		if( mergeoralign[l] == '2' )
		{
//			if( localkeeplength ) ndeleted += deletenewinsertions( clus1, clus2, mseq1, mseq2, NULL );
//			for( i=0; i<clus1; i++ ) reporterr(       ">STEP0 mseq1[%d] = \n%s\n", i, mseq1[i] );
//			for( i=0; i<clus2; i++ ) reporterr(       ">STEP0 mseq2[%d] = \n%s\n", i, mseq2[i] );
			gapmaplen = strlen( mseq1[0] )-len1nocommongap+len1;
			adjustgapmap( gapmaplen, gapmap, mseq1[0] ); //defined in addfunctions.c. update gapmap values based on mseq1[0] chars
#if 0
			reporterr( "\n" );
			for( i=0; i<clus1; i++ ) reporterr(       ">STEP1 mseq1[%d] = \n%s\n", i, mseq1[i] );
			for( i=0; i<clus2; i++ ) reporterr(       ">STEP1 mseq2[%d] = \n%s\n", i, mseq2[i] );
#endif
//			if( clus1 + clus2 < njob ) restorecommongaps( njob, aseq, topol[l][0], topol[l][1], gapmap, *alloclen, '-' );
			if( smoothing )
			{
				restorecommongapssmoothly( njob, njob-(clus1+clus2), aseq, localmem[0], localmem[1], gapmap, *alloclen, '-' ); //defined in addfunctions.c
				findnewgaps( clus1, 0, mseq1, gaplen ); //defined in addfunctions.c. update gaplen based on mseq1 and 0 values
				insertnewgaps_bothorders( njob, alreadyaligned, aseq, localmem[0], localmem[1], gaplen, gapmap, gapmaplen, *alloclen, alg, '-' );
			}
			else
			{
				//update aseq chars - gaps - based on gapmap values
				restorecommongaps( njob, njob-(clus1+clus2), aseq, localmem[0], localmem[1], gapmap, *alloclen, '-' ); //defined in addfunctions.c
				findnewgaps( clus1, 0, mseq1, gaplen ); //defined in addfunctions.c. update gaplen based on mseq1 and 0 values
				insertnewgaps( njob, alreadyaligned, aseq, localmem[0], localmem[1], gaplen, gapmap, *alloclen, alg, '-' );
			}

#if 0
			reporterr( "\n" );
			for( i=0; i<clus1; i++ ) reporterr(       ">STEP3 mseq1[%d] = \n%s\n", i, mseq1[i] );
			for( i=0; i<clus2; i++ ) reporterr(       ">STEP3 mseq2[%d] = \n%s\n", i, mseq2[i] );
#endif

#if 0
			for( i=0; i<njob; i++ ) eq2dash( aseq[i] );
			for( i=0; i<clus1; i++ ) 
			{
				reporterr( "mseq1[%d] bef change = %s\n", i, mseq1[i] );
				eq2dash( mseq1[i] );
				reporterr( "mseq1[%d] aft change = %s\n", i, mseq1[i] );
			}
			for( i=0; i<clus2; i++ ) 
			{
				reporterr( "mseq2[%d] bef change = %s\n", i, mseq2[i] );
				eq2dash( mseq2[i] );
				reporterr( "mseq2[%d] aft change = %s\n", i, mseq2[i] );
			}
			for( i=0; i<clus1; i++ ) eq2dash( mseq1[i] );
			for( i=0; i<clus2; i++ ) eq2dash( mseq2[i] );
#endif


			eq2dashmatometehayaku( mseq1, clus1 ); //defined in addfunctions.c.
			eq2dashmatometehayaku( mseq2, clus2 );

			for( i=0; (m=localmem[1][i])>-1; i++ ) alreadyaligned[m] = 1;
		}

		if( newdistmtx ) // tsukawanai
		{
#if 0
			reporterr( "group1 = " );
			for( i=0; i<clus1; i++ ) reporterr( "%d ", topol[l][0][i] );
			reporterr( "\n" );
			reporterr( "group2 = " );
			for( m=0; m<clus2; m++ ) reporterr( "%d ", topol[l][1][m] );
			reporterr( "\n" );
#endif
#if SKIP
			skiptable1 = AllocateIntMtx( clus1, 0 );
			skiptable2 = AllocateIntMtx( clus2, 0 );
			makeskiptable( clus1, skiptable1, mseq1 ); // allocate suru. //defined in mltaln9.c. fill skiptable1 with values based on gaps in mseq1
			makeskiptable( clus2, skiptable2, mseq2 ); // allocate suru. //fill skiptable2 with values based on gaps in mseq2
#endif
			for( i=0; i<clus1; i++ ) 
			{
#if SKIP
//				makeskiptable( 1, skiptable1, mseq1+i ); // allocate suru.
#endif
				ti = localmem[0][i];
				ssi = selfscore[localmem[0][i]];
				for( m=0; m<clus2; m++ )
				{
					ssm = selfscore[localmem[1][m]];
					tm = localmem[1][m];
					if( ti<tm )
					{
						immin = ti;
						immax = tm;
					}
					else
					{
						immin = tm;
						immax = ti;
					}
					bunbo = MIN( ssi, ssm );
					if( bunbo == 0.0 )
						newdistmtx[immin][immax-immin] = 2.0; // 2013/Oct/17
					else
#if SKIP
						newdistmtx[immin][immax-immin] = ( 1.0 - naivepairscorefast( mseq1[i], mseq2[m], skiptable1[i], skiptable2[m], penalty_dist ) / bunbo ) * 2.0;
#else
						newdistmtx[immin][immax-immin] = ( 1.0 - naivepairscore11( mseq1[i], mseq2[m], penalty_dist ) / bunbo ) * 2.0;
#endif
				}
			}
#if SKIP
			FreeIntMtx( skiptable1 ); skiptable1 = NULL;
			FreeIntMtx( skiptable2 ); skiptable2 = NULL;
#endif
		}

//		free( topol[l][0] ); topol[l][0] = NULL;
//		free( topol[l][1] ); topol[l][1] = NULL;
//		free( topol[l] ); topol[l] = NULL;


//		reporterr(       ">514\n%s\n", aseq[514] );
	}
#if SCOREOUT
	reporterr(       "totalscore = %10.2f\n\n", tscore );
#endif
	Falign( NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, 0, 0, 0, NULL, NULL, 0, NULL );
	Falign_udpari_long( NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, 0, 0, 0, NULL );
	A__align( NULL,  NULL, NULL, NULL, NULL, 0, 0, 0, NULL, NULL, NULL, NULL, NULL, NULL, NULL, 0, NULL, 0, 0, -1, -1 );
	D__align( NULL,  NULL, NULL, NULL, NULL, 0, 0, 0, NULL, NULL, NULL, NULL, NULL, NULL, NULL, 0, NULL, 0, 0 );
	G__align11( NULL, NULL, NULL, 0, 0, 0 ); // iru?
	free( effarr1 );
	free( effarr2 );
	free( indication1 );
	free( indication2 );
	free( fftlog );
	free( gaplen );
	free( gapmap );
	FreeDoubleMtx( dynamicmtx );
	free( alreadyaligned );
	FreeIntMtx( localmem );
	effarr1 = NULL;
	return( 0 );

	chudan_tbfast:

	Falign( NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, 0, 0, 0, NULL, NULL, 0, NULL );
	Falign_udpari_long( NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, 0, 0, 0, NULL );
	A__align( NULL,  NULL, NULL, NULL, NULL, 0, 0, 0, NULL, NULL, NULL, NULL, NULL, NULL, NULL, 0, NULL, 0, 0, -1, -1 );
	D__align( NULL,  NULL, NULL, NULL, NULL, 0, 0, 0, NULL, NULL, NULL, NULL, NULL, NULL, NULL, 0, NULL, 0, 0 );
	G__align11( NULL, NULL, NULL, 0, 0, 0 ); // iru?
	if( effarr1 ) free( effarr1 ); effarr1 = NULL;
	if( effarr2 ) free( effarr2 ); effarr2 = NULL;
	if( indication1 ) free( indication1 ); indication1 = NULL;
	if( indication2 ) free( indication2 ); indication2 = NULL;
	if( fftlog ) free( fftlog ); fftlog = NULL;
	if( gaplen ) free( gaplen ); gaplen = NULL;
	if( gapmap ) free( gapmap ); gapmap = NULL;
	if( alreadyaligned ) free( alreadyaligned ); alreadyaligned = NULL;
	if( dynamicmtx ) FreeDoubleMtx( dynamicmtx ); dynamicmtx = NULL;
	if( localmem ) FreeIntMtx( localmem ); localmem = NULL;
#if SKIP
	if( skiptable1 ) FreeIntMtx( skiptable1 ); skiptable1 = NULL;
	if( skiptable2 ) FreeIntMtx( skiptable2 ); skiptable2 = NULL;
#endif

	return( 1 );
}

//write all options to fp: dorp, penalty, FFT, tree method, algorithm and related values
static void WriteOptions( FILE *fp )
{

	if( dorp == 'd' ) fprintf( fp, "DNA\n" );
	else
	{
		if     ( scoremtx ==  0 ) fprintf( fp, "JTT %dPAM\n", pamN );
		else if( scoremtx ==  1 ) fprintf( fp, "BLOSUM %d\n", nblosum );
		else if( scoremtx ==  2 ) fprintf( fp, "M-Y\n" );
	}
    reporterr(       "Gap Penalty = %+5.2f, %+5.2f, %+5.2f\n", (double)ppenalty/1000, (double)ppenalty_ex/1000, (double)poffset/1000 );
    if( use_fft ) fprintf( fp, "FFT on\n" );

	fprintf( fp, "tree-base method\n" );
	if( tbrweight == 0 ) fprintf( fp, "unweighted\n" );
	else if( tbrweight == 3 ) fprintf( fp, "clustalw-like weighting\n" );
	if( tbitr || tbweight ) 
	{
		fprintf( fp, "iterate at each step\n" );
		if( tbitr && tbrweight == 0 ) fprintf( fp, "  unweighted\n" ); 
		if( tbitr && tbrweight == 3 ) fprintf( fp, "  reversely weighted\n" ); 
		if( tbweight ) fprintf( fp, "  weighted\n" ); 
		fprintf( fp, "\n" );
	}

   	 fprintf( fp, "Gap Penalty = %+5.2f, %+5.2f, %+5.2f\n", (double)ppenalty/1000, (double)ppenalty_ex/1000, (double)poffset/1000 );

	if( alg == 'a' )
		fprintf( fp, "Algorithm A\n" );
	else if( alg == 'A' ) 
		fprintf( fp, "Algorithm A+\n" );
	else
		fprintf( fp, "Unknown algorithm\n" );

	if( treemethod == 'X' )
		fprintf( fp, "Tree = UPGMA (mix).\n" );
	else if( treemethod == 'E' )
		fprintf( fp, "Tree = UPGMA (average).\n" );
	else if( treemethod == 'q' )
		fprintf( fp, "Tree = Minimum linkage.\n" );
	else
		fprintf( fp, "Unknown tree.\n" );

    if( use_fft )
    {
        fprintf( fp, "FFT on\n" );
        if( dorp == 'd' )
            fprintf( fp, "Basis : 4 nucleotides\n" );
        else
        {
            if( fftscore )
                fprintf( fp, "Basis : Polarity and Volume\n" );
            else
                fprintf( fp, "Basis : 20 amino acids\n" );
        }
        fprintf( fp, "Threshold   of anchors = %d%%\n", fftThreshold );
        fprintf( fp, "window size of anchors = %dsites\n", fftWinSize );
    }
	else
        fprintf( fp, "FFT off\n" );
	fflush( fp );
}

static double **preparepartmtx( int nseq )
{
	int i;
	double **val;
	double size;

	val = (double **)calloc( nseq, sizeof( double *) );;
	size = 0;

	if( compacttree == 1 )
	{
		for( i=0; i<nseq; i++ )
		{
			size += (double)sizeof( double ) * nseq;
			if( size > maxdistmtxsize ) 
			{
				reporterr( "\n\nThe size of full distance matrix is estimated to exceed %.2fGB.\n", maxdistmtxsize / 1000 / 1000 /1000 );
				reporterr( "Will try the calculation using a %d x %d matrix.\n", nseq, i );
				reporterr( "This calculation will be slow due to the limited RAM space.\n", i, nseq );
				reporterr( "To avoid the slowdown, please try '--initialramusage xGB' (x>>%.2f),\n", maxdistmtxsize / 1000 / 1000 /1000 );
				reporterr( "if larger RAM space is available.\n" );
				reporterr( "Note that xGB is NOT the upper limit of RAM usage.\n" );
				reporterr( "Two to three times larger space may be used for building a guide tree.\n" );
				reporterr( "Memory usage of the MSA stage depends on similarity of input sequences.\n\n" );
//				reporterr( "If the RAM is small, try '--initialramusage xGB' with a smaller x value.\n" );
				reporterr( "The '--memsavetree' option uses smaller RAM space.\n" );
				reporterr( "If tree-like relationship can be ignored, try '--pileup' or '--randomchain'.\n\n" );
				reporterr( "The result of --initialramusage xGB is almost identical to the default, except for rounding differences.\n" );

				reporterr( "In the cases of --memsavetree, --pileup and --randomchain, the result differs from the default.\n\n" );
				break;
			}
			val[i] = (double *)calloc( nseq, sizeof( double ) );
		}
		if( i == nseq ) reporterr( "The full matrix will be used.\n" );

		for( ;i<nseq; i++ ) val[i] = NULL; // nen no tame
	}
	else
	{
		for( i=0; i<nseq; i++ ) val[i] = NULL; // nen no tame
	}
	return( val );
}


// 1. Parse arguments.
// 2. open inputfile for reading and find sequences count, max length and dna or protein from infp file.
// 3. Check sequences number max and min and also max length of sequence.
// 4. set stack size to 200 * njob
// 5. Read sub alignment file and fill matrix with it.
// 6. Read sequences and their names from infp file into seq, name and nlen arrays.
// 7. Fill arrays from constants call.
// 8. Init prep_g and trap_g files which are used for tracing.
// 9. Check seq characters and report error if unusual character is found.
// 10. read first line from '_guidetree' file and determine the type used and memory size required, then returns char to represent type found
// 11. if treein -> fill randomseed, npickx and maxdistmtxsize based on values read from '_guidetree’ file
// 12. If treein == t -> align sequences from seq using G__align11_noalign based on npickx and other args
// 13. if !treein -> Make distance matrix.
// 14. For each guide tree -> load tree / create chain -> write result to hat2 / construct UPGMA tree - with different parameters -> write order to order file -> check original and additional sequences lengths and make sure they are equal -> if subalignment: check sub alignments sequences -> then progressive alignment -> Making a distance matrix from msa / Making a compact tree from msa
// 15. if keep length -> write deletelist and its indices in '_deletelist’ file -> print deleted letters from sequences to '_deletemap’ file
// 16. Print score
// 17. Write sequences and their names to stdout
// 18. Free allocated memory
int disttbfast( int ngui, int lgui, char **namegui, char **seqgui, int argc, char **argv, int (*callback)(int, int, char*))
{
	int  *nlen = NULL;	
	int  *nogaplen = NULL;	
	char **name = NULL, **seq = NULL;
	char **mseq1 = NULL, **mseq2 = NULL;
	char **bseq = NULL;
	double *eff = NULL;
	int i, j;
	int ***topol = NULL;
	int *addmem = NULL;
	Treedep *dep = NULL; //Treedep is a structure defined in mltaln.h.
	double **len = NULL;
	FILE *infp = NULL;
//	FILE *adfp;
	char c;
	int alloclen;
	double longer, shorter;
	double lenfac;
	double bunbo;

	FILE *orderfp = NULL, *hat2p = NULL;
	int *grpseq = NULL;
	char *tmpseq = NULL;
	int  **pointt = NULL;
	double **mtx = NULL; // by D. Mathog
	int *table1 = NULL;
	char b[B];
	int ien, nlim;
	int includememberres0, includememberres1;
	double unweightedspscore;
	int alignmentlength;
	char *mergeoralign = NULL;
	int foundthebranch;
	int nsubalignments = 0, maxmem;
	int **subtable = NULL;
	int *insubtable = NULL;
	int *preservegaps = NULL;
	char ***subalnpt = NULL;
	int val;
	char **tmpargv = NULL;
	int iguidetree;
	int *selfscore = NULL;
	int calcpairdists;
	int	**skiptable = NULL;
	char algbackup;
	char *originalgaps = NULL;
	char **addbk = NULL;
	int **deletelist = NULL;
	FILE *dlf = NULL;
	int randomseed;
	int **localmem = NULL;
	int posinmem;
// for compacttree
	int *mindistfrom = NULL;
	double *mindist = NULL;
	double **partmtx = NULL;
// for compacttree

	//in this call, ngui = 0
	if( ngui )
	{
		initglobalvariables();
		njob = ngui;
		nlenmax = 0;
		for( i=0; i<njob; i++ )
		{
			ien = strlen( seqgui[i] );
			if( ien > nlenmax ) nlenmax = ien;
		}
		infp = NULL;
//		stderr = fopen( "/dev/null", "a" ); // Windows????
		tmpargv = AllocateCharMtx( argc, 0 );
		for( i=0; i<argc; i++ ) tmpargv[i] = argv[i];
		gmsg = 1;
	}
	else
		gmsg = 0; // iranai  //defined in defs.c.

	dist_tb_fast_arguments( argc, argv ); //parse arguments
	algbackup = alg; // tbfast wo disttbfast ni ketsugou shitatame.
#ifndef enablemultithread
	nthread = 0;
#endif

	//in this call, ngui = 0
	if( ngui )
	{
		for( i=0; i<argc; i++ ) 
		{
//			free( tmpargv[i] );
			argv[i] = tmpargv[i];
		}
		free( tmpargv );
	}
	else
	{
		if( inputfile )
		{
			infp = fopen( inputfile, "r" ); //open inputfile for reading
			if( !infp )
			{
				reporterr(       "Cannot open %s\n", inputfile );
				exit( 1 );
			}
		}
		else
			infp = stdin;
	
		getnumlen( infp ); //defined in io.c. Finds sequences count, max length and dna or protein from infp file.
		rewind( infp );
	}
	
	if( njob > 1000000 ) //if number of sequences > 1,000,000
	{
		reporterr(       "The number of sequences must be < %d\n", 1000000 );
		reporterr(       "Please try the --parttree option for such large data.\n" );
		exit( 1 );
	}

	if( njob < 2 )
	{
		reporterr(       "At least 2 sequences should be input!\n"
						 "Only %d sequence found.\n", njob ); 
		exit( 1 );
	}

	if( specificityconsideration != 0.0 && nlenmax)
	{
		if( nlenmax > 100000 ) //if max length of sequences > 100,000
		{
			reporterr( "\n" );
			reporterr( "Too long to apply --allowshift or --unalignlevel>0\n" );
			reporterr( "Please use the normal mode.\n" );
			reporterr( "Please also note that MAFFT does not assume genomic rearrangements.\n" );
			reporterr( "\n" );
			exit( 1 );
		}
	}

#ifndef mingw
	setstacksize( 200 * njob ); // topolorder() de ookime no stack wo shiyou. //defined in io.c.
	//set stack size to 200 * njob
#endif

	if( subalignment )
	{
		//First call of this method with table = NULL reads number of subalignments in subalignments file and assign to nsubpt
		//also reads max number of spaces in all sequences into maxmem
		//Second call reads data from the file and fills table with it
		readsubalignmentstable( njob, NULL, NULL, &nsubalignments, &maxmem ); //defined in io.c.
		reporterr(       "nsubalignments = %d\n", nsubalignments );
		reporterr(       "maxmem = %d\n", maxmem );
		subtable = AllocateIntMtx( nsubalignments, maxmem+1 );
		insubtable = AllocateIntVec( njob );
		preservegaps = AllocateIntVec( njob );
		for( i=0; i<njob; i++ ) insubtable[i] = 0;
		for( i=0; i<njob; i++ ) preservegaps[i] = 0;
		subalnpt = AllocateCharCub( nsubalignments, maxmem, 0 );
		readsubalignmentstable( njob, subtable, preservegaps, NULL, NULL );
	}


	seq = AllocateCharMtx( njob, nlenmax*1+1 );
	mseq1 = AllocateCharMtx( njob, 0 );
	mseq2 = AllocateCharMtx( njob, 0 );

	eff = AllocateDoubleVec( njob );
	mergeoralign = AllocateCharVec( njob );

	dep = (Treedep *)calloc( njob, sizeof( Treedep ) );

	if( nadd ) addmem = AllocateIntVec( nadd+1 );

	localmem = AllocateIntMtx( 2, njob+1 );

#if 0
	Read( name, nlen, seq );
	readData( infp, name, nlen, seq );
#else
    name = AllocateCharMtx( njob, B+1 );
    nlen = AllocateIntVec( njob ); 
    nogaplen = AllocateIntVec( njob ); 
	if( ngui ) //in this call, ngui = 0
	{
		if( copydatafromgui( namegui, seqgui, name, nlen, seq ) )
			exit( 1 );
	}
	else
	{
		readData_pointer( infp, name, nlen, seq ); //defined in io.c. reads sequences and their names from infp file into seq, name and nlen arrays.
		fclose( infp );
	}
#endif


	constants( njob, seq ); //defined in constants.c.
	//after all this method, n_dis, ribosumdis, amino_dis, amino_dis_consweight_multi, n_dis_consweight_multi,
	//n_disLN, n_disFFT, polarity, volume arrays are initialized and some constants are set.


#if 0
	reporterr(       "params = %d, %d, %d\n", penalty, penalty_ex, offset );
#endif

	initSignalSM(); //defined in io.c. inits signalSM value which is defined in defs.h.

	initFiles(); //defined in io.c. init prep_g and trap_g files. I think these files are for tracing

	WriteOptions( trap_g ); //defined here.  //trap_g defined in defs.h.
	//write all options to trap_g: dorp, penalty, FFT, tree method, algorithm and related values

	c = seqcheck( seq ); //defined in mltaln9.c. check seq characters and report error if unusual character is found
	if( c )
	{
		reporterr(       "Illegal character %c\n", c );
		exit( 1 );
	}

	reporterr(       "\n" );

//	reporterr(       "tuplesize = %d, dorp = %c\n", tuplesize, dorp );
	if( dorp == 'p' && tuplesize != 6 )
	{
		reporterr(       "tuplesize must be 6 for aa sequence\n" );
		exit( 1 );
	}
	if( dorp == 'd' && tuplesize != 6 && tuplesize != 10 )
	{
		reporterr(       "tuplesize must be 6 or 10 for dna sequence\n" );
		exit( 1 );
	}

	if( treein )
	{
		int npickx;
		//read first line from '_guidetree' file and determine the type used and memory size required, then returns char to represent type found
		//fill randomseed, npickx and maxdistmtxsize based on values read from '_guidetree' file
		treein = check_guidetreefile( &randomseed, &npickx, &maxdistmtxsize ); //defined in mltaln9.c.
		if( treein == 't' )
		{
			varpairscore( njob, npickx, nlenmax, seq, randomseed ); //defined here.
			//align sequences from seq using G__align11_noalign based on npickx and other args
			exit( 1 );
		}
		else if( treein == 'c' )
		{
			compacttree = 1;
			treein = 0;
			use_fft = 0; // kankeinai?
//			maxdistmtxsize = 5 * 1000 * 1000; // 5GB. ato de kahen ni suru.
//			maxdistmtxsize =  1.0 * 1000 * 1000 * 1000; // 5GB. ato de kahen ni suru.
		}
		else if( treein == 'C' )
		{
			compacttree = 2;
			treein = 0;
			use_fft = 0; // kankeinai?
		}
		else if( treein == 'a' )
		{
//			reporterr( "Compute pairwise scores\n" );
			if( njob > 200000 )
			{
				reporterr( "Chain?\n" );
				treein = 's';
				nguidetree = 1;
			}
			else if( njob < 100 || 't' == varpairscore( njob, npickx, nlenmax, seq, randomseed ) ) //defined here. align sequences from seq using G__align11_noalign based on npickx and other args
			{
				if( treein == 'c' ) exit( 1 );
				reporterr( "Tree!\n" );
				treein = 0;
				nguidetree = 2;
			}
			else
			{
				reporterr( "Chain!\n" );
				treein = 's';
				nguidetree = 1;
			}
		}
		else if ( treein != 0 ) // auto no toki arieru
			nguidetree = 1;
	}

	if( compacttree == 1 )
	{
		if( maxdistmtxsize > (double)njob * (njob-1) * sizeof( double ) / 2 ) 
		{
			reporterr( "Use conventional tree.\n" );
			compacttree = 0;
		}
	}

	if( !treein )
	{
		reporterr(       "\n\nMaking a distance matrix ..\n" );
		if( callback && callback( 0, 0, "Distance matrix" ) ) goto chudan;

		tmpseq = AllocateCharVec( nlenmax+1 );
		grpseq = AllocateIntVec( nlenmax+1 );
		pointt = AllocateIntMtx( njob, nlenmax+1 );
		if( !compacttree ) mtx = AllocateFloatHalfMtx( njob ); 
		if( dorp == 'd' ) tsize = (int)pow( 4, tuplesize );
		else              tsize = (int)pow( 6, 6 );

		if( dorp == 'd' && tuplesize == 6 )
		{
			lenfaca = D6LENFACA;
			lenfacb = D6LENFACB;
			lenfacc = D6LENFACC;
			lenfacd = D6LENFACD;
		}
		else if( dorp == 'd' && tuplesize == 10 )
		{
			lenfaca = D10LENFACA; //these constants are defined here
			lenfacb = D10LENFACB;
			lenfacc = D10LENFACC;
			lenfacd = D10LENFACD;
		}
		else    
		{
			lenfaca = PLENFACA; //these constants are defined here
			lenfacb = PLENFACB;
			lenfacc = PLENFACC;
			lenfacd = PLENFACD;
		}

		maxl = 0;
		for( i=0; i<njob; i++ ) 
		{
			gappick0( tmpseq, seq[i] ); //defined in mltaln9.c. copy 'seq[i]' chars to 'tmpseq' without gaps chars
			nogaplen[i] = strlen( tmpseq );
			if( nogaplen[i] < 6 )
			{
//				reporterr(       "Seq %d, too short, %d characters\n", i+1, nogaplen[i] );
//				reporterr(       "Please use mafft-ginsi, mafft-linsi or mafft-ginsi\n\n\n" );
//				exit( 1 );
			}
			if( nogaplen[i] > maxl ) maxl = nogaplen[i];
			if( dorp == 'd' ) /* nuc */
			{
				disttbfast_seq_grp_nuc( grpseq, tmpseq ); //defined here. fill grpseq with values from 'amino_grp' based on tmpseq chars
//				makepointtable_nuc( pointt[i], grpseq );
//				makepointtable_nuc_octet( pointt[i], grpseq );
				if( tuplesize == 10 )
					//fill pointt[i] with values based on grpseq and some calculations on it
					disttbfast_makepointtable_nuc_dectet( pointt[i], grpseq ); //defined here.
				else if( tuplesize == 6 )
					//fill pointt[i] with values based on grpseq and some calculations on it
					disttbfast_makepointtable_nuc( pointt[i], grpseq ); //defined here.
				else
				{
					reporterr(       "tuplesize=%d: not supported\n", tuplesize );
					exit( 1 );
				}
			}
			else                 /* amino */
			{
				disttbfast_seq_grp( grpseq, tmpseq ); //defined here. fill grpseq with values from 'amino_grp' based on tmpseq chars
				disttbfast_makepointtable( pointt[i], grpseq ); //defined here. fill pointt[i] with values based on grpseq and some calculations on it
			}
		}
		if( nunknown ) reporterr(       "\nThere are %d ambiguous characters.\n", nunknown );


		if( compacttree )
		{

			reporterr( "Compact tree, step 1\n" );
			mindistfrom = (int *)calloc( njob, sizeof( int ) );
			mindist = (double *)calloc( njob, sizeof( double ) );
			selfscore = (int *)calloc( njob, sizeof( int ) );
			partmtx = preparepartmtx( njob );


			for( i=0; i<njob; i++ )
			{
				table1 = (int *)calloc( tsize, sizeof( int ) );
				if( !table1 ) ErrorExit( "Cannot allocate table1\n" );
				disttbfast_makecompositiontable_p( table1, pointt[i] ); //defined here. update table1 values based on pointt[i] values
				selfscore[i] = commonsextet_p( table1, pointt[i] ); //defined in mltaln9.c. return value based on table1 and pointt[i] values comparison
				free( table1 );
				table1 = NULL;
			}
			commonsextet_p( NULL, NULL ); //defined in mltaln9.c. free allocated memory from previous call.

#ifdef enablemultithread
			if( nthread > 0 )
			{
				compactdistmtxthread_arg_t *targ;
				int jobpos;
				pthread_t *handle;
				pthread_mutex_t mutex;
				double **mindistthread;
				int **mindistfromthread;
				jobpos = 0;
				targ = calloc( nthread, sizeof( compactdistmtxthread_arg_t ) );
				handle = calloc( nthread, sizeof( pthread_t ) );
				mindistthread = AllocateDoubleMtx( nthread, njob );
				mindistfromthread = AllocateIntMtx( nthread, njob );
				pthread_mutex_init( &mutex, NULL );

		
				for( j=0; j<nthread; j++ )
				{
					for( i=0; i<njob; i++ )
					{
						mindistthread[j][i] = 999.9;
						mindistfromthread[j][i] = -1;
					}
					targ[j].thread_no = j;
					targ[j].nogaplen = nogaplen;
					targ[j].pointt = pointt;
					targ[j].selfscore = selfscore;
					targ[j].partmtx = partmtx;
					targ[j].njob = njob;
					targ[j].mindist = mindistthread[j];
					targ[j].mindistfrom = mindistfromthread[j];
					targ[j].jobpospt = &jobpos;
					targ[j].mutex = &mutex;
	
					pthread_create( handle+j, NULL, compactdisthalfmtxthread, (void *)(targ+j) );
				}
		
				for( j=0; j<nthread; j++ ) pthread_join( handle[j], NULL );

				for( i=0; i<njob; i++ )
				{
					mindist[i] = 999.9;
					mindistfrom[i] = -1;
					for( j=0; j<nthread; j++ )
					{
						if( mindistthread[j][i] < mindist[i] )
						{
							mindist[i] = mindistthread[j][i];
							mindistfrom[i] = mindistfromthread[j][i];
						}
					}
				}
				for( i=0; i<njob; i++ ) mindist[i] -= preferenceval( i, mindistfrom[i], njob ); // for debug


				pthread_mutex_destroy( &mutex );
				FreeDoubleMtx( mindistthread );
				FreeIntMtx( mindistfromthread );
				free( handle );
				free( targ );
	
			}
			else
#endif
			{
				compactdistmtxthread_arg_t *targ; //compact dist mtx thread arg t
				int jobpos;
		
				jobpos = 0;
				targ = calloc( 1, sizeof( compactdistmtxthread_arg_t ) );
		
				{
					for( i=0; i<njob; i++ )
					{
						mindist[i] = 999.9;
						mindistfrom[i] = -1;
					}
					targ[0].thread_no = 0;
					targ[0].nogaplen = nogaplen;
					targ[0].pointt = pointt;
					targ[0].selfscore = selfscore;
					targ[0].partmtx = partmtx;
					targ[0].njob = njob;
					targ[0].mindist = mindist;
					targ[0].mindistfrom = mindistfrom;
					targ[0].jobpospt = &jobpos;
	
					compactdisthalfmtxthread( targ ); //defined here.
					//update targ->partmtx, targ->jobpospt, targ->mindist, targ->mindistfrom values based on calculations on other args.
				}

				free( targ );

				for( i=0; i<njob; i++ ) mindist[i] -= preferenceval( i, mindistfrom[i], njob ); // for debug
				//defined here. return value based on i and mindistfrom[i] subtraction
			}

//			for( i=0; i<njob; i++ ) printf( "mindist[%d] = %f, mindistfrom[%d] = %d\n", i, mindist[i], i, mindistfrom[i] );
			reporterr( "\ndone.\n" );

#if 0
			reporterr( "\npartmtx = .\n" );
			for( i=0; i<njob; i++ )
			{
				reporterr( "i=%d\n", i );
				if( partmtx[i] ) for( j=0; j<njob; j++ ) reporterr( "%f ", partmtx[i][j]);
				else reporterr( "nil" );
				reporterr( "\n", i );
			}
#endif
		}
		else
		{
#ifdef enablemultithread
			if( nthread > 0 )
			{
				distancematrixthread_arg_t *targ; 
				int jobpos;
				pthread_t *handle;
				pthread_mutex_t mutex;
	
				jobpos = 0; 
				targ = calloc( nthread, sizeof( distancematrixthread_arg_t ) ); 
				handle = calloc( nthread, sizeof( pthread_t ) ); 
				pthread_mutex_init( &mutex, NULL );
	
				for( i=0; i<nthread; i++ )
				{
					targ[i].thread_no = i;
					targ[i].njob = njob;
					targ[i].jobpospt = &jobpos;
					targ[i].pointt = pointt;
					targ[i].mtx = mtx;
					targ[i].mutex = &mutex;
	
					pthread_create( handle+i, NULL, distancematrixthread, (void *)(targ+i) );
				}
			
				for( i=0; i<nthread; i++ )
				{
					pthread_join( handle[i], NULL );
				}
				pthread_mutex_destroy( &mutex );
				free( handle );
				free( targ );
			}
			else
#endif
			{
				for( i=0; i<njob; i++ )
				{
					table1 = (int *)calloc( tsize, sizeof( int ) );
					if( !table1 ) ErrorExit( "Cannot allocate table1\n" );
					if( i % 100 == 0 )
					{
						reporterr(       "\r% 5d / %d", i+1, njob );
						if( callback && callback( 0, i*25/njob, "Distance matrix" ) ) goto chudan;
					}
					disttbfast_makecompositiontable_p( table1, pointt[i] ); //defined here.
					//update table1 values based on pointt[i] values
			
					for( j=i; j<njob; j++ ) 
					{
						mtx[i][j-i] = (double)commonsextet_p( table1, pointt[j] ); //defined in mltaln9.c. return value based on table1 and pointt[j] values comparison
					} 
					free( table1 ); table1 = NULL;
				}
			}
			reporterr(       "\ndone.\n\n" );
			ien = njob-1;
	
			for( i=0; i<ien; i++ )
			{
				for( j=i+1; j<njob; j++ ) 
				{
					if( nogaplen[i] > nogaplen[j] )
					{
						longer=(double)nogaplen[i];
						shorter=(double)nogaplen[j];
					}
					else
					{
						longer=(double)nogaplen[j];
						shorter=(double)nogaplen[i];
					}
//					if( tuplesize == 6 )
					lenfac = 1.0 / ( shorter / longer * lenfacd + lenfacb / ( longer + lenfacc ) + lenfaca );
//					else
//						lenfac = 1.0;
//					reporterr(       "lenfac = %f (%.0f,%.0f)\n", lenfac, longer, shorter );
					bunbo = MIN( mtx[i][0], mtx[j][0] );
					if( bunbo == 0.0 )
						mtx[i][j-i] = 2.0; // 2013/Oct/17 -> 2bai
					else
						mtx[i][j-i] = ( 1.0 - mtx[i][j-i] / bunbo ) * lenfac * 2.0; // 2013/Oct/17 -> 2bai
//					reporterr(       "##### mtx = %f, mtx[i][0]=%f, mtx[j][0]=%f, bunbo=%f\n", mtx[i][j-i], mtx[i][0], mtx[j][0], bunbo );
				}
			}
			if( disopt )
			{
				for( i=0; i<njob; i++ ) 
				{
					sprintf( b, "=lgth = %04d", nogaplen[i] );
					strins( b, name[i] );
				}
			}
			FreeIntMtx( pointt ); pointt = NULL;
			commonsextet_p( NULL, NULL ); //defined in mltaln9.c. free memory allocated in previous call.
		}
		free( grpseq ); grpseq = NULL;
		free( tmpseq ); tmpseq = NULL;

#if 0 // writehat2 wo kakinaosu -> iguidetree loop nai ni idou
		if( distout )
		{
			hat2p = fopen( "hat2", "w" );
			WriteFloatHat2_pointer_halfmtx( hat2p, njob, name, mtx );
			fclose( hat2p );
		}
#endif

	}
#if 0 
	else 
	{
		reporterr(       "Loading 'hat2' ... " );
		prep = fopen( "hat2", "r" );
		if( prep == NULL ) ErrorExit( "Make hat2." );
		readhat2_double( prep, njob, name, mtx ); // name chuui
		fclose( prep );
		reporterr(       "done.\n" );
	}
#endif

//	reporterr( "after computing distance matrix," );
//	use_getrusage();

	if( nadd && keeplength )
	{
		originalgaps = (char *)calloc( nlenmax+1, sizeof( char) );
		recordoriginalgaps( originalgaps, njob-nadd, seq ); //defined in addfunctions.c. fill originalgaps array with indicators based on gaps found in seq

		if( mapout )
		{
			addbk = (char **)calloc( nadd+1, sizeof( char * ) );
			for( i=0; i<nadd; i++ )
			{
				ien = strlen( seq[njob-nadd+i] );
				addbk[i] = (char *)calloc( ien + 1, sizeof( char ) );
				gappick0( addbk[i], seq[njob-nadd+i] ); //defined in mltaln9.c. copy 'seq[njob-nadd+i]' chars to 'addbk[i]' without gaps chars
			}
			addbk[nadd] = NULL;
		}
		else
			addbk = NULL;
	}
	else
	{
		originalgaps = NULL;
		addbk = NULL;
	}
	


	for( iguidetree=0; iguidetree<nguidetree; iguidetree++ )
//	for( iguidetree=0; ; iguidetree++ )
	{

		alg = algbackup; // tbfast wo disttbfast ni ketsugou shitatame.



		topol = AllocateIntCub( njob, 2, 0 );
		len = AllocateFloatMtx( njob, 2 );

		if( iguidetree == nguidetree - 1 ) calcpairdists = 0;
		else                               calcpairdists = 1;
	
		if( treein )
		{
			nguidetree = 1; //  iranai
			calcpairdists = 0; // iranai
			if( treein == (int)'l' )
			{
				loadtree( njob, topol, len, name, nogaplen, dep, treeout ); //defined in mltaln9.c.
				//read tree from '_guidetree' file and fill nogaplen and dep with values from it. i think if treeout != 0, it writes the tree read to infile.tree
			}
			else if( treein == (int)'s' )
			{
				createchain( njob, topol, len, name, nogaplen, dep, treeout, 1, randomseed ); //defined in mltaln9.c.
				//update topol values based on shuffled integers of njob and write them to _guidetree
				//update len values based on calculations on njob and write them to _guidetree
				//update dep. treeout and shuffle (1) determines some decisions like writing in 'infile.tree' and shuffling order array
				nthread = 0;
				weight = 0; // mafft.tmpl kara idou
				tbrweight = 0; // mafft.tmpl kara idou
			}
			else if( treein == (int)'p' )
			{
				createchain( njob, topol, len, name, nogaplen, dep, treeout, 0, randomseed ); //defined in mltaln9.c.
				//update topol values based on shuffled integers of njob and write them to _guidetree
				//update len values based on calculations on njob and write them to _guidetree
				//update dep. treeout and shuffle (1) determines some decisions like writing in 'infile.tree' and shuffling order array
				nthread = 0;
				weight = 0; // mafft.tmpl kara idou
				tbrweight = 0; // mafft.tmpl kara idou
			}
			else
			{
				reporterr( "Error. treein = %d or %c\n", treein, treein );
				exit( 1 );
			}
		}
		else if( topin )
		{
			reporterr(       "Loading a topology ... " );
			reporterr(       "--topin has been disabled\n" );
			exit( 1 );
//			loadtop( njob, mtx, topol, len );
//			FreeFloatHalfMtx( mtx, njob );
		}
		else
		{
			if( distout )
			{
				hat2p = fopen( "hat2", "w" ); //open 'hat2' file for writing
				WriteFloatHat2_pointer_halfmtx( hat2p, njob, name, mtx ); //defined in io.c.
				//write name and mtx content to hat2p file in specific format
				// writehat2 wo kakinaosu
				fclose( hat2p );

				if( !treeout && noalign )  // 2016Jul31
				{
					writeData_pointer( stdout, njob, name, nlen, seq ); //defined in io.c. write sequences and their names to stdout
					reporterr(       "\n" );
					SHOWVERSION;
					goto chudan; //label defined later to free all allocated memory and return to the caller
//					return( 0 );
				}
			}

			if( subalignment ) // merge error no tame
			{
				reporterr(       "Constructing a UPGMA tree ... " );
				fixed_supg_double_realloc_nobk_halfmtx_treeout_constrained( njob, mtx, topol, len, name, nlen, dep, nsubalignments, subtable, !calcpairdists );
				//defined in mltaln9.c.
				//update mtx, topol, len, dep, subtable values based on other args values and calculations to determine tree shape and values - to be studied in details later -.
				if( !calcpairdists ) 
				{
					FreeFloatHalfMtx( mtx, njob ); mtx = NULL;
				}
			}
			else if( compacttree ) // merge error no tame
			{
				reporterr(       "Constructing a tree ... " );
				compacttree_memsaveselectable( njob, partmtx, mindistfrom, mindist, pointt, selfscore, bseq, skiptable, topol, len, name, nogaplen, dep, treeout, compacttree, 1 );
				//defined in mltaln9.c.
				//update mindistfrom, mindist, topol, len, dep values based on tree construction - needs to be studied in details later -
				if( mindistfrom ) free( mindistfrom ); mindistfrom = NULL;
				if( mindist ) free( mindist );; mindist = NULL;
				if( selfscore ) free( selfscore ); selfscore = NULL;
				if( bseq ) FreeCharMtx( bseq ); bseq = NULL; // nikaime dake
				if( skiptable) FreeIntMtx( skiptable ); skiptable = NULL; // nikaime dake
				if( pointt ) FreeIntMtx( pointt ); pointt = NULL; // ikkaime dake.
				free( partmtx );
			}
			else if( treeout ) // merge error no tame
			{
				reporterr(       "Constructing a UPGMA tree (treeout, efffree=%d) ... ", !calcpairdists );
				fixed_musclesupg_double_realloc_nobk_halfmtx_treeout_memsave( njob, mtx, topol, len, name, nogaplen, dep, !calcpairdists );
				//defined in mltaln9.c.
				//update topol, len, dep values based on tree construction - needs to be studied in details later -
				if( !calcpairdists )
				{
					FreeFloatHalfMtx( mtx, njob ); mtx = NULL;
				}
			}
			else
			{
				reporterr(       "Constructing a UPGMA tree (efffree=%d) ... ", !calcpairdists );
				fixed_musclesupg_double_realloc_nobk_halfmtx_memsave( njob, mtx, topol, len, dep, 1, !calcpairdists );
				//defined in mltaln9.c.
				//update topol, len, dep values based on tree construction - needs to be studied in details later -
				if( !calcpairdists )
				{
					FreeFloatHalfMtx( mtx, njob ); mtx = NULL;
				}
			}
		}
//		else 
//			ErrorExit( "Unknown tree method\n" );




		if( calcpairdists ) selfscore = AllocateIntVec( njob );


		if( callback && callback( 0, 25, "Guide tree" ) ) goto chudan; //label defined later to free all allocated memory and return to the caller
		reporterr(       "\ndone.\n\n" );
		if( callback && callback( 0, 50, "Guide tree" ) ) goto chudan;

		if( sparsepickup && iguidetree == nguidetree-1 )
		{
			reporterr(       "Sparsepickup! \n" );
			pickup( njob, nogaplen, topol, name, seq ); //defined here.
			//pickup part of sequences based on sparsepickup value, and write selected sequences and their names to stdout, and write not used ones to notused file
			reporterr(       "done. \n" );
			SHOWVERSION;
			goto chudan;
		}
//		reporterr( "after tree building" );
//		use_getrusage();


		if( treein == 's' || treein == 'p' )
		{
			localmem[0][0] = topol[0][0][0];
			for( i=1; i<njob; i++ )
				localmem[0][i] = topol[i-1][1][0];
		}
		else
		{
			localmem[0][0] = -1;
			posinmem = 0;
			topolorder( njob, localmem[0], &posinmem, topol, dep, njob-2, 2 ); //defined in mltaln9.c.
			//recursion method. update localmem[0] and posinmem values based on other arguments
		}
	
		orderfp = fopen( "order", "w" ); //open 'order' file for writing
		if( !orderfp )
		{
			reporterr(       "Cannot open 'order'\n" );
			exit( 1 );
		}
#if 0
		for( i=0; (j=topol[njob-2][0][i])!=-1; i++ )
		{
			fprintf( orderfp, "%d\n", j );
		}
		for( i=0; (j=topol[njob-2][1][i])!=-1; i++ )
		{
			fprintf( orderfp, "%d\n", j );
		}
#else
		for( i=0; i<njob; i++ )
			fprintf( orderfp, "%d\n", localmem[0][i] ); //write localmem - which contains topol order - to 'order' file
#endif
		fclose( orderfp );



	
		if( ( treeout || distout )  && noalign ) 
		{
			writeData_pointer( stdout, njob, name, nlen, seq ); //defined in io.c. write sequences and their names to stdout
			reporterr(       "\n" );
			SHOWVERSION;
			goto chudan;
//			return( 0 );
		}
		
		if( tbrweight )
		{
			weight = 3; 
#if 0
			utree = 0; counteff( njob, topol, len, eff ); utree = 1;
#else
			counteff_simple_double_nostatic_memsave( njob, topol, len, dep, eff ); //defined in mltaln9.c.
			//update eff values based on len and internal eff values and other values calculations
//			counteff_simple_double_nostatic( njob, topol, len, eff );
#endif
		}
		else
		{
			for( i=0; i<njob; i++ ) eff[i] = 1.0;
		}
	
#if 0
		for( i=0; i<njob; i++ )
			reporterr(       "eff[%d] = %20.16f\n", i, eff[i] );
		exit( 1 );
#endif
	
	
		FreeFloatMtx( len ); len = NULL;
	
		bseq = AllocateCharMtx( njob, nlenmax*2+1 );
		alloclen = nlenmax*2+1;


	
		if( nadd )
		{
			alignmentlength = strlen( seq[0] );
			for( i=0; i<njob-nadd; i++ ) //check original sequences lengths and make sure they are equal
			{
				if( alignmentlength != strlen( seq[i] ) )
				{
					reporterr(       "#################################################################################\n" );
					reporterr(       "# ERROR!\n" );
					reporterr(       "# The original %d sequences must be aligned\n", njob-nadd );
					reporterr(       "# alignmentlength = %d, but strlen(seq[%d])=%d\n", alignmentlength, i, (int)strlen( seq[i] ) );
					reporterr(       "#################################################################################\n" );
					goto chudan; // TEST!!
					//exit( 1 );
				}
			}
			if( addprofile )
			{
				alignmentlength = strlen( seq[njob-nadd] );
				for( i=njob-nadd; i<njob; i++ ) //check additional sequences lengths and make sure they are equal
				{
					if( alignmentlength != strlen( seq[i] ) )
					{
						reporterr(       "###############################################################################\n" );
						reporterr(       "# ERROR!\n" );
						reporterr(       "# The %d additional sequences must be aligned\n", nadd );
						reporterr(       "# Otherwise, try the '--add' option, instead of '--addprofile' option.\n" );
						reporterr(       "###############################################################################\n" );
						exit( 1 );
					}
				}
				for( i=0; i<nadd; i++ ) addmem[i] = njob-nadd+i;
				addmem[nadd] = -1;
				foundthebranch = 0;
				for( i=0; i<njob-1; i++ )
				{
					localmem[0][0] = -1;
					posinmem = 0;
					topolorder( njob, localmem[0], &posinmem, topol, dep, i, 0 ); //defined in mltaln9.c.
					//recursion method. update localmem[0] and posinmem values based on other arguments
					localmem[1][0] = -1;
					posinmem = 0;
					topolorder( njob, localmem[1], &posinmem, topol, dep, i, 1 );
					//recursion method. update localmem[1] and posinmem values based on other arguments

					//I think this method returns 1 if localmem[0] and addmem contains the same values, 0 otherwise.
					if( samemember( localmem[0], addmem ) ) // jissainiha nai   //defined in mltaln9.c.
					{
						mergeoralign[i] = '1';
						foundthebranch = 1;
					}
					else if( samemember( localmem[1], addmem ) ) // samemembern ni henkou kanou
					{
						mergeoralign[i] = '2';
						foundthebranch = 1;
					}
					else
					{
						mergeoralign[i] = 'n';
					}
				}
				if( !foundthebranch )
				{
					reporterr(       "###############################################################################\n" );
					reporterr(       "# ERROR!\n" );
					reporterr(       "# There is no appropriate position to add the %d sequences in the guide tree.\n", nadd );
					reporterr(       "# Check whether the %d sequences form a monophyletic cluster.\n", nadd );
					reporterr(       "# If not, try the '--add' option, instead of the '--addprofile' option.\n" );
					reporterr(       "############################################################################### \n" );
					exit( 1 );
				}
				commongappick( nadd, seq+njob-nadd ); //defined in mltaln9.c.
				//update seq+njob-nadd values based on gaps positions in it and nadd value
				for( i=njob-nadd; i<njob; i++ ) strcpy( bseq[i], seq[i] ); //copy seq to bseq
			}
			else
			{
				for( i=0; i<njob-1; i++ ) mergeoralign[i] = 'n';
#if 0
				for( j=njob-nadd; j<njob; j++ )
				{
					addmem[0] = j;
					addmem[1] = -1;
					for( i=0; i<njob-1; i++ )
					{
						reporterr( "Looking for samemember, %d-%d/%d\n", j, i, njob );
						localmem[0][0] = -1;
						posinmem = 0;
						topolorder( njob, localmem[0], &posinmem, topol, dep, i, 0 );
						localmem[1][0] = -1;
						posinmem = 0;
						topolorder( njob, localmem[1], &posinmem, topol, dep, i, 1 );

						if( samemembern( localmem[0], addmem, 1 ) ) // arieru
						{
//							reporterr(       "HIT!\n" );
							if( mergeoralign[i] != 'n' ) mergeoralign[i] = 'w';
							else mergeoralign[i] = '1';
						}
						else if( samemembern( localmem[1], addmem, 1 ) )
						{
//							reporterr(       "HIT!\n" );
							if( mergeoralign[i] != 'n' ) mergeoralign[i] = 'w';
							else mergeoralign[i] = '2';
						}
					}
				}
#else
				for( i=0; i<njob-1; i++ )
				{
//					reporterr( "Looking for samemember, %d-%d/%d\n", j, i, njob );
					localmem[0][0] = -1;
					posinmem = 0;
					topolorder( njob, localmem[0], &posinmem, topol, dep, i, 0 ); //defined in mltaln9.c.
					//recursion method. update localmem[0] and posinmem values based on other arguments
					localmem[1][0] = -1;
					posinmem = 0;
					topolorder( njob, localmem[1], &posinmem, topol, dep, i, 1 );
					//recursion method. update localmem[1] and posinmem values based on other arguments

					for( j=njob-nadd; j<njob; j++ )
					{
						addmem[0] = j;
						addmem[1] = -1;

						//I think this method returns 1 if localmem[0] and addmem contains the same values, 0 otherwise.
						if( samemembern( localmem[0], addmem, 1 ) ) // arieru   //defined in mltaln9.c.
						{
//							reporterr(       "HIT!\n" );
							if( mergeoralign[i] != 'n' ) mergeoralign[i] = 'w';
							else mergeoralign[i] = '1';
						}
						//I think this method returns 1 if localmem[1] and addmem contains the same values, 0 otherwise.
						else if( samemembern( localmem[1], addmem, 1 ) )
						{
//							reporterr(       "HIT!\n" );
							if( mergeoralign[i] != 'n' ) mergeoralign[i] = 'w';
							else mergeoralign[i] = '2';
						}
					}
				}
#endif
		
				for( i=0; i<nadd; i++ ) addmem[i] = njob-nadd+i;
				addmem[nadd] = -1;
				nlim = njob-1;
//				for( i=0; i<njob-1; i++ )
				for( i=0; i<nlim; i++ )
				{
					localmem[0][0] = -1;
					posinmem = 0;
					topolorder( njob, localmem[0], &posinmem, topol, dep, i, 0 ); //defined in mltaln9.c.
					//recursion method. update localmem[0] and posinmem values based on other arguments
					localmem[1][0] = -1;
					posinmem = 0;
					topolorder( njob, localmem[1], &posinmem, topol, dep, i, 1 );
					//recursion method. update localmem[1] and posinmem values based on other arguments

					includememberres0 = includemember( localmem[0], addmem ); //defined in mltaln9.c. return 0 if localmem[0] not included in addmem, 1 otherwise
					includememberres1 = includemember( localmem[1], addmem ); //return 0 if localmem[1] not included in addmem, 1 otherwise
//					if( includemember( topol[i][0], addmem ) && includemember( topol[i][1], addmem ) )
					if( includememberres0 && includememberres1 )
					{
						mergeoralign[i] = 'w';
					}
					else if( includememberres0 )
					{
						mergeoralign[i] = '1';
					}
					else if( includememberres1 )
					{
						mergeoralign[i] = '2';
					}
				}
#if 0
				for( i=0; i<njob-1; i++ )
				{
					reporterr(       "mem0 = " );
					for( j=0; topol[i][0][j]>-1; j++ )	reporterr(       "%d ", topol[i][0][j] );
					reporterr(       "\n" );
					reporterr(       "mem1 = " );
					for( j=0; topol[i][1][j]>-1; j++ )	reporterr(       "%d ", topol[i][1][j] );
					reporterr(       "\n" );
					reporterr(       "i=%d, mergeoralign[] = %c\n", i, mergeoralign[i] );
				}
#endif
				for( i=njob-nadd; i<njob; i++ ) gappick0( bseq[i], seq[i] ); //defined in mltaln9.c.
				//copy 'seq' chars to 'bseq' without gaps chars
			}
	
//			if( !keeplength ) commongappick( njob-nadd, seq );
			commongappick( njob-nadd, seq ); //defined in mltaln9.c. update seq values based on gaps positions in it and njob-nadd value

			for( i=0; i<njob-nadd; i++ ) strcpy( bseq[i], seq[i] ); //copy seq to bseq

		}
//--------------- kokokara ----
		else if( subalignment )
		{
			for( i=0; i<njob-1; i++ ) mergeoralign[i] = 'a';
			for( i=0; i<nsubalignments; i++ )
			{
				reporterr(       "Checking subalignment %d:\n", i+1 );
				alignmentlength = strlen( seq[subtable[i][0]] );
//				for( j=0; subtable[i][j]!=-1; j++ )
//					reporterr(       " %d. %-30.30s\n", subtable[i][j]+1, name[subtable[i][j]]+1 );
				for( j=0; subtable[i][j]!=-1; j++ ) //check subalignments lengths and make sure they are aligned
				{
					if( subtable[i][j] >= njob ) // check sumi
					{
						reporterr(       "No such sequence, %d.\n", subtable[i][j]+1 );
						exit( 1 );
					}
					if( alignmentlength != strlen( seq[subtable[i][j]] ) )
					{
						reporterr(       "\n" );
						reporterr(       "###############################################################################\n" );
						reporterr(       "# ERROR!\n" );
						reporterr(       "# Subalignment %d must be aligned.\n", i+1 );
						reporterr(       "# Please check the alignment lengths of following sequences.\n" );
						reporterr(       "#\n" );
						reporterr(       "# %d. %-10.10s -> %d letters (including gaps)\n", subtable[i][0]+1, name[subtable[i][0]]+1, alignmentlength );
						reporterr(       "# %d. %-10.10s -> %d letters (including gaps)\n", subtable[i][j]+1, name[subtable[i][j]]+1, (int)strlen( seq[subtable[i][j]] ) );
						reporterr(       "#\n" );
						reporterr(       "# See http://mafft.cbrc.jp/alignment/software/merge.html for details.\n" );
						if( subalignmentoffset )
						{
							reporterr(       "#\n" );
							reporterr(       "# You specified seed alignment(s) consisting of %d sequences.\n", subalignmentoffset );
							reporterr(       "# In this case, the rule of numbering is:\n" );
							reporterr(       "#   The aligned seed sequences are numbered as 1 .. %d\n", subalignmentoffset );
							reporterr(       "#   The input sequences to be aligned are numbered as %d .. %d\n", subalignmentoffset+1, subalignmentoffset+njob );
						}
						reporterr(       "###############################################################################\n" );
						reporterr(       "\n" );
						goto chudan; // TEST!!
						//exit( 1 );
					}
					insubtable[subtable[i][j]] = 1;
				}
				for( j=0; j<njob-1; j++ )
				{
					//return 0 if topol[j][0] not included in subtable[i], 1 otherwise
					//return 0 if topol[j][1] not included in subtable[i], 1 otherwise
					if( includemember( topol[j][0], subtable[i] ) && includemember( topol[j][1], subtable[i] ) ) //defined in mltaln9.c.
					{
						mergeoralign[j] = 'n';
					}
				}
				foundthebranch = 0;
				for( j=0; j<njob-1; j++ )
				{
					//I think this method returns 1 if topol[j][0] and subtable[i] contains the same values, 0 otherwise.
					//I think this method returns 1 if topol[j][1] and subtable[i] contains the same values, 0 otherwise.
					if( samemember( topol[j][0], subtable[i] ) || samemember( topol[j][1], subtable[i] ) ) //defined in mltaln9.c.
					{
						foundthebranch = 1;
						reporterr(       " -> OK\n" );
						break;
					}
				}
				if( !foundthebranch )
				{
					system( "cp infile.tree GuideTree" ); // tekitou   //copy GuideTree to infile.tree
					reporterr(       "\n" );
					reporterr(       "###############################################################################\n" );
					reporterr(       "# ERROR!\n" );
					reporterr(       "# Subalignment %d does not seem to form a monophyletic cluster\n", i+1 );
					reporterr(       "# in the guide tree ('GuideTree' in this directory) internally computed.\n" );
					reporterr(       "# If you really want to use this subalignment, pelase give a tree with --treein \n" );
					reporterr(       "# http://mafft.cbrc.jp/alignment/software/treein.html\n" );
					reporterr(       "# http://mafft.cbrc.jp/alignment/software/merge.html\n" );
					if( subalignmentoffset )
					{
						reporterr(       "#\n" );
						reporterr(       "# You specified seed alignment(s) consisting of %d sequences.\n", subalignmentoffset );
						reporterr(       "# In this case, the rule of numbering is:\n" );
						reporterr(       "#   The aligned seed sequences are numbered as 1 .. %d\n", subalignmentoffset );
						reporterr(       "#   The input sequences to be aligned are numbered as %d .. %d\n", subalignmentoffset+1, subalignmentoffset+njob );
					}
					reporterr(       "############################################################################### \n" );
					reporterr(       "\n" );
					goto chudan; // TEST!!
					//exit( 1 );
				}
//				commongappick( seq[subtable[i]], subalignment[i] ); // irukamo
			}
#if 0
			for( i=0; i<njob-1; i++ )
			{
				reporterr(       "STEP %d\n", i+1 );
				reporterr(       "group1 = " );
				for( j=0; topol[i][0][j] != -1; j++ )
					reporterr(       "%d ", topol[i][0][j]+1 );
				reporterr(       "\n" );
				reporterr(       "group2 = " );
				for( j=0; topol[i][1][j] != -1; j++ )
					reporterr(       "%d ", topol[i][1][j]+1 );
				reporterr(       "\n" );
				reporterr(       "%d -> %c\n\n", i, mergeoralign[i] );
			}
#endif
	
			for( i=0; i<njob; i++ ) 
			{
				if( insubtable[i] ) strcpy( bseq[i], seq[i] ); //copy seq to bseq
				else gappick0( bseq[i], seq[i] ); //defined in mltaln9.c. copy 'seq' chars to 'bseq' without gaps chars
			}
	
			for( i=0; i<nsubalignments; i++ ) 
			{
				for( j=0; subtable[i][j]!=-1; j++ ) subalnpt[i][j] = bseq[subtable[i][j]];
				if( !preservegaps[i] ) commongappick( j, subalnpt[i] ); //defined in mltaln9.c. update subalnpt[i] values based on gaps positions in it and j value
			}
	
#if 0 // --> iguidetree loop no soto he
			FreeIntMtx( subtable );
			free( insubtable );
			for( i=0; i<nsubalignments; i++ ) free( subalnpt[i] );
			free( subalnpt );
			free( preservegaps );
#endif
		}
//--------------- kokomade ----
		else
		{
			for( i=0; i<njob; i++ ) gappick0( bseq[i], seq[i] ); //defined in mltaln9.c. copy 'seq' chars to 'bseq' without gaps chars
			for( i=0; i<njob-1; i++ ) mergeoralign[i] = 'a';
		}

		//defined in mltaln9.c. calculates score between seq[i] and seq[i] based on penalty_dist and amino_dis values
		//so this loop fills selfscore with naive score values between each seq char and itself
		if( calcpairdists ) for( i=0; i<njob; i++ ) selfscore[i] = naivepairscore11( seq[i], seq[i], penalty_dist ); // (int)?
	
		reporterr(       "Progressive alignment %d/%d... \n", iguidetree+1, nguidetree );

//		reporterr( "\nbefore treebase" );
//		use_getrusage();
	
#ifdef enablemultithread
		if( nthread > 0 && nadd == 0 )
		{
			treebasethread_arg_t *targ; 
			int jobpos;
			pthread_t *handle;
			pthread_mutex_t mutex;
			pthread_cond_t treecond;
			int *fftlog;
			int nrun;
			int nthread_yoyu;
	
			nthread_yoyu = nthread * 1;
			nrun = 0;
			jobpos = 0; 
			targ = calloc( nthread_yoyu, sizeof( treebasethread_arg_t ) ); 
			fftlog = AllocateIntVec( njob );
			handle = calloc( nthread_yoyu, sizeof( pthread_t ) ); 
			pthread_mutex_init( &mutex, NULL );
			pthread_cond_init( &treecond, NULL );
	
			for( i=0; i<njob; i++ ) dep[i].done = 0; 
			for( i=0; i<njob; i++ ) fftlog[i] = 1; 
	
			for( i=0; i<nthread_yoyu; i++ )
			{
				targ[i].thread_no = i;
				targ[i].njob = njob;
				targ[i].nrunpt = &nrun;
				targ[i].nlen = nlen;
				targ[i].jobpospt = &jobpos;
				targ[i].topol = topol;
				targ[i].dep = dep;
				targ[i].aseq = bseq;
				targ[i].effarr = eff;
				targ[i].alloclenpt = &alloclen;
				targ[i].fftlog = fftlog;
				targ[i].mergeoralign = mergeoralign;
#if 1 // tsuneni SEPARATELYCALCPAIRDISTS
				targ[i].newdistmtx = NULL;
				targ[i].selfscore = NULL;
#else
				if( calcpairdists ) // except for last cycle
				{
					targ[i].newdistmtx = mtx;
					targ[i].selfscore = selfscore;
				}
				else
				{
					targ[i].newdistmtx = NULL;
					targ[i].selfscore = NULL;
				}
#endif
				targ[i].mutex = &mutex;
				targ[i].treecond = &treecond;
	
				pthread_create( handle+i, NULL, treebasethread, (void *)(targ+i) );
			}

			for( i=0; i<nthread_yoyu; i++ )
			{
				pthread_join( handle[i], NULL );
			}
			pthread_mutex_destroy( &mutex );
			pthread_cond_destroy( &treecond );
			free( handle );
			free( targ );
			free( fftlog );

//			reporterr( "after treebasethread, " );
//			use_getrusage();
		}
		else
#endif
		{
#if 0
			if( calcpairdists ) // except for last
			{
				if( disttbfast_treebase( nlen, bseq, nadd, mergeoralign, mseq1, mseq2, topol, dep, eff, mtx, selfscore, &alloclen, callback ) ) goto chudan;
			}
			else
#endif
			{
//				if( disttbfast_treebase( keeplength && (iguidetree==nguidetree-1), nlen, bseq, nadd, mergeoralign, mseq1, mseq2, topol, dep, eff, NULL, NULL, deletemap, deletelag, &alloclen, callback ) ) goto chudan;

				//defined here. I think it is for progressive alignment
				if( disttbfast_treebase( nlen, bseq, nadd, mergeoralign, mseq1, mseq2, topol, dep, eff, NULL, NULL, &alloclen, callback ) ) goto chudan;
			}
		}
//		reporterr( "after treebase, " );
//		use_getrusage();
		reporterr(       "\ndone.\n\n" );
		if( callback && callback( 0, 100, "Progressive alignment" ) ) goto chudan;
//		free( topol[njob-1][0] ); topol[njob-1][0]=NULL;
//		free( topol[njob-1][1] ); topol[njob-1][1]=NULL;
//		free( topol[njob-1] ); topol[njob-1]=NULL;
//		free( topol ); topol=NULL;
		FreeIntCub( topol ); topol = NULL;
//		reporterr( "after freeing topol, " );
//		use_getrusage();


//		reporterr( "compacttree = %d, calcpairdist = %d\n", compacttree, calcpairdists );


//		reporterr( "\nbseq[njob-3] = %s\n", bseq[njob-3] );
//		reporterr( "bseq[njob-2] = %s\n", bseq[njob-2] );
//		reporterr( "bseq[njob-1] = %s\n", bseq[njob-1] );



// Distance matrix from MSA SEPARATELYCALCPAIRDISTS
//		if( iguidetree < nguidetree-1 )
#ifdef enablemultithread
//		if( nthread>0 && nadd==0 ) if( calcpairdists )
		if( calcpairdists && !compacttree )
#else
//		if( 0 && nadd==0 ) if( calcpairdists ) // zettai nai
		if( calcpairdists && !compacttree )
#endif
		{
			reporterr( "Making a distance matrix from msa.. \n" );
			skiptable = AllocateIntMtx( njob, 0 );
			makeskiptable( njob, skiptable, bseq ); // allocate suru.
#ifdef enablemultithread
			if( nthread > 0 )
			{
				msadistmtxthread_arg_t *targ;
				Jobtable jobpos;
				pthread_t *handle;
				pthread_mutex_t mutex;
	
				jobpos.i = 0;
				jobpos.j = 0;
	
				targ = calloc( nthread, sizeof( msadistmtxthread_arg_t ) );
				handle = calloc( nthread, sizeof( pthread_t ) );
				pthread_mutex_init( &mutex, NULL );
	
				for( i=0; i<nthread; i++ )
				{
					targ[i].thread_no = i;
					targ[i].njob = njob;
					targ[i].selfscore = selfscore;
					targ[i].iscore = mtx;
					targ[i].seq = bseq;
					targ[i].skiptable = skiptable;
					targ[i].jobpospt = &jobpos;
					targ[i].mutex = &mutex;
	
					pthread_create( handle+i, NULL, msadistmtxthread, (void *)(targ+i) );
				}
	
				for( i=0; i<nthread; i++ )
				{
					pthread_join( handle[i], NULL );
				}
				pthread_mutex_destroy( &mutex );
				free( handle );
				free( targ );
			}
			else
#endif
			{
//				reporterr( "Check source!\n" );
//				exit( 1 );

#if 1
				msadistmtxthread_arg_t *targ;
				Jobtable jobpos;

				jobpos.i = 0;
				jobpos.j = 0;
	
				targ = calloc( 1, sizeof( msadistmtxthread_arg_t ) );
	
				{
					targ[0].thread_no = 0;
					targ[0].njob = njob;
					targ[0].selfscore = selfscore;
					targ[0].iscore = mtx;
					targ[0].seq = bseq;
					targ[0].skiptable = skiptable;
					targ[0].jobpospt = &jobpos;
	
					msadistmtxthread( targ );
				}
	
				free( targ );
#endif
			}
			if( skiptable) FreeIntMtx( skiptable ); skiptable = NULL;
			reporterr(       "\ndone.\n\n" );
			free( selfscore ); selfscore = NULL;
			FreeCharMtx( bseq ); bseq = NULL;
		}
		else if( calcpairdists && compacttree )
		{
			reporterr( "Making a compact tree from msa, step 1.. \n" );
			skiptable = AllocateIntMtx( njob, 0 );
			makeskiptable( njob, skiptable, bseq ); // allocate suru.   //defined in mltaln9.c. fill skiptable matrix with values based on gaps in bseq
			mindistfrom = (int *)calloc( njob, sizeof( int ) );
			mindist = (double *)calloc( njob, sizeof( double ) );
			partmtx = preparepartmtx( njob ); //defined here.
#ifdef enablemultithread
			if( nthread > 0 )
			{
				msacompactdistmtxthread_arg_t *targ;
				int jobpos;
				pthread_t *handle;
				pthread_mutex_t mutex;
				double **mindistthread;
				int **mindistfromthread;

				mindistthread = AllocateDoubleMtx( nthread, njob );
				mindistfromthread = AllocateIntMtx( nthread, njob );
				targ = calloc( nthread, sizeof( msacompactdistmtxthread_arg_t ) );
				handle = calloc( nthread, sizeof( pthread_t ) );
				pthread_mutex_init( &mutex, NULL );
				jobpos = 0;

				for( i=0; i<nthread; i++ )
				{
					for( j=0; j<njob; j++ )
					{
						mindistthread[i][j] = 999.9;
						mindistfromthread[i][j] = -1;
					}
					targ[i].thread_no = i;
					targ[i].njob = njob;
					targ[i].selfscore = selfscore;
					targ[i].partmtx = partmtx;
					targ[i].seq = bseq;
					targ[i].skiptable = skiptable;
					targ[i].jobpospt = &jobpos;
					targ[i].mindistfrom = mindistfromthread[i];
					targ[i].mindist = mindistthread[i];
					targ[i].mutex = &mutex;
	
					pthread_create( handle+i, NULL, msacompactdisthalfmtxthread, (void *)(targ+i) );
				}
	
				for( i=0; i<nthread; i++ ) pthread_join( handle[i], NULL );
				pthread_mutex_destroy( &mutex );

				for( i=0; i<njob; i++ )
				{
					mindist[i] = 999.9;
					mindistfrom[i] = -1;
					for( j=0; j<nthread; j++ )
					{
						if( mindistthread[j][i] < mindist[i] )
						{
							mindist[i] = mindistthread[j][i];
							mindistfrom[i] = mindistfromthread[j][i];
						}
					}
				}
				for( i=0; i<njob; i++ ) mindist[i] -= preferenceval( i, mindistfrom[i], njob ); // for debug

				free( handle );
				free( targ );
				FreeDoubleMtx( mindistthread );
				FreeIntMtx( mindistfromthread );
			}
			else
#endif
			{
				msacompactdistmtxthread_arg_t *targ;
				int jobpos;
				jobpos = 0;
				targ = calloc( 1, sizeof( msacompactdistmtxthread_arg_t ) );

				{
					for( j=0; j<njob; j++ )
					{
						mindist[j] = 999.9;
						mindistfrom[j] = -1;
					}
					targ[0].thread_no = 0;
					targ[0].njob = njob;
					targ[0].selfscore = selfscore;
					targ[0].partmtx = partmtx;
					targ[0].seq = bseq;
					targ[0].skiptable = skiptable;
					targ[0].jobpospt = &jobpos;
					targ[0].mindistfrom = mindistfrom;
					targ[0].mindist = mindist;
	
					msacompactdisthalfmtxthread( targ );
//					msacompactdistmtxthread( targ );
				}
				free( targ );
				for( i=0; i<njob; i++ ) mindist[i] -= preferenceval( i, mindistfrom[i], njob ); // for debug  //defined here. return value based on i and mindistfrom[i] subtraction
			}
//			free( selfscore ); selfscore = NULL; // mada tsukau
//			FreeCharMtx( bseq ); bseq = NULL; // mada tsukau
//			if( skiptable) FreeIntMtx( skiptable ); skiptable = NULL;

//			for( i=0; i<njob; i++ ) printf( "mindist[%d] = %f\n", i, mindist[i] );
//			exit( 1 );
		}
// Distance matrix from MSA end
//		reporterr( "at the end of guidetree loop, " );
//		use_getrusage();

	}
#if DEBUG
	reporterr(       "closing trap_g\n" );
#endif
//	fclose( trap_g );
//	reporterr( "after guidetree loop, " );
//	use_getrusage();

	if( keeplength )
	{

		dlf = fopen( "_deletelist", "w" ); //open '_deletelist' file for writing
		deletelist = (int **)calloc( nadd+1, sizeof( int * ) );
		for( i=0; i<nadd; i++ )
		{
			deletelist[i] = calloc( 1, sizeof( int ) );
			deletelist[i][0] = -1;
		}
		deletelist[nadd] = NULL;
		//update bseq, bseq+njob-nadd and deletelist based on conditions on bseq and bseq+njob-nadd values
		ndeleted = deletenewinsertions_whole( njob-nadd, nadd, bseq, bseq+njob-nadd, deletelist ); //defined in addfunctions.c.

		for( i=0; i<nadd; i++ )  //write deletelist and its indices in '_deletelist' file
		{
			if( deletelist[i] )
				for( j=0; deletelist[i][j]!=-1; j++ )
					fprintf( dlf, "%d %d\n", njob-nadd+i, deletelist[i][j] ); // 0origin
		}
		fclose( dlf );

		 //defined in addfunctions.c. restore gaps from originalgaps to bseq.
		restoreoriginalgaps( njob, bseq, originalgaps );

		if( mapout ) //print deleted letters from sequences to '_deletemap' file
		{
			dlf = fopen( "_deletemap", "w" ); //open '_deletemap' file for writing
			reconstructdeletemap( nadd, addbk, deletelist, bseq+njob-nadd, dlf, name+njob-nadd ); //defined in addfunctions.c. print addbk and nadd content to dlf file
			FreeCharMtx( addbk );
			addbk = NULL;
			fclose( dlf );
		}

		FreeIntMtx( deletelist );
		deletelist = NULL;
	}

	if( scoreout )
	{
		unweightedspscore = plainscore( njob, bseq ); //defined in mltaln9.c. calculates score between all sequences in bseq based on naive score method
		reporterr(       "\nSCORE %s = %.0f, ", "(treebase)", unweightedspscore );
		reporterr(       "SCORE / residue = %f", unweightedspscore / ( njob * strlen( bseq[0] ) ) );
		reporterr(       "\n\n" );
	}

#if DEBUG
	reporterr(       "writing alignment to stdout\n" );
#endif


	val = 0;
	if( ngui ) 
	{
		ien = strlen( bseq[0] );
		if( ien > lgui )
		{
			reporterr( "alignmentlength = %d, gui allocated %d", ien, lgui );
			val = GUI_LENGTHOVER;
		}
		else
		{
			for( i=0; i<njob; i++ ) //copy bseq to seqgui
			{
#if 1
				strcpy( seqgui[i], bseq[i] );
#else
				free( seqgui[i] );
				seqgui[i] =  bseq[i];
#endif
			}
		}
	}
	else
	{
		writeData_pointer( stdout, njob, name, nlen, bseq ); //defined in io.c. write sequences and their names to stdout
	} 

	 //defined in mltaln9.c. return score of all sequences in bseq based on naive pair score method
	if( spscoreout ) reporterr( "Unweighted sum-of-pairs score = %10.5f\n", sumofpairsscore( njob, bseq ) );
	SHOWVERSION;
	if( ndeleted > 0 )
	{
		reporterr( "\nTo keep the alignment length, %d letters were DELETED.\n", ndeleted );
		if( mapout )
			reporterr( "The deleted letters are shown in the (filename).map file.\n" );
		else
			reporterr( "To know the positions of deleted letters, rerun the same command with the --mapout option.\n" );
	}



	if( subalignment )
	{
		FreeIntMtx( subtable );
		free( insubtable );
		for( i=0; i<nsubalignments; i++ ) free( subalnpt[i] );
		free( subalnpt );
		free( preservegaps );
	}


#if 1 // seqgui[i] =  bseq[i] no toki bseq ha free shinai
	FreeCharMtx( bseq );
#endif
	FreeCharMtx( name );
    free( nlen );

	free( mergeoralign );
	FreeCharMtx( seq );
    free( nogaplen );

	free( mseq1 );
	free( mseq2 );
//	FreeIntCub( topol ); // 
//	FreeFloatMtx( len ); //
//	free( mergeoralign ); //
	free( dep );

	if( nadd ) free( addmem );
	FreeIntMtx( localmem );
	free( eff );
	freeconstants(); //defined in constants.c. free all constants allocated memory.
	closeFiles(); //in io.c. close prep_g and trap_g files.
	FreeCommonIP();
	if( originalgaps ) free( originalgaps ); originalgaps = NULL;
	if( deletelist ) FreeIntMtx( deletelist ); deletelist = NULL;

//	use_getrusage();

	return( val );

chudan: //free all allocated memory in the file

	if( nlen ) free( nlen ); nlen = NULL;
	if( seq ) FreeCharMtx( seq ); seq = NULL;
	if( mseq1 ) free( mseq1 ); mseq1 = NULL;
	if( mseq2 ) free( mseq2 ); mseq2 = NULL;
//	if( topol ) 
//	{
//		for( i=0; i<njob; i++ )
//		{
//			if( topol[i] && topol[i][0] ) 
//			{
//				free( topol[i][0] ); topol[i][0] = NULL;
//			}
//			if( topol[i] && topol[i][1] ) 
//			{
//				free( topol[i][1] ); topol[i][1] = NULL;
//			}
//			if( topol[i] ) free( topol[i] ); topol[i] = NULL;
//		}
//		free( topol ); topol = NULL;
//	}
	if( topol ) FreeIntCub( topol ); topol = NULL;
	if( len ) FreeFloatMtx( len ); len = NULL;
	if( eff ) free( eff ); eff = NULL;
	if( mergeoralign ) free( mergeoralign ); mergeoralign = NULL;
	if( dep ) free( dep ); dep = NULL;
	if( addmem ) free( addmem ); addmem = NULL;
	if( localmem ) FreeIntMtx( localmem ); localmem = NULL;
	if( name ) FreeCharMtx( name ); name = NULL;
	if( nogaplen ) free( nogaplen ); nogaplen = NULL;

	if( tmpseq ) free( tmpseq ); tmpseq = NULL;
	if( grpseq ) free( grpseq ); grpseq = NULL;
	if( pointt ) FreeIntMtx( pointt ); pointt = NULL;
	if( mtx ) FreeFloatHalfMtx( mtx, njob ); mtx = NULL;
	if( table1 ) free( table1 ); table1 = NULL;

	if( bseq ) FreeCharMtx( bseq ); bseq = NULL;
	if( selfscore ) free( selfscore ); selfscore = NULL;
	if( skiptable ) FreeIntMtx( skiptable ); skiptable = NULL;
	if( originalgaps ) free( originalgaps ); originalgaps = NULL;
	if( deletelist ) FreeIntMtx( deletelist ); deletelist = NULL;


	if( subtable ) FreeIntMtx( subtable ); subtable = NULL;
	if( insubtable ) free( insubtable ); insubtable = NULL;
	for( i=0; i<nsubalignments; i++ ) 
	{
		if( subalnpt[i] ) free( subalnpt[i] ); subalnpt[i] = NULL;
	}
	if( subalnpt ) free( subalnpt ); subalnpt = NULL;
	if( preservegaps ) free( preservegaps ); preservegaps = NULL;


	if( mindistfrom ) free( mindistfrom ); mindistfrom = NULL;
	if( mindist ) free( mindist ); mindist = NULL;

	freeconstants();
	closeFiles();
	FreeCommonIP();

	return( GUI_CANCEL ); //defined in mafft.h and = 3.
}

int disttbfast_main( int argc, char **argv )
{
	int res = disttbfast( 0, 0, NULL, NULL, argc, argv, NULL );
	if( res == GUI_CANCEL ) res = 0; // treeout de goto chudan wo riyousuru
	return res;
}
