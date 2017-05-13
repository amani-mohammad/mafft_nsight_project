//what is direction list?
#include "mltaln.h"

#define DEBUG 0
#define IODEBUG 0
#define SCOREOUT 0

#define GLOBAL 0

#define END_OF_VEC -1

int nadd;
double thresholdtorev;
int dodp;
int addfragment;
int mode = '2';
int reflim = 1000;
int contrastsort = 1;

typedef struct _thread_arg
{
	int iend; 
	char **seq;
	int *map;
	char *tmpseq;
	int *res;
	int **spointt;
	short *table1;
	int iq;
#ifdef enablemultithread
	int *jshare;
	int thread_no;
	pthread_mutex_t *mutex_counter;
#endif
} thread_arg_t;

typedef struct _selfdpthread_arg
{
	int iend;
	char **seq;
	double *res;
#ifdef enablemultithread
	int *jshare;
	int thread_no;
	pthread_mutex_t *mutex_counter;
#endif
} selfdpthread_arg_t;

typedef struct _contrast
{
	int pos; 
	double dif;
} contrastarr;

static void	*selfdpthread( void *arg )
{
	selfdpthread_arg_t *targ = (selfdpthread_arg_t *)arg;
	int iend = targ->iend;
	char **seq = targ->seq;
	double *res = targ->res;
#ifdef enablemultithread
	int thread_no = targ->thread_no;
	int *jshare = targ->jshare; 
#endif
	int j;
	char **revseq;

	revseq = AllocateCharMtx( 1, nlenmax+1 );

	j = -1;
	while( 1 )
	{
#ifdef enablemultithread
		if( nthread )
		{
			pthread_mutex_lock( targ->mutex_counter );
			j = *jshare;
			if( j%100 == 0 ) reporterr( "%d / %d (thread %d)   \r", j, iend, thread_no );
			if( j == iend )
			{
				pthread_mutex_unlock( targ->mutex_counter );
				break;
			}
			++(*jshare);
			pthread_mutex_unlock( targ->mutex_counter );
		}
		else
#endif
		{
			j++;
			if( j%100 == 0 ) reporterr( "%d / %d      \r", j, iend );
			if( j == iend ) 
			{
				break;
			}
		}

		sreverse( revseq[0], seq[j] );
#if GLOBAL
		res[j] =  G__align11_noalign( n_dis_consweight_multi, penalty, penalty_ex, seq+j, seq+j, 0 );
		res[j] -= G__align11_noalign( n_dis_consweight_multi, penalty, penalty_ex, seq+j, revseq, 0 );
#else
		res[j] =  L__align11_noalign( n_dis_consweight_multi, seq+j, seq+j );
		res[j] -= L__align11_noalign( n_dis_consweight_multi, seq+j, revseq );
#endif
	}

	creverse( 0 );
	FreeCharMtx( revseq );
#if GLOBAL
	G__align11_noalign( NULL, 0, 0, NULL, NULL, 0 );
#else
	L__align11_noalign( NULL, NULL, NULL );
#endif
	return( NULL );
}

#if 0
static void partshuffle( int size, int outsize, int *ary )
{
	int i;

//	reporterr( "ary before shuffle = \n" );
 //   for(i=0;i<size;i++) reporterr( "%d ", ary[i] );
//	reporterr( "\n" );

    for(i=0;i<outsize;i++)
    {
        int j = rand()%size;
        int t = ary[i];
        ary[i] = ary[j];
        ary[j] = t;
    }

//	reporterr( "ary after shuffle = \n" );
 //   for(i=0;i<outsize;i++) reporterr( "%d ", ary[i] );
//	reporterr( "|" );
 //   for(i=outsize;i<size;i++) reporterr( "%d ", ary[i] );
//	reporterr( "\n" );
}
#endif

void make_direction_list_arguments( int argc, char *argv[] )
{
    int c;

    //set default values before parsing input arguments
	nthread = 1; //defined in defs.c
	inputfile = NULL; //defined in defs.c
	nadd = 0;
	dodp = 0;
	alg = 'a'; //defined in defs.h
	alg = 'm'; //why set to 'a' then 'm'
	dorp = NOTSPECIFIED; //dna or protein
	fmodel = 0; //defined in defs.h
//	ppenalty = (int)( -2.0 * 1000 - 0.5 );
//	ppenalty_ex = (int)( -0.1 * 1000 - 0.5 );
//	poffset = (int)( 0.1 * 1000 - 0.5 ); 
	ppenalty = NOTSPECIFIED; //in defs.h
	ppenalty_ex = NOTSPECIFIED; //in defs.h
	poffset = NOTSPECIFIED; //in defs.h
	kimuraR = 2; //in defs.h
	pamN = 200; //in defs.h
	thresholdtorev = 0.0;
	addfragment = 0;


    while( --argc > 0 && (*++argv)[0] == '-' )
	{
        while ( (c = *++argv[0]) )
		{
            switch( c )
            {
				case 'i':
					inputfile = *++argv;
					fprintf( stderr, "inputfile = %s\n", inputfile );
					--argc;
					goto nextoption;
				case 'I':
					nadd = myatoi( *++argv );
					fprintf( stderr, "nadd = %d\n", nadd );
					--argc;
					goto nextoption;
				case 'C':
					nthread = myatoi( *++argv );
					fprintf( stderr, "nthread = %d\n", nthread );
					--argc; 
					goto nextoption;
				case 'f':
					ppenalty = (int)( atof( *++argv ) * 1000 - 0.5 );
//					fprintf( stderr, "ppenalty = %d\n", ppenalty );
					--argc;
					goto nextoption;
				case 'g':
					ppenalty_ex = (int)( atof( *++argv ) * 1000 - 0.5 );
					fprintf( stderr, "ppenalty_ex = %d\n", ppenalty_ex );
					--argc;
					goto nextoption;
				case 'h':
					poffset = (int)( atof( *++argv ) * 1000 - 0.5 );
//					fprintf( stderr, "poffset = %d\n", poffset );
					--argc;
					goto nextoption;
				case 'k':
					kimuraR = myatoi( *++argv );
					fprintf( stderr, "kappa = %d\n", kimuraR );
					--argc;
					goto nextoption;
				case 'j':
					pamN = myatoi( *++argv );
					scoremtx = 0;
					TMorJTT = JTT; //defined in mltaln.h and = 201
					fprintf( stderr, "jtt/kimura %d\n", pamN );
					--argc;
					goto nextoption;
				case 't':
					thresholdtorev = atof( *++argv );
					fprintf( stderr, "thresholdtorev = %f\n", thresholdtorev );
					--argc; 
					goto nextoption;
				case 'o':
					mode = *(*++argv);
					fprintf( stderr, "mode = %c\n", mode );
					--argc; 
					goto nextoption;
				case 'r':
					reflim = myatoi(*++argv);
					fprintf( stderr, "reflim = %d\n", reflim );
					--argc; 
					goto nextoption;
				case 'c':
					contrastsort = 0;
					break;
				case 'd':
					dodp = 1;
					break;
				case 'F':
					addfragment = 1;
					break;
#if 1
				case 'a':
					fmodel = 1;
					break;
#endif
				case 'S':
					alg = 'S';
					break;
				case 'M':
					alg = 'M';
					break;
				case 'm':
					alg = 'm';
					break;
				case 'G':
					alg = 'G';
					break;
				case 'D':
					dorp = 'd';
					break;
				case 'P':
					dorp = 'p';
					break;
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
        fprintf( stderr, "options: Check source file !\n" );
        exit( 1 );
    }
	if( tbitr == 1 && outgap == 0 )
	{
		fprintf( stderr, "conflicting options : o, m or u\n" );
		exit( 1 );
	}
}




//set values from amino_grp array to 'grp' based on 'seq' chars
void make_direction_list_seq_grp_nuc( int *grp, char *seq )
{
	int tmp;
	int *grpbk = grp;
	while( *seq )
	{
		tmp = amino_grp[(int)*seq++]; //amino_grp is defined in defs.h and set in constants.c
		if( tmp < 4 )
			*grp++ = tmp;
		else
//			fprintf( stderr, "WARNING : Unknown character %c\r", *(seq-1) );
			;
	}
	*grp = END_OF_VEC;
	if( grp - grpbk < 6 )
	{
//		fprintf( stderr, "\n\nWARNING: Too short.\nPlease also consider use mafft-ginsi, mafft-linsi or mafft-ginsi.\n\n\n" );
//		exit( 1 );
		*grpbk = -1;
	}
}

void make_direction_list_seq_grp( int *grp, char *seq )
{
	int tmp;
	int *grpbk = grp;
	while( *seq )
	{
		tmp = amino_grp[(int)*seq++];
		if( tmp < 6 )
			*grp++ = tmp;
		else
//			fprintf( stderr, "WARNING : Unknown character %c\r", *(seq-1) );
			;
	}
	*grp = END_OF_VEC;
	if( grp - grpbk < 6 )
	{
//		fprintf( stderr, "\n\nWARNING: Too short.\nPlease also consider use mafft-ginsi, mafft-linsi or mafft-ginsi.\n\n\n" );
//		exit( 1 );
		*grpbk = -1;
	}
}

void make_direction_list_makecompositiontable_p( short *table, int *pointt )
{
	int point;

	while( ( point = *pointt++ ) != END_OF_VEC )
		table[point]++;
}


void make_direction_list_makepointtable_nuc( int *pointt, int *n )
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

void make_direction_list_makepointtable( int *pointt, int *n )
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

static int localcommonsextet_p2( short *table, int *pointt )
{
	int value = 0;
	short tmp;
	int point;
	short *memo;
	int *ct;
	int *cp;

	if( *pointt == -1 )
		return( 0 );

	memo = (short *)calloc( tsize, sizeof( short ) );
	if( !memo ) ErrorExit( "Cannot allocate memo\n" );
	ct = (int *)calloc( MIN( maxl, tsize )+1, sizeof( int ) ); // chuui!!
	if( !ct ) ErrorExit( "Cannot allocate memo\n" );

	cp = ct;
	while( ( point = *pointt++ ) != END_OF_VEC )
	{
		tmp = memo[point]++;
		if( tmp < table[point] )
			value++;
		if( tmp == 0 ) *cp++ = point;
	}
	*cp = END_OF_VEC;
	
	cp =  ct;
	while( *cp != END_OF_VEC )
		memo[*cp++] = 0;

	free( memo );
	free( ct );
	return( value );
}

static int compfunc( const void *a, const void *b )
{
	return ((contrastarr *)b)->dif - ((contrastarr *)a)->dif; // correct
//	return ((contrastarr *)a)->dif - ((contrastarr *)b)->dif; // incorrect!
} 

static void makecontrastorder6mer( int *order, int **pointt, int **pointt_rev, char **seq, int iend, int shift )
{
	int i;
	double *res;
	contrastarr *arr;
	short *table1, *table1_rev;


	arr = calloc( iend, sizeof( contrastarr ) );
	res = calloc( iend, sizeof( double ) );

	for( i=0; i<iend; i++ )
	{
		if( i % 100 == 1 ) reporterr( "%d   \r", i );
		table1 = (short *)calloc( tsize, sizeof( short ) );
		if( !table1 ) ErrorExit( "Cannot allocate table1\n" );
		make_distance_makecompositiontable_p( table1, pointt[i] ); //defined in mafft-distance.c
		res[i] = localcommonsextet_p2( table1, pointt[i] );
		free( table1 );

		table1_rev = (short *)calloc( tsize, sizeof( short ) );
		if( !table1_rev ) ErrorExit( "Cannot allocate table1\n" );
		make_distance_makecompositiontable_p( table1_rev, pointt_rev[i] );
		res[i] -= localcommonsextet_p2( table1_rev, pointt[i] );
		free( table1_rev );

	}

	for( i=0; i<iend; i++ )
	{
		arr[i].pos = i;
		arr[i].dif = res[i];
	}

	//compfunc compares arr elements based on 'diff' attribute
	qsort( arr, iend, sizeof( contrastarr ), compfunc ); //qsort sorts array. it sorts 'arr' with compfunc comparing function

	for( i=0; i<iend; i++ )
		order[i] = arr[i].pos + shift;

//	for( i=0; i<iend; i++ ) reporterr( "%f\n", arr[i].dif );
//	reporterr( "highest contrast, %s\n", seq[order[0]] );
//	reporterr( "lowest contrast, %s\n", seq[order[iend-1]] );

	free( arr );
	free( res );

}
static void makecontrastorder( int *order, char **seq, int iend, int shift )
{
	int i;
	double *res;
	contrastarr *arr;

	arr = calloc( iend, sizeof( contrastarr ) );
	res = calloc( iend, sizeof( double ) );

#ifdef enablemultithread
	if( nthread )
	{
		int j;
		pthread_t *handle;
		pthread_mutex_t mutex_counter;
		selfdpthread_arg_t *targ;
		int *jsharept;
		
		targ = calloc( nthread, sizeof( selfdpthread_arg_t ) );
		handle = calloc( nthread, sizeof( pthread_t ) );
		pthread_mutex_init( &mutex_counter, NULL );
		jsharept = calloc( 1, sizeof(int) );
		*jsharept = 0;
		
		for( j=0; j<nthread; j++ )
		{
			targ[j].iend = iend;
			targ[j].seq = seq;
			targ[j].res = res; 
			targ[j].jshare = jsharept;
			targ[j].mutex_counter = &mutex_counter;
			targ[j].thread_no = j;
			pthread_create( handle+j, NULL, selfdpthread, (void *)(targ+j) );
		}
		for( j=0; j<nthread; j++ ) pthread_join( handle[j], NULL );
		pthread_mutex_destroy( &mutex_counter );
		free( handle );
		free( targ );
		free( jsharept );
	}
	else
#endif
	{
		selfdpthread_arg_t *targ;
		targ = calloc( 1, sizeof( selfdpthread_arg_t ) );
		targ[0].iend = iend;
		targ[0].seq = seq;
		targ[0].res = res; 
		selfdpthread( targ );
		free( targ );
	}

	for( i=0; i<iend; i++ )
	{
		arr[i].pos = i;
		arr[i].dif = res[i];
	}

	//compfunc compares arr elements based on 'diff' attribute
	qsort( arr, iend, sizeof( contrastarr ), compfunc ); //qsort sorts array. it sorts 'arr' with compfunc comparing function

	for( i=0; i<iend; i++ )
		order[i] = arr[i].pos + shift;

//	for( i=0; i<iend; i++ ) reporterr( "%f\n", arr[i].dif );
//	reporterr( "highest contrast, %s\n", seq[order[0]] );
//	reporterr( "lowest contrast, %s\n", seq[order[iend-1]] );

	free( arr );
	free( res );

}


static void	*directionthread( void *arg )
{
	thread_arg_t *targ = (thread_arg_t *)arg;
	int iend = targ->iend;
	char **seq = targ->seq;
	int *map = targ->map;
	char *tmpseq = targ->tmpseq;
	int *res = targ->res;
	int **spointt = targ->spointt;
	short *table1 = targ->table1;
//	int iq = targ->iq;
#ifdef enablemultithread
//	int thread_no = targ->thread_no;
	int *jshare = targ->jshare; 
#endif
	int j;
	char **mseq1, **mseq2;


	if( dodp ) // nakuserukamo
	{
		mseq1 = AllocateCharMtx( 1, 0 );
		mseq2 = AllocateCharMtx( 1, 0 );
	}

	j = -1;
	while( 1 )
	{
#ifdef enablemultithread
		if( nthread )
		{
			pthread_mutex_lock( targ->mutex_counter );
			j = *jshare;
			if( j == iend )
			{
				pthread_mutex_unlock( targ->mutex_counter );
				break;
			}
			++(*jshare);
			pthread_mutex_unlock( targ->mutex_counter );
		}
		else
#endif
		{
			j++;
			if( j == iend ) 
			{
//				if( iq%100==1 ) fprintf( stderr, "\r %d / %d  \r", iq, njob );
				break;
			}
		}


		if( dodp )
		{
//			strcpy( mseq1[0], tmpseq );
//			strcpy( mseq2[0], seq[j] );
			mseq1[0] = tmpseq;
			mseq2[0] = seq[map[j]];
#if GLOBAL
			res[j] = G__align11_noalign( n_dis_consweight_multi, penalty, penalty_ex, mseq1, mseq2, 0 );
#else
			res[j] = L__align11_noalign( n_dis_consweight_multi, mseq1, mseq2 );
#endif
		}
		else
		{
//			reporterr( "\n\nj=%d, map[j]=%d\n\n", j, map[j] );
			res[j] = localcommonsextet_p2( table1, spointt[map[j]] );
		}
	}
	if( dodp ) // nakuserukamo
	{
		free( mseq1 );
		free( mseq2 );
#if GLOBAL
		G__align11_noalign( NULL, 0, 0, NULL, NULL, 0 );
#else
		L__align11_noalign( NULL, NULL, NULL );
#endif
	}
//	else
//		if( nthread )  // inthread == 0 no toki free suru to, error. nazeda
//			localcommonsextet_p( NULL, NULL );
	return( NULL );
}

//int make_direction_list_main( int argc, char *argv[] )
int main( int argc, char *argv[] )
{
	static int  *nlen;	
	static int  *nogaplen;	
	static char **name, **seq;
	int i, j, istart, iend, ic;
	FILE *infp;
//	FILE *adfp;
	char c;

	int *grpseq;
	char *tmpseq, *revseq;
	int  **pointt, **pointt_rev, **spointt;
	double res_forward, res_reverse, res_max;
	int ires, mres, mres2;
	int *res, *resr, *resf;
	int *map;
	static short *table1, *table1_rev;
	static char **mseq1f, **mseq1r, **mseq2;
	int *contrastorder;

	make_direction_list_arguments( argc, argv );
#ifndef enablemultithread
	nthread = 0;
#endif

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

	getnumlen( infp ); //defined in io.c. finds sequences count, max length and dna or protein
	rewind( infp );

	if( alg == 'a' ) //I need to know what a, G, S stands for ??
	{
		if( nlenmax < 10000 )
			alg = 'G';
		else
			alg = 'S';
	}

	seq = AllocateCharMtx( njob, nlenmax*1+1 ); //allocate 2d char matrix for sequences

#if 0
	Read( name, nlen, seq );
	readData( infp, name, nlen, seq );
#else
    name = AllocateCharMtx( njob, B+1 ); //B is defined in mltaln.h and = 256
    nlen = AllocateIntVec( njob ); //allocate int vector with num of sequences cells
    nogaplen = AllocateIntVec( njob ); 
	readData_pointer( infp, name, nlen, seq ); //defined in io.c. fill matrices of sequences, sequences names and lengths
	fclose( infp ); //close input stream

	if( dorp != 'd' ) //if input sequences are not DNA, i. e. proteins
	{
		fprintf( stderr, "Not necessary!\n" );
		for( i=0; i<njob; i++ ) 
			fprintf( stdout, "_F_%-10.10s\n", name[i]+1 );
		exit( 1 );
	}
#endif

	constants( njob, seq ); //this method is defined in constants.c
	//after all this method, n_dis, ribosumdis, amino_dis, amino_dis_consweight_multi, n_dis_consweight_multi,
	//n_disLN, n_disFFT, polarity, volume arrays are initialized and some constants are set.

#if 0
	fprintf( stderr, "params = %d, %d, %d\n", penalty, penalty_ex, offset );
#endif

	initSignalSM(); //defined in io.c - inits signalSM value

	initFiles(); //defined in io.c - inits prep_g and trap_g files. I think these files are for tracing

	c = seqcheck( seq ); //seqcheck is defined in mltaln9.c. It checks 'seq' characters for unusual characters
	if( c )
	{
		fprintf( stderr, "Illegal character %c\n", c );
		exit( 1 );
	}

	fprintf( stderr, "\n" );
	//I need to know what is 'G' ?
	if( alg == 'G' ) // dp to the first sequence
	{
		mseq1f = AllocateCharMtx( 1, nlenmax+nlenmax );
		mseq1r = AllocateCharMtx( 1, nlenmax+nlenmax );
		mseq2 = AllocateCharMtx( 1, nlenmax+nlenmax );
	    tmpseq = AllocateCharVec( MAX( nlenmax, B ) +1 ); //B is defined in mltaln.h and = 256

		gappick0( mseq1f[0], seq[0] ); //defined in mltaln9.c. copies seq[0] chars to mseq1f[0] without gaps
		sreverse( mseq1r[0], mseq1f[0] ); //defined in io.c. fills mseq1r[0] with reversed chars from mseq1f[0]
		strcpy( seq[0], mseq1f[0] ); //copy seq[0] without gaps to seq[0]. i.e. remove gaps from seq[0]

		if( nadd )
			istart = njob - nadd;
		else
			istart = 1;

		fprintf( stderr, "\n" );

		for( i=0; i<istart; i++ )
		{
			gappick0( tmpseq, seq[i] ); //copies seq[i] chars to tmpseq without gaps
			strcpy( seq[i], tmpseq ); //copy tmpseq to seq[i]
			strcpy( tmpseq, name[i] ); //copy seq name to tmpseq
			strcpy( name[i], "_F_" ); //fill name[i] with "_F_"
			strncpy( name[i]+3, tmpseq+1, 10 ); //then copy first 10 characters from tmpseq to name
			name[i][13] = 0; //and end 14-chars name with 0 character
		}
		//by end of for loop above, 'seq' are replaced by 'seq' without gaps and 'name' are replaced by short names

		for( i=istart; i<njob; i++ ) 
		{
			gappick0( mseq2[0], seq[i] ); //copy seq[i] to mseq2[0]

#if GLOBAL
			res_forward = G__align11_noalign( n_dis_consweight_multi, penalty, penalty_ex, mseq1f, mseq2, 0 );
			res_reverse = G__align11_noalign( n_dis_consweight_multi, penalty, penalty_ex, mseq1r, mseq2, 0 );
#else
			//this method returns max weight between mseq1f[0] and mseq2[0]
			res_forward = L__align11_noalign( n_dis_consweight_multi, mseq1f, mseq2 ); //defined in Lalign11.c
			//this method returns max weight between mseq1r[0] and mseq2[0]
			res_reverse = L__align11_noalign( n_dis_consweight_multi, mseq1r, mseq2 ); //n_dis_consweight_multi defined in defs.c. it is filled in in constants.c
#endif
#if 0

			strcpy( mseq2[0], seq[i] );
			strcpy( mseq1f[0], seq[0] );
			res_forward = G__align11( n_dis_consweight_multi, mseq1f, mseq2, nlenmax*2, 0, 0 );
			fprintf( stdout, "%s\n", mseq1f[0] );
			fprintf( stdout, "%s\n", mseq2[0] );

			strcpy( mseq2[0], seq[i] );
			sreverse( mseq1r[0], seq[0] );
			res_reverse = G__align11( n_dis_consweight_multi, mseq1r, mseq2, nlenmax*2, 0, 0 );
			fprintf( stdout, "%s\n", mseq1r[0] );
			fprintf( stdout, "%s\n", mseq2[0] );
#endif

//			fprintf( stdout, "\nscore_for(%d,%d) = %f\n", 0, i, res_forward );
//			fprintf( stdout, "score_rev(%d,%d) = %f\n", 0, i, res_reverse );
			res_max = MAX(res_reverse,res_forward); //MAX defined in fft.h
			if( (res_reverse-res_forward)/res_max > thresholdtorev ) // tekitou
			{
//				fprintf( stderr, "REVERSE!!!\n" );
				sreverse( seq[i], mseq2[0] ); //fill seq[i] with reversed chars from mseq2[0]
				//and set name to _R_shortname
				strcpy( tmpseq, name[i] );
				strcpy( name[i], "_R_" );
				strncpy( name[i]+3, tmpseq+1, 10 );
				name[i][13] = 0;
			}
			else
			{
				strcpy( seq[i], mseq2[0] ); //fill seq[i] with chars from mseq2[0]

				strcpy( tmpseq, name[i] );
				strcpy( name[i], "_F_" );
				strncpy( name[i]+3, tmpseq+1, 10 );
				name[i][13] = 0;
			}
		}
		FreeCharMtx( mseq1f );
		FreeCharMtx( mseq1r );
		FreeCharMtx( mseq2 );
		free( tmpseq );
		//by end of this 'G' algorithm, first sequence of all sequences read from input file is matched with all other sequences.
		//each sequence is replaced either by the original sequence or the reversed one - based on best weight -.
	}
	else if( alg == 'm' ) //what is 'm' algorithm ?
	{

		if( dodp ) // nakuserukamo  //moni: dodp defined here and default val = 0, else = 1
		{
			mseq1f = AllocateCharMtx( 1, nlenmax+1);
			mseq1r = AllocateCharMtx( 1, nlenmax+1 );
			mseq2 = AllocateCharMtx( 1, nlenmax+1 );
		}
		else
		{
//			nthread = 0; // heiretsu keisan no kouritsu ha warui node
			spointt = AllocateIntMtx( njob, 0 ); 
			pointt = AllocateIntMtx( njob, nlenmax+1 );
			pointt_rev = AllocateIntMtx( njob, nlenmax+1 );
		}
	    tmpseq = AllocateCharVec( MAX( nlenmax, B ) +1 );
	    revseq = AllocateCharVec( nlenmax+1 );
		grpseq = AllocateIntVec( nlenmax+1 );
		res = AllocateIntVec( njob );
		resr = AllocateIntVec( njob );
		resf = AllocateIntVec( njob );
		map = AllocateIntVec( njob );
		contrastorder = AllocateIntVec( njob );
		if( dorp == 'd' ) tsize = (int)pow( 4, 6 ); //tsize defined in both defs.h and mltaln.h
		else              tsize = (int)pow( 6, 6 ); // iranai

		maxl = 0; //maxl defined in both defs.h and mltaln.h
		for( i=0; i<njob; i++ )
		{
			gappick0( tmpseq, seq[i] ); //copy seq[i] chars to tmpseq without gaps
			nogaplen[i] = strlen( tmpseq ); //get length of tmpseq, i.e. seq[i] without gaps
			if( nogaplen[i] > maxl ) maxl = nogaplen[i]; //set maxl to max seq length without gaps
		}

		reporterr( "Step 1/2\n" );

		if( !dodp ) //if dodp is not set
		{
			if( nadd )
				iend = njob - nadd;
			else
				iend = 0; // keisan shinai
	
			for( i=0; i<iend; i++ )
			{
				gappick0( tmpseq, seq[i] ); //copy seq[i] chars to tmpseq without gaps
				strcpy( seq[i], tmpseq ); //copy tmpseq to seq[i]
				make_direction_list_seq_grp_nuc( grpseq, tmpseq );
				make_direction_list_makepointtable_nuc( pointt[i], grpseq ); //set pointt[i] to some value based on calculations for values in grpseq
				spointt[i] = pointt[i];
			}
	
			if( nadd )
				istart = njob - nadd;
			else
				istart = 0;
			for( i=istart; i<njob; i++ ) 
			{

				gappick0( tmpseq, seq[i] ); //copy seq[i] chars to tmpseq without gaps
				strcpy( seq[i], tmpseq ); //copy tmpseq to seq[i]
				sreverse( revseq, tmpseq ); //fill revseq with reversed chars from tmpseq

				make_direction_list_seq_grp_nuc( grpseq, tmpseq );
				make_direction_list_makepointtable_nuc( pointt[i], grpseq );
//				makecompositiontable_p( table1, pointt[i] ); -> moto no basho ni modosu
				make_direction_list_seq_grp_nuc( grpseq, revseq );
				make_direction_list_makepointtable_nuc( pointt_rev[i], grpseq );
//				makecompositiontable_p( table1_rev, pointt_rev[i] ); -> moto no basho ni modosu
				spointt[i] = pointt[i];


//				reporterr( "pointt[i] = %p\n", pointt[i] );
//				reporterr( "pointt[i][0] = %p\n", pointt[i][0] );

			}	
		}


		if( contrastsort ) // sukoshi chuui //default = 1, and if set from parameters = 0
		{


			if( nadd )
			{
				iend = njob-nadd;
				for( i=0; i<iend; i++ ) contrastorder[i] = i;
				istart = njob-nadd;
				iend = nadd;
			}
			else
			{
				istart = 0;
				iend = njob;
			}

			if( dodp )
				makecontrastorder( contrastorder+istart, seq+istart, iend, istart ); //orders something and save in contrastorder
			else
				makecontrastorder6mer( contrastorder+istart, pointt+istart, pointt_rev+istart, seq+istart, iend, istart ); //orders something and save in contrastorder
		}
		else
		{
			for( i=0; i<njob; i++ ) contrastorder[i] = i;
		}


//		reporterr( "contrastorder = \n" );
//		for( i=0; i<njob; i++ )
//			reporterr( "%d ", contrastorder[i] );
//		reporterr( "\n" );



		if( nadd )
			iend = njob - nadd;
		else
			iend = 1;
		for( i=0; i<iend; i++ )
		{
			ic = contrastorder[i];
//			fprintf( stdout, "%d, SKIP\n", i );
			gappick0( tmpseq, seq[ic] ); //copy seq[ic] chars to tmpseq without gaps
			strcpy( seq[ic], tmpseq ); //copy tmpseq to seq[ic]
//			if( !nadd ) strcpy( seq[i], tmpseq ); // seq ha tsukawanaikara ii.

#if 0 // -> makecontrastorder() no mae ni idou
			if( !dodp )
			{
				make_direction_list_seq_grp_nuc( grpseq, tmpseq );
				make_direction_list_makepointtable_nuc( pointt[ic], grpseq );
				spointt[ic] = pointt[ic];
			}
#endif

			strcpy( tmpseq, name[ic] ); //copy also seq name and trim it to 14 chars
			strcpy( name[ic], "_F_" );
			strncpy( name[ic]+3, tmpseq+1, 10 );
			name[ic][13] = 0;
		}

		reporterr( "\n\nStep 2/2\n" );

		if( nadd )
			istart = njob - nadd;
		else
			istart = 1;
		for( i=istart; i<njob; i++ ) 
		{
//			fprintf( stderr, "\r %d / %d ", i, njob );
			ic = contrastorder[i];
			gappick0( tmpseq, seq[ic] ); //copy seq[ic] chars to tmpseq without gaps
			strcpy( seq[ic], tmpseq ); //copy tmpseq to seq[ic]
			sreverse( revseq, tmpseq ); //fill revseq with reversed chars from tmpseq

#if 0 // -> makecontrastorder() no mae ni idou
			if( !dodp )
			{
				table1 = (short *)calloc( tsize, sizeof( short ) );
				if( !table1 ) ErrorExit( "Cannot allocate table1\n" );
				table1_rev = (short *)calloc( tsize, sizeof( short ) );
				if( !table1_rev ) ErrorExit( "Cannot allocate table1_rev\n" );
				make_direction_list_seq_grp_nuc( grpseq, tmpseq );
				make_direction_list_makepointtable_nuc( pointt[ic], grpseq );
				make_distance_makecompositiontable_p( table1, pointt[ic] );
				make_direction_list_seq_grp_nuc( grpseq, revseq );
				make_direction_list_makepointtable_nuc( pointt_rev[ic], grpseq );
				make_distance_makecompositiontable_p( table1_rev, pointt_rev[ic] );
			}
#else
			if( !dodp )
			{
				table1 = (short *)calloc( tsize, sizeof( short ) );
				if( !table1 ) ErrorExit( "Cannot allocate table1\n" );
				table1_rev = (short *)calloc( tsize, sizeof( short ) );
				if( !table1_rev ) ErrorExit( "Cannot allocate table1_rev\n" );
				make_distance_makecompositiontable_p( table1, pointt[ic] ); //increment table1 values in pointt[ic] values
				make_distance_makecompositiontable_p( table1_rev, pointt_rev[ic] ); //increment table1_rev values in pointt_rev[ic] values
			}
#endif

			if( nadd && addfragment )
				iend = njob-nadd;
			else
				iend = i;


			if( iend > reflim ) 
			{
//				reporterr( "iend = %d -> %d\n", iend, reflim );
#if 0
				for( j=0; j<iend; j++ ) map[j] = j;
				partshuffle( iend, reflim, map );
#else
				for( j=0; j<iend; j++ ) map[j] = contrastorder[j];
#endif
				iend = reflim; // approximation
			}
			else
			{
#if 0
				for( j=0; j<iend; j++ ) map[j] = j;
#else
				for( j=0; j<iend; j++ ) map[j] = contrastorder[j];
#endif
			}

//			reporterr( "reflim = %d, seq[%d] = %s\n", reflim, contrastorder[0], seq[contrastorder[0]] );

#ifdef enablemultithread
			if( nthread )
			{
				pthread_t *handle;
				pthread_mutex_t mutex_counter;
				thread_arg_t *targ;
				int *jsharept;
		
				targ = calloc( nthread, sizeof( thread_arg_t ) );
				handle = calloc( nthread, sizeof( pthread_t ) );
				pthread_mutex_init( &mutex_counter, NULL );
				jsharept = calloc( 1, sizeof(int) );
				*jsharept = 0;
		
				if( i%100==1 ) fprintf( stderr, " %d / %d (%d threads)   \r", i, njob, nthread );
				for( j=0; j<nthread; j++ )
				{
					targ[j].iend = iend;
					targ[j].map = map;
					targ[j].seq = seq;
					targ[j].tmpseq = tmpseq; 
					targ[j].res = resf; 
					targ[j].spointt = spointt; 
					targ[j].table1 = table1; 
					targ[j].jshare = jsharept;
					targ[j].iq = i; // iranai
					targ[j].mutex_counter = &mutex_counter;
					targ[j].thread_no = j;
					pthread_create( handle+j, NULL, directionthread, (void *)(targ+j) );
				}
				for( j=0; j<nthread; j++ ) pthread_join( handle[j], NULL );
				pthread_mutex_destroy( &mutex_counter );
				free( handle );
				free( targ );
				free( jsharept );
			}
			else
#endif
			{
				thread_arg_t *targ;

				if( i%100==1 ) fprintf( stderr, " %d / %d   \r", i, njob );
				targ = calloc( 1, sizeof( thread_arg_t ) );
				targ[0].iend = iend;
				targ[0].map = map;
				targ[0].seq = seq;
				targ[0].tmpseq = tmpseq; 
				targ[0].res = resf; 
				targ[0].spointt = spointt; 
				targ[0].table1 = table1; 
				targ[0].iq = i;  // iranai
				directionthread( targ );
				free( targ );
			}



#ifdef enablemultithread
			if( nthread )
			{
				pthread_t *handle;
				pthread_mutex_t mutex_counter;
				thread_arg_t *targ;
				int *jsharept;
		
				targ = calloc( nthread, sizeof( thread_arg_t ) );
				handle = calloc( nthread, sizeof( pthread_t ) );
				pthread_mutex_init( &mutex_counter, NULL );
				jsharept = calloc( 1, sizeof(int) );
				*jsharept = 0;
		
				for( j=0; j<nthread; j++ )
				{
					targ[j].iend = iend;
					targ[j].seq = seq;
					targ[j].map = map;
					targ[j].tmpseq = revseq; 
					targ[j].res = resr; 
					targ[j].spointt = spointt; 
					targ[j].table1 = table1_rev; 
					targ[j].jshare = jsharept;
					targ[j].iq = i; // iranai
					targ[j].mutex_counter = &mutex_counter;
					targ[j].thread_no = j;
					pthread_create( handle+j, NULL, directionthread, (void *)(targ+j) );
				}
				for( j=0; j<nthread; j++ ) pthread_join( handle[j], NULL );
				pthread_mutex_destroy( &mutex_counter );
				free( handle );
				free( targ );
				free( jsharept );
			}
			else
#endif
			{
				thread_arg_t *targ;
				targ = calloc( 1, sizeof( thread_arg_t ) );
				targ[0].iend = iend;
				targ[0].seq = seq;
				targ[0].map = map;
				targ[0].tmpseq = revseq; 
				targ[0].res = resr; 
				targ[0].spointt = spointt;
				targ[0].table1 = table1_rev; 
				targ[0].iq = i;  // iranai
				directionthread( targ );
				free( targ );
			}

			if( mode == '2' )
			{
				mres = mres2 = 0;
				for( j=0; j<iend; j++ )
				{
					ires = resf[j];
//					fprintf( stdout, "ires (%d,%d) = %d\n", i, j, ires );
//					fflush( stdout );
					if( ires>mres2 ) 
					{
						if( ires>mres ) 
						{
							mres2 = mres;
							mres = ires;
						}
						else
							mres2 = ires;
					}
				}
				res_forward = (double)( mres + mres2 ) / 2;
				mres = mres2 = 0;
				for( j=0; j<iend; j++ )
				{
					ires = resr[j];
					if( ires>mres2 )
					{
						if( ires>mres ) 
						{
							mres2 = mres;
							mres = ires;
						}
						else
							mres2 = ires;
					}
				}
				res_reverse = (double)( mres + mres2 ) / 2;
				res_max = MAX(res_reverse,res_forward);
			}
//			reporterr( "i=%d, res_reverse = %f\n", i, res_reverse );
			else if( mode == '1' )
			{
				res_reverse = 0.0;
				for( j=0; j<iend; j++ ) if( res_reverse < (double)resr[j] ) res_reverse = (double)resr[j];
				res_forward = 0.0;
				for( j=0; j<iend; j++ ) if( res_forward < (double)resf[j] ) res_forward = (double)resf[j];
				res_max = 1.0;
			}

			else if( mode == 'd' )
			{
				res_reverse = 0.0;
				for( j=0; j<iend; j++ ) if( res_reverse < (double)(resr[j]-resf[j]) ) res_reverse = (double)(resr[j]-resf[j]);
				res_forward = 0.0;
				for( j=0; j<iend; j++ ) if( res_forward < (double)(resf[j]-resr[j]) ) res_forward = (double)(resf[j]-resr[j]);
				res_max = 1.0;
			}

			else if( mode == 'a' )
			{
				res_reverse = 0.0;
				for( j=0; j<iend; j++ ) res_reverse += (double)resr[j];
				res_reverse /= (double)iend;
				res_forward = 0.0;
				for( j=0; j<iend; j++ ) res_forward += (double)resf[j];
				res_forward /= (double)iend;
				res_max = 1.0;
			}
			else
			{
				reporterr( "Unknown mode!\n" );
				exit( 1 );
			}


			if( (res_reverse>res_forward) ) // tekitou
//			if( (res_reverse-res_forward)/res_max > thresholdtorev ) // tekitou
			{
				strcpy( seq[ic], revseq );

				strcpy( tmpseq, name[ic] );
				strcpy( name[ic], "_R_" );
				strncpy( name[ic]+3, tmpseq+1, 10 );
				name[ic][13] = 0;
				if( !dodp ) spointt[ic] = pointt_rev[ic];
			}
			else
			{
				strcpy( tmpseq, name[ic] );
				strcpy( name[ic], "_F_" );
				strncpy( name[ic]+3, tmpseq+1, 10 );
				name[ic][13] = 0;
				if( !dodp ) spointt[ic] = pointt[ic];
			}

			if( !dodp )
			{
				free( table1 );
				free( table1_rev );
			}
		}

		if( name[0][1] == 'R' )
		{
			for( j=0; j<njob; j++ ) 
			{
				if( name[j][1] == 'R' ) 
					name[j][1] = 'F';
				else
					name[j][1] = 'R';
			}
		}

		creverse( 0 ); //free table in io.c
		free( tmpseq );
		free( revseq );
		free( grpseq );
		free( res );
		free( resr );
		free( resf );
		free( map );
		free( nlen );
		free( nogaplen );
		free( contrastorder );
		if( dodp )
		{
			FreeCharMtx( mseq1f );
			FreeCharMtx( mseq1r );
			FreeCharMtx( mseq2 );
		}
		else
		{
			FreeIntMtx( pointt );
			FreeIntMtx( pointt_rev );
			free( spointt );
		}
	}
	else
	{
		fprintf( stderr, "Unknown alg %c\n", alg );
		exit( 1 );
	}
//	writeData_pointer( stdout, njob, name, nlen, seq );
	for( i=0; i<njob; i++ ) 
	{
//		fprintf( stdout, ">%s\n", name[i] );
//		fprintf( stdout, "%s\n", seq[i] );
		fprintf( stdout, "%s\n", name[i] );
	}

	FreeCharMtx( seq );
	FreeCharMtx( name );
	freeconstants(); //free memory allocated in constants.c
	closeFiles(); //close prep_g and trap_g files in io.c

	fprintf( stderr, "\n" );
	SHOWVERSION; //i think this file finds score of sequences and make them forward or reverse
	return( 0 );
}

