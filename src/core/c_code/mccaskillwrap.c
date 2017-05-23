//what is mccaskill ?!!!
#include "mltaln.h"

#define DEBUG 0

static char *whereismccaskillmea;

#ifdef enablemultithread
typedef struct _thread_arg
{
	int thread_no;
	int njob;
	int *jobpospt;
	int **gapmap;
	char **nogap;
	int nlenmax;
	RNApair ***pairprob;
	pthread_mutex_t *mutex;
} thread_arg_t;
#endif

void outmccaskill( FILE *fp, RNApair **pairprob, int length )
{
	int i;
	RNApair *pt;
	for( i=0; i<length; i++ ) for( pt=pairprob[i]; pt->bestpos!=-1; pt++ )
	{
		if( pt->bestpos > i ) 
			fprintf( fp, "%d %d %50.40f\n", i, pt->bestpos, pt->bestscore );
	}
}

#if 1
static void readrawmccaskill( FILE *fp, RNApair **pairprob, int length )
{
	char gett[1000];
	int *pairnum;
	int i;
	int left, right;
	double prob;

	pairnum = (int *)calloc( length, sizeof( int ) );
	for( i=0; i<length; i++ ) pairnum[i] = 0; //length = max length of sequences

	while( 1 )
	{
		fgets( gett, 999, fp ); //read line from fp to gett with max 999 chars
		if( feof( fp ) ) break;
		if( gett[0] == '>' ) continue;
		sscanf( gett, "%d %d %lf", &left, &right, &prob ); //read formatted input from gett string
		if( prob < 0.01 ) continue; // mxscarna to mafft ryoho ni eikyou
//fprintf( stderr, "gett = %s\n", gett );

		if( left != right && prob > 0.0 )
		{
			//realloc resizes the memory block pointed to by pairprob[left] that was previously allocated
			pairprob[left] = (RNApair *)realloc( pairprob[left], (pairnum[left]+2) * sizeof( RNApair ) );
			pairprob[left][pairnum[left]].bestscore = prob;
			pairprob[left][pairnum[left]].bestpos = right;
			pairnum[left]++;
			pairprob[left][pairnum[left]].bestscore = -1.0;
			pairprob[left][pairnum[left]].bestpos = -1;
//			fprintf( stderr, "%d-%d, %f\n", left, right, prob );

			pairprob[right] = (RNApair *)realloc( pairprob[right], (pairnum[right]+2) * sizeof( RNApair ) );
			pairprob[right][pairnum[right]].bestscore = prob;
			pairprob[right][pairnum[right]].bestpos = left;
			pairnum[right]++;
			pairprob[right][pairnum[right]].bestscore = -1.0;
			pairprob[right][pairnum[right]].bestpos = -1;
//			fprintf( stderr, "%d-%d, %f\n", right, left, prob );
		}
	}
	free( pairnum );
}
#endif

#ifdef enablemultithread
static void *athread( void *arg )
{
	thread_arg_t *targ = (thread_arg_t *)arg;
	int thread_no = targ->thread_no;
	int njob = targ->njob;
	int *jobpospt = targ->jobpospt;
	int **gapmap = targ->gapmap;
	char **nogap = targ->nogap;
	int nlenmax = targ->nlenmax;
	RNApair ***pairprob = targ->pairprob;

	int i, res;
	FILE *infp;
	char *com;
	char *dirname;

	dirname = calloc( 100, sizeof( char ) );
	com = calloc( 1000, sizeof( char ) );
	

	while( 1 )
	{
		pthread_mutex_lock( targ->mutex );
		i = *jobpospt;
		if( i == njob )
		{
			pthread_mutex_unlock( targ->mutex );
//			return( NULL );
			break;
		}
		*jobpospt = i+1;
		pthread_mutex_unlock( targ->mutex );

		commongappick_record( 1, nogap+i, gapmap[i] );
		if( strlen( nogap[i] ) == 0 ) continue;

		sprintf( dirname, "_%d", i );
		sprintf( com, "rm -rf %s", dirname );
		system( com );
		sprintf( com, "mkdir %s", dirname );
		system( com );

		fprintf( stderr, "%d / %d (by thread %4d)\n", i+1, njob, thread_no );
		sprintf( com, "%s/_mccaskillinorg", dirname );
		infp = fopen( com, "w" );
//		fprintf( infp, ">in\n%s\n", nogap[i] );
		fprintf( infp, ">in\n" );
		write1seq( infp, nogap[i] );
		fclose( infp );

		sprintf( com, "tr -d '\\r' < %s/_mccaskillinorg > %s/_mccaskillin", dirname, dirname );
		system( com ); // for cygwin, wakaran
		if( alg == 'G' )
			sprintf( com, "cd %s; %s/dafs --mafft-out _mccaskillout _mccaskillin > _dum1 2>_dum", dirname, whereismccaskillmea );
		else
			sprintf( com, "cd %s; %s/mxscarnamod -m -writebpp  _mccaskillin > _mccaskillout 2>_dum", dirname, whereismccaskillmea );
		res = system( com );

		if( res )
		{
			fprintf( stderr, "ERROR IN mccaskill_mea\n" );
			exit( 1 );
		}

		sprintf( com, "%s/_mccaskillout", dirname );
		infp = fopen( com, "r" );
		readrawmccaskill( infp, pairprob[i], nlenmax );
		fclose( infp );

		sprintf( com, "rm -rf %s > /dev/null 2>&1", dirname );
		if( system( com ) )
		{
			fprintf( stderr, "retrying to rmdir\n" );
//			nanosleep( 100000 );
			sleep( 1 );
			system( com );
		}
	}
	free( dirname );
	free( com );
	return( NULL );
}
#endif

void kill_wrap_arguments( int argc, char *argv[] )
{
    int c;
	nthread = 1; //defined in defs.c
	inputfile = NULL; //defined in defs.h
	dorp = NOTSPECIFIED; //defined in defs.c
	kimuraR = NOTSPECIFIED; //defined in defs.h
	pamN = NOTSPECIFIED; //defined in defs.h
	whereismccaskillmea = NULL;
	alg = 's'; //defined in defs.h

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
				case 'd':
					whereismccaskillmea = *++argv;
					fprintf( stderr, "whereismccaskillmea = %s\n", whereismccaskillmea );
					--argc;
					goto nextoption;
				case 'C':
					nthread = myatoi( *++argv );
					fprintf( stderr, "nthread = %d\n", nthread );
					--argc; 
					goto nextoption;
				case 's':
					alg = 's'; // use scarna; default
					break;
				case 'G':
					alg = 'G'; // use dafs, instead of scarna
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
    if( argc != 0 ) 
    {
        fprintf( stderr, "options: Check source file !\n" );
        exit( 1 );
    }
}


int mccaskillwrap_main( int argc, char *argv[] )
{
	static char com[10000];
	static int  *nlen;	
	int left, right;
	int res;
	static char **name, **seq, **nogap;
	static int **gapmap;
	static int *order;
	int i, j;
	FILE *infp;
	RNApair ***pairprob; //RNApair is a structure defined in mltaln.h
	RNApair **alnpairprob;
	RNApair *pairprobpt;
	RNApair *pt;
	int *alnpairnum;
	double prob;
	int adpos;

	kill_wrap_arguments( argc, argv );
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

	if( !whereismccaskillmea )
		whereismccaskillmea = "";

	getnumlen( infp ); //defined in io.c. finds sequences count, max length and dna or protein from input file
	rewind( infp );

	if( dorp != 'd' ) //protein, so exit. this means that this file is executed for only nucleotides
	{
		fprintf( stderr, "nuc only\n" );
		exit( 1 );
	}

	seq = AllocateCharMtx( njob, nlenmax*2+1 );
	nogap = AllocateCharMtx( njob, nlenmax*2+1 );
	gapmap = AllocateIntMtx( njob, nlenmax*2+1 );
	order = AllocateIntVec( njob );
	name = AllocateCharMtx( njob, B+1 );
	nlen = AllocateIntVec( njob );
	pairprob = (RNApair ***)calloc( njob, sizeof( RNApair ** ) );
	alnpairprob = (RNApair **)calloc( nlenmax, sizeof( RNApair * ) );
	alnpairnum = AllocateIntVec( nlenmax );

	for( i=0; i<nlenmax; i++ ) alnpairnum[i] = 0;

	readData_pointer( infp, name, nlen, seq ); //defined in io.c. it reads sequences and their names in seq, name and nlen arrays.
	fclose( infp );

	//initialize pairprob, alnpairprob and nogap matrices
	for( i=0; i<njob; i++ )
	{
		pairprob[i] = (RNApair **)calloc( nlenmax, sizeof( RNApair * ) );
		for( j=0; j<nlenmax; j++ )
		{
			pairprob[i][j] = (RNApair *)calloc( 1, sizeof( RNApair ) );
			pairprob[i][j][0].bestpos = -1;
			pairprob[i][j][0].bestscore = -1.0;
		}
		strcpy( nogap[i], seq[i] );
		order[i] = i;
	}
	for( j=0; j<nlenmax; j++ )
	{
		alnpairprob[j] = (RNApair *)calloc( 1, sizeof( RNApair ) );
		alnpairprob[j][0].bestpos = -1;
		alnpairprob[j][0].bestscore = -1.0;
	}


	constants( njob, seq ); //defined in constants.c.
	//after all this method, n_dis, ribosumdis, amino_dis, amino_dis_consweight_multi, n_dis_consweight_multi,
	//n_disLN, n_disFFT, polarity, volume arrays are initialized and some constants are set.

	if( alg == 'G' )
		fprintf( stderr, "Running DAFS (Sato et al. 2012; http://www.ncrna.org/).\n" );
	else
		fprintf( stderr, "Running mxscarna with the mccaskill_mea mode.\n" );
#ifdef enablemultithread
	if( nthread > 0 )
	{
		int jobpos;
		pthread_t *handle;
		pthread_mutex_t mutex;
		thread_arg_t *targ;
		jobpos = 0;

		targ = calloc( nthread, sizeof( thread_arg_t ) );
		handle = calloc( nthread, sizeof( pthread_t ) );
		pthread_mutex_init( &mutex, NULL );

		for( i=0; i<nthread; i++ )
		{
			targ[i].thread_no = i;
			targ[i].njob = njob;
			targ[i].jobpospt = &jobpos;
			targ[i].gapmap = gapmap;
			targ[i].nogap = nogap;
			targ[i].nlenmax = nlenmax;
			targ[i].pairprob = pairprob;
			targ[i].mutex = &mutex;

//			athread( targ );
			pthread_create( handle+i, NULL, athread, (void *)(targ+i) );
			
		}

		for( i=0; i<nthread; i++ )
		{
			pthread_join( handle[i], NULL );
		}
		pthread_mutex_destroy( &mutex );

		free( handle );
		free( targ );


		for( i=0; i<njob; i++ )
		{
			fprintf( stdout, ">%d\n", i );
			outmccaskill( stdout, pairprob[i], nlenmax );
		}
	}
	else
#endif
	{
		for( i=0; i<njob; i++ )
		{
			fprintf( stderr, "%d / %d\n", i+1, njob );
			commongappick_record( 1, nogap+i, gapmap[i] ); //defined in mltaln9.c.
			//at end of this method, nogap[i] contains all chars in nogap filled in sequentially without gaps
			//and other remaining chars at the end. Also, gapmap contains indices of chars in nogap[i]
			if( strlen( nogap[i] ) == 0 ) 
			{
				fprintf( stdout, ">%d\n", i );
				continue; //jump directly to next loop iteration
			}

			infp = fopen( "_mccaskillinorg", "w" );
//			fprintf( infp, ">in\n%s\n", nogap[i] );
			fprintf( infp, ">in\n" );
			write1seq( infp, nogap[i] ); //defined in io.c. write nogap[i] in infp
			fclose( infp );
	
			//this function passes the command parameter to the host environment to be executed by the command-processor
			//and returns after the command has been executed.
			system( "tr -d '\\r' < _mccaskillinorg > _mccaskillin" ); // for cygwin, wakaran
			if( alg == 'G' )
				sprintf( com, "env PATH=%s dafs --mafft-out _mccaskillout _mccaskillin > _dum1 2>_dum", whereismccaskillmea );
			else
				sprintf( com, "env PATH=%s mxscarnamod -m -writebpp  _mccaskillin > _mccaskillout 2>_dum", whereismccaskillmea );
			res = system( com ); //execute the last given command in 'com'
	
			if( res ) //res == -1, then error
			{
				fprintf( stderr, "ERROR IN mccaskill_mea\n" );
				exit( 1 );
			}
	
			infp = fopen( "_mccaskillout", "r" );
			readrawmccaskill( infp, pairprob[i], nlenmax ); //defined here. fill pairprob[i] with values based on '_mccaskillout' file content
			fclose( infp );
			fprintf( stdout, ">%d\n", i );
			outmccaskill( stdout, pairprob[i], nlenmax ); //defined here. print pairprob[i] content to stdout
		}
	}

	for( i=0; i<njob; i++ )
	{
		for( j=0; j<nlen[i]; j++ ) for( pairprobpt=pairprob[i][j]; pairprobpt->bestpos!=-1; pairprobpt++ )
		{
			left = gapmap[i][j];
			right = gapmap[i][pairprobpt->bestpos];
			prob = pairprobpt->bestscore;

			for( pt=alnpairprob[left]; pt->bestpos!=-1; pt++ )
				if( pt->bestpos == right ) break;

			if( pt->bestpos == -1 )
			{
				alnpairprob[left] = (RNApair *)realloc( alnpairprob[left], (alnpairnum[left]+2) * sizeof( RNApair ) );
				adpos = alnpairnum[left];
				alnpairnum[left]++;
				alnpairprob[left][adpos].bestscore = 0.0;
				alnpairprob[left][adpos].bestpos = right;
				alnpairprob[left][adpos+1].bestscore = -1.0;
				alnpairprob[left][adpos+1].bestpos = -1;
				pt = alnpairprob[left]+adpos;
			}
			else
				adpos = pt-alnpairprob[left];

			pt->bestscore += prob;
			if( pt->bestpos != right )
			{
				fprintf( stderr, "okashii!\n" ); //okashii = funny
				exit( 1 );
			}
//			fprintf( stderr, "adding %d-%d, %f\n", left, right, prob );
		}
	}

	//free all allocated memory
	for( i=0; i<njob; i++ )
	{
		for( j=0; j<nlenmax; j++ ) free( pairprob[i][j] );
		free( pairprob[i] );
	}
	free( pairprob );
	for( j=0; j<nlenmax; j++ ) free( alnpairprob[j] );
	free( alnpairprob );
	free( alnpairnum );
	free( order );
	free( nlen );
	FreeCharMtx( seq );
	FreeCharMtx( nogap );
	FreeCharMtx( name );
	FreeIntMtx( gapmap );
	freeconstants();
	fprintf( stderr, "%d thread(s)\n", nthread );
	return( 0 );

	//so, after this file, I think it reads sequences from input file then performs some operations on them
	//based on data read from _mccaskillinorg and _mccaskillout files. All output are stored to output file

#if 0
	fprintf( stdout, "result=\n" );

	for( i=0; i<nlenmax; i++ ) for( pairprobpt=alnpairprob[i]; pairprobpt->bestpos!=-1; pairprobpt++ )
	{
		pairprobpt->bestscore /= (double)njob;
		left = i;
		right = pairprobpt->bestpos;
		prob = pairprobpt->bestscore;
		fprintf( stdout, "%d-%d, %f\n", left, right, prob );
	}

	return( 0 );
#endif
}
