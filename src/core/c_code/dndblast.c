//I think this file uses BLAST to make alignment
#include "mltaln.h"
#include <sys/types.h>
#include <unistd.h>
#define DEBUG 0
#define TEST 0


int dndblast_howmanyx( char *s )
{
	int val = 0;
	if( scoremtx == -1 )
	{
		do
		{
			if( !strchr( "atgcuATGCU", *s ) ) val++;
		} while( *++s );
	}
	else
	{
		do
		{
			if( !strchr( "ARNDCQEGHILKMFPSTWYV", *s ) ) val++;
		} while( *++s );
	}
	return( val );
}

void dndblastArguments( int argc, char *argv[] )
{
	int c;

	inputfile = NULL; //defined in defs.h
	disopt = 0; //defined in defs.h
	divpairscore = 0; //defined in defs.h

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
					disopt = 1;
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
        fprintf( stderr, "options: -i\n" );
        exit( 1 );
    }
}

int dndblast_main( int argc, char *argv[] )
{
	int ktuple;
	int i, j;
	FILE *infp;
	FILE *hat2p;
	FILE *hat3p;
	char **seq = NULL; // by D.Mathog
	char **seq1;
	static char **name;
	static char **name1;
	static int nlen1[M]; //M is constant in mltaln.h and = 500,000
	double **mtx;
	double **mtx2;
	static int nlen[M];
	char b[B]; //B is constant in mltaln.h and = 256
	double max;
	char com[1000];
	int opt[M];
	int res;
	char *home;
	char queryfile[B];
	char datafile[B];
	char fastafile[B];
	char hat2file[B];
	int pid = (int)getpid(); //returns the process ID of the current process. (This is often used by routines that generate unique temporary filenames.)
	LocalHom **localhomtable, *tmpptr; //LocalHom is a structure defined in mltaln.
#if 1
	home = getenv( "HOME" );
#else /* $HOME wo tsukau to fasta ni watasu hikisuu ga afureru */ 
	home = NULL;
#endif

#if DEBUG
	if( home ) fprintf( stderr, "home = %s\n", home );
#endif
	if( !home ) home = "";
	sprintf( queryfile, "%s/tmp/query-%d", home, pid ); //sprintf sends formatted output to first parameter string
	sprintf( datafile, "%s/tmp/data-%d", home, pid );
	sprintf( fastafile, "%s/tmp/fasta-%d", home, pid );
	sprintf( hat2file, "hat2-%d", pid );


	dndblastArguments( argc, argv );

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
#if 0
	PreRead( infp, &njob, &nlenmax );
#else
	dorp = NOTSPECIFIED;
	getnumlen( infp ); //defined in io.c. Finds sequences count, max length and dna or protein from input file
#endif

	if( dorp == 'd' )
	{
		scoremtx = -1; //defined in defs.h
		pamN = NOTSPECIFIED; //defined in defs.h
	}
	else
	{
		nblosum = 62; //defined in defs.h
		scoremtx = 1; //defined in defs.h
	}
	constants( njob, seq ); //defined in constants.c
	//after all this method, n_dis, ribosumdis, amino_dis, amino_dis_consweight_multi, n_dis_consweight_multi,
	//n_disLN, n_disFFT, polarity, volume arrays are initialized and some constants are set.

	rewind( infp ); //set the file position to the beginning of the infp

	name = AllocateCharMtx( njob, B+1 );
	name1 = AllocateCharMtx( njob, B+1 );
	seq = AllocateCharMtx( njob, nlenmax+1 );
	seq1 = AllocateCharMtx( 2, nlenmax+1 );
	mtx = AllocateDoubleMtx( njob, njob );
	mtx2 = AllocateDoubleMtx( njob, njob );
	localhomtable = (LocalHom **)calloc( njob, sizeof( LocalHom *) );
	for( i=0; i<njob; i++) //initialize localhomtable
	{
		localhomtable[i] = (LocalHom *)calloc( njob, sizeof( LocalHom ) );
		for( j=0; j<njob; j++) 
		{
			localhomtable[i][j].start1 = -1;
			localhomtable[i][j].end1 = -1;
			localhomtable[i][j].start2 = -1;
			localhomtable[i][j].end2 = -1;
			localhomtable[i][j].opt = -1.0; 
			localhomtable[i][j].next = NULL; 

		}
    }

#if 0
	FRead( infp, name, nlen, seq );
#else
	readData_pointer( infp, name, nlen, seq ); //defined in io.c. It reads sequences and their names in seq, name and nlen arrays.
#endif
	fclose( infp );
	
	if( scoremtx == -1 ) ktuple = 6; //if DNA
	else                 ktuple = 1; //if protein

	for( i=0; i<njob; i++ )
	{
		gappick0( seq1[0], seq[i] ); //copy 'seq[i]' chars to 'seq1[0]' without gaps chars
		strcpy( seq[i], seq1[0] ); //copy seq1[0] to seq[i]
	}
	for( j=0; j<njob; j++ )
	{
		sprintf( name1[j], "+==========+%d                      ", j );
		nlen1[j] = nlen[j];
	}

	for( i=0; i<njob; i++ ) //here
	{
//		fprintf( stderr, "###  i = %d\n", i );

		if( i % 10 == 0 )
		{
			fprintf( stderr, "formatting .. " );
			hat2p = fopen( datafile, "w" ); //open data file
			if( !hat2p ) ErrorExit( "Cannot open datafile." ); //ErrorExit defined in io.c. Print error message to stderr then exit.
			WriteForFasta( hat2p, njob-i, name1+i, nlen1+i, seq+i ); //defined in io.c. Write njob-i sequences - names and chars - to data file
			fclose( hat2p );
		
			if( scoremtx == -1 ) //DNA
				sprintf( com, "formatdb  -p f -i %s -o F", datafile );
			else //protein
				sprintf( com, "formatdb  -i %s -o F", datafile );
			system( com ); //execute 'com' command.
			fprintf( stderr, "done.\n" );
		}

		hat2p = fopen( queryfile, "w" ); //open query file
		if( !hat2p ) ErrorExit( "Cannot open queryfile." );
		WriteForFasta( hat2p, 1, name1+i, nlen+i, seq+i ); //defined in io.c. Write one sequence - name and char - to query file
		fclose( hat2p );


		if( scoremtx == -1 ) //DNA
			sprintf( com, "blastall -e 1e10 -p blastn -m 7  -i %s -d %s >  %s", queryfile, datafile, fastafile );
		else //protein
			sprintf( com, "blastall -G 10 -E 1 -e 1e10 -p blastp -m 7  -i %s -d %s >  %s", queryfile, datafile, fastafile );
		res = system( com ); //execute 'com' command, which I think executes blast command.
		if( res ) ErrorExit( "error in fasta" );


		hat2p = fopen( fastafile, "r" ); //open fasta file
		if( hat2p == NULL ) 
			ErrorExit( "file 'fasta.$$' does not exist\n" );
		res = ReadBlastm7( hat2p, mtx[i], i, name1, localhomtable[i] ); //defined in io.c.
		//this method reads data from fasta file and get some values from it, then fill in localhomtable with values based on read data.
		fclose( hat2p );

#if 0
		for( j=0; j<njob; j++ )
		{
			for( tmpptr=localhomtable[i]+j; tmpptr; tmpptr=tmpptr->next )
			{
				if( tmpptr->opt == -1.0 ) continue;
//				fprintf( stderr, "%d %d %d %6.3f %d %d %d %d %p\n", i, j, tmpptr->overlapaa, tmpptr->opt, tmpptr->start1, tmpptr->end1, tmpptr->start2, tmpptr->end2, tmpptr->next );
			}
		}
#endif

		if( res < njob-i+i%10 ) //there is an error - missing sequences -
		{
			fprintf( stderr, "WARNING: count (blast) = %d < %d\n", res, njob-i+i%10 );
		}


#if 0
		{
			int ii, jj;
			if( i < njob-1 ) for( jj=i; jj<i+5; jj++ ) 
				fprintf( stdout, "mtx[%d][%d] = %f\n", i+1, jj+1, mtx[i][jj] );
		}
#endif
		fprintf( stderr, "query : %4d / %d\n", i+1, njob );
	}

#if 1
	fprintf( stderr, "##### writing hat3\n" );
	hat3p = fopen( "hat3", "w" );
	if( !hat3p ) ErrorExit( "Cannot open hat3." );
	for( i=0; i<njob; i++ ) for( j=0; j<njob; j++ )
	{
//		fprintf( stderr, "mtx[%d][%d] = %f, mtx[%d][%d] = %f\n", i, j, mtx[i][j], j, i, mtx[j][i] );
		if( i == j ) continue;
		if( mtx[i][j] == mtx[j][i] && i > j ) continue;
		if( mtx[j][i] > mtx[i][j] ) continue;
		for( tmpptr=localhomtable[i]+j; tmpptr; tmpptr=tmpptr->next )
		{
			if( tmpptr->opt == -1.0 ) continue; //go to next iteration
			//write localhomtable data to hat3 file
			fprintf( hat3p, "%d %d %d %6.3f %d %d %d %d %p\n", i, j, tmpptr->overlapaa, tmpptr->opt, tmpptr->start1, tmpptr->end1, tmpptr->start2, tmpptr->end2, (void *)tmpptr->next );
		}
	}
	fclose( hat3p );
#endif

	for( i=0; i<njob; i++ )
	{
//		fprintf( stderr, "###  i = %d\n", i );
		hat2p = fopen( datafile, "w" ); //open datafile for write
		if( !hat2p ) ErrorExit( "Cannot open datafile." );
		WriteForFasta( hat2p, njob-i, name1+i, nlen1+i, seq+i ); //defined in io.c. Write njob-i sequences - name and char - to data file
		fclose( hat2p );

//		seq1[0] = seq[i];
//		nlen1[0] = nlen[i];

		hat2p = fopen( queryfile, "w" );
		if( !hat2p ) ErrorExit( "Cannot open queryfile." );
		WriteForFasta( hat2p, 1, name1+i, nlen+i, seq+i ); //defined in io.c. Write one sequence - name and char - to data file
		fclose( hat2p );


		if( scoremtx == -1 )
			sprintf( com, "fasta34 -z3 -m10  -n -Q  -b%d -E%d -d%d %s %s %d > %s", M, M, 0, queryfile, datafile, ktuple, fastafile );
		else
			sprintf( com, "fasta34 -z3 -m10  -Q  -b%d -E%d -d%d %s %s %d > %s", M, M, 0, queryfile, datafile, ktuple, fastafile );
		res = system( com ); //execute fasta34 command on query file, data file, ktuple and fasta file
		if( res ) ErrorExit( "error in fasta" );


		hat2p = fopen( fastafile, "r" );
		if( hat2p == NULL ) 
			ErrorExit( "file 'fasta.$$' does not exist\n" );
		res = ReadFasta34noalign( hat2p, mtx[i], i, name1, localhomtable[i] ); //defined in io.c.
		//fill mtx[i] with values read from fasta file and returns count of sequences - i think - from the file
		fclose( hat2p );
		if( res < njob - i )
		{
			fprintf( stderr, "count (fasta34 -z 3) = %d\n", res );
			exit( 1 );
		}


		if( i == 0 )
			for( j=0; j<njob; j++ ) opt[j] = (int)mtx[0][j];


#if 0
		{
			int ii, jj;
			if( i < njob-1 ) for( jj=i; jj<i+5; jj++ ) 
				fprintf( stdout, "mtx[%d][%d] = %f\n", i+1, jj+1, mtx[i][jj] );
		}
#endif
		fprintf( stderr, "query : %4d\r", i+1 );
	}




	for( i=0; i<njob; i++ )
	{
		max = mtx[i][i];
		if( max == 0.0 )
		{
			for( j=0; j<njob; j++ )
				mtx2[i][j] = 2.0;
		}
		else
		{
			for( j=0; j<njob; j++ )
			{
				mtx2[i][j] = ( max - mtx[MIN(i,j)][MAX(i,j)] ) / max * 2.0;
//				fprintf( stdout, "max = %f, mtx[%d][%d] = %f -> %f\n", max, i+1, j+1, mtx[i][j], mtx2[i][j] );
			}
		}
	}
	for( i=0; i<njob-1; i++ ) for( j=i+1; j<njob; j++ )
	{
//		fprintf( stdout, "mtx2[%d][%d] = %f, %f\n", i+1, j+1, mtx2[i][j], mtx2[j][i] );
		mtx2[i][j] = MIN( mtx2[i][j], mtx2[j][i] );
	}

#if 0
	{
		int ii, jj;
		if( i < njob-1 ) for( jj=i+1; jj<njob; jj++ ) 
			fprintf( stderr, "mtx2[][] = %f\n", mtx2[i][jj] );
	}
#endif

	for( i=0; i<njob; i++ ) name[i][0] = '=';

	if( disopt ) //update name array content with some additional info
	{
		strcpy( b, name[0] );
		sprintf( name[0], "=query====lgth=%04d-%04d %.*s", nlen[0], dndblast_howmanyx( seq[0] ), B-30, b );
#if 0
		strins(  b, name[0] );
#endif
		for( i=1; i<njob; i++ ) 
		{	
			strcpy( b, name[i] );
			sprintf( name[i], "=opt=%04d=lgth=%04d-%04d %.*s", opt[i], nlen[i], dndblast_howmanyx( seq[i] ), B-30, b );
#if 0
			strins( b, name[i] );
#endif
		}
	}

	hat2p = fopen( hat2file, "w" ); //open hat2file for write
	if( !hat2p ) ErrorExit( "Cannot open hat2." );
	WriteHat2_pointer( hat2p, njob, name, mtx2 ); //defined in io.c. write number of sequences, max distance value and mtx values to hat2 file
	fclose( hat2p );


	sprintf( com, "/bin/rm %s %s %s", queryfile, datafile, fastafile ); //remove temporary files
	system( com );

#if 0
	sprintf( com, ALNDIR "/supgsdl < %s", hat2file );
	res = system( com );
	if( res ) ErrorExit( "error in spgsdl" );
#endif

	sprintf( com, "mv %s hat2", hat2file ); //copy hat2 file
	res = system( com );
	if( res ) ErrorExit( "error in mv" );

	//so this file reads input file sequences, then apply blast and fasta programs on them and calculate distance between each sequence and the others.
	//then write this data in hat2 and hat3 files

	SHOWVERSION;
	exit( 0 );
}
