#include "mltaln.h"

#define DEBUG 0 //i need to know what does this mean ?

static int seedoffset;

static void replace_unusual( int n, char **seq, char *usual, char unknown, int (*uporlow)( int ) )
{
	int i;
	char *pt;
	for( i=0; i<n; i++ )
	{
		pt = seq[i];
		while( *pt )
		{
			if( !strchr( usual, *pt ) ) *pt = unknown;
			else *pt = uporlow( *pt );
			pt++;
		}
	}
}


void replaceu_arguments( int argc, char *argv[] )
{
    int c;

	seedoffset = 0;
	inputfile = NULL;
	dorp = NOTSPECIFIED;

    while( --argc > 0 && (*++argv)[0] == '-' )
	{
        while ( (c = *++argv[0]) )
		{
            switch( c )
            {
				case 'o':
					seedoffset = myatoi( *++argv ); //what is myatoi? and where is it defined?
					fprintf( stderr, "seedoffset = %d\n", seedoffset );
					--argc;
					goto nextoption;
				case 'i':
					inputfile = *++argv;
					fprintf( stderr, "inputfile = %s\n", inputfile );
					--argc;
					goto nextoption;
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
    if( argc != 0 ) 
    {
        fprintf( stderr, "options: Check source file !\n" );
        exit( 1 );
    }
}


//finallyyyyyyy, the main method :D
int replaceu_main( int argc, char *argv[] )
{
	FILE *infp;
	int nlenmin;
	char **name;
	char **seq;
	int *nlen;
	int i;
	char *usual;

	replaceu_arguments( argc, argv );

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
		infp = stdin; //what is stdin?


//	dorp = NOTSPECIFIED;
	getnumlen_casepreserve( infp, &nlenmin );
    //what are njob, nlenmax and dorp? and where are they defined?
	fprintf( stderr, "%d x %d - %d %c\n", njob, nlenmax, nlenmin, dorp );

	seq = AllocateCharMtx( njob, nlenmax+1 ); //where are those methods and what they do ?
	name = AllocateCharMtx( njob, B+1 );
	nlen = AllocateIntVec( njob );

	readData_pointer_casepreserve( infp, name, nlen, seq );

//	for( i=0; i<njob; i++ ) gappick_samestring( seq[i] );

#if 0 //this code is not compiled
	FILE *origfp;
	origfp = fopen( "_original", "w" );
	if( !origfp )
	{
		fprintf( stderr, "Cannot open _original\n" );
		exit( 1 );
	}
	for( i=0; i<njob; i++ )
	{
		nlen[i] = strlen( seq[i] );
		fprintf( origfp, ">%s\n", name[i]+1 );
		if( seq[i][nlen[i]-1] == '\n' ) seq[i][nlen[i]-1] = 0;
		fprintf( origfp, "%s\n", seq[i] );
	}
	fclose( origfp );
#endif

	if( dorp == 'p' )
	{
		usual = "ARNDCQEGHILKMFPSTWYVarndcqeghilkmfpstwyv-.";
		replace_unusual( njob, seq, usual, 'X', toupper );
	}
	else
	{
		usual = "ATGCUatgcuBDHKMNRSVWYXbdhkmnrsvwyx-";
		replace_unusual( njob, seq, usual, 'n', tolower );
	}
	


	for( i=0; i<njob; i++ )
	{
		fprintf( stdout, ">_os_%d_oe_%s\n", i+seedoffset, name[i]+1 );
		fprintf( stdout, "%s\n", seq[i] );
	}

	free( nlen );
	FreeCharMtx( seq );
	FreeCharMtx( name );

	return( 0 );
}
