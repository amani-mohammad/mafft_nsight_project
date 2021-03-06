//Iterative refinement method

/* Tree-dependent-iteration */
 /* Devide to segments       */ 

#include "mltaln.h"

extern char **seq_g;
extern char **res_g;
static int subalignment;
static int subalignmentoffset;
static int specifictarget;

static int intop;
static int intree;
static double autosubalignment;

//calculate maxdistclass based on some calcs on ndistclass
static void calcmaxdistclass( void )
{
	int c;
	double rep;
	for( c=0; c<ndistclass; c++ ) //ndistclass defined in defs.c. and = 10
	{
		rep = (double) 2 * c / ndistclass; // dist:0-2 for dist2offset 
//		fprintf( stderr, "c=%d, rep=%f, offset=%f\n", c, rep, dist2offset( rep )  );
		if( dist2offset( rep ) == 0.0 )
			break;
	}
	fprintf( stderr, "ndistclass = %d, maxdistclass = %d\n", ndistclass, c+1 );
	maxdistclass = c + 1; //maxdistclass defined in defs.c.
//	maxdistclass = ndistclass; // CHUUI!!!!
	return;
}

void dvtditrArguments( int argc, char *argv[] )
{
	int c;
	char *argkey;

	outnumber = 0; //defined in defs.c.
	nthread = 1; //defined in defs.c.
	randomseed = 0; //defined in defs.c.
	scoreout = 0; //defined in defs.c.
	spscoreout = 0; //defined in defs.c.
	parallelizationstrategy = BAATARI1;  //defined in defs.c. BAATARI1 defined in mltaln.h.
	intop = 0; //defined here.
	intree = 0; //defined here.
	inputfile = NULL; //defined in defs.h.
	rnakozo = 0; //defined in defs.h.
	rnaprediction = 'm'; //defined in defs.h.
	nevermemsave = 0; //defined in defs.h.
	score_check = 1; //defined in defs.h.
	fftkeika = 1; //defined in defs.h.
	constraint = 0; //defined in defs.h.
	fmodel = 0; //defined in defs.h.
	kobetsubunkatsu = 1; //defined in defs.h.
	bunkatsu = 1; //defined in defs.h.
	nblosum = 62; //defined in defs.h.
	niter = 100; //defined in defs.h.
	calledByXced = 0; //defined in defs.h.
	devide = 1; //defined in defs.h.
	divWinSize = 20; /* 70 */    //defined in defs.h.
	divThreshold = 65; //defined in defs.h.
	fftscore = 1; //defined in defs.h.
	fftRepeatStop = 0; //defined in defs.h.
	fftNoAnchStop = 0; //defined in defs.h.
    scmtd = 5; //defined in defs.h.
	cooling = 1; //defined in defs.h.
    weight = 4; //defined in defs.h.
    utree = 1; //defined in defs.h.
    refine = 1; //defined in defs.h.
    check = 1; //defined in defs.h.
    cut = 0.0; //defined in defs.h.
	disp = 0; //defined in defs.h.
	outgap = 1; //defined in defs.c.
	use_fft = 0; // CHUUI dochira demo mafft.tmpl deha F    //defined in defs.h.
	force_fft = 0; //defined in defs.h.
	alg = 'A';  /* chuui */    //defined in defs.h.
	mix = 0; //defined in defs.h.
	checkC = 0; //defined in defs.h.
	tbitr = 0; //defined in defs.h.
	treemethod = 'X'; //defined in defs.h.
	sueff_global = 0.1; //defined in defs.c.
	scoremtx = 1; //defined in defs.h.
	dorp = NOTSPECIFIED; //defined in defs.c.
	ppenalty = NOTSPECIFIED; //defined in defs.h.
	penalty_shift_factor = 1000.0; //defined in defs.c.
	ppenalty_ex = NOTSPECIFIED; //defined in defs.h.
	poffset = NOTSPECIFIED; //defined in defs.h.
	kimuraR = NOTSPECIFIED; //defined in defs.h.
	pamN = NOTSPECIFIED; //defined in defs.h.
	geta2 = GETA2; //defined in defs.h. GETA2 defined in mltaln.h.
	fftWinSize = NOTSPECIFIED; //defined in defs.h.
	fftThreshold = NOTSPECIFIED; //defined in defs.h.
	RNAppenalty = NOTSPECIFIED; //defined in defs.h.
	RNAppenalty_ex = NOTSPECIFIED; //defined in defs.h.
	RNApthr = NOTSPECIFIED; //defined in defs.h.
	TMorJTT = JTT; //defined in defs.h. JTT defined in mltaln.h.
	consweight_multi = 1.0; //defined in defs.c.
	consweight_rna = 0.0; //defined in defs.c.
	subalignment = 0; //defined here.
	subalignmentoffset = 0; //defined here.
	legacygapcost = 0; //defined in defs.c.
	specificityconsideration = 0.0; //defined in defs.c.
	autosubalignment = 0.0; //defined here.
	specifictarget = 0; //defined here.
	nwildcard = 0; //defined in defs.c.

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
					niter = myatoi( *++argv );
					fprintf( stderr, "niter = %d\n", niter );
					--argc;
					goto nextoption;
				case 'e':
					RNApthr = (int)( atof( *++argv ) * 1000 - 0.5 );
					--argc;
					goto nextoption;
				case 'o':
					RNAppenalty = (int)( atof( *++argv ) * 1000 - 0.5 );
					--argc;
					goto nextoption;
				case 'f':
					ppenalty = (int)( atof( *++argv ) * 1000 - 0.5 );
//					fprintf( stderr, "ppenalty = %d\n", ppenalty );
					--argc;
					goto nextoption;
				case 'Q':
					penalty_shift_factor = atof( *++argv );
					if( penalty_shift_factor < 100.0 && penalty_shift_factor != 2.0 )
					{
						fprintf( stderr, "%f, penalty_shift is fixed to penalty x 2 in the iterative refinement phase.\n", penalty_shift_factor );
						penalty_shift_factor = 2.0;
					}
					--argc;
					goto nextoption;
				case 'g':
					ppenalty_ex = (int)( atof( *++argv ) * 1000 - 0.5 );
//					fprintf( stderr, "ppenalty_ex = %d\n", ppenalty_ex );
					--argc;
					goto nextoption;
				case 'h':
					poffset = (int)( atof( *++argv ) * 1000 - 0.5 );
					fprintf( stderr, "poffset = %d\n", poffset );
					--argc;
					goto nextoption;
				case 'k':
					kimuraR = myatoi( *++argv );
					fprintf( stderr, "kappa = %d\n", kimuraR );
					--argc; 
					goto nextoption;
				case 'b':
					nblosum = myatoi( *++argv );
					scoremtx = 1;
					fprintf( stderr, "blosum %d / kimura 200\n", nblosum );
					--argc; 
					goto nextoption;
				case 'j':
					pamN = myatoi( *++argv );
					scoremtx = 0;
					TMorJTT = JTT;
					fprintf( stderr, "jtt/kimura %d\n", pamN );
					--argc; 
					goto nextoption;
				case 'm':
					pamN = myatoi( *++argv );
					scoremtx = 0;
					TMorJTT = TM;
					fprintf( stderr, "tm %d\n", pamN );
					--argc; 
					goto nextoption;
				case 'l':
					fastathreshold = atof( *++argv );
					constraint = 2;
					--argc;
					goto nextoption;
				case 'r':
					consweight_rna = atof( *++argv );
					rnakozo = 1;
					--argc;
					goto nextoption;
				case 'c':
					consweight_multi = atof( *++argv );
					--argc;
					goto nextoption;
				case 'C':
					nthread = myatoi( *++argv );
					fprintf( stderr, "nthread = %d\n", nthread );
					--argc; 
					goto nextoption;
				case 'H':
					subalignment = 1;
					subalignmentoffset = myatoi( *++argv );
					--argc;
					goto nextoption;
				case 't':
					randomseed = myatoi( *++argv );
					fprintf( stderr, "randomseed = %d\n", randomseed );
					--argc; 
					goto nextoption;
				case 'p':
					argkey = *++argv;
					if( !strcmp( argkey, "BESTFIRST" ) ) parallelizationstrategy = BESTFIRST;
					else if( !strcmp( argkey, "BAATARI0" ) ) parallelizationstrategy = BAATARI0;
					else if( !strcmp( argkey, "BAATARI1" ) ) parallelizationstrategy = BAATARI1;
					else if( !strcmp( argkey, "BAATARI2" ) ) parallelizationstrategy = BAATARI2;
					else
					{
						fprintf( stderr, "Unknown parallelization strategy, %s\n", argkey );
						exit( 1 );
					}
//					exit( 1 );
					--argc; 
					goto nextoption;
				case 's':
					specificityconsideration = (double)myatof( *++argv );
//					fprintf( stderr, "specificityconsideration = %f\n", specificityconsideration );
					--argc; 
					goto nextoption;
#if 0
				case 'S' :
					scoreout = 1; // for checking parallel calculation
					break;
#else
				case 'S' :
					spscoreout = 1; // 2014/Dec/30, sp score
					break;
#endif
#if 0
				case 's' :
					RNAscoremtx = 'r';
					break;
#endif
#if 1
				case 'a':
					fmodel = 1;
					break;
#endif
				case 'N':
					nevermemsave = 1;
					break;
				case 'D':
					dorp = 'd';
					break;
				case 'P':
					dorp = 'p';
					break;
#if 0
				case 'Q':
					alg = 'Q';
					break;
#endif
				case 'R':
					rnaprediction = 'r';
					break;
				case 'O':
					fftNoAnchStop = 1;
					break;
#if 0
				case 'e':
					fftscore = 0;
					break;
				case 'r':
					fmodel = -1;
					break;
				case 'R':
					fftRepeatStop = 1;
					break;
#endif
				case 'T':
					kobetsubunkatsu = 0;
					break;
				case 'B':
					bunkatsu = 0;
					break;
#if 0
				case 'c':
					cooling = 1;
					break;
				case 'a':
					alg = 'a';
					break;
				case 's' :
					treemethod = 's';
					break;
				case 'H':
					alg = 'H';
					break;
#endif
				case 'A':
					alg = 'A';
					break;
				case 'M':
					alg = 'M';
					break;
				case '@':
					alg = 'd';
					break;
				case 'F':
					use_fft = 1;
					break;
#if 0
				case 't':
					weight = 4;
					break;
#endif
				case 'u':
					weight = 0;
					break;
				case 'U':
					intree = 1;
					break;
				case 'V':
					intop = 1;
					break;
				case 'J':
					utree = 0;
					break;
#if 0
				case 'd':
					disp = 1;
					break;
#endif
				case 'Z':
					score_check = 0;
					break;
				case 'Y':
					score_check = 2;
					break;
				case 'L':
					legacygapcost = 1;
					break;
#if 0
				case 'n' :
					treemethod = 'n';
					break;
#endif
				case 'n' :
					outnumber = 1;
					break;
				case 'X':
					treemethod = 'X';
					sueff_global = atof( *++argv );
					fprintf( stderr, "sueff_global = %f\n", sueff_global );
					--argc;
					goto nextoption;
#if 0
				case 'E' :
					treemethod = 'E';
					break;
				case 'q' :
					treemethod = 'q';
					break;
#endif
				case 'E':
					autosubalignment = atof( *++argv );
					fprintf( stderr, "autosubalignment = %f\n", autosubalignment );
					--argc;
					goto nextoption;
				case 'W':
					minimumweight = atof( *++argv );
					fprintf( stderr, "minimumweight = %f\n", minimumweight );
					--argc;
					goto nextoption;
				case 'z':
					fftThreshold = myatoi( *++argv );
					--argc;
					goto nextoption;
				case 'w':
					fftWinSize = myatoi( *++argv );
					--argc;
					goto nextoption;
				case '=':
					specifictarget = 1;
					break;
				case ':':
					nwildcard = 1;
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
		fprintf( stderr, "options : Check source file!\n" );
		exit( 1 );
	}
#if 0
	if( alg == 'A' && weight == 0 ) 
		ErrorExit( "ERROR : Algorithm A+ and un-weighted\n" ); 
#endif
}


int dvtditr_main( int argc, char *argv[] )
{
    int identity;
	static int nlen[M];
	static char **name, **seq, **aseq, **bseq;
	static Segment *segment = NULL; //structure defined in mltaln.h.
	static int anchors[MAXSEG];
	int i, j;
	int iseg, nseg;
	int ***topol;
	double **len;
	double **eff;
	FILE *prep;
	FILE *infp;
	FILE *orderfp;
	int alloclen;
	int returnvalue;
	char c;
	int ocut;
	char **seq_g_bk;
	LocalHom **localhomtable = NULL; // by D.Mathog   //structure defined in mltaln.h.
	RNApair ***singlerna; //structure defined in mltaln.h.
	int nogaplen;
	static char **nogap1seq;
	static char *kozoarivec;
	int nkozo;
	int alignmentlength;
	int **skipthisbranch;
	int foundthebranch;
	int nsubalignments, maxmem;
	int **subtable;
	int *insubtable;
	int *preservegaps;
	char ***subalnpt;
	int ntarget, *targetmap, *targetmapr;
	int ilim;

	dvtditrArguments( argc, argv ); //parse arguments
#ifndef enablemultithread
	nthread = 0;
#endif
	if( fastathreshold < 0.0001 ) constraint = 0;

	if( inputfile )
	{
		infp = fopen( inputfile, "r" ); //open inputfile for reading
		if( !infp ) 
		{
			fprintf( stderr, "Cannot open %s\n", inputfile );
			exit( 1 );
		}
	}
	else    
		infp = stdin;

#if 0
	PreRead( stdin, &njob, &nlenmax );
#else
	getnumlen( infp ); //defined in io.c. finds sequences count, max length and dna or protein from infp file
#endif
	rewind( infp );

	nkozo = 0;

	if( njob < 2 ) 
	{
		fprintf( stderr, "At least 2 sequences should be input!\n"
						 "Only %d sequence found.\n", njob ); 
		exit( 1 );
	}


	if( subalignment )
	{
		//First call of this method with table = NULL reads number of subalignments in subalignments file and assign to nsubpt
		//also reads max number of spaces in all sequences into maxmem
		//Second call reads data from the file and fills table with it
		readsubalignmentstable( njob, NULL, NULL, &nsubalignments, &maxmem ); //defined in io.c.
		fprintf( stderr, "nsubalignments = %d\n", nsubalignments );
		fprintf( stderr, "maxmem = %d\n", maxmem );
		subtable = AllocateIntMtx( nsubalignments, maxmem+1 );
		insubtable = AllocateIntVec( njob );
		preservegaps = AllocateIntVec( njob );
		for( i=0; i<njob; i++ ) insubtable[i] = 0;
		for( i=0; i<njob; i++ ) preservegaps[i] = 0;
		subalnpt = AllocateCharCub( nsubalignments, maxmem, 0 );
		readsubalignmentstable( njob, subtable, preservegaps, NULL, NULL );
		for( i=0; i<nsubalignments; i++ ) for( j=0; j<insubtable[i]; j++ )
		{
			if( subtable[i][j] < 0 )
			{
				fprintf( stderr, "Not supported in the iterative refinmenment mode.\n" );
				fprintf( stderr, "Please use a positive number to represent a sequence.\n" );
			}
		}
	}

	ocut = cut;

	segment = (Segment *)calloc( MAXSEG, sizeof( Segment ) ); //MAXSEG defined in mltaln.h.
//	if( treemethod == 'X' || treemethod == 'E' || treemethod == 'q' )
	topol = AllocateIntCub( njob, 2, njob );
	len = AllocateDoubleMtx( njob, 2 );
	eff = AllocateDoubleMtx( njob, njob );
	seq = AllocateCharMtx( njob, nlenmax*9+1 );
	name = AllocateCharMtx( njob, B+1 );
	seq_g = AllocateCharMtx( njob, nlenmax*9+1 ); //seq_g defined in defs.h.
	res_g = AllocateCharMtx( njob, nlenmax*9+1 );
	aseq = AllocateCharMtx( njob, nlenmax*9+1 );
	bseq = AllocateCharMtx( njob, nlenmax*9+1 );
	nogap1seq = AllocateCharMtx( 1, nlenmax*1+1 );
	alloclen = nlenmax * 9;
	seq_g_bk = AllocateCharMtx( njob, 0 );
	for( i=0; i<njob; i++ ) seq_g_bk[i] = seq_g[i];
	kozoarivec = AllocateCharVec( njob );
	skipthisbranch = AllocateIntMtx( njob, 2 );
	for( i=0; i<njob; i++ ) skipthisbranch[i][0] = skipthisbranch[i][1] = 0;


#if 0
	Read( name, nlen, seq_g );
#else
	readData_pointer( infp, name, nlen, seq_g ); //defined in io.c. reads sequences and their names from infp file into seq_g, name and nlen arrays.
#endif
	fclose( infp );

	if( specifictarget ) //fill targetmap and targetmapr with values based on 'focus' parameters in sequences names
	{
		targetmap = calloc( njob, sizeof( int ) );
		ntarget = 0;
		for( i=0; i<njob; i++ )
		{
			targetmap[i] = -1;
			if( !strncmp( name[i]+1, "_focus_", 7 ) )
				targetmap[i] = ntarget++;
		}
		targetmapr = calloc( ntarget, sizeof( int ) );
		for( i=0; i<njob; i++ )
			if( targetmap[i] != -1 ) targetmapr[targetmap[i]] = i;
	}
	else //target all sequences
	{
		ntarget = njob;
		targetmap = calloc( njob, sizeof( int ) );
		targetmapr = calloc( njob, sizeof( int ) );
		for( i=0; i<njob; i++ )
			targetmap[i] = targetmapr[i] = i;
	}

	if( constraint )
	{
		ilim = njob;
		localhomtable = (LocalHom **)calloc( ntarget, sizeof( LocalHom *) );
		for( i=0; i<ntarget; i++) //initialize localhomtable
		{
			localhomtable[i] = (LocalHom *)calloc( ilim, sizeof( LocalHom ) );
			for( j=0; j<ilim; j++)
			{
				localhomtable[i][j].start1 = -1;
				localhomtable[i][j].end1 = -1;
				localhomtable[i][j].start2 = -1;
				localhomtable[i][j].end2 = -1;
				localhomtable[i][j].overlapaa = -1.0; 
				localhomtable[i][j].opt = -1.0;
				localhomtable[i][j].importance = -1.0; 
				localhomtable[i][j].next = NULL; 
				localhomtable[i][j].nokori = 0;
				localhomtable[i][j].extended = -1;
				localhomtable[i][j].last = localhomtable[i]+j;
				localhomtable[i][j].korh = 'h';
			}
			if( !specifictarget ) ilim--;
		}
		fprintf( stderr, "Loading 'hat3' ... " );
		fflush( stderr );
		prep = fopen( "hat3", "r" ); //open hat3 file for reading
		if( prep == NULL ) ErrorExit( "Make hat3." );
		if( specifictarget ) readlocalhomtable2_target( prep, njob, localhomtable, kozoarivec, targetmap ); //defined in io.c. updates localhomtable values based on prep content. also set kozoarivec value based on it.
		else readlocalhomtable2_half( prep, njob, localhomtable, kozoarivec ); //defined in io.c. updates localhomtable values based on prep content. also set kozoarivec value based on it.
		//the difference between this and the previous one is the absence of targetmap here
		fclose( prep ); 

//		for( i=0; i<njob-1; i++ ) for( j=i+1; j<njob; j++ )
//			fprintf( stdout, "%d %d %d %d %d %d %d\n", i, j, localhomtable[i][j].opt, localhomtable[i][j].start1, localhomtable[i][j].end1, localhomtable[i][j].start2, localhomtable[i][j].end2 );
		fprintf( stderr, "done.\n" );
		fflush( stderr );
#if 0
		fprintf( stderr, "Extending localhom ... " );
		extendlocalhom( njob, localhomtable );
		fprintf( stderr, "done.\n" );
#endif
		nkozo = 0;
		for( i=0; i<njob; i++ ) if( kozoarivec[i] ) nkozo++;
	}


//		if( nlenmax > 30000 )
		if( nlenmax > 50000 ) // version >= 6.823
		{
#if 0
			if( constraint )
			{
				fprintf( stderr, "\nnlenmax=%d, nagasugi!\n", nlenmax ); 
				exit( 1 );
			}
			if( nevermemsave )
			{
				fprintf( stderr, "\nnevermemsave=1, nlenmax=%d, nagasugi!\n", nlenmax ); 
				exit( 1 );
			}
#endif
			if( !constraint && !nevermemsave && alg != 'M' )
			{
				fprintf( stderr, "\nnlenmax=%d, Switching to the memsave mode\n", nlenmax ); 
				alg = 'M';
			}
		}


	if( specificityconsideration ) calcmaxdistclass(); //defined here. calculate maxdistclass based on some calcs on ndistclass

	for( i=0; i<njob; i++ )
	{
		res_g[i][0] = 0;
	}

	identity = 1;
	for( i=0; i<njob; i++ ) 
	{
		identity *= ( nlen[i] == nlen[0] );
	}
	if( !identity ) //check equal length of all sequences
	{
		fprintf( stderr, "Input pre-aligned data\n" );
		exit( 1 );
	}
	constants( njob, seq_g ); //defined in constants.c.
	//after all this method, n_dis, ribosumdis, amino_dis, amino_dis_consweight_multi, n_dis_consweight_multi,
	//n_disLN, n_disFFT, polarity, volume arrays are initialized and some constants are set.

#if 0
	fprintf( stderr, "penalty = %d\n", penalty ); 
	fprintf( stderr, "penalty_ex = %d\n", penalty_ex ); 
	fprintf( stderr, "offset = %d\n", offset ); 
#endif

	initSignalSM(); //defined in io.c. inits signalSM value which is defined in defs.h.

	initFiles(); //defined in io.c. init prep_g and trap_g files. I think these files are for tracing

#if 0
	if( njob == 2 )
	{
		writePre( njob, name, nlen, seq_g, 1 );
		SHOWVERSION;
		return( 0 );
	}
#else
	if( njob == 2 )
	{
		weight = 0;
		niter = 1;
	}
#endif

	c = seqcheck( seq_g ); //defined in mltaln9.c. check sequence characters and report error if unusual character is found
	if( c )
	{
		fprintf( stderr, "Illegal character %c\n", c );
		exit( 1 );
	}
	commongappick( njob, seq_g ); //defined in mltaln9.c. update seq_g values based on gaps positions in it and njob value

	if( rnakozo && rnaprediction == 'm' )
	{
		singlerna = (RNApair ***)calloc( njob, sizeof( RNApair ** ) );
		prep = fopen( "hat4", "r" ); //what is hat4 ?  //open hat4 file for reading
		if( prep == NULL ) ErrorExit( "Make hat4 using mccaskill." );
		fprintf( stderr, "Loading 'hat4' ... " );
		fflush( stderr );
		for( i=0; i<njob; i++ )
		{
			gappick0( nogap1seq[0], seq_g[i] ); //defined in mltaln9.c. copy 'seq_g' chars to 'nogap1seq' without gaps chars.
			nogaplen = strlen( nogap1seq[0] );
			singlerna[i] = (RNApair **)calloc( nogaplen+1, sizeof( RNApair * ) );
			for( j=0; j<nogaplen; j++ ) //initialize singlerna[i]
			{
				singlerna[i][j] = (RNApair *)calloc( 1, sizeof( RNApair ) );
				singlerna[i][j][0].bestpos = -1;
				singlerna[i][j][0].bestscore = -1.0;
			}
			singlerna[i][nogaplen] = NULL;
			readmccaskill( prep, singlerna[i], nogaplen ); //defined in io.c. fill singlerna[i] with values from prep file
		}
		fclose( prep );
		fprintf( stderr, "\ndone.\n" );
		fflush( stderr );
	}
	else
		singlerna = NULL;



	if( utree )
	{
		prep = fopen( "hat2", "r" ); //what is hat2 ?  //open hat2 file for reading
		if( !prep ) ErrorExit( "Make hat2." );
		readhat2_pointer( prep, njob, name, eff ); //defined in io.c. fill eff with values read from prep
		fclose( prep );
#if 0
		fprintf( stderr, "eff = \n" );
		for( i=0; i<njob-1; i++ ) 
		{
			for( j=i+1; j<njob; j++ ) 
			{
				fprintf( stderr, "%d-%d,  %f\n", i, j, eff[i][j] );
			}
			fprintf( stderr, "\n" );
		}
#endif
		if( intree )
		{
			veryfastsupg_double_loadtree( njob, eff, topol, len, name ); //defined in mltaln9.c.
			//update eff and len values based on values and computations on data from guidetree file and write results to infile.tree file
#if 0
			fprintf( stderr, "eff = \n" );
			for( i=0; i<njob-1; i++ ) 
			{
				for( j=i+1; j<njob; j++ ) 
				{
					fprintf( stderr, "%d-%d,  %f\n", i, j, eff[i][j] );
				}
				fprintf( stderr, "\n" );
			}
exit( 1 );
#endif
		}
		else if( intop ) // v6.528 deha if( intop ) dattanode intree ga mukou datta.
		{
			fprintf( stderr, "--topin has been disabled\n" );
			exit( 1 );
//			veryfastsupg_double_loadtop( njob, eff, topol, len );
		}
		else if( subalignment )
			fixed_supg_double_treeout_constrained( njob, eff, topol, len, name, nsubalignments, subtable ); //defined in mltaln9.c.
			//update eff and write tree to infile.tree. I think this method is related to tree construction/update and groups
		else if( treemethod == 'X' || treemethod == 'E' || treemethod == 'q' ) 
//			veryfastsupg_double_outtree( njob, eff, topol, len, name );
			fixed_musclesupg_double_treeout( njob, eff, topol, len, name ); //defined in mltaln9.c.
			//update eff and write tree to infile.tree. I think this method is related to tree construction/update and groups
		else if( treemethod == 'n' ) 
			nj( njob, eff, topol, len ); //defined in nj.c. may be meaning neighbor joining ?!!
		else if( treemethod == 's' )
			spg( njob, eff, topol, len ); //defined in mltaln9.c.
		else if( treemethod == 'p' )
			upg2( njob, eff, topol, len ); //defined in mltaln9.c.
		else ErrorExit( "Incorrect treemethod.\n" );
	}
#if DEBUG
	printf( "utree = %d\n", utree );
#endif

	if( autosubalignment > 0.0 && subalignment == 0 )
	{
//		reporterr( "Computing skipthisbranch..\n" );
		insubtable = AllocateIntVec( njob );
		preservegaps = AllocateIntVec( njob );
		subtable = calloc( 1, sizeof( char * ) );
		subtable[0] = NULL; // for FreeIntMtx
		for( i=0; i<njob; i++ ) insubtable[i] = 0;
		for( i=0; i<njob; i++ ) preservegaps[i] = 0; // tsukawanaikamo
		if( generatesubalignmentstable( njob, &subtable, &nsubalignments, &maxmem, topol, len, autosubalignment ) ) // subtable ha allocate sareru.
		{
			reporterr( "################################################################################################ \n" );
			reporterr( "#\n" );
			reporterr( "# WARNING: Iterative refinment was not done because you gave a large --fix value (%6.3f).\n", autosubalignment );
			reporterr( "#\n" );
			reporterr( "################################################################################################ \n" );
			writePre( njob, name, nlen, seq_g, 1 );

			FreeCharMtx( seq_g_bk );
			FreeIntCub( topol );
			FreeDoubleMtx( len );
			FreeDoubleMtx( eff );
			FreeIntMtx( skipthisbranch );
			FreeIntMtx( subtable );
			free( preservegaps );
			free( insubtable );
			SHOWVERSION;
			return( 0 );
		}
//		subtable = AllocateIntMtx( nsubalignments, maxmem+1 );
		fprintf( stderr, "nsubalignments = %d, maxmem = %d\n", nsubalignments, maxmem );
		subalnpt = AllocateCharCub( nsubalignments, maxmem, 0 );
#if 0
		for( i=0; i<nsubalignments; i++ )
		{
			reporterr( "subalignment %d\n", i );
			for( j=0; subtable[i][j]!=-1; j++ )
			{
				reporterr( "%5d", subtable[i][j] );
			}
			reporterr( "\n" );
		}
#endif
#if 0 // wakaran
		for( i=0; i<nsubalignments; i++ ) for( j=0; j<insubtable[i]; j++ )
		{
			if( subtable[i][j] < 0 )
			{
				fprintf( stderr, "Not supported in the iterative refinmenment mode.\n" );
				fprintf( stderr, "Please use a positive number to represent a sequence.\n" );
			}
		}
#endif
//		reporterr( "done.\n" );
	}


	orderfp = fopen( "order", "w" ); //open order file for writing
	if( !orderfp )
	{
		fprintf( stderr, "Cannot open 'order'\n" );
		exit( 1 );
	}
	for( i=0; (j=topol[njob-2][0][i])!=-1; i++ ) //write values from topol to order file
	{
		fprintf( orderfp, "%d\n", j );
	}
	for( i=0; (j=topol[njob-2][1][i])!=-1; i++ ) //write values from topol to order file
	{
		fprintf( orderfp, "%d\n", j );
	}
	fclose( orderfp );



	fprintf( stderr, "\n" );
	if( ( !utree && kobetsubunkatsu ) || constraint || !bunkatsu )
	{
		nseg = 0;
		anchors[0] =0;
		anchors[1] =strlen( seq_g[0] );
		nseg += 2;
	}
	else
	{
		nseg = searchAnchors( njob, seq_g, segment ); //defined in mltaln9.c.
#if 0
		fprintf( stderr, "### nseg = %d\n", nseg );
		fprintf( stderr, "### seq_g[0] = %s\n", seq_g[0] );
		fprintf( stderr, "### nlenmax = %d\n", nlenmax );
		fprintf( stderr, "### len = %d\n", strlen( seq_g[0] ) );

#endif

		anchors[0] = 0;
		for( i=0; i<nseg; i++ ) anchors[i+1] = segment[i].center;
		anchors[nseg+1] = strlen( seq_g[0] );
		nseg +=2;

#if 0
		for( i=0; i<nseg; i++ )
			fprintf( stderr, "anchor[%d] = %d\n", i, anchors[i] );
#endif
	}


//--------------- kokokara ----
	if( subalignment || autosubalignment )
	{
		for( i=0; i<nsubalignments; i++ )
		{
			fprintf( stderr, "\nChecking subalignment %d:\n", i+1 );
			alignmentlength = strlen( seq[subtable[i][0]] );
			for( j=0; subtable[i][j]!=-1; j++ )
				fprintf( stderr, " %d ", subtable[i][j]+1 );
//				fprintf( stderr, " ##### %d-%d. %-30.30s\n", i, subtable[i][j]+1, name[subtable[i][j]]+1 );
			fprintf( stderr, "\n" );
			for( j=0; subtable[i][j]!=-1; j++ )
			{
				if( subtable[i][j] >= njob )
				{
					fprintf( stderr, "No such sequence, %d.\n", subtable[i][j]+1 );
					exit( 1 );
				}
				if( alignmentlength != strlen( seq[subtable[i][j]] ) )
				{
					fprintf( stderr, "\n" );
					fprintf( stderr, "###############################################################################\n" );
					fprintf( stderr, "# ERROR!\n" );
					fprintf( stderr, "# Subalignment %d must be aligned.\n", i+1 );
					fprintf( stderr, "# Please check the alignment lengths of following sequences.\n" );
					fprintf( stderr, "#\n" );
					fprintf( stderr, "# %d. %-10.10s -> %d letters (including gaps)\n", subtable[i][0]+1, name[subtable[i][0]]+1, alignmentlength );
					fprintf( stderr, "# %d. %-10.10s -> %d letters (including gaps)\n", subtable[i][j]+1, name[subtable[i][j]]+1, (int)strlen( seq[subtable[i][j]] ) );
					fprintf( stderr, "#\n" );
					fprintf( stderr, "# See http://mafft.cbrc.jp/alignment/software/merge.html for details.\n" );
					if( subalignmentoffset )
					{
						fprintf( stderr, "#\n" );
						fprintf( stderr, "# You specified seed alignment(s) consisting of %d sequences.\n", subalignmentoffset );
						fprintf( stderr, "# In this case, the rule of numbering is:\n" );
						fprintf( stderr, "#   The aligned seed sequences are numbered as 1 .. %d\n", subalignmentoffset );
						fprintf( stderr, "#   The input sequences to be aligned are numbered as %d .. %d\n", subalignmentoffset+1, subalignmentoffset+njob );
					}
					fprintf( stderr, "###############################################################################\n" );
					fprintf( stderr, "\n" );
					exit( 1 );
				}
				insubtable[subtable[i][j]] = 1;
			}
			for( j=0; j<njob-1; j++ )
			{
#if 0
				int k;
				reporterr( "####  STEP%d\n", j );
				for( k=0; topol[j][0][k]!=-1; k++ ) reporterr( "%3d ", topol[j][0][k] );
				reporterr( "len=%f\n", len[j][0] );
				for( k=0; topol[j][1][k]!=-1; k++ ) reporterr( "%3d ", topol[j][1][k] );
				reporterr( "len=%f\n", len[j][1] );
				reporterr( "\n" );
#endif
				if( includemember( topol[j][0], subtable[i] ) && !samemember( topol[j][0], subtable[i] ) )
				{
					skipthisbranch[j][0] = 1;
//					reporterr( "SKIP 0 !!!!!!\n" );
				}
				if( includemember( topol[j][1], subtable[i] ) && !samemember( topol[j][1], subtable[i] ) )
				{
					skipthisbranch[j][1] = 1;
//					reporterr( "SKIP 1 !!!!!!\n" );
				}
			}
			foundthebranch = 0;
			for( j=0; j<njob-1; j++ )
			{
				if( samemember( topol[j][0], subtable[i] ) || samemember( topol[j][1], subtable[i] ) )
				{
					foundthebranch = 1;
					fprintf( stderr, " -> OK\n" );
					break;
				}
			}
			if( !foundthebranch )
			{
				system( "cp infile.tree GuideTree" ); // tekitou
				fprintf( stderr, "\n" );
				fprintf( stderr, "###############################################################################\n" );
				fprintf( stderr, "# ERROR!\n" );
				fprintf( stderr, "# Subalignment %d does not seem to form a monophyletic cluster\n", i+1 );
				fprintf( stderr, "# in the guide tree ('GuideTree' in this directory) internally computed.\n" );
				fprintf( stderr, "# If you really want to use this subalignment, pelase give a tree with --treein \n" );
				fprintf( stderr, "# http://mafft.cbrc.jp/alignment/software/treein.html\n" );
				fprintf( stderr, "# http://mafft.cbrc.jp/alignment/software/merge.html\n" );
				if( subalignmentoffset )
				{
					fprintf( stderr, "#\n" );
					fprintf( stderr, "# You specified seed alignment(s) consisting of %d sequences.\n", subalignmentoffset );
					fprintf( stderr, "# In this case, the rule of numbering is:\n" );
					fprintf( stderr, "#   The aligned seed sequences are numbered as 1 .. %d\n", subalignmentoffset );
					fprintf( stderr, "#   The input sequences to be aligned are numbered as %d .. %d\n", subalignmentoffset+1, subalignmentoffset+njob );
				}
				fprintf( stderr, "############################################################################### \n" );
				fprintf( stderr, "\n" );
				exit( 1 );
			}
//			commongappick( seq[subtable[i]], subalignment[i] ); // irukamo
		}
#if 0
		for( i=0; i<njob-1; i++ )
		{
			fprintf( stderr, "STEP %d\n", i+1 );
			fprintf( stderr, "group1 = " );
			for( j=0; topol[i][0][j] != -1; j++ )
				fprintf( stderr, "%d ", topol[i][0][j]+1 );
			fprintf( stderr, "\n" );
			fprintf( stderr, "SKIP -> %d\n\n", skipthisbranch[i][0] );
			fprintf( stderr, "group2 = " );
			for( j=0; topol[i][1][j] != -1; j++ )
				fprintf( stderr, "%d ", topol[i][1][j]+1 );
			fprintf( stderr, "\n" );
			fprintf( stderr, "SKIP -> %d\n\n", skipthisbranch[i][1] );
		}
#endif

		for( i=0; i<njob; i++ ) 
		{
			if( insubtable[i] ) strcpy( bseq[i], seq[i] );
			else gappick0( bseq[i], seq[i] );
		}

		for( i=0; i<nsubalignments; i++ ) 
		{
			for( j=0; subtable[i][j]!=-1; j++ ) subalnpt[i][j] = bseq[subtable[i][j]];
			commongappick( j, subalnpt[i] );
		}

		FreeIntMtx( subtable );
		free( insubtable );
		for( i=0; i<nsubalignments; i++ ) free( subalnpt[i] );
		free( subalnpt );
		free( preservegaps );
	}
//--------------- kokomade ----




	for( i=0; i<njob; i++ ) res_g[i][0] = 0;

	for( iseg=0; iseg<nseg-1; iseg++ )
	{
		int tmplen = anchors[iseg+1]-anchors[iseg];
		int pos = strlen( res_g[0] );
		for( j=0; j<njob; j++ )
		{
			strncpy( seq[j], seq_g[j], tmplen );
			seq[j][tmplen]= 0;
			seq_g[j] += tmplen;	

		}
		fprintf( stderr, "Segment %3d/%3d %4d-%4d\n", iseg+1, nseg-1, pos+1, pos+1+tmplen );
		fflush( stderr );
		fprintf( trap_g, "Segment %3d/%3d %4d-%4d\n", iseg+1, nseg-1, pos+1, pos+1+tmplen );
	
		cut = ocut;
		returnvalue = TreeDependentIteration( njob, name, nlen, seq, bseq, topol, len, eff, skipthisbranch, alloclen, localhomtable, singlerna, nkozo, kozoarivec, ntarget, targetmap, targetmapr );
		//defined in tditeration.c.

		for( i=0; i<njob; i++ )
			strcat( res_g[i], bseq[i] );
	}
	FreeCharMtx( seq_g_bk );
	FreeIntCub( topol );
	FreeDoubleMtx( len );
	FreeDoubleMtx( eff );
	FreeIntMtx( skipthisbranch );
	free( kozoarivec );
	if( constraint ) 
	{
		if( specifictarget ) FreeLocalHomTable_part( localhomtable, ntarget, njob );
		else FreeLocalHomTable_half( localhomtable, njob );
	}
	free( targetmap );
	free( targetmapr );
	if( rnakozo && rnaprediction == 'm' ) 
	{
		if( singlerna ) // nen no tame
		{
			for( i=0; i<njob; i++ ) 
			{
				for( j=0; singlerna[i][j]!=NULL; j++ )
				{
					if( singlerna[i][j] ) free( singlerna[i][j] );
				}
				if( singlerna[i] ) free( singlerna[i] );
			}
			free( singlerna );
			singlerna = NULL;
		}
	}

#if 0
	Write( stdout, njob, name, nlen, bseq );
#endif

	fprintf( stderr, "done\n" );
	fprintf( trap_g, "done\n" );
	fclose( trap_g );
	freeconstants();


	devide = 0; 
	writePre( njob, name, nlen, res_g, 1 ); //defined in io.c. write sequences and their names to prep_g file
#if 0
	writeData( stdout, njob, name, nlen, res_g, 1 );
#endif


	if( spscoreout ) reporterr( "Unweighted sum-of-pairs score = %10.5f\n", sumofpairsscore( njob, res_g ) );

	SHOWVERSION;
	return( 0 );
}

#if 0
signed int main( int argc, char *argv[] )
{
	int i, nlen[M];
	char b[B];
	char a[] = "=";
	int value;

	gets( b ); njob = atoi( b );

/*
	scoremtx = 0;
	if( strstr( b, "ayhoff" ) ) scoremtx = 1;
	else if( strstr( b, "dna" ) || strstr( b, "DNA" ) ) scoremtx = -1;
	else if( strstr( b, "M-Y" ) || strstr( b, "iyata" ) ) scoremtx = 2;
	else scoremtx = 0;
*/
	if( strstr( b, "constraint" ) ) cnst = 1;

	nlenmax = 0;
	i = 0;
	while( i<njob )
	{
		gets( b );
		if( !strncmp( b, a, 1 ) ) 
		{
			gets( b ); nlen[i] = atoi( b );
			if( nlen[i] > nlenmax ) nlenmax = nlen[i];
			i++;
		}
	}
	if( nlenmax > N || njob > M ) 
	{
		fprintf( stderr, "ERROR in main\n" );
		exit( 1 );
	}
	/*
	nlenmax = Na;
	*/
	rewind( stdin );
	value = main1( nlen, argc, argv );
	exit( 0 );
}
#endif
