#! /bin/sh

er=0;
myself=`dirname "$0"`/`basename "$0"`; export myself
version="v7.305b (2016/Aug/16)"; export version
LANG=C; export LANG
os=`uname`
progname=`basename "$0"`
if [ `echo $os | grep -i cygwin` ]; then
	os="cygwin"
elif [ `echo $os | grep -i mingw` ]; then
	os="mingw"
elif [ `echo $os | grep -i darwin` ]; then
	os="darwin"
elif [ `echo $os | grep -i sunos` ]; then
	os="sunos"
elif [ `echo $os | grep -i linux` ]; then
	os="linux"
else
	os="unix"
fi
export os

if [ "$MAFFT_BINARIES" ]; then
	prefix="$MAFFT_BINARIES"
else
	prefix=_LIBDIR
fi
export prefix

if [ $# -gt 0 ]; then
	if [ "$1" = "--man" ]; then 
		man "$prefix/mafft.1"
		exit 0;
	fi
fi

if [ -x "$prefix/version" ]; then
		versionbin=`"$prefix/version" | awk '{print $1}'` # for cygwin
	else
		versionbin="0.000"
fi

if ! expr "$version" : v"$versionbin" > /dev/null ; then
	echo "" 1>&2
	echo "v$versionbin != $version" 1>&2
	echo "" 1>&2
	echo "There is a problem in the configuration of your shell." 1>&2
	echo "Check the MAFFT_BINARIES environmental variable by" 1>&2
	echo "$ echo \$MAFFT_BINARIES" 1>&2
	echo "" 1>&2
	echo "This variable must be *unset*, unless you have installed MAFFT" 1>&2
	echo "with a special configuration.  To unset this variable, type" 1>&2
	echo "$ unset MAFFT_BINARIES" 1>&2
	echo "or" 1>&2
	echo "% unsetenv MAFFT_BINARIES" 1>&2
	echo "Then retry" 1>&2
	echo "$ mafft input > output" 1>&2
	echo "" 1>&2
	echo "To keep this change permanently, edit setting files" 1>&2
	echo "(.bash_profile, .profile, .cshrc, etc) in your home directory" 1>&2
	echo "to delete the MAFFT_BINARIES line." 1>&2
	echo "On MacOSX, also edit or remove the .MacOSX/environment.plist file" 1>&2
	echo "and then re-login (MacOSX 10.6) or reboot (MacOSX 10.7)." 1>&2
	echo "" 1>&2
	echo "Please send a problem report to kazutaka.katoh@aist.go.jp," 1>&2
	echo "if this problem remains." 1>&2
	echo "" 1>&2
	exit 1
	er=1
fi

defaultiterate=0
defaultcycle=2
defaultgop="1.53"
#defaultaof="0.123"
defaultaof="0.000"
defaultlaof="0.100"
defaultlgop="-2.00"
defaultfft=1
defaultrough=0
defaultdistance="ktuples"
#defaultdistance="local"
defaultweighti="2.7"
defaultweightr="0.0"
defaultweightm="1.0"
defaultdafs=0
defaultmccaskill=0
defaultcontrafold=0
defaultalgopt="  "
defaultalgoptit="  "
defaultsbstmodel=" -b 62 "
defaultfmodel=" "
defaultkappa=" "
if [ $progname = "xinsi" -o $progname = "mafft-xinsi" ]; then
	defaultfft=1
	defaultcycle=1
	defaultiterate=1000
	defaultdistance="scarna"
	defaultweighti="3.2"
	defaultweightr="8.0"
	defaultweightm="2.0"
	defaultmccaskill=1
	defaultcontrafold=0
	defaultdafs=0
	defaultalgopt=" -A "
	defaultalgoptit=" -AB " ## chui
	defaultaof="0.0"
	defaultsbstmodel=" -b 62 "
	defaultkappa=" "
	defaultfmodel="    "  # 2013/06/18
elif [ $progname = "qinsi" -o $progname = "mafft-qinsi" ]; then
	defaultfft=1
	defaultcycle=1
	defaultiterate=1000
	defaultdistance="global"
	defaultweighti="3.2"
	defaultweightr="8.0"
	defaultweightm="2.0"
	defaultmccaskill=1
	defaultcontrafold=0
	defaultdafs=0
	defaultalgopt=" -A "
	defaultalgoptit=" -AB " ## chui
	defaultaof="0.0"
	defaultsbstmodel=" -b 62 "
	defaultkappa=" "
	defaultfmodel="    "  # 2013/06/18
elif [ $progname = "linsi" -o $progname = "mafft-linsi" ]; then
	defaultfft=0
	defaultcycle=1
	defaultiterate=1000
	defaultdistance="local"
elif [ $progname = "ginsi" -o $progname = "mafft-ginsi" ]; then
	defaultfft=1
	defaultcycle=1
	defaultiterate=1000
	defaultdistance="global"
elif [ $progname = "einsi" -o $progname = "mafft-einsi" ]; then
	defaultfft=0
	defaultcycle=1
	defaultiterate=1000
	defaultdistance="localgenaf"
elif [ $progname = "fftns" -o $progname = "mafft-fftns" ]; then
	defaultfft=1
	defaultcycle=2
	defaultdistance="ktuples"
elif [ $progname = "fftnsi" -o $progname = "mafft-fftnsi" ]; then
	defaultfft=1
	defaultcycle=2
	defaultiterate=2
	defaultdistance="ktuples"
elif [ $progname = "nwns" -o $progname = "mafft-nwns" ]; then
	defaultfft=0
	defaultcycle=2
	defaultdistance="ktuples"
elif [ $progname = "nwnsi" -o $progname = "mafft-nwnsi" ]; then
	defaultfft=0
	defaultcycle=2
	defaultiterate=2
	defaultdistance="ktuples"
fi
outputfile=""
namelength=-1
anysymbol=0
parallelizationstrategy="BAATARI2"
kappa=$defaultkappa
sbstmodel=$defaultsbstmodel
fmodel=$defaultfmodel
nmodel=" "
gop=$defaultgop
gopdist=$defaultgop
aof=$defaultaof
cycle=$defaultcycle
iterate=$defaultiterate
fft=$defaultfft
rough=$defaultrough
distance=$defaultdistance
forcefft=0
memopt=" "
weightopt=" "
GGOP="-6.00"
LGOP="-6.00"
LEXP="-0.000"
GEXP="-0.000"
lgop=$defaultlgop
lexp="-0.100"
laof=$defaultlaof
pggop="-2.00"
pgexp="-0.10"
pgaof="0.10"
rgop="-1.530"
rgep="-0.000"
seqtype="  "
weighti=$defaultweighti
weightr=$defaultweightr
weightm=$defaultweightm
rnaalifold=0
dafs=$defaultdafs
mccaskill=$defaultmccaskill
contrafold=$defaultcontrafold
progressfile="/dev/stderr"
debug=0
sw=0
algopt=$defaultalgopt
algoptit=$defaultalgoptit
#algspecified=0
pairspecified=0
scorecalcopt=" "
coreout=0
corethr="0.5"
corewin="100"
coreext=" "
outputformat="pir"
f2clext="-N"
outorder="input"
seed="x"
seedtable="x"
auto=0
groupsize=-1
partsize=50
partdist="ktuples"
partorderopt=" -x "
treeout=0
distout=0
treein=0
topin=0
treeinopt="  "
seedfiles="/dev/null"
seedtablefile="/dev/null"
pdblist="/dev/null"
ownlist="/dev/null"
strdir="$PWD"
aamatrix="/dev/null"
treeinfile="/dev/null"
rnascoremtx=" "
laraparams="/dev/null"
foldalignopt=" "
treealg=" -X 0.1 "
sueff="1.0"
scoreoutarg=" "
numthreads=0
numthreadsit=-1
numthreadstb=-1
randomseed=0
addfile="/dev/null"
addarg0=" "
addarg=" "
addsinglearg=" "
add2ndhalfarg=" "
mapoutfile="/dev/null"
fragment=0
legacygapopt=" "
mergetable="/dev/null"
mergearg=" "
seedoffset=0
outnum=" "
last_e=5000
last_m=3
last_subopt=" "
last_once=" "
adjustdirection=0
tuplesize=6
termgapopt=" -O "
#termgapopt=" " # gap/gap ga kakenai node
similarityoffset="0.0"
unalignlevel="0.0"
unalignspecified=0
spfactor="100.0"
shiftpenaltyspecified=0
opdistspecified=0
allowshift=0
enrich=0
enrichseq=0
enrichstr=0
seektarget=""
fixthreshold="0.0"
bunkatsuopt=" "
npickup=0
minimumweight="0.00001" # 2016/Mar
usenaivepairscore=" "
oldgenafparam=0
sprigorous=0
pileuporshuffle="l"
initialramusage="20GB"
focusarg=" "
if [ $# -gt 0 ]; then
	if [ "$1" = "--version" ]; then
		echo "$version" 1>&2
		exit 0;
	elif [ "$1" = "--help" -o "$1" = "--info" ]; then 
		shift
		er=1;
	fi
	while [ $# -gt 1 ];
	do
		if [ "$1" = "--auto" ]; then 
			auto=1
		elif [ "$1" = "--anysymbol" ]; then 
			anysymbol=1
		elif [ "$1" = "--preservecase" ]; then 
			anysymbol=1
		elif [ "$1" = "--clustalout" ]; then 
			outputformat="clustal"
		elif [ "$1" = "--phylipout" ]; then 
			outputformat="phylip"
		elif [ "$1" = "--reorder" ]; then 
			outorder="aligned"
			partorderopt=" "
		elif [ "$1" = "--inputorder" ]; then 
			outorder="input"
			partorderopt=" -x "
		elif [ "$1" = "--unweight" ]; then 
			weightopt=" -u "
		elif [ "$1" = "--termgappenalty" ]; then 
			termgapopt=" "
		elif [ "$1" = "--alga" ]; then 
			algopt=" "
			algoptit=" "
#			algspecified=1
		elif [ "$1" = "--algq" ]; then 
			algopt=" -Q "
			algoptit=" "
			echo "" 1>&2
			echo "--algq is no longer supported!" 1>&2
			echo "" 1>&2
			exit 1;
#			algspecified=1
		elif [ "$1" = "--namelength" ]; then 
			shift   
			namelength=`expr "$1" - 0`
			if ! expr "$1" : "[0-9]" > /dev/null ; then
				echo "Specify the length of name in clustal format output!" 1>&2
				exit
			fi
		elif [ "$1" = "--groupsize" ]; then 
			shift   
			groupsize=`expr "$1" - 0`
			if ! expr "$1" : "[0-9]" > /dev/null ; then
				echo "Specify groupsize!" 1>&2
				exit
			fi
		elif [ "$1" = "--partsize" ]; then 
			shift   
			partsize=`expr "$1" - 0`
			if ! expr "$1" : "[0-9]" > /dev/null ; then
				echo "Specify partsize!" 1>&2
				exit
			fi
		elif [ "$1" = "--parttree" ]; then 
			distance="parttree"
			partdist="ktuples"
		elif [ "$1" = "--dpparttree" ]; then 
			distance="parttree"
			partdist="localalign"
		elif [ "$1" = "--fastaparttree" ]; then 
			distance="parttree"
			partdist="fasta"
		elif [ "$1" = "--treeout" ]; then 
			treeout=1
		elif [ "$1" = "--distout" ]; then 
			distout=1
		elif [ "$1" = "--fastswpair" ]; then
			distance="fasta"
			pairspecified=1
			sw=1
		elif [ "$1" = "--fastapair" ]; then
			distance="fasta"
			pairspecified=1
			sw=0
		elif [ "$1" = "--averagelinkage" ]; then
			treealg=" -X 1.0 "
			sueff="1.0"
		elif [ "$1" = "--minimumlinkage" ]; then
			treealg=" -X 0.0 "
			sueff="0.0"
		elif [ "$1" = "--mixedlinkage" ]; then
			shift   
			sueff="$1"
			treealg=" -X $1"
		elif [ "$1" = "--noscore" ]; then
			scorecalcopt=" -Z "
		elif [ "$1" = "--6mermultipair" ]; then
			distance="ktuplesmulti"
			tuplesize=6
			pairspecified=1
		elif [ "$1" = "--10mermultipair" ]; then
			distance="ktuplesmulti"
			tuplesize=10
			pairspecified=1
		elif [ "$1" = "--6merpair" ]; then
			distance="ktuples"
			tuplesize=6
			pairspecified=1
		elif [ "$1" = "--10merpair" ]; then
			distance="ktuples"
			tuplesize=10
			pairspecified=1
		elif [ "$1" = "--blastpair" ]; then
			distance="blast"
			pairspecified=1
		elif [ "$1" = "--lastmultipair" ]; then
			distance="lastmulti"
			pairspecified=1
		elif [ "$1" = "--globalpair" ]; then
			distance="global"
			pairspecified=1
		elif [ "$1" = "--shortlongpair" ]; then
			distance="local"
			usenaivepairscore="-Z"
			laof=0.0  # addfull no tokini tsukawareru.
			lexp=0.0  # addfull no tokini tsukawareru.
			pgaof=0.0 # local nara iranai
			pgexp=0.0 # local nara iranai
			pairspecified=1
		elif [ "$1" = "--longshortpair" ]; then
			distance="local"
			usenaivepairscore="-Z"
			laof=0.0  # addfull no tokini tsukawareru.
			lexp=0.0  # addfull no tokini tsukawareru.
			pgaof=0.0 # local nara iranai
			pgexp=0.0 # local nara iranai
			pairspecified=1
		elif [ "$1" = "--localpair" ]; then
			distance="local"
			pairspecified=1
		elif [ "$1" = "--lastpair" ]; then
			distance="last"
			pairspecified=1
		elif [ "$1" = "--multipair" ]; then
			distance="multi"
			pairspecified=1
		elif [ "$1" = "--hybridpair" ]; then
			distance="hybrid"
			pairspecified=1
		elif [ "$1" = "--scarnapair" ]; then
			distance="scarna"
			pairspecified=1
		elif [ "$1" = "--dafspair" ]; then
			distance="dafs"
			pairspecified=1
		elif [ "$1" = "--larapair" ]; then
			distance="lara"
			pairspecified=1
		elif [ "$1" = "--slarapair" ]; then
			distance="slara"
			pairspecified=1
		elif [ "$1" = "--foldalignpair" ]; then
			distance="foldalignlocal"
			pairspecified=1
		elif [ "$1" = "--foldalignlocalpair" ]; then
			distance="foldalignlocal"
			pairspecified=1
		elif [ "$1" = "--foldalignglobalpair" ]; then
			distance="foldalignglobal"
			pairspecified=1
		elif [ "$1" = "--globalgenafpair" ]; then
			distance="globalgenaf"
			pairspecified=1
			echo "" 1>&2
			echo "--globalgenaf is no longer supported!" 1>&2
			echo "" 1>&2
			exit 1;
		elif [ "$1" = "--localgenafpair" ]; then
			distance="localgenaf"
			pairspecified=1
		elif [ "$1" = "--genafpair" ]; then
			distance="localgenaf"
			pairspecified=1
		elif [ "$1" = "--oldgenafpair" ]; then
			distance="localgenaf"
			pairspecified=1
			oldgenafparam=1
		elif [ "$1" = "--memsave" ]; then
			memopt=" -M -B "         # -B (bunkatsunashi no riyu ga omoidasenai)
		elif [ "$1" = "--nomemsave" ]; then
			memopt=" -N "
		elif [ "$1" = "--nuc" ]; then 
			seqtype=" -D "
		elif [ "$1" = "--amino" ]; then 
			seqtype=" -P "
		elif [ "$1" = "--fft" ]; then 
			fft=1
			forcefft=1
		elif [ "$1" = "--nofft" ]; then 
			fft=0
		elif [ "$1" = "--quiet" ]; then 
			if [ $os = "mingw" ]; then
				progressfile="nul"
			else
				progressfile="/dev/null"
			fi
		elif [ "$1" = "--debug" ]; then 
			debug=1
		elif [ "$1" = "--coreext" ]; then 
			coreext=" -c "
		elif [ "$1" = "--core" ]; then 
			coreout=1
		elif [ "$1" = "--adjustdirection" ]; then 
			adjustdirection=1
		elif [ "$1" = "--adjustdirectionaccurately" ]; then 
			adjustdirection=2
		elif [ "$1" = "--progress" ]; then 
			shift   
			progressfile="$1"
			if ! ( expr "$progressfile" : "\/" > /dev/null || expr "$progressfile" : "[A-Za-z]\:" > /dev/null ) ; then
				echo "Specify a progress file name with the absolute path!" 1>&2
				exit
			fi
		elif [ "$1" = "--out" ]; then 
			shift   
			outputfile="$1"
		elif [ "$1" = "--thread" ]; then 
			shift
			if ! expr "$1" : "[0-9\-]" > /dev/null ; then
				echo "Specify the number of threads.  Or, use --thread -1" 1>&2
				exit
			fi
			numthreads=`expr "$1" - 0` 
		elif [ "$1" = "--threadtb" ]; then 
			shift
			if ! expr "$1" : "[0-9\-]" > /dev/null ; then
				echo "Specify the number of threads for the iterative step!" 1>&2
				exit
			fi
			numthreadstb=`expr "$1" - 0` 
		elif [ "$1" = "--threadit" ]; then 
			shift
			if ! expr "$1" : "[0-9\-]" > /dev/null ; then
				echo "Specify the number of threads for the iterative step!" 1>&2
				exit
			fi
			numthreadsit=`expr "$1" - 0` 
		elif [ "$1" = "--last_subopt" ]; then 
			last_subopt="-S"
		elif [ "$1" = "--last_once" ]; then 
			last_once="-U"
		elif [ "$1" = "--last_m" ]; then 
			shift
			last_m=`expr "$1" - 0` 
		elif [ "$1" = "--last_e" ]; then 
			shift
			last_e=`expr "$1" - 0` 
		elif [ "$1" = "--randomseed" ]; then 
			shift
			randomseed=`expr "$1" - 0` 
		elif [ "$1" = "--bestfirst" ]; then 
			parallelizationstrategy="BESTFIRST"
		elif [ "$1" = "--adhoc0" ]; then 
			parallelizationstrategy="BAATARI0"
		elif [ "$1" = "--adhoc1" ]; then 
			parallelizationstrategy="BAATARI1"
		elif [ "$1" = "--adhoc2" ]; then 
			parallelizationstrategy="BAATARI2"
		elif [ "$1" = "--simplehillclimbing" ]; then 
			parallelizationstrategy="BAATARI2"
		elif [ "$1" = "--scoreout" ]; then 
			scoreoutarg="-S -B"
		elif [ "$1" = "--outnum" ]; then 
			scoreoutarg="-n"
		elif [ "$1" = "--leavegappyregion" ]; then 
			legacygapopt="-L"
		elif [ "$1" = "--legacygappenalty" ]; then 
			legacygapopt="-L"
		elif [ "$1" = "--merge" ]; then
			shift
			mergetable="$1"
			if [ ! -e "$mergetable" ]; then
				echo "Cannot open $mergetable" 1>&2
				echo "" 1>&2
				exit
			fi
		elif [ "$1" = "--addprofile" ]; then 
			shift   
			addarg0="-I"
			addfile="$1"
		elif [ "$1" = "--add" ]; then 
			shift   
			addarg0="-K -I"
			addfile="$1"
		elif [ "$1" = "--addfragments" ]; then 
			shift   
			addarg0="-K -I"
			addfile="$1"
			fragment=1
		elif [ "$1" = "--addfull" ]; then 
			shift   
			addarg0="-K -I"
			addfile="$1"
			fragment=-1
		elif [ "$1" = "--addlong" ]; then 
			shift   
			addarg0="-K -I"
			addfile="$1"
			fragment=-2
		elif [ "$1" = "--smoothing" ]; then 
			add2ndhalfarg=$add2ndhalfarg" -p "
		elif [ "$1" = "--keeplength" ]; then 
			add2ndhalfarg=$add2ndhalfarg" -Y "
		elif [ "$1" = "--mapout" ]; then 
			add2ndhalfarg=$add2ndhalfarg" -Z -Y "
		elif [ "$1" = "--mapoutfile" ]; then 
			shift
			add2ndhalfarg=$add2ndhalfarg" -Z -Y "
			mapoutfile="$1"
		elif [ "$1" = "--maxiterate" ]; then 
			shift   
			iterate=`expr "$1" - 0` 
			if ! expr "$1" : "[0-9]" > /dev/null ; then
				echo "Specify the number of iterations!" 1>&2
				exit
			fi
		elif [ "$1" = "--retree" ]; then 
			shift   
			cycle=`expr "$1" - 0`
			if ! expr "$1" : "[0-9]" > /dev/null ; then
				echo "Specify the number of tree rebuilding!" 1>&2
				exit
			fi
		elif [ "$1" = "--text" ]; then 
			sbstmodel=" -b -2 -a "
			f2clext="-E"
			seqtype="-P"
			fft=0
		elif [ "$1" = "--aamatrix" ]; then 
			shift   
			sbstmodel=" -b -1 "
			aamatrix="$1"
			if [ ! -e "$aamatrix" ]; then
				echo "Cannot open $aamatrix" 1>&2
				echo "" 1>&2
				exit
			fi
		elif [ "$1" = "--treein" ]; then 
			shift   
			treeinopt=" -U "
			treein=1
			treeinfile="$1"
			if [ ! -e "$treeinfile" ]; then
				echo "Cannot open $treeinfile" 1>&2
				echo "" 1>&2
				exit
			fi
		elif [ "$1" = "--pileup" ]; then 
			treeinopt=" -U "
			treein=1
			pileuporshuffle="p"
		elif [ "$1" = "--randomchain" ]; then 
			treeinopt=" -U "
			treein=1
			pileuporshuffle="s"
		elif [ "$1" = "--topin" ]; then 
			shift   
			treeinopt=" -V "
			treein=1
			treeinfile="$1"
			echo "The --topin option has been disabled." 1>&2
			echo "There was a bug in version < 6.530."   1>&2
			echo "This bug has not yet been fixed."      1>&2
			exit 1
		elif [ "$1" = "--memsavetree" ]; then
			treeinopt=" -U "
			treein=1
			pileuporshuffle="C"
		elif [ "$1" = "--memsavetreex" ]; then
			treeinopt=" -U "
			treein=1
			pileuporshuffle="c"
		elif [ "$1" = "--initialramusage" ]; then 
			shift   
			treeinopt=" -U "
			treein=1
			initialramusage="$1"
			pileuporshuffle="c"
		elif [ "$1" = "--kappa" ]; then 
			shift   
			kappa=" -k $1 "
			if ! expr "$1" : "[0-9]" > /dev/null ; then
				echo "Specify kappa value!" 1>&2
				exit
			fi
		elif [ "$1" = "--fmodel" ]; then 
			fmodel=" -a "
		elif [ "$1" = "--nwildcard" ]; then 
			nmodel=" -: "
		elif [ "$1" = "--nzero" ]; then 
			nmodel="  "
		elif [ "$1" = "--jtt" ]; then 
			shift   
			sbstmodel=" -j $1"
#			if ! expr "$1" : "[0-9]" > /dev/null ; then
#				echo "Specify pam value!" 1>&2
#				exit
#			fi
		elif [ "$1" = "--kimura" ]; then 
			shift   
			sbstmodel=" -j $1"
#			if ! expr "$1" : "[0-9]" > /dev/null ; then
#				echo "Specify pam value!" 1>&2
#				exit
#			fi
		elif [ "$1" = "--tm" ]; then 
			shift   
			sbstmodel=" -m $1"
#			if ! expr "$1" : "[0-9]" > /dev/null ; then
#				echo "Specify pam value!" 1>&2
#				exit
#			fi
		elif [ "$1" = "--bl" ]; then 
			shift   
			sbstmodel=" -b $1"
			if ! expr "$1" : "[0-9]" > /dev/null ; then
				echo "blosum $1?" 1>&2
				exit
			fi
		elif [ "$1" = "--weighti" ]; then
			shift   
			weighti="$1"
			if ! expr "$1" : "[0-9]" > /dev/null ; then
				echo "Specify weighti value!" 1>&2
				exit
			fi
		elif [ "$1" = "--weightr" ]; then
			shift   
			weightr="$1"
			if ! expr "$1" : "[0-9]" > /dev/null ; then
				echo "Specify weightr value!" 1>&2
				exit
			fi
		elif [ "$1" = "--weightm" ]; then
			shift   
			weightm="$1"
			if ! expr "$1" : "[0-9]" > /dev/null ; then
				echo "Specify weightm value!" 1>&2
				exit
			fi
		elif [ "$1" = "--rnaalifold" ]; then
			rnaalifold=1
		elif [ "$1" = "--mccaskill" ]; then
			mccaskill=1
			contrafold=0
			dafs=0
		elif [ "$1" = "--contrafold" ]; then
			mccaskill=0
			contrafold=1
			dafs=0
		elif [ "$1" = "--dafs" ]; then
			mccaskill=0
			contrafold=0
			dafs=1
		elif [ "$1" = "--ribosum" ]; then
			rnascoremtx=" -s "
		elif [ "$1" = "--op" ]; then 
			shift   
			gop="$1"
			if ! expr "$1" : "[0-9]" > /dev/null ; then
				echo "Specify op!" 1>&2
				exit
			fi
		elif [ "$1" = "--opdist" ]; then 
			shift   
			gopdist="$1"
			if ! expr "$1" : "[0-9]" > /dev/null ; then
				echo "Specify opdist!" 1>&2
				exit
			fi
			opdistspecified=1
		elif [ "$1" = "--allowshift" ]; then 
			allowshift=1
		elif [ "$1" = "--shiftpenalty" ]; then 
			shift   
			spfactor="$1"
			if ! expr "$1" : "[0-9]" > /dev/null ; then
				echo "Specify sf!" 1>&2
				exit
			fi
			shiftpenaltyspecified=1
		elif [ "$1" = "--ep" ]; then 
			shift   
#			aof="$1"
			tmpval="$1"
			aof=`awk "BEGIN{ print -1.0 * \"$tmpval\"}"`
			if ! expr "$aof" : "[0-9\-]" > /dev/null ; then
				printf "\nSpecify a number for ep, like --ep 0.1\n" 1>&2
				printf "'$1' cannot be interpreted as a number..\n\n" 1>&2
				exit
			fi
		elif [ "$1" = "--rop" ]; then 
			shift   
			rgop="$1"
# Atode check
		elif [ "$1" = "--rep" ]; then 
			shift   
			rgep="$1"
		elif [ "$1" = "--lop" ]; then 
			shift   
			lgop="$1"
		elif [ "$1" = "--LOP" ]; then 
			shift   
			LGOP="$1"
		elif [ "$1" = "--lep" ]; then 
			shift   
			laof="$1"
		elif [ "$1" = "--lexp" ]; then 
			shift   
			lexp="$1"
		elif [ "$1" = "--LEXP" ]; then 
			shift   
			LEXP="$1"
		elif [ "$1" = "--GEXP" ]; then 
			shift   
			GEXP="$1"
		elif [ "$1" = "--GOP" ]; then 
			shift   
			GGOP="$1"
		elif [ "$1" = "--gop" ]; then 
			shift   
			pggop="$1"
		elif [ "$1" = "--gep" ]; then 
			shift   
			pgaof="$1"
		elif [ "$1" = "--gexp" ]; then 
			shift   
			pgexp="$1"
		elif [ "$1" = "--laraparams" ]; then 
			shift   
			laraparams="$1"
		elif [ "$1" = "--corethr" ]; then 
			shift   
			corethr="$1"
		elif [ "$1" = "--corewin" ]; then 
			shift   
			corewin="$1"
		elif [ "$1" = "--strdir" ]; then
			shift
			strdir="$1"
		elif [ "$1" = "--pdbidlist" ]; then
			shift
			pdblist="$1"
			if [ ! -e "$pdblist" ]; then
				echo "Cannot open $pdblist" 1>&2
				echo "" 1>&2
				exit
			fi
		elif [ "$1" = "--pdbfilelist" ]; then
			shift
			ownlist="$1"
			if [ ! -e "$ownlist" ]; then
				echo "Cannot open $ownlist" 1>&2
				echo "" 1>&2
				exit
			fi
		elif [ "$1" = "--enrich" ]; then
			enrich=1
			enrichseq=1
			enrichstr=1
			seektarget=""
		elif [ "$1" = "--enrichseq" ]; then
			enrich=1
			enrichseq=1
			enrichstr=0
			seektarget="-seq"
		elif [ "$1" = "--enrichstr" ]; then
			enrich=1
			enrichseq=0
			enrichstr=1
			seektarget="-str"
		elif [ "$1" = "--seedtable" ]; then
			shift
			seedtable="y"
			seedtablefile="$1"
		elif [ "$1" = "--seed" ]; then
			shift
			seed="m"
			seedfiles="$seedfiles $1"
		elif [ "$1" = "--minimumweight" ]; then
			shift
			minimumweight="$1"
		elif [ "$1" = "--similaritylevel" ]; then
			shift
			similarityoffset="$1"
		elif [ "$1" = "--unalignlevel" ]; then
			shift
			unalignlevel="$1"
			unalignspecified=1
		elif [ "$1" = "--skipiterate" ]; then
			shift
			fixthreshold="$1"
		elif [ "$1" = "--bunkatsunashi" ]; then
			bunkatsuopt=" -B "
		elif [ "$1" = "--sp" ]; then
			sprigorous=1
		elif [ "$1" = "--focus" ]; then
			focusarg=" -= "
		elif [ "$1" = "--sparsepickup" ]; then
			shift
			npickup="$1"
		elif [ $progname = "fftns" -o  $progname = "nwns" ]; then
			if [ "$1" -gt 0 ]; then
				cycle=`expr "$1" - 0`
			fi
		else
			echo "Unknown option:  $1" 1>&2
			er=1;
#			exit 1;
		fi
		shift   
	done;


	echo "" 1>"$progressfile"

#	TMPFILE=/tmp/$progname.$$
	TMPFILE=`mktemp -dt $progname.XXXXXXXXXX`
	if [ $? -ne 0 ]; then
		echo "mktemp seems to be obsolete. Re-trying without -t" 1>&2
		TMPFILE=`mktemp -d /tmp/$progname.XXXXXXXXXX`
	fi	

	if [ $os = "cygwin" ]; then
		TMPFILE=`cygpath -w $TMPFILE`
	fi

	umask 077
#	mkdir  $TMPFILE  || er=1
	if [ $debug -eq 1 ]; then
		trap "tar cfvz debuginfo.tgz $TMPFILE; rm -rf $TMPFILE " 0
	else
		trap "rm -rf $TMPFILE" 0
	fi
	if [ $# -eq 1 ]; then
		if [ -r "$1" -o "$1" = - ]; then

			if [ -r "$addfile" ]; then
				printf '';
			else
				echo "$0": Cannot open "$addfile". 1>&2
				echo "" 1>&2
				exit 1;
			fi

			cat "$1"              | tr "\r" "\n" > $TMPFILE/infile 
			echo ""                             >> $TMPFILE/infile
			cat "$addfile"        | tr "\r" "\n" | grep -v "^$" >> $TMPFILE/infile
			cat "$addfile"        | tr "\r" "\n" | grep -v "^$" > $TMPFILE/_addfile
			cat "$aamatrix"       | tr "\r" "\n" | grep -v "^$" > $TMPFILE/_aamtx
			cat "$mergetable"     | tr "\r" "\n" | grep -v "^$" > $TMPFILE/_subalignmentstable
			cat "$treeinfile"     | tr "\r" "\n" | grep -v "^$" > $TMPFILE/_guidetree
			cat "$seedtablefile"  | tr "\r" "\n" | grep -v "^$" > $TMPFILE/_seedtablefile
			cat "$laraparams"     | tr "\r" "\n" | grep -v "^$" > $TMPFILE/_lara.params
			cat "$pdblist"        | tr "\r" "\n" | grep -v "^$" > $TMPFILE/pdblist
			cat "$ownlist"        | tr "\r" "\n" | grep -v "^$" > $TMPFILE/ownlist

#			echo $seedfiles
			infilename="$1"
			seedfilesintmp="/dev/null"
			seednseq="0"
			set $seedfiles > /dev/null
			while [ $# -gt 1 ];
			do
				shift
				if [ -r "$1" ]; then
					cat "$1" | tr "\r" "\n" >  $TMPFILE/seed$#
				else
					echo "$0": Cannot open "$1". 1>&2
					echo "" 1>&2
					exit 1;
				fi
				seednseq=$seednseq" "`grep -c '^[>|=]' $TMPFILE/seed$#`
				seedfilesintmp=$seedfilesintmp" "seed$#
			done
#			ls $TMPFILE
#			echo $seedfilesintmp
#			echo $seednseq


		else
			echo "$0": Cannot open "$1". 1>&2
			echo "" 1>&2
			er=1
#			exit 1;
		fi
	else
#		echo '$#'"=$#" 1>&2
		er=1
	fi



	if [ $numthreads -lt 0 ]; then
		if [ $os = "linux" ]; then
			nlogicalcore=`cat /proc/cpuinfo | grep "^processor" | uniq | wc -l`
			ncoresinacpu=`cat /proc/cpuinfo | grep 'cpu cores' | uniq | awk '{print $4}'`
			nphysicalcpu=`cat /proc/cpuinfo | grep 'physical id' | sort | uniq | wc -l`
			if [ $nlogicalcore -eq 0 ]; then
				echo "Cannot get the number of processors from /proc/cpuinfo" 1>>"$progressfile"
				exit 1
			fi
			if [ ${#ncoresinacpu} -gt 0 -a $nphysicalcpu -gt 0 ]; then
				numthreads=`expr $ncoresinacpu '*' $nphysicalcpu`
#				if [ $nlogicalcore -gt $numthreads ]; then    # Hyperthreading
#					numthreads=`expr $numthreads '+' 1`
#				fi
			else
				numthreads=$nlogicalcore
			fi
		elif [ $os = "darwin" ]; then
			numthreads=`sysctl -n hw.physicalcpu`
			if [ -z $numthreads ]; then
				echo "Cannot get the number of physical cores from sysctl" 1>>"$progressfile"
				exit 1
			fi
#			nlogicalcore=`sysctl -n hw.logicalcpu`
#			if [ $nlogicalcore -gt $numthreads ]; then    # Hyperthreading
#				numthreads=`expr $numthreads '+' 1`
#			fi
		elif [ $os = "mingw" -o $os = "cygwin" ]; then
			numthreads=`wmic cpu get NumberOfCores | head -2 | tail -1 | awk '{print $1}'`
		else
			echo "Cannot count the number of physical cores." 1>>"$progressfile"
			exit 1
		fi
		echo "OS = "$os 1>>"$progressfile"
		echo "The number of physical cores = " $numthreads 1>>"$progressfile"
	fi

	if [ $numthreadstb -lt 0 ]; then
		numthreadstb=$numthreads
	fi

	if [ $numthreadsit -lt 0 ]; then
		if [ $numthreads -lt 11 ]; then
			numthreadsit=$numthreads
		else
			numthreadsit=10
		fi
	fi

	if [ $numthreadsit -eq 0 -a $parallelizationstrategy = "BESTFIRST" ]; then
		echo 'Impossible' 1>&2;
		exit 1;
	fi

	if [ "$addarg0" != " " ]; then
		iterate=0  # 2013/03/23
		"$prefix/countlen" < $TMPFILE/_addfile > $TMPFILE/addsize 2>>"$progressfile"
		nadd=`awk '{print $1}' $TMPFILE/addsize`
		if [ $nadd -eq "0" ]; then
			echo Check $addfile 1>&2
			exit 1;
		fi
		if [ $seed != "x" -o $seedtable != "x" ]; then
			echo 'Impossible' 1>&2;
			echo 'Use either ONE of --seed, --seedtable, --addprofile and --add.' 1>&2
			exit 1;
		fi
	else
		nadd="0"
	fi

	if [ $auto -eq 1 ]; then
		"$prefix/countlen" < $TMPFILE/infile > $TMPFILE/size 2>>"$progressfile"
		nseq=`awk '{print $1}' $TMPFILE/size`
		nlen=`awk '{print $3}' $TMPFILE/size`

		if [ $nlen -lt 3000 -a $nseq -lt 100 ]; then
			distance="local"
			iterate=1000
			cycle=1
		elif [ $nlen -lt 1000 -a $nseq -lt 200 ]; then
			distance="local"
			iterate=2
			cycle=1
		elif [ $nlen -lt 10000 -a $nseq -lt 500 ]; then
			distance="ktuples"
			iterate=2
			cycle=2
		elif [ $nseq -lt 50000 ]; then  # changed from 10000 2014/Oct/4
			distance="ktuples"
			iterate=0
			cycle=2
		elif [ $nseq -lt 90000 ]; then  # changed from 30000 2014/Oct/4
			distance="ktuples"
			iterate=0
			cycle=1
		elif [ $nlen -lt 3000 ]; then
			distance="parttree"
			partdist="localalign"
			algopt=" "
			algoptit=" "
#			algspecified=1
			cycle=1
		else
			distance="parttree"
			partdist="ktuples"
			algopt=" "
			algoptit=" "
#			algspecified=1
			cycle=1
		fi


#		if [ $nlen -lt 3000 -a $nseq -lt 100 ]; then
#			distance="local"
#			iterate=1000
#			cycle=1
#		elif [ $nlen -lt 1000 -a $nseq -lt 200 ]; then
#			distance="local"
#			iterate=2
#			cycle=1
#		elif [ $nlen -lt 10000 -a $nseq -lt 500 ]; then
#			distance="ktuples"
#			iterate=2
#			cycle=2
#		elif [ $nseq -lt 200000 ]; then
#			distance="ktuples"
#			iterate=0
#			treeinopt=" -U "
#			treein=1
#			pileuporshuffle="a"
#		elif [ $nlen -lt 3000 ]; then
#			distance="parttree"
#			partdist="localalign"
#			algopt=" "
#			algoptit=" "
##			algspecified=1
#			cycle=1
#		else
#			distance="parttree"
#			partdist="ktuples"
#			algopt=" "
#			algoptit=" "
##			algspecified=1
#			cycle=1
#		fi


		if [ $fragment -ne 0 ]; then
			norg=`expr $nseq '-' $nadd`
			npair=`expr $norg '*' $nadd`
			echo "nadd = " $nadd               1>>"$progressfile"
			echo "npair = " $npair             1>>"$progressfile"
			echo "nseq = " $nseq               1>>"$progressfile"
			echo "nlen = " $nlen               1>>"$progressfile"
# nagasa check!
#
			if [ $npair -gt 10000000 -o $nlen -gt 500000 ]; then  # 2015/Jun
				distance="ktuples"
				echo "use ktuples, size=$tuplesize!" 1>>"$progressfile"
			elif [ $npair -gt 3000000 -o $nlen -gt 100000 ]; then  # 2015/Jun
				distance="multi"
				weighti="0.0"
				echo "use multipair, weighti=0.0!" 1>>"$progressfile"
			else
				distance="multi"
				echo "use multipair, weighti=$weighti!" 1>>"$progressfile"
			fi
			pairspecified=1
		fi
	fi

	if [ `awk "BEGIN {print( 0.0+\"$sueff\" < 0.0 || 0.0+\"$sueff\" > 1.0 )}"` -gt 0 ]; then
		printf "\n%s\n\n" "The argument of --mixedlinkage must be between 0.0 and 1.0" 1>>"$progressfile"
		exit 1;
	fi

	if [ $allowshift -eq 1 ]; then
		if [ $unalignspecified -ne 1 ]; then
			unalignlevel="0.8"
		fi
		if [ $shiftpenaltyspecified -ne 1 ]; then
			spfactor="2.00"
		fi
	fi

	if [ $opdistspecified -ne 1 ]; then
		gopdist=$gop
	fi

	if [ $unalignlevel != "0.0" -o `awk "BEGIN {print( 0.0+\"$spfactor\" < 100.0 )}"` -gt 0 ]; then
		termgapopt=" "
		if [ $distance = "localgenaf" ]; then
			printf "\n%s\n" "The combination of --allowshift and --genafpair (E-INS-i/-1) is not supported." 1>>"$progressfile"
			printf "%s\n" "Instead, please try --allowshift --globalpair (G-INS-i/-1 in the web version)," 1>>"$progressfile"
			printf "%s\n\n" "which covers the situation for --genafpair (E-INS-i/-1), too." 1>>"$progressfile"
			exit 1;
		fi
		if [ $distance != "global" -o `awk "BEGIN {print( 0.0+\"$weighti\" < 1.0 )}"` -gt 0 ]; then
			printf "\n%s\n\n" "At present, --unalignlevel # or --allowshift is supported only with the --globalpair option." 1>>"$progressfile"
			exit 1;
		fi
		if [ $fragment -ne 0 ]; then
			printf "\n%s\n\n" "At present, --unalignlevel # or --allowshift is not supported with the --addfragments option." 1>>"$progressfile"
			exit 1;
		fi
	fi

	if [ `awk "BEGIN {print( 0.0+\"$spfactor\" < 1.0 )}"` -gt 0 ]; then
			printf "\n%s\n" "shiftpenalty must be >1." 1>>"$progressfile"
			exit 1;
	fi

	if [ `awk "BEGIN {print( 0.0+\"$fixthreshold\" < 0.0 )}"` -gt 0 ]; then
		printf "\n%s\n\n" "The 'fix' parameter must be >= 0.0" 1>>"$progressfile"
		exit 1;
	fi

	if [ `awk "BEGIN {print( 0.0+\"$unalignlevel\" < 0.0 || 0.0+\"$unalignlevel\" > 1.0 )}"` -gt 0 ]; then
		printf "\n%s\n\n" "The 'unalignlevel' parameter must be between 0.0 and 1.0" 1>>"$progressfile"
		exit 1;
	fi
	if [ `awk "BEGIN {print( 0.0+\"$unalignlevel\" > 0.0 )}"` -gt 0 ]; then
		laof="0"
		lexp="0"
		pgaof="0"
		pgexp="0"
		LEXP="0"
		GEXP="0"
		termgapopt=" "
#		if [ $auto -eq 1 -o $fragment -ne 0 -o $iterate -gt 0 ]; then
		if [ $fragment -ne 0 ]; then
			printf "\n%s\n\n" "At present, the 'unalignlevel > 0' mode is not supported with the --addfragments option." 1>>"$progressfile"
			exit 1;
		fi
		if [ $distance = "parttree" ]; then
			printf "\n%s\n\n" "At present, the 'unalignlevel > 0' mode is not supported in the (dp)parttree option." 1>>"$progressfile"
			exit 1;
		fi
		if [ $distance = "localgenaf" ]; then
			printf "\n%s\n" "The --genafpair is not supported in the 'unalignlevel > 0' mode." 1>>"$progressfile"
			printf "%s\n" "Instead, please try --unalignlevel xx --globalpair," 1>>"$progressfile"
			printf "%s\n\n" "which covers the situation for --genafpair (E-INS-i), too." 1>>"$progressfile"
			exit 1;
		fi
#		if [ $distance != "ktuples" -a `awk "BEGIN {print( 0.0+\"$weighti\" > 0.0 )}"` -gt 0 -a $iterate -gt 0 ]; then
#			printf "\n%s\n\n" "Please add --weighti 0.0, for now." 1>>"$progressfile"
#			exit 1;
#		fi
	fi

	if [ `awk "BEGIN {print( 0.0+\"$similarityoffset\" != 0.0 && 0.0+\"$unalignlevel\" != 0.0 )}"` -gt 0 ]; then
		printf "\n%s\n\n" "Do not simultaneously specify --similaritylevel and --unalignlevel" 1>>"$progressfile"
		exit 1;
	fi

	if [ `awk "BEGIN {print( 0.0+\"$similarityoffset\" < -1.0 || 0.0+\"$similarityoffset\" > 1.0 )}"` -gt 0 ]; then
		printf "\n%s\n\n" "Similarity must be between -1.0 and +1.0" 1>>"$progressfile"
		exit 1;
	fi
	aof=`awk "BEGIN{print 0.0 + \"$similarityoffset\" + $aof}"`
	laof=`awk "BEGIN{print 0.0 + \"$similarityoffset\" + $laof}"`
	pgaof=`awk "BEGIN{print 0.0 + \"$similarityoffset\" + $pgaof}"`


	if [ $parallelizationstrategy = "BESTFIRST" -o  $parallelizationstrategy = "BAATARI0" ]; then
		iteratelimit=254
	else
		iteratelimit=16
	fi
	if [ $iterate -gt $iteratelimit ]; then    #??
		iterate=$iteratelimit
	fi

	if [ $rnaalifold -eq 1 ]; then
		rnaopt=" -e $rgep -o $rgop -c $weightm -r $weightr -R $rnascoremtx "
#		rnaoptit=" -o $rgop -BT -c $weightm -r $weightr -R "
		rnaoptit=" -o $rgop -F -c $weightm -r $weightr -R "
	elif [ $mccaskill -eq 1 -o $dafs -eq 1 -o $contrafold -eq 1 ]; then
		rnaopt=" -o $rgop -c $weightm -r $weightr "
#		rnaoptit=" -e $rgep -o $rgop -BT -c $weightm -r $weightr $rnascoremtx "
		rnaoptit=" -e $rgep -o $rgop -F -c $weightm -r $weightr $rnascoremtx "
	else
		rnaopt="  "
		rnaoptit=" -F "
	fi

#	if [ $algspecified -eq 0 ]; then
#		if [ $distance = "parttree" ]; then 
#			algopt=" -Q "
#			algoptit=" "
#		else
#			algopt=" "
#			algoptit=" "
#		fi
#	fi

	if [ $sprigorous -eq 1 ]; then
		algopt=" -@ "
		if [ $iterate -gt 0 ]; then
			if [ $numthreadsit -eq 0 ]; then
				algoptit=" -@ -B -Z -z 1000 "
			else
				echo "" 1>>"$progressfile"
				echo "At present, the combination of --sp and iterative refinement is supported only in a single thread." 1>>"$progressfile"
				echo "Please try \"--thread -1 --threadit 0\", which runs the iterative refinment calculation on a single thread." 1>>"$progressfile"
				echo "" 1>>"$progressfile"
				exit 1;
#				algoptit=" -@ -B -z 1000 "
			fi
		fi
		termgapopt=" "
		fft=0
		memopt=" -N "
	fi

	model="$sbstmodel $kappa $fmodel $nmodel"

	if [ $er -eq 1 ]; then
		echo "------------------------------------------------------------------------------" 1>&2
		echo "  MAFFT" $version 1>&2
#		echo "" 1>&2
#		echo "  Input format: fasta" 1>&2
#		echo ""  1>&2
#		echo "  Usage: `basename $0` [options] inputfile > outputfile" 1>&2
	    echo "  http://mafft.cbrc.jp/alignment/software/" 1>&2
		echo "  MBE 30:772-780 (2013), NAR 30:3059-3066 (2002)"        1>&2
#		echo "------------------------------------------------------------------------------" 1>&2
#		echo "  % mafft in > out" 1>&2
		echo "------------------------------------------------------------------------------" 1>&2
#		echo ""  1>&2
		echo "High speed:" 1>&2
		echo "  % mafft in > out" 1>&2
		echo "  % mafft --retree 1 in > out (fast)" 1>&2
		echo "" 1>&2
		echo "High accuracy (for <~200 sequences x <~2,000 aa/nt):" 1>&2
		echo "  % mafft --maxiterate 1000 --localpair  in > out (% linsi in > out is also ok)" 1>&2
		echo "  % mafft --maxiterate 1000 --genafpair  in > out (% einsi in > out)" 1>&2
		echo "  % mafft --maxiterate 1000 --globalpair in > out (% ginsi in > out)" 1>&2
		echo "" 1>&2
		echo "If unsure which option to use:" 1>&2
		echo "  % mafft --auto in > out" 1>&2
		echo "" 1>&2
#		echo "Other options:" 1>&2
		echo "--op # :         Gap opening penalty, default: 1.53" 1>&2
		echo "--ep # :         Offset (works like gap extension penalty), default: 0.0" 1>&2
		echo "--maxiterate # : Maximum number of iterative refinement, default: 0" 1>&2
		echo "--clustalout :   Output: clustal format, default: fasta" 1>&2
		echo "--reorder :      Outorder: aligned, default: input order" 1>&2
		echo "--quiet :        Do not report progress" 1>&2
		echo "--thread # :     Number of threads (if unsure, --thread -1)" 1>&2
#		echo "" 1>&2
#		echo " % mafft --maxiterate 1000 --localpair in > out (L-INS-i)" 1>&2
#		echo " most accurate in many cases, assumes only one alignable domain" 1>&2 
#		echo "" 1>&2
#		echo " % mafft --maxiterate 1000 --genafpair in > out (E-INS-i)" 1>&2
#		echo " works well if many unalignable residues exist between alignable domains" 1>&2
#		echo "" 1>&2
#		echo " % mafft --maxiterate 1000 --globalpair in > out (G-INS-i)" 1>&2
#		echo " suitable for globally alignable sequences            " 1>&2
#		echo "" 1>&2
#		echo " % mafft --maxiterate 1000 in > out (FFT-NS-i)" 1>&2
#		echo " accurate and slow, iterative refinement method      " 1>&2
#		echo "" 1>&2
#		echo "If the input sequences are long (~1,000,000nt)," 1>&2
#		echo " % mafft --retree 1 --memsave --fft in > out (FFT-NS-1-memsave, new in v5.8)" 1>&2
#		echo "" 1>&2
#		echo "If many (~5,000) sequences are to be aligned," 1>&2
#		echo "" 1>&2
#		echo " % mafft --retree 1 [--memsave] --nofft in > out (NW-NS-1, new in v5.8)" 1>&2
#		echo "" 1>&2
#		echo " --localpair :      All pairwise local alignment information is included"  1>&2
#		echo "                    to the objective function, default: off"  1>&2
#		echo " --globalpair :     All pairwise global alignment information is included"  1>&2
#		echo "                    to the objective function, default: off"  1>&2
#		echo " --op # :           Gap opening penalty, default: $defaultgop " 1>&2
#		echo " --ep # :           Offset (works like gap extension penalty), default: $defaultaof " 1>&2
#		echo " --bl #, --jtt # :  Scoring matrix, default: BLOSUM62" 1>&2
#		echo "                    Alternatives are BLOSUM (--bl) 30, 45, 62, 80, " 1>&2
#		echo "                    or JTT (--jtt) # PAM. " 1>&2
#		echo " --nuc or --amino : Sequence type, default: auto" 1>&2
#		echo " --retree # :       The number of tree building in progressive method " 1>&2
#		echo "                    (see the paper for detail), default: $defaultcycle " 1>&2
#		echo " --maxiterate # :   Maximum number of iterative refinement, default: $defaultiterate " 1>&2
#		if [ $defaultfft -eq 1 ]; then
#			echo " --fft or --nofft:  FFT is enabled or disabled, default: enabled" 1>&2
#		else
#			echo " --fft or --nofft:  FFT is enabled or disabled, default: disabled" 1>&2
#		fi
#		echo " --memsave:         Memory saving mode" 1>&2
#		echo "                    (for long genomic sequences), default: off" 1>&2
#		echo " --clustalout :     Output: clustal format, default: fasta" 1>&2
#		echo " --reorder :        Outorder: aligned, default: input order" 1>&2
#		echo " --quiet :          Do not report progress" 1>&2
#		echo "-----------------------------------------------------------------------------" 1>&2
		exit 1; 
	fi
	if [ $sw -eq 1 ]; then
		swopt=" -A "
	else
		swopt=" "
	fi

	if [ $distance = "fasta" -o $partdist = "fasta" ]; then
		if [ ! "$FASTA_4_MAFFT" ]; then
			FASTA_4_MAFFT=`which fasta34`
		fi

		if [ ! -x "$FASTA_4_MAFFT" ]; then
			echo ""       1>&2
			echo "== Install FASTA ========================================================" 1>&2
			echo "This option requires the fasta34 program (FASTA version x.xx or higher)"   1>&2
			echo "installed in your PATH.  If you have the fasta34 program but have renamed" 1>&2
			echo "(like /usr/local/bin/myfasta), set the FASTA_4_MAFFT environment variable" 1>&2
			echo "to point your fasta34 (like setenv FASTA_4_MAFFT /usr/local/bin/myfasta)." 1>&2
			echo "=========================================================================" 1>&2
			echo "" 1>&2
			exit 1
		fi
	fi
	if [ $distance = "last" -o $distance = "lastmulti" ]; then
		if [ ! -x "$prefix/lastal" -o ! -x "$prefix/lastdb" ]; then
			echo ""       1>&2
			echo "== Install LAST ============================================================" 1>&2
			echo "LAST (Kielbasa, Wan, Sato, Horton, Frith 2011 Genome Res. 21:487) is required." 1>&2
			echo "http://last.cbrc.jp/"                                                       1>&2
			echo "http://mafft.cbrc.jp/alignment/software/xxxxxxx.html "                      1>&2
			echo "============================================================================" 1>&2
			echo "" 1>&2
			exit 1
		fi
	fi
	if [ $distance = "lara" -o $distance = "slara" ]; then
		if [ ! -x "$prefix/mafft_lara" ]; then
			echo ""       1>&2
			echo "== Install LaRA =========================================================" 1>&2
			echo "This option requires LaRA (Bauer et al. http://www.planet-lisa.net/)."     1>&2
			echo "The executable have to be renamed to 'mafft_lara' and installed into "     1>&2
			echo "the $prefix directory. "                                                   1>&2
			echo "A configuration file of LaRA also have to be given"                        1>&2
			echo "mafft-xinsi --larapair --laraparams parameter_file"                        1>&2
			echo "mafft-xinsi --slarapair --laraparams parameter_file"                       1>&2
			echo "=========================================================================" 1>&2
			echo "" 1>&2
			exit 1
		fi
		if [ ! -s "$laraparams" ]; then
			echo ""       1>&2
			echo "== Configure LaRA =======================================================" 1>&2
			echo "A configuration file of LaRA have to be given"                             1>&2
			echo "mafft-xinsi --larapair --laraparams parameter_file"                        1>&2
			echo "mafft-xinsi --slarapair --laraparams parameter_file"                       1>&2
			echo "=========================================================================" 1>&2
			echo "" 1>&2
			exit 1
		fi
	fi
	if [ $distance = "foldalignlocal" -o $distance = "foldalignglobal" ]; then
	if [ ! -x "$prefix/foldalign210" ]; then
			echo ""       1>&2
			echo "== Install FOLDALIGN ====================================================" 1>&2
			echo "This option requires FOLDALIGN (Havgaard et al. http://foldalign.ku.dk/)." 1>&2
			echo "The executable have to be renamed to 'foldalign210' and installed into "   1>&2
			echo "the $prefix directory. "                                                   1>&2
			echo "=========================================================================" 1>&2
			echo "" 1>&2
			exit 1
		fi
	fi
	if [ $distance = "scarna" -o $mccaskill -eq 1 ]; then
		if [ ! -x "$prefix/mxscarnamod" ]; then
			echo ""       1>&2
			echo "== Install MXSCARNA ======================================================" 1>&2
			echo "MXSCARNA (Tabei et al. BMC Bioinformatics 2008 9:33) is required."          1>&2
			echo "Please 'make' at the 'extensions' directory of the MAFFT source package,"   1>&2
			echo "which contains the modified version of MXSCARNA."                           1>&2
			echo "http://mafft.cbrc.jp/alignment/software/source.html "                1>&2
			echo "==========================================================================" 1>&2
			echo "" 1>&2
			exit 1
		fi
	fi
	if [ $distance = "dafs" -o $dafs -eq 1 ]; then
		if [ ! -x "$prefix/dafs" ]; then
			echo ""       1>&2
			echo "== Install DAFS===========================================================" 1>&2
			echo "DAFS (Sato et al. Journal 2012 issue:page) is required."                    1>&2
			echo "http://www.ncrna.org/ "                                                     1>&2
			echo "==========================================================================" 1>&2
			echo "" 1>&2
			exit 1
		fi
	fi
	if [ $contrafold -eq 1 ]; then
		if [ ! -x "$prefix/contrafold" ]; then
			echo ""       1>&2
			echo "== Install CONTRAfold ===================================================" 1>&2
			echo "This option requires CONTRAfold"                                           1>&2
			echo "(Do et al. http://contra.stanford.edu/contrafold/)."                       1>&2
			echo "The executable 'contrafold' have to be installed into "                    1>&2
			echo "the $prefix directory. "                                                   1>&2
			echo "=========================================================================" 1>&2
			echo "" 1>&2
			exit 1
		fi
	fi

#old
#	if [ $treeout -eq 1 ]; then
#		parttreeoutopt="-t"
#		if [ $cycle -eq 0 ]; then
#			treeoutopt="-t -T"
#			groupsize=1
#			iterate=0 
#			if [ $distance = "global" -o $distance = "local" -o $distance = "localgenaf" -o $distance = "globalgenaf" ]; then
#				distance="distonly"
#			fi
#		else
#			treeoutopt="-t"
#		fi
#	else
#		parttreeoutopt=" "
#		if [ $cycle -eq 0 ]; then
#			treeoutopt="-t -T"
#			iterate=0 
#			if [ $distance = "global" -o $distance = "local" -o $distance = "localgenaf" -o $distance = "globalgenaf" ]; then
#				distance="distonly"
#			fi
#		else
#			treeoutopt=" "
#		fi
#	fi

#new
	if [ $cycle -eq 0 ]; then
		treeoutopt="-t -T"
		iterate=0 
		weighti="0.0"  # 2016Jul31, tbfast.c kara idou
#		if [ $distance = "global" -o $distance = "local" -o $distance = "localgenaf" -o $distance = "globalgenaf" ]; then # 2012/04, localpair --> local alignment distance
#		if [ $distance = "global" ]; then
#			distance="distonly"
#		fi
		if [ $treeout -eq 1 ]; then
			parttreeoutopt="-t"
			groupsize=1
		else
			parttreeoutopt=" "
		fi
		if [ $distout -eq 1 ]; then
			distoutopt="-y -T"
			if [ $treeout -eq 0 ]; then
				treeoutopt=""
			fi
		fi
	else
		if [ $treeout -eq 1 ]; then
			parttreeoutopt="-t"
			treeoutopt="-t"
		else
			parttreeoutopt=" "
			treeoutopt=" "
		fi
		if [ $distout -eq 1 ]; then
			distoutopt="-y"
		fi
	fi
#

	formatcheck=`grep -c '^[[:blank:]]\+>' $TMPFILE/infile | head -1 `
	if [ $formatcheck -gt 0 ]; then
		echo "The first character of a description line must be " 1>&2
		echo "the greater-than (>) symbol, not a blank."           1>&2
		echo "Please check the format around the following line(s):"  1>&2
		grep -n '^[[:blank:]]\+>' $TMPFILE/infile  1>&2
		exit 1
	fi

	nseq=`grep -c '^[>|=]' $TMPFILE/infile | head -1 ` 
	if [ $nseq -eq 2 ]; then
		cycle=1
	fi
	if [ $cycle -gt 3 ]; then
		cycle=3
	fi

	if [ $nseq -gt 60000 -a $iterate -gt 1 ]; then   # 2014/Oct/22, test
		echo "Too many sequences to perform iterative refinement!" 1>&2
		echo "Please use a progressive method." 1>&2
		exit 1
	fi
	if [ $distance = "lastmulti" -o $distance = "multi" ]; then
		if [ $fragment -eq 0 ]; then
			echo 'Specify --addfragments too' 1>&2
			exit 1
		fi
	fi

	if [ $fragment -ne 0 ]; then
		if [ $pairspecified -eq 0 ]; then
			distance="multi"
		fi
		if [ $distance != "multi" -a $distance != "hybrid" -a $distance != "lastmulti" -a $distance != "local" -a $distance != "last" -a  $distance != "ktuples" -a $distance != "ktuplesmulti" ]; then
			echo 'Specify --multipair, --lastmultipair, --lastpair, --localpair, --6merpair, --6mermultipair or --hybridpair' 1>&2
			exit 1
		fi
	fi

	if [ "$memopt" = " -M -B " -a "$distance" != "ktuples" ]; then
		echo "Impossible" 1>&2
		exit 1
	fi

	if [ $distance = "parttree" ]; then
		if [ $mergetable != "/dev/null" ]; then
			echo "The combination of (dp)parttree and merge is Impossible.  " 1>&2
			exit 1
		fi
		if [ $addfile != "/dev/null" ]; then
			echo "The combination of (dp)parttree and add(fragments) is Impossible.  " 1>&2
			exit 1
		fi
		if [ $seed != "x" -o $seedtable != "x" ]; then
			echo "Impossible" 1>&2
			exit 1
		fi
		if [ $iterate -gt 1 ]; then
			echo "Impossible" 1>&2
			exit 1
		fi
		if [ $outorder = "aligned" ]; then
			outorder="input"
		fi
		outorder="input"   # partorder ga kiku
		if [ $partdist = "localalign" ]; then
			splitopt=" -U "    # -U -l -> fast 
			cycle=1
		elif [ $partdist = "fasta" ]; then
			splitopt=" -S "
			cycle=1
		else
			splitopt="  "
		fi
	fi


	if [ \( $distance = "ktuples" -o $distance = "ktuplesmulti" \) -a \( $seed = "x" -a $seedtable = "x" -a $ownlist = "/dev/null" -a $pdblist = "/dev/null" -a $enrichstr -eq 0 \) ]; then
		localparam=""
		weighti="0.0"
	elif [ \( $distance = "ktuples" -o $distance = "ktuplesmulti" \) -a \( $seed != "x" -o $seedtable != "x" -o $ownlist != "/dev/null" -o $pdblist != "/dev/null" -o $enrichstr -eq 1 \) ]; then
		if [ $cycle -lt 2 ]; then
			cycle=2                # disttbfast ha seed hi-taiou #  chuui 2014Aug21
		fi
		if [ $iterate -lt 2 ]; then
			echo "############################################################################" 1>&2
			echo "# Warning:" 1>&2
			echo "#   Progressive alignment method is incompatible with the --seed option." 1>&2
			echo "#   Automatically switched to the iterative refinement method." 1>&2
			echo "#   " 1>&2
			echo "# Also consider using the '--add' option, which is compatible with" 1>&2
			echo "#   the progressive method and FASTER than the '--seed' option." 1>&2
			echo "#   Usage is:" 1>&2
			echo "#   % mafft --add newSequences existingAlignment > output" 1>&2
			echo "############################################################################" 1>&2
			iterate=2
		fi
		localparam="-l "$weighti
	elif [ $distance = "parttree" ]; then
		localparam=""
		weighti="0.0"
		if [ $groupsize -gt -1 ]; then
			cycle=1
		fi
	else
		localparam="-B -l "$weighti   # weighti=0 demo bunkatsu nashi
		if [ $cycle -gt 1 ]; then  # 09/01/08
			cycle=1
		fi
	fi


	if [ $distance = "localgenaf" -o $distance = "globalgenaf" ]; then
		aof="0.000"
		if [ $oldgenafparam -ne 1 ]; then
			laof="0.0"
			lexp="0.0"
#			LEXP="0.0" # default = 0.0
			usenaivepairscore="-Z"
		fi
	fi


#	if [ $nseq -gt 5000 ]; then
#		fft=0
#	fi
	if [ $forcefft -eq 1 ]; then
		param_fft=" -G "
		fft=1
	elif [ $fft -eq 1 ]; then
		param_fft=" -F "
	else
		param_fft=" "
	fi

	if [ $seed != "x" -a $seedtable != "x" ]; then
			echo 'Use either one of seedtable and seed.  Not both.' 1>&2
			exit 1
	fi
	if [ $f2clext = "-E" -a $anysymbol -gt 0 ]; then
			echo '' 1>&2
			echo 'At present, the combination of --text and ( --anysymbol or --preservecase ) is impossible.' 1>&2
			echo '' 1>&2
			exit 1
	fi

	if [ $f2clext = "-E" -a $aamatrix != "/dev/null" ]; then
			echo '' 1>&2
			echo 'At present, the combination of --text and (--aamatrix) is impossible.' 1>&2
			echo '' 1>&2
			exit 1
	fi

	if [ $treein -eq 1 ]; then
#		if [ $iterate -gt 0 ]; then
#			echo 'Not supported yet.' 1>&2
#			exit 1
#		fi
		if [ ! -s $TMPFILE/_guidetree ]; then
			if [ $distance != "ktuples" ]; then
				echo "Not supported yet"  1>>"$progressfile"
				exit 1
			fi
			if [ $pileuporshuffle = "p" ]; then
				echo "pileup" > $TMPFILE/_guidetree
#				weightopt=" -u " -> disttbfast.c?
#				numthreadstb=0 -> disttbfast.c
				cycle=1       #  disttbfast. shitei
			elif [ $pileuporshuffle = "s" ]; then
				echo "shuffle $randomseed" > $TMPFILE/_guidetree
#				numthreadstb=0 -> disttbfast.c
#				weightopt=" -u " -> disttbfast.c?
				cycle=1       #  disttbfast.c dem shitei
			elif [ $pileuporshuffle = "C" ]; then
				echo "very compact" > $TMPFILE/_guidetree
			elif [ $pileuporshuffle = "c" ]; then
				echo "compact " "$initialramusage" > $TMPFILE/_guidetree
			elif [ $pileuporshuffle = "a" ]; then
				echo "auto $randomseed 200" > $TMPFILE/_guidetree
			fi
		fi
	fi

	if [ $nadd -gt "0" ]; then
		if [ $fragment -eq "1" ]; then
			addarg="$addarg0 $nadd -g -0.01"
			addsinglearg=""
			cycle=1 # chuui 2014Aug25
		elif [ $fragment -eq "-1" ]; then
			addarg="$addarg0 $nadd"
			addsinglearg="-V"       # allowlongadds, 2014/04/02
			cycle=1 # chuui 2014Aug25
		elif [ $fragment -eq "-2" ]; then
			addarg="$addarg0 $nadd"
			addsinglearg="-V"       # allowlongadds + smoothing
			add2ndhalfarg=$add2ndhalfarg" -p "
			cycle=1 # chuui 2014Aug25
			usenaivepairscore="-Z" # 2015Jun01
			laof=0.0               # 2015Jun01
			lexp=0.0               # 2015Jun01
		else
			addarg="$addarg0 $nadd"
			addsinglearg=""
		fi

#		cycle=1 # chuui 2014Aug19
		iterate=0
#		treealg=" -q "  ## 2012/01/24  ## removed 2012/02/06
	fi


	if [ -z "$localparam" -a $fragment -eq 0 -a $distance != "parttree" ]; then
#		echo "use disttbfast"
#		echo cycle = $cycle
		cycletbfast=1            # tbfast wo jikkou shinai
		cycledisttbfast=$cycle   # disttbfast ni -E cycle wo watasu
		if [ $cycledisttbfast -eq 0 ]; then # --treeout de tsukau
			cycledisttbfast=1
		fi
	else
#		echo "use tbfast"
#		echo cycle = $cycle
		cycletbfast=$cycle       # 1 ijou nara jikkou
		cycledisttbfast=1        # disttbfast ha ikkai dake
	fi

#	echo localparam=
#	echo $localparam
#	echo cycletbfast=
#	echo $cycletbfast
#	echo cycledisttbfast=
#	echo $cycledisttbfast

#exit

	if [ $adjustdirection -gt 0 -a $seed != "x" ]; then
			echo '' 1>&2
			echo 'The combination of --adjustdirection(accurately) and --seed is not supported.' 1>&2
			echo '' 1>&2
			exit 1
	fi


	if [ $mccaskill -eq 1 -o $dafs -eq 1 -o $rnaalifold -eq 1 -o $contrafold -eq 1 ]; then
		if [ $distance = "ktuples" ]; then
			echo 'Not supported.' 1>&2
			echo 'Please add --globalpair, --localpair, --scarnapair, --dafspair' 1>&2
			echo '--larapair, --slarapair, --foldalignlocalpair or --foldalignglobalpair' 1>&2
			exit 1
		fi
		if [ $f2clext = "-E" ]; then
				echo '' 1>&2
				echo 'For RNA alignment, the --text mode is impossible.' 1>&2
				echo '' 1>&2
				exit 1
		fi
	fi

# cycle ga atode henkou sareru node koko de strategy no namae wo kimeru.
# kokokara
	if [ $pileuporshuffle = "p" ]; then
		strategy="Pileup-"
	elif [ $pileuporshuffle = "s" ]; then
		strategy="Randomchain-"
	elif [ $mccaskill -eq 1 -o $dafs -eq 1 -o $rnaalifold -eq 1 -o $contrafold -eq 1 ]; then
		if [ $distance = "scarna" -o $distance = "dafs" -o $distance = "lara" -o $distance = "slara" -o $distance = "foldalignlocal" -o $distance = "foldalignglobal" ]; then
			strategy="X-"
		elif [ $distance = "global" -o $distance = "local" -o $distance = "localgenaf" -o "globalgenaf" ]; then
			strategy="Q-"
		fi
	elif [ $distance = "fasta" -a $sw -eq 0 ]; then
		strategy="F-"
	elif [ $distance = "fasta" -a $sw -eq 1 ]; then
		strategy="H-"
	elif [ $distance = "blast" ]; then
		strategy="B-"
	elif [ $distance = "global" -o $distance = "distonly" ]; then
		strategy="G-"
	elif [ $distance = "local" ]; then
		strategy="L-"
	elif [ $distance = "last" ]; then
		strategy="Last-"
	elif [ $distance = "hybrid" ]; then
		strategy="Hybrid-"
	elif [ $distance = "multi" ]; then
		strategy="Multi-"
	elif [ $distance = "lastmulti" ]; then
		strategy="LastMulti-"
	elif [ $distance = "localgenaf" ]; then
		strategy="E-"
	elif [ $distance = "globalgenaf" ]; then
		strategy="K-"
	elif [ $fft -eq 1 ]; then
		strategy="FFT-"
	else
		strategy="NW-"
	fi
#	if [ `echo "$weighti>0.0" | bc` -gt 0 ]; then
	if [ `awk "BEGIN {print(0.0+\"$weighti\">0.0)}"` -gt 0 ]; then
		strategy=$strategy"I"
	fi
	strategy=$strategy"NS-"
	if [ $iterate -gt 0 ]; then
		strategy=$strategy"i"
	elif [ $distance = "parttree" ]; then
		if [ $partdist = "fasta" ]; then
			strategy=$strategy"FastaPartTree-"$cycle
		elif [ $partdist = "localalign" ]; then
			strategy=$strategy"DPPartTree-"$cycle
		else
			strategy=$strategy"PartTree-"$cycle
		fi
	elif [ $fragment -eq 1 ]; then
		strategy=$strategy"fragment"
	elif [ $fragment -eq -1 ]; then
		strategy=$strategy"full"
	elif [ $fragment -eq -2 ]; then
		strategy=$strategy"long"
	else
		strategy=$strategy$cycle
	fi

	explanation='?'
	performance='Not tested.'
	if [ $strategy = "F-INS-i" ]; then
		explanation='Iterative refinement method (<'$iterate') with LOCAL pairwise alignment information'
		performance='Most accurate, but very slow'
	elif [ $strategy = "L-INS-i" ]; then
		explanation='Iterative refinement method (<'$iterate') with LOCAL pairwise alignment information'
		performance='Probably most accurate, very slow'
	elif [ $strategy = "E-INS-i" ]; then
		explanation='Iterative refinement method (<'$iterate') with LOCAL pairwise alignment with generalized affine gap costs (Altschul 1998)'
		performance='Suitable for sequences with long unalignable regions, very slow'
	elif [ $strategy = "G-INS-i" ]; then
		explanation='Iterative refinement method (<'$iterate') with GLOBAL pairwise alignment information'
		performance='Suitable for sequences of similar lengths, very slow'
	elif [ $strategy = "X-INS-i" ]; then
		explanation='RNA secondary structure information is taken into account.'
		performance='For short RNA sequences only, extremely slow'
	elif [ $strategy = "F-INS-1" ]; then
		explanation='Progressive method incorporating LOCAL pairwise alignment information'
	elif [ $strategy = "L-INS-1" ]; then
		explanation='Progressive method incorporating LOCAL pairwise alignment information'
	elif [ $strategy = "G-INS-1" ]; then
		explanation='Progressive method incorporating GLOBAL pairwise alignment information'
	elif [ $strategy = "FFT-NS-i" -o $strategy = "NW-NS-i" ]; then
		explanation='Iterative refinement method (max. '$iterate' iterations)'
		if [ $iterate -gt 2 ]; then
			performance='Accurate but slow'
		else
			performance='Standard'
		fi
	elif [ $strategy = "FFT-NS-2" -o $strategy = "NW-NS-2" ]; then
		explanation='Progressive method (guide trees were built '$cycle' times.)'
		performance='Fast but rough'
	elif [ $strategy = "FFT-NS-1" -o $strategy = "NW-NS-1" ]; then
		explanation='Progressive method (rough guide tree was used.)'
		performance='Very fast but very rough'
	fi

	if [ $outputformat = "clustal" -a $outorder = "aligned" ]; then
		outputopt=" -c $strategy -r $TMPFILE/order $f2clext "
	elif [ $outputformat = "clustal" -a $outorder = "input" ]; then
		outputopt=" -c $strategy  $f2clext "
	elif [ $outputformat = "phylip" -a $outorder = "aligned" ]; then
		outputopt=" -y -r $TMPFILE/order "
	elif [ $outputformat = "phylip" -a $outorder = "input" ]; then
		outputopt=" -y "
	elif [ $outputformat = "pir" -a $outorder = "aligned" ]; then
		outputopt=" -f -r $TMPFILE/order "
	else
		outputopt="null"
	fi
# kokomade

	
	
	(
		cd $TMPFILE;

		cat /dev/null > pre

#		echo "nseq = " $nseq              1>>"$progressfile"
#		echo "distance = " $distance      1>>"$progressfile"
#		echo "iterate = " $iterate        1>>"$progressfile"
#		echo "cycle = " $cycle            1>>"$progressfile"

		if [ $anysymbol -eq 1 ]; then
			mv infile orig
			"$prefix/replaceu" $seqtype -i orig > infile 2>>"$progressfile" || exit 1
		fi

		if [ $mergetable != "/dev/null" ]; then
			if [ $nadd -gt "0" ]; then
				echo "Impossible" 1>&2
				exit 1
			fi
#			if [ $seed != "x" -o $seedtable != "x" ]; then
#				echo "This version does not support the combination of merge and seed." 1>&2
#				exit 1
#			fi
#			iterate=0 # 2013/04/16
			mergearg="-H $seedoffset"
		fi

		if [ $adjustdirection -gt 0 ]; then
			if [ $fragment -ne 0 ]; then
				fragarg="-F" #
			else
				fragarg="-F" # 2014/02/06, do not consider other additional sequences, even in the case of --add
			fi
			if [ $adjustdirection -eq 1 ]; then
				"$prefix/makedirectionlist" $fragarg -C $numthreads -m -I $nadd -i infile -t 0.00 -r 5000 -o a > _direction
			elif [ $adjustdirection -eq 2 ]; then
				"$prefix/makedirectionlist" $fragarg -C $numthreads -m -I $nadd -i infile -t 0.00 -r 100 -o a -d > _direction
			fi
			"$prefix/setdirection" $mergearg -d _direction -i infile > infiled || exit
			mv infiled infile
			if [ $anysymbol -eq 1 ]; then
				"$prefix/setdirection" $mergearg -d _direction -i orig -r  > origd || exit
				mv origd orig
			fi
		fi

		if [ $seed != "x" -o $seedtable != "x" ]; then
			if [ $pdblist != "/dev/null" -o $ownlist != "/dev/null" ]; then
				echo "The combination of --seed and (--pdbidlist or --pdbfilelist) is impossible."  1>>"$progressfile"
				exit 1
			fi
			if [ $enrich -eq 1 ]; then
				echo "The combination of --seed and (--enrich, --enrichseq or --enrichstr) is impossible at present."  1>>"$progressfile"
				exit 1
			fi
		fi

		if [ $enrich -eq 1 ]; then
			if [ $ownlist != "/dev/null" ]; then
				echo "Warning: Sequence homologs of the structures given with the --pdbfilelist option cannot be collected.\n" 1>>"$progressfile"
			fi
			echo "SEEKQUENCER (http://sysimm.ifrec.osaka-u.ac.jp/seekquencer/) is" 1>>"$progressfile"
			if [ $pdblist != "/dev/null" ]; then
				echo "collecting homoplogs of the input sequences and the structures given with the --pdbidlist option." 1>>"$progressfile"
				perl "$prefix/seekquencer_premafft.pl" $seektarget -run thread -trd 2 -seqd uniref90 -blim 1000 -noin -seqf infile -idf pdblist -out seekout -mod mafftash-split 2>>"seekerr"
				seekres="$?"
			else
				echo "collecting homologs of the input sequences." 1>>"$progressfile"
				perl "$prefix/seekquencer_premafft.pl" $seektarget -run thread -trd 2 -seqd uniref90 -blim 1000 -noin -seqf infile -out seekout -mod mafftash-split 2>>"seekerr"
				seekres="$?"
			fi
			cat seekerr  1>>"$progressfile"

			if [ $seekres -ne "0" ]; then
				echo "Error in SEEKQUENCER" 1>>"$progressfile"
				exit 1;
			fi
			echo "Done." 1>>"$progressfile"

			if [ $enrichseq -eq 1 ]; then
#				cat seekout.seq >> infile
				if [ $anysymbol -eq 1 ]; then
					"$prefix/replaceu" $seqtype -i seekout.seq -o $nseq >> infile
					cat seekout.seq >> orig
				else
					"$prefix/replaceu" $seqtype -i seekout.seq | sed 's/_os_[0-9]*_oe_//' >> infile
				fi

			fi
			if [ $enrichstr -eq 1 ]; then
				nseekstr=`wc -l < seekout.str`
				if [ $nseekstr -gt 1 ]; then
					cat seekout.str >> pdblist
					pdblist="tsukaimasu"
				fi
			fi
		fi

		if [ $seed != "x" ]; then
			mv infile infile2
			if [ $anysymbol -eq 1 ]; then
				mv orig orig2
				cat /dev/null > orig
			fi
			cat /dev/null > infile
			cat /dev/null > hat3.seed
			seedoffset=0
#			echo "seednseq="$seednseq
#			echo "seedoffset="$seedoffset
			set $seednseq >> "$progressfile"
#			echo $#
			while [ $# -gt 1 ]
			do
				shift
#				echo "num="$#

				if [ $anysymbol -eq 1 ]; then
					cat seed$# >> orig
					"$prefix/replaceu" $seqtype -i seed$# -o $seedoffset > clean 2>>"$progressfile" || exit 1
					mv clean seed$#
				fi
				"$prefix/multi2hat3s" -t $nseq -o $seedoffset -i seed$# >> infile 2>>"$progressfile" || exit 1
				cat hat3 >> hat3.seed
#				echo "$1"
				seedoffset=`expr $seedoffset + $1`
#				echo "$1"
#				echo "seedoffset="$seedoffset
			done;
#			echo "seedoffset="$seedoffset
			if [ $anysymbol -eq 1 ]; then
				"$prefix/replaceu" $seqtype -i orig2 -o $seedoffset >> infile 2>>"$progressfile" || exit 1  # yarinaoshi
				cat orig2 >> orig
			else
				cat infile2 >> infile
			fi
		elif [ $seedtable != "x" ]; then
			cat _seedtablefile > hat3.seed
		elif [ $pdblist != "/dev/null" -o $ownlist != "/dev/null" ]; then
			mv infile infile2
			if [ $anysymbol -eq 1 ]; then
				mv orig orig2
				cat /dev/null > orig
			fi
			cat /dev/null > infile

			echo "strdir = " 1>>"$progressfile"
			echo $strdir 1>>"$progressfile"

			echo "Calling DASH (http://sysimm.ifrec.osaka-u.ac.jp/dash/)" 1>>"$progressfile"
			perl "$prefix/mafftash_premafft.pl" -p pdblist -o ownlist -d "$strdir" 2>>"dasherr"
			dashres="$?"
			cat dasherr  1>>"$progressfile"

			if [ $dashres -ne "0" ]; then
				echo "Error in DASH" 1>>"$progressfile"
				exit 1;
			fi
			echo "Done." 1>>"$progressfile"

			seedoffset=`grep -c '^[>|=]' instr | head -1 ` 

			echo "# of structures = " 1>>"$progressfile"
			echo $seedoffset 1>>"$progressfile"
			mv hat3 hat3.seed

			if [ $anysymbol -eq 1 ]; then
				cat instr >> orig
				"$prefix/replaceu" $seqtype -i instr -o 0 > clean 2>>"$progressfile" || exit 1
				mv clean infile

				"$prefix/replaceu" $seqtype -i orig2 -o $seedoffset >> infile 2>>"$progressfile" || exit 1  # yarinaoshi
				cat orig2 >> orig
			else
				cat instr > infile
				cat infile2 >> infile
			fi
		else
			cat /dev/null > hat3.seed
		fi
#		cat hat3.seed




		if [ $mccaskill -eq 1 ]; then
			"$prefix/mccaskillwrap" -s -C $numthreads -d "$prefix" -i infile > hat4 2>>"$progressfile" || exit 1
		elif [ $dafs -eq 1 ]; then
			"$prefix/mccaskillwrap" -G -C $numthreads -d "$prefix" -i infile > hat4 2>>"$progressfile" || exit 1
		elif [ $contrafold -eq 1 ]; then
			"$prefix/contrafoldwrap" -d "$prefix" -i infile > hat4 2>>"$progressfile" || exit 1
		fi
		if [ $distance = "fasta" ]; then
			"$prefix/dndfast7" $swopt < infile > /dev/null  2>>"$progressfile"    || exit 1
			cat hat3.seed hat3 > hatx
			mv hatx hat3
			"$prefix/tbfast" -W $minimumweight -V "-"$gopdist -s $unalignlevel $legacygapopt $mergearg $outnum $addarg $add2ndhalfarg -C $numthreadstb $rnaopt $weightopt $treeinopt $treeoutopt $distoutopt $seqtype $model -f "-"$gop -Q $spfactor -h $aof  $param_fft $localparam   $algopt $treealg $scoreoutarg < infile   > /dev/null 2>>"$progressfile" || exit 1
		elif [ $distance = "blast" ]; then
			"$prefix/dndblast" < infile > /dev/null  2>>"$progressfile"      || exit 1
			cat hat3.seed hat3 > hatx
			mv hatx hat3
			"$prefix/tbfast" -W $minimumweight -V "-"$gopdist -s $unalignlevel $legacygapopt $mergearg $outnum $addarg $add2ndhalfarg -C $numthreadstb $rnaopt $weightopt $treeinopt $treeoutopt $distoutopt $seqtype $model -f "-"$gop -Q $spfactor -h $aof  $param_fft $localparam   $algopt $treealg $scoreoutarg < infile   > /dev/null 2>>"$progressfile" || exit 1
		elif [ $distance = "foldalignlocal" ]; then
			"$prefix/pairlocalalign" -C $numthreads $seqtype $foldalignopt $model -g $lexp -f $lgop -Q $spfactor -h $laof -H -d "$prefix" < infile > /dev/null  2>>"$progressfile"      || exit 1
			cat hat3.seed hat3 > hatx
			mv hatx hat3
			"$prefix/tbfast" -W $minimumweight -V "-"$gopdist -s $unalignlevel $legacygapopt $mergearg $outnum $addarg $add2ndhalfarg -C $numthreadstb $rnaopt $weightopt $treeinopt $treeoutopt $distoutopt $seqtype $model -f "-"$gop -Q $spfactor -h $aof  $param_fft $localparam   $algopt $treealg $scoreoutarg < infile   > /dev/null 2>>"$progressfile" || exit 1
		elif [ $distance = "foldalignglobal" ]; then
			"$prefix/pairlocalalign" -C $numthreads $seqtype $foldalignopt $model -g $pgexp -f $pggop -Q $spfactor -h $pgaof -H -o -global -d "$prefix" < infile > /dev/null  2>>"$progressfile"      || exit 1
			cat hat3.seed hat3 > hatx
			mv hatx hat3
			"$prefix/tbfast" -W $minimumweight -V "-"$gopdist -s $unalignlevel $legacygapopt $mergearg $outnum $addarg $add2ndhalfarg -C $numthreadstb $rnaopt $weightopt $treeinopt $treeoutopt $distoutopt $seqtype $model -f "-"$gop -Q $spfactor -h $aof  $param_fft $localparam   $algopt $treealg $scoreoutarg < infile   > /dev/null 2>>"$progressfile" || exit 1
		elif [ $distance = "slara" ]; then
			"$prefix/pairlocalalign" -C $numthreads -p $laraparams  $seqtype $model  -f $lgop -Q $spfactor -T -d "$prefix" < infile > /dev/null  2>>"$progressfile"      || exit 1
			cat hat3.seed hat3 > hatx
			mv hatx hat3
			"$prefix/tbfast" -W $minimumweight -V "-"$gopdist -s $unalignlevel $legacygapopt $mergearg $outnum $addarg $add2ndhalfarg -C $numthreadstb $rnaopt $weightopt $treeinopt $treeoutopt $distoutopt $seqtype $model -f "-"$gop -Q $spfactor -h $aof  $param_fft $localparam   $algopt $treealg $scoreoutarg < infile   > /dev/null 2>>"$progressfile" || exit 1
		elif [ $distance = "lara" ]; then
			"$prefix/pairlocalalign" -C $numthreads -p $laraparams  $seqtype $model  -f $lgop -Q $spfactor -B -d "$prefix" < infile > /dev/null  2>>"$progressfile"      || exit 1
			cat hat3.seed hat3 > hatx
			mv hatx hat3
			"$prefix/tbfast" -W $minimumweight -V "-"$gopdist -s $unalignlevel $legacygapopt $mergearg $outnum $addarg $add2ndhalfarg -C $numthreadstb $rnaopt $weightopt $treeinopt $treeoutopt $distoutopt $seqtype $model -f "-"$gop -Q $spfactor -h $aof  $param_fft $localparam   $algopt $treealg $scoreoutarg < infile   > /dev/null 2>>"$progressfile" || exit 1
		elif [ $distance = "scarna" ]; then
#			"$prefix/pairlocalalign"   -C $numthreads $seqtype $model  -f $pggop -Q $spfactor -s -d "$prefix" < infile > /dev/null  2>>"$progressfile"      || exit 1
#			cat hat3.seed hat3 > hatx
#			mv hatx hat3
#			"$prefix/tbfast" -W $minimumweight -V "-"$gopdist -s $unalignlevel $legacygapopt $mergearg $outnum $addarg $add2ndhalfarg -C $numthreadstb $rnaopt $weightopt $treeinopt $treeoutopt $distoutopt $seqtype $model -f "-"$gop -Q $spfactor -h $aof  $param_fft $localparam   $algopt $treealg $scoreoutarg < infile   > /dev/null 2>>"$progressfile" || exit 1
			"$prefix/tbfast" _  -C $numthreads $seqtype $model  -f $pggop -Q $spfactor -s -d "$prefix" _ -+ $iterate -W $minimumweight -V "-"$gopdist -s $unalignlevel $legacygapopt $mergearg $outnum $addarg $add2ndhalfarg -C $numthreadstb $rnaopt $weightopt $treeinopt $treeoutopt $distoutopt $seqtype $model -f "-"$gop -Q $spfactor -h $aof  $param_fft $localparam   $algopt $treealg $scoreoutarg < infile   > /dev/null 2>>"$progressfile" || exit 1
		elif [ $distance = "dafs" ]; then
			"$prefix/pairlocalalign"  -C $numthreads $seqtype $model  -f $pggop -Q $spfactor -G -d "$prefix" < infile > /dev/null  2>>"$progressfile"      || exit 1
			cat hat3.seed hat3 > hatx
			mv hatx hat3
			"$prefix/tbfast" -W $minimumweight -V "-"$gopdist -s $unalignlevel $legacygapopt $mergearg $outnum $addarg $add2ndhalfarg -C $numthreadstb $rnaopt $weightopt $treeinopt $treeoutopt $distoutopt $seqtype $model -f "-"$gop -Q $spfactor -h $aof  $param_fft $localparam   $algopt $treealg $scoreoutarg < infile   > /dev/null 2>>"$progressfile" || exit 1
		elif [ $distance = "global" ]; then
#			"$prefix/pairlocalalign" -u $unalignlevel $localparam  -C $numthreads $seqtype $model -g $pgexp -f $pggop -Q $spfactor -h $pgaof  -A  $usenaivepairscore $focusarg < infile > /dev/null  2>>"$progressfile"      || exit 1
#			cat hat3.seed hat3 > hatx
#			mv hatx hat3
#			"$prefix/tbfast" -W $minimumweight -V "-"$gopdist -s $unalignlevel $legacygapopt $mergearg $outnum $addarg $add2ndhalfarg -C $numthreadstb $rnaopt $weightopt $treeinopt $treeoutopt $distoutopt $seqtype $model -f "-"$gop -Q $spfactor -h $aof  $param_fft $localparam   $algopt $treealg $scoreoutarg $focusarg < infile   > /dev/null 2>>"$progressfile" || exit 1
			"$prefix/tbfast" _  -u $unalignlevel $localparam  -C $numthreads $seqtype $model -g $pgexp -f $pggop -Q $spfactor -h $pgaof  -A  $usenaivepairscore $focusarg  _ -+ $iterate -W $minimumweight -V "-"$gopdist -s $unalignlevel $legacygapopt $mergearg $outnum $addarg $add2ndhalfarg -C $numthreadstb $rnaopt $weightopt $treeinopt $treeoutopt $distoutopt $seqtype $model -f "-"$gop -Q $spfactor -h $aof  $param_fft $localparam   $algopt $treealg $scoreoutarg $focusarg < infile   > /dev/null 2>>"$progressfile" || exit 1
			
		elif [ $distance = "local" ]; then
			if [ $fragment -ne 0 ]; then 
				"$prefix/pairlocalalign" $localparam $addarg   -C $numthreads $seqtype $model  -g $lexp -f $lgop -Q $spfactor -h $laof -L $usenaivepairscore < infile > /dev/null  2>>"$progressfile"      || exit 1
				cat hat3.seed hat3 > hatx
				mv hatx hat3
				"$prefix/addsingle" -Q 100 $legacygapopt -O $outnum $addsinglearg $addarg $add2ndhalfarg -C $numthreads $memopt $weightopt $treeinopt $treeoutopt $distoutopt $seqtype $model -f "-"$gop  -h $aof  $param_fft $localparam   $algopt $treealg $scoreoutarg < infile   > /dev/null 2>>"$progressfile" || exit 1
			else
#				"$prefix/pairlocalalign" -u $unalignlevel $localparam -C $numthreads $seqtype $model  -g $lexp -f $lgop -Q $spfactor -h $laof -L $usenaivepairscore $focusarg < infile > /dev/null  2>>"$progressfile"      || exit 1
#				addarg wo watasanai
#				cat hat3.seed hat3 > hatx
#				mv hatx hat3
#				"$prefix/tbfast" -W $minimumweight -V "-"$gopdist -s $unalignlevel $legacygapopt $mergearg $termgapopt $outnum $addarg $add2ndhalfarg -C $numthreadstb $rnaopt $weightopt $treeinopt $treeoutopt $distoutopt $seqtype $model -f "-"$gop -Q $spfactor -h $aof  $param_fft $localparam   $algopt $treealg $scoreoutarg $focusarg < infile   > /dev/null 2>>"$progressfile" || exit 1
				"$prefix/tbfast" _  -u $unalignlevel $localparam -C $numthreads $seqtype $model  -g $lexp -f $lgop -Q $spfactor -h $laof -L $usenaivepairscore $focusarg _ -+ $iterate -W $minimumweight -V "-"$gopdist -s $unalignlevel $legacygapopt $mergearg $termgapopt $outnum $addarg $add2ndhalfarg -C $numthreadstb $rnaopt $weightopt $treeinopt $treeoutopt $distoutopt $seqtype $model -f "-"$gop -Q $spfactor -h $aof  $param_fft $localparam   $algopt $treealg $scoreoutarg $focusarg < infile   > /dev/null 2>>"$progressfile" || exit 1
			fi
		elif [ $distance = "globalgenaf" ]; then
			"$prefix/pairlocalalign"  -u $unalignlevel $localparam -C $numthreads $seqtype $model  -g $pgexp -f $pggop -Q $spfactor -h $pgaof -O $GGOP -E $GEXP -K $usenaivepairscore < infile > /dev/null 2>>"$progressfile"    || exit 1
			cat hat3.seed hat3 > hatx
			mv hatx hat3
			"$prefix/tbfast" -W $minimumweight -V "-"$gopdist -s $unalignlevel $legacygapopt $mergearg $outnum $addarg $add2ndhalfarg -C $numthreadstb $rnaopt $weightopt $treeinopt $treeoutopt $distoutopt $seqtype $model -f "-"$gop -Q $spfactor -h $aof  $param_fft $localparam   $algopt $treealg $scoreoutarg < infile   > /dev/null 2>>"$progressfile" || exit 1
		elif [ $distance = "localgenaf" ]; then
#			"$prefix/pairlocalalign"  -u $unalignlevel $localparam -C $numthreads $seqtype $model -g $lexp -f $lgop -Q $spfactor -h $laof -O $LGOP -E $LEXP -N $usenaivepairscore $focusarg < infile > /dev/null  2>>"$progressfile"      || exit 1
#			cat hat3.seed hat3 > hatx
#			mv hatx hat3
#			"$prefix/tbfast" -W $minimumweight -V "-"$gopdist -s $unalignlevel $legacygapopt $mergearg $termgapopt $outnum $addarg $add2ndhalfarg -C $numthreadstb $rnaopt $weightopt $treeinopt $treeoutopt $distoutopt $seqtype $model -f "-"$gop -Q $spfactor -h $aof  $param_fft $localparam   $algopt $treealg $scoreoutarg $focusarg < infile   > /dev/null 2>>"$progressfile" || exit 1
			"$prefix/tbfast" _ -u $unalignlevel $localparam -C $numthreads $seqtype $model -g $lexp -f $lgop -Q $spfactor -h $laof -O $LGOP -E $LEXP -N $usenaivepairscore $focusarg _ -+ $iterate -W $minimumweight -V "-"$gopdist -s $unalignlevel $legacygapopt $mergearg $termgapopt $outnum $addarg $add2ndhalfarg -C $numthreadstb $rnaopt $weightopt $treeinopt $treeoutopt $distoutopt $seqtype $model -f "-"$gop -Q $spfactor -h $aof  $param_fft $localparam   $algopt $treealg $scoreoutarg $focusarg < infile > /dev/null  2>>"$progressfile"      || exit 1
		elif [ $distance = "last" ]; then
			if [ $fragment -ne 0 ]; then 
				"$prefix/pairlocalalign" $addarg   -C $numthreads $seqtype $model -e $last_e -w $last_m -g $lexp -f $lgop -Q $spfactor -h $laof -R $last_subopt $last_once -d "$prefix" < infile > /dev/null  2>>"$progressfile"      || exit 1
				cat hat3.seed hat3 > hatx
				mv hatx hat3
				"$prefix/addsingle" -Q 100 $legacygapopt -O $outnum $addsinglearg $addarg $add2ndhalfarg -C $numthreads $memopt $weightopt $treeinopt $treeoutopt $distoutopt $seqtype $model -f "-"$gop  -h $aof  $param_fft $localparam   $algopt $treealg $scoreoutarg < infile   > /dev/null 2>>"$progressfile" || exit 1
			else
				"$prefix/pairlocalalign" -C $numthreads $seqtype $model -e $last_e -w $last_m -g $lexp -f $lgop -Q $spfactor -h $laof -R $last_subopt $last_once -d "$prefix" < infile > /dev/null  2>>"$progressfile"      || exit 1
#				addarg wo watasanai
				cat hat3.seed hat3 > hatx
				mv hatx hat3
				"$prefix/tbfast" -W $minimumweight -V "-"$gopdist -s $unalignlevel $legacygapopt $mergearg $termgapopt $outnum $addarg $add2ndhalfarg -C $numthreadstb $rnaopt $weightopt $treeinopt $treeoutopt $distoutopt $seqtype $model -f "-"$gop -Q $spfactor -h $aof  $param_fft $localparam   $algopt $treealg $scoreoutarg < infile   > /dev/null 2>>"$progressfile" || exit 1
			fi
		elif [ $distance = "lastmulti" ]; then
			"$prefix/dndpre" $model -M 2 $addarg -C $numthreads $seqtype $model -g $lexp -f $lgop -Q $spfactor -h $laof < infile > /dev/null 2>>"$progressfile"      || exit 1
			mv hat2 hat2i
			"$prefix/pairlocalalign" $addarg   -C $numthreads $seqtype $model -e $last_e -w $last_m -g $lexp -f $lgop -Q $spfactor -h $laof -r $last_subopt $last_once -d "$prefix" < infile > /dev/null  2>>"$progressfile"      || exit 1
			cat hat3.seed hat3 > hatx
			mv hat2 hat2n
			mv hatx hat3
			if [ $fragment -ne 0 ]; then 
				"$prefix/addsingle" -Q 100 $legacygapopt -d -O $outnum $addsinglearg $addarg $add2ndhalfarg -C $numthreads $memopt $weightopt $treeinopt $treeoutopt $distoutopt $seqtype $model -f "-"$gop  -h $aof  $param_fft $localparam   $algopt $treealg $scoreoutarg < infile   > /dev/null 2>>"$progressfile" || exit 1
			else
				echo "Impossible" 1>&2
				exit 1
			fi
		elif [ $distance = "multi" ]; then
			"$prefix/dndpre" $model -M 2 $addarg -C $numthreads $seqtype $model -g $lexp -f $lgop -h $laof  $usenaivepairscore < infile > /dev/null 2>>"$progressfile"      || exit 1
			mv hat2 hat2i
			"$prefix/pairlocalalign" $localparam $addarg   -C $numthreads $seqtype $model  -g $lexp -f $lgop -Q $spfactor -h $laof -Y $usenaivepairscore < infile > /dev/null  2>>"$progressfile"      || exit 1
			cat hat3.seed hat3 > hatx
			mv hat2 hat2n
			mv hatx hat3
			if [ $fragment -ne 0 ]; then 
				"$prefix/addsingle" -Q 100 $legacygapopt -d -O $outnum $addsinglearg $addarg $add2ndhalfarg -C $numthreads $memopt $weightopt $treeinopt $treeoutopt $distoutopt $seqtype $model -f "-"$gop  -h $aof  $param_fft $localparam   $algopt $treealg $scoreoutarg < infile   > /dev/null 2>>"$progressfile" || exit 1
			else
				echo "Impossible" 1>&2
				exit 1
			fi
		elif [ $distance = "hybrid" ]; then
			"$prefix/pairlocalalign" $addarg   -C $numthreads $seqtype $model  -g $lexp -f $lgop -Q $spfactor -h $laof -Y < infile > /dev/null  2>>"$progressfile"      || exit 1
			cat hat3.seed hat3 > hatx
			mv hatx hat3
			"$prefix/disttbfast" -E 1 -s $unalignlevel $legacygapopt -W $tuplesize $termgapopt $outnum $addarg $add2ndhalfarg -C $numthreadstb $memopt $weightopt $treeinopt $treeoutopt -T -y $seqtype $model -f "-"$gop -Q $spfactor -h $aof  $param_fft $algopt $treealg $scoreoutarg < infile   > /dev/null 2>>"$progressfile" || exit 1
			if [ $fragment -ne 0 ]; then 
				"$prefix/addsingle" -Q 100 $legacygapopt -O $outnum $addsinglearg $addarg $add2ndhalfarg -C $numthreads $memopt $weightopt $treeinopt $treeoutopt $distoutopt $seqtype $model -f "-"$gop  -h $aof  $param_fft $localparam   $algopt $treealg $scoreoutarg < infile   > /dev/null 2>>"$progressfile" || exit 1
			else
				"$prefix/tbfast" -W $minimumweight -V "-"$gopdist -s $unalignlevel $legacygapopt $mergearg $termgapopt $outnum $addarg $add2ndhalfarg -C $numthreadstb $rnaopt $weightopt $treeinopt $treeoutopt $distoutopt $seqtype $model -f "-"$gop -Q $spfactor -h $aof  $param_fft $localparam   $algopt $treealg $scoreoutarg < infile   > /dev/null 2>>"$progressfile" || exit 1
			fi
#		elif [ $distance = "distonly" ]; then
#			"$prefix/pairlocalalign"   -C $numthreads $seqtype $model -g $pgexp -f $pggop -Q $spfactor -h $pgaof  -t < infile > /dev/null  2>>"$progressfile"      || exit 1
#			"$prefix/tbfast" -W $minimumweight -V "-"$gopdist -s $unalignlevel $legacygapopt $outnum $addarg $add2ndhalfarg -C $numthreadstb $rnaopt $weightopt $treeinopt $treeoutopt $distoutopt $seqtype $model -f "-"$gop -Q $spfactor -h $aof  $param_fft $localparam   $algopt $treealg $scoreoutarg < infile   > /dev/null 2>>"$progressfile" || exit 1
		elif [ $distance = "parttree" ]; then
			"$prefix/splittbfast" $legacygapopt $algopt $splitopt $partorderopt $parttreeoutopt $memopt $seqtype $model -f "-"$gop -Q $spfactor -h $aof -p $partsize -s $groupsize $treealg -i infile   > pre 2>>"$progressfile" || exit 1
			mv hat3.seed hat3
		elif [ $distance = "ktuplesmulti" ]; then
#			"$prefix/dndpre" $model -M 1 $addarg -C $numthreads $seqtype $model -g $lexp -f $lgop -h $laof < infile > /dev/null 2>>"$progressfile"      || exit 1
#			mv hat2 hat2i
#			"$prefix/disttbfast" -E 1 -s $unalignlevel $legacygapopt -W $tuplesize $termgapopt $outnum $addarg $add2ndhalfarg -C $numthreadstb $memopt $weightopt $treeinopt $treeoutopt -T -y $seqtype $model -f "-"$gop -Q $spfactor -h $aof  $param_fft $algopt $treealg $scoreoutarg < infile   > /dev/null 2>>"$progressfile" || exit 1
#			mv hat2 hat2n
			if [ $fragment -ne 0 ]; then 
				"$prefix/addsingle" -Q 100 $legacygapopt -d -W $tuplesize -O $outnum $addsinglearg $addarg $add2ndhalfarg -C $numthreads $memopt $weightopt $treeinopt $treeoutopt $distoutopt $seqtype $model -f "-"$gop  -h $aof  $param_fft $localparam   $algopt $treealg $scoreoutarg < infile   > /dev/null 2>>"$progressfile" || exit 1
#				"$prefix/addsingle" -Q 100 $legacygapopt -d -O $outnum $addsinglearg $addarg $add2ndhalfarg -C $numthreads $memopt $weightopt $treeinopt $treeoutopt $distoutopt $seqtype $model -f "-"$gop  -h $aof  $param_fft $localparam   $algopt $treealg $scoreoutarg < infile   > /dev/null 2>>"$progressfile" || exit 1
			else
				echo "Impossible" 1>&2
				exit 1
			fi
		else
			if [ $fragment -ne 0 ]; then 
				"$prefix/addsingle" -Q 100 $legacygapopt -W $tuplesize -O $outnum $addsinglearg $addarg $add2ndhalfarg -C $numthreads $memopt $weightopt $treeinopt $treeoutopt $distoutopt $seqtype $model -f "-"$gop  -h $aof  $param_fft $localparam   $algopt $treealg $scoreoutarg < infile   > /dev/null 2>>"$progressfile" || exit 1
			else
				"$prefix/disttbfast" -q $npickup -E $cycledisttbfast -V "-"$gopdist  -s $unalignlevel $legacygapopt $mergearg -W $tuplesize $termgapopt $outnum $addarg $add2ndhalfarg -C $numthreadstb $memopt $weightopt $treeinopt $treeoutopt $distoutopt $seqtype $model -f "-"$gop -Q $spfactor -h $aof  $param_fft $algopt $treealg $scoreoutarg < infile   > pre 2>>"$progressfile" || exit 1
				mv hat3.seed hat3
			fi
		fi
		while [ $cycletbfast -gt 1 ]
		do
			if [ $distance = "parttree" ]; then
				mv pre infile
				"$prefix/splittbfast" $legacygapopt -Z $algopt $splitopt $partorderopt $parttreeoutopt $memopt $seqtype $model -f "-"$gop -Q $spfactor -h $aof  -p $partsize -s $groupsize $treealg -i infile   > pre 2>>"$progressfile" || exit 1
			else
				"$prefix/tbfast" -W $minimumweight -V "-"$gopdist -s $unalignlevel $legacygapopt $mergearg $termgapopt $outnum -C $numthreadstb $rnaopt $weightopt $treeoutopt $distoutopt $memopt $seqtype $model  -f "-"$gop -Q $spfactor -h $aof $param_fft  $localparam $algopt -J $treealg $scoreoutarg < pre > /dev/null 2>>"$progressfile" || exit 1
# fragment>0 no baai, nanimoshinai
# seed youchuui!!
			fi
			cycletbfast=`expr $cycletbfast - 1`
		done
		if [ $iterate -gt 0 ]; then
			if [ $distance = "ktuples" ]; then
			    "$prefix/dndpre" $model -M 2 -C $numthreads < pre     > /dev/null 2>>"$progressfile" || exit 1
			fi
			"$prefix/dvtditr" -W $minimumweight $bunkatsuopt -E $fixthreshold -s $unalignlevel  $legacygapopt $mergearg -C $numthreadsit -t $randomseed $rnaoptit $memopt $scorecalcopt $localparam -z 50 $seqtype $model -f "-"$gop -Q $spfactor -h $aof  -I $iterate $weightopt $treeinopt $algoptit $treealg -p $parallelizationstrategy  $scoreoutarg < pre     > /dev/null 2>>"$progressfile" || exit 1
		fi
		if [ $coreout -eq 1 ]; then
			"$prefix/setcore" -w $corewin -i $corethr $coreext < pre > pre2
			mv pre2 pre
		elif [ $anysymbol -eq 1 ]; then
			"$prefix/restoreu" $add2ndhalfarg -a pre -i orig > restored || exit 1
			mv restored pre
		fi




		echo '' 1>>"$progressfile"
		if [ $mccaskill -eq 1 ]; then
			echo "RNA base pairing probaility was calculated by the McCaskill algorithm (1)" 1>>"$progressfile"
			echo "implemented in Vienna RNA package (2) and MXSCARNA (3), and then" 1>>"$progressfile"
			echo "incorporated in the iterative alignment process (4)." 1>>"$progressfile"
			echo "(1) McCaskill, 1990, Biopolymers 29:1105-1119" 1>>"$progressfile"
			echo "(2) Hofacker et al., 2002, J. Mol. Biol. 319:3724-3732" 1>>"$progressfile"
			echo "(3) Tabei et al., 2008, BMC Bioinformatics 9:33" 1>>"$progressfile"
			echo "(4) Katoh and Toh, 2008, BMC Bioinformatics 9:212" 1>>"$progressfile"
			echo "" 1>>"$progressfile"
		elif [ $contrafold -eq 1 ]; then
			echo "RNA base pairing probaility was calculated by the CONTRAfold algorithm (1)" 1>>"$progressfile"
			echo "and then incorporated in the iterative alignment process (4)." 1>>"$progressfile"
			echo "(1) Do et al., 2006, Bioinformatics 22:e90-98" 1>>"$progressfile"
			echo "(2) Katoh and Toh, 2008, BMC Bioinformatics 9:212" 1>>"$progressfile"
			echo "" 1>>"$progressfile"
		fi
		if [ $pdblist != "/dev/null" -o $ownlist != "/dev/null" ]; then
			echo "Input structures are decomposed into structural domains using" 1>>"$progressfile"
			echo "Protein Domain Parser (Alexandrov & Shindyalov 2003)."         1>>"$progressfile"
			echo "Domain pairs are aligned using the rash function in"           1>>"$progressfile"
			echo "the ASH structural alignment package (Standley et al. 2007)."  1>>"$progressfile"
		fi
		if [ $pdblist != "/dev/null" ]; then
			echo "Pre-computed alignments stored in "                            1>>"$progressfile"
			echo "DASH (http://sysimm.ifrec.osaka-u.ac.jp/dash/) are used. "     1>>"$progressfile"
		fi
		if [ $distance = "fasta" -o $partdist = "fasta" ]; then
			echo "Pairwise alignments were computed by FASTA" 1>>"$progressfile"
			echo "(Pearson & Lipman, 1988, PNAS 85:2444-2448)" 1>>"$progressfile"
		fi
		if [ $distance = "blast" ]; then
			echo "Pairwise alignments were computed by BLAST" 1>>"$progressfile"
			echo "(Altschul et al., 1997, NAR 25:3389-3402)" 1>>"$progressfile"
		fi
		if [ $distance = "last" -o $distance = "lastmulti" ]; then
			echo "Pairwise alignments were computed by LAST" 1>>"$progressfile"
			echo "http://last.cbrc.jp/" 1>>"$progressfile"
			echo "Kielbasa, Wan, Sato, Horton, Frith 2011 Genome Res. 21:487" 1>>"$progressfile"
		fi
		if [ $distance = "scarna" ]; then
			echo "Pairwise alignments were computed by MXSCARNA" 1>>"$progressfile"
			echo "(Tabei et al., 2008, BMC Bioinformatics 9:33)." 1>>"$progressfile"
		fi
		if [ $distance = "dafs" ]; then
			echo "Pairwise alignments were computed by DAFS" 1>>"$progressfile"
			echo "(Sato et al., 2012,,,,)." 1>>"$progressfile"
		fi
		if [ $distance = "lara" -o $distance = "slara" ]; then
			echo "Pairwise alignments were computed by LaRA" 1>>"$progressfile"
			echo "(Bauer et al., 2007, BMC Bioinformatics 8:271)." 1>>"$progressfile"
		fi
		if [ $distance = "foldalignlocal" ]; then
			echo "Pairwise alignments were computed by FOLDALIGN (local)" 1>>"$progressfile"
			echo "(Havgaard et al., 2007, PLoS Computational Biology 3:e193)." 1>>"$progressfile"
		fi
		if [ $distance = "foldalignglobal" ]; then
			echo "Pairwise alignments were computed by FOLDALIGN (global)" 1>>"$progressfile"
			echo "(Havgaard et al., 2007, PLoS Computational Biology 3:e193)." 1>>"$progressfile"
		fi
#		printf "\n" 1>>"$progressfile"
		echo 'Strategy:' 1>>"$progressfile"
		printf ' '$strategy 1>>"$progressfile"
		echo ' ('$performance')' 1>>"$progressfile"
		echo ' '$explanation 1>>"$progressfile"
		echo '' 1>>"$progressfile"
		echo "If unsure which option to use, try 'mafft --auto input > output'." 1>>"$progressfile"
		echo "For more information, see 'mafft --help', 'mafft --man' and the mafft page." 1>>"$progressfile"
		echo "" 1>>"$progressfile"
		echo "The default gap scoring scheme has been changed in version 7.110 (2013 Oct)." 1>>"$progressfile"
		echo "It tends to insert more gaps into gap-rich regions than previous versions." 1>>"$progressfile"
		echo "To disable this change, add the --leavegappyregion option." 1>>"$progressfile"
#		echo "If long gaps are expected, try 'mafft --ep 0.0 --auto input > output'." 1>>"$progressfile"
#		echo "If the possibility of long gaps can be excluded, add '--ep 0.123'." 1>>"$progressfile"
		if [ $distance = "localgenaf" -o $distance = "globalgenaf" ]; then
			echo "" 1>>"$progressfile"
			if [ $oldgenafparam -eq 1 ]; then
				echo "Obsolete parameters used for this calculation." 1>>"$progressfile"
				echo "Also try the new parameters for E-INS-i, by not specifying --oldgenafpair." 1>>"$progressfile"
			else
				echo "Parameters for the E-INS-i option have been changed in version 7.243 (2015 Jun)." 1>>"$progressfile"
				echo "To switch to the old parameters, use --oldgenafpair, instead of --genafpair." 1>>"$progressfile"
			fi
		fi
		echo '' 1>>"$progressfile"

		
		if [ $pdblist != "/dev/null" -o $ownlist != "/dev/null" ]; then
			cat dasherr >>"$progressfile"
		fi

	)


	if [ "$outputfile" = "" ]; then
		if [ "$outputopt" = "null" ]; then
			cat < $TMPFILE/pre || exit 1
		else
			"$prefix/f2cl" -n $namelength $outputopt < $TMPFILE/pre 2>>/dev/null || exit 1
		fi
	else
		if [ "$outputopt" = "null" ]; then
			cat < $TMPFILE/pre > "$outputfile" || exit 1
		else
			"$prefix/f2cl" -n $namelength $outputopt < $TMPFILE/pre > "$outputfile" 2>>/dev/null || exit 1
		fi
	fi

	if [ $treeout -eq 1 ]; then
		cp $TMPFILE/infile.tree "$infilename.tree"
	fi

	if [ -s $TMPFILE/GuideTree ]; then # --merge no toki dake
		cp $TMPFILE/GuideTree .
	fi

	if [ $distout -eq 1 ]; then
		cp $TMPFILE/hat2 "$infilename.hat2"
	fi

	if [ $npickup -ne 0 ]; then
		cp $TMPFILE/notused "$infilename.notused"
	fi

	if [ -s $TMPFILE/_deletemap ]; then
		if [ "$mapoutfile" = "/dev/null" ]; then
			cp $TMPFILE/_deletemap "$addfile.map"
		else
			cp $TMPFILE/_deletemap "$mapoutfile"
		fi
	fi

	exit 0;
fi

prog="awk"

tmpawk=`which nawk 2>/dev/null | awk '{print $1}'`
if [ -x "$tmpawk" ]; then
	prog="$tmpawk"
fi

tmpawk=`which gawk 2>/dev/null | awk '{print $1}'`
if [ -x "$tmpawk" ]; then
	prog="$tmpawk"
fi

#echo "prog="$prog 1>&2

umask 077
(
$prog '
BEGIN {
	prefix = ENVIRON["prefix"];
	version = ENVIRON["version"];
	myself = ENVIRON["myself"];
	while( 1 )
	{
		options = ""
		printf( "\n" ) > "/dev/tty";
		printf( "---------------------------------------------------------------------\n" )      > "/dev/tty";
		printf( "\n" )                                                                           > "/dev/tty";
		printf( "   MAFFT %s\n", version )                                                       > "/dev/tty";
		printf( "\n" )                                                                           > "/dev/tty";
		printf( "        Copyright (c) 2016 Kazutaka Katoh\n" )                                  > "/dev/tty";
		printf( "        MBE 30:772-780 (2013), NAR 30:3059-3066 (2002)\n" )                     > "/dev/tty";
		printf( "        http://mafft.cbrc.jp/alignment/software/\n" )                           > "/dev/tty";
		printf( "---------------------------------------------------------------------\n" )      > "/dev/tty";
		printf( "\n" ) > "/dev/tty";
	
		while( 1 )
		{
			printf( "\n" ) > "/dev/tty";
			printf( "Input file? (fasta format)\n@ " ) > "/dev/tty";
			res = getline < "/dev/tty";
			close( "/dev/tty" )
			if( res == 0 || NF == 0 )
				continue;
			infile = sprintf( "%s", $0 );
	
			res = getline < infile;
			close( infile );
			if( res == -1 )
			{
				printf( "%s: No such file.\n\n", infile ) > "/dev/tty";
				printf( "Filename extension (eg., .txt) must be typed, if any.\n\n" ) > "/dev/tty";
			}
			else if( res == 0 )
				printf( "%s: Empty.\n", infile ) > "/dev/tty";
			else
			{
				printf( "OK. infile = %s\n\n", infile ) > "/dev/tty";
				break;
			}
		}
		nseq = 0;
	
		while( 1 )
		{
			printf( "\n" ) > "/dev/tty";
			printf( "Output file?\n" ) > "/dev/tty";
			printf( "@ " ) > "/dev/tty";
			res = getline < "/dev/tty";
			close( "/dev/tty" );
			if( res == 0 || NF == 0 )
				continue;
			else
			{
				outfile = sprintf( "%s", $0 );
				printf( "OK. outfile = %s\n\n", outfile ) > "/dev/tty";
				break;
			}
		}
	
		while( 1 )
		{
			outargs = "";
			printf( "\n" ) > "/dev/tty";
			printf( "Output format?\n" ) > "/dev/tty";
			printf( "  1. Clustal format / Sorted\n" ) > "/dev/tty";
			printf( "  2. Clustal format / Input order\n" ) > "/dev/tty";
			printf( "  3. Fasta format   / Sorted\n" ) > "/dev/tty";
			printf( "  4. Fasta format   / Input order\n" ) > "/dev/tty";
			printf( "  5. Phylip format  / Sorted\n" ) > "/dev/tty";
			printf( "  6. Phylip format  / Input order\n" ) > "/dev/tty";
			printf( "@ " ) > "/dev/tty";
			res = getline < "/dev/tty";
			close( "/dev/tty" );
#			printf( "res=%d, NF=%d\n", res, NF );

			resnum = 0 + $1;
#			printf( "resnum=%d\n", resnum );

			if( resnum < 1 || 6 < resnum )
				continue;
			else
			{
				if( resnum == 1 )
					outargs = "--clustalout --reorder";
				else if( resnum == 2 )
					outargs = "--clustalout --inputorder";
				else if( resnum == 3 )
					outargs = "--reorder";
				else if( resnum == 4 )
					outargs = "--inputorder";
				else if( resnum == 5 )
					outargs = "--phylipout --reorder";
				else if( resnum == 6 )
					outargs = "--phylipout --inputorder";
				else
					continue;
				printf( "OK. arguments = %s\n\n", outargs ) > "/dev/tty";
				break;
			}
		}
	
		while( 1 )
		{
			arguments = "";
			printf( "\n" ) > "/dev/tty";
			printf( "Strategy?\n" ) > "/dev/tty";
			printf( "  1. --auto\n" ) > "/dev/tty";
			printf( "  2. FFT-NS-1 (fast)\n" ) > "/dev/tty";
			printf( "  3. FFT-NS-2 (default)\n" ) > "/dev/tty";
			printf( "  4. G-INS-i  (accurate)\n" ) > "/dev/tty";
			printf( "  5. L-INS-i  (accurate)\n" ) > "/dev/tty";
			printf( "  6. E-INS-i  (accurate)\n" ) > "/dev/tty";
			printf( "@ " ) > "/dev/tty";
			res = getline < "/dev/tty";
			close( "/dev/tty" );
#			printf( "res=%d, NF=%d\n", res, NF );

			resnum = 0 + $1;
#			printf( "resnum=%d\n", resnum );

			if( resnum < 1 || 6 < resnum )
				continue;
			else
			{
				if( resnum == 1 )
					arguments = "--auto";
				else if( resnum == 2 )
					arguments = "--retree 1";
				else if( resnum == 3 )
					arguments = "--retree 2";
				else if( resnum == 4 )
					arguments = "--globalpair --maxiterate 16";
				else if( resnum == 5 )
					arguments = "--localpair  --maxiterate 16";
				else if( resnum == 6 )
					arguments = "--genafpair  --maxiterate 16";
				else
					arguments = sprintf( "%s", $0 );
				printf( "OK. arguments = %s %s\n\n", arguments, outargs ) > "/dev/tty";
				break;
			}
		}


		while( 1 )
		{
			printf( "\n" ) > "/dev/tty";
			printf( "Additional arguments? (--ep # --op # --kappa # etc)\n" ) > "/dev/tty";
			printf( "@ " ) > "/dev/tty";
			res = getline < "/dev/tty";
			close( "/dev/tty" );
			if( res == 0 || NF == 0 )
			{
				break;
			}
			else
			{
				addargs = sprintf( "%s", $0 );
				printf( "OK. arguments = %s %s %s\n\n", addargs, arguments, outargs ) > "/dev/tty";
				break;
			}
		}

		arguments = sprintf( "%s %s %s", addargs, arguments, outargs );

		print ""
		command = sprintf( "\"%s\" %s \"%s\" > \"%s\"", myself, arguments, infile, outfile );
		gsub( /\\/, "/", command );


		printf( "command=\n%s\n", command ) > "/dev/tty";
	
	
		while( 1 )
		{
			go = 0;
			printf( "OK?\n" ) > "/dev/tty";
			printf( "@ [Y] " ) > "/dev/tty";
			res = getline < "/dev/tty";
			close( "/dev/tty" );
			if( res == 0 )
				continue;
			else if( NF == 0 || $0 ~ /^[Yy]/ )
			{
				go=1;
				break;
			}
			else
				break;
		}
		if( go ) break;
		printf( "\n" ) > "/dev/tty";
		printf( "\n" ) > "/dev/tty";
	}
	system( command );
	command = sprintf( "less \"%s\"", outfile );
	system( command );
	printf( "Press Enter to exit." ) > "/dev/tty";
	res = getline < "/dev/tty";
}
'
)
exit 0;
