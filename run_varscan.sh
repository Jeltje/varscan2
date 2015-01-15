#! /bin/bash

print_usage(){
cat <<EOF
$0 -c <control bam file> -t <tumor bam file> -i <genome index> -a <arm sizes> -s <scratch output>
	Wrapper script for Varscan2
	Runs the following steps:
	1. samtools flagstat on each bam file
	2. samtools mpileup on both bam files
	3. determine unique mapped read ratio
	4. Varscan2.3.7 copynumber
	5. Varscan2.3.7 copyCaller
	6. DNAcopy 
	7. MergeSegments.pl
	8. create bed output

OPTIONS:
   -h      Show this message
   -t      tumor bam file
   -c      control bam file
   -i      path to samtools genome index
   -a      ref_arm_sizes file for mergeSegments.pl
           Format: 
		1       0       124300000       p
		1       124300000       247249719       q
   -s      directory for temporary files
   -n      <varscan.N> instead of creating a new temporary directory, use this one
   -d      debug (print commands but do not run)

EOF
}

cBam=
tBam=
idx=
arms=
scratch=
prevDir=
debug=

while getopts "ht:c:i:s:n:a:d" OPTION
do
     case $OPTION in
         h)
             usage
             exit 1
             ;;
         t)
             tBam=$OPTARG
             ;;
         c)
             cBam=$OPTARG
             ;;
         i)
             idx=$OPTARG
             ;;
         a)
             arms=$OPTARG
             ;;
         s)
             scratch=$OPTARG
             ;;
         d)
             DEBUG=1
             ;;
         n)
             prevDir=$OPTARG
             ;;
         ?)
             usage
             exit
             ;;
     esac
done


if [[ -z $cBam ]] || [[ -z $tBam ]] || [[ -z $idx ]] || [[ -z $arms ]]
then
     print_usage
     echo "ERROR: Please give in two bam files, samtools index file, arm sizes file, and a directory to put large temporary files"
     exit 1
fi

graceful_death() {
	echo "ERROR: Cannot finish $0 because $1";
	exit;
}

if [ ! -e "$cBam" ]; then
	graceful_death "cannot find control bam file $cBam"
fi 
if [ ! -e "$tBam" ]; then
	graceful_death "cannot find tumor bam file $tBam"
fi 
if [ ! -e "$idx" ]; then
	graceful_death "cannot find samtools index $idx"
fi 
if [ ! -e "$arms" ]; then
	graceful_death "cannot find mergeSegments input file $arms"
fi 

tmpdir=
if [[ -z "$prevDir" ]] && [[ -z $scratch ]]; then
	graceful_death "Please give either the -n OR -s option"
fi

# select correct temp dir
if  [[ -z "$prevDir" ]]; then
	if [ ! -d "$scratch" ]; then
		graceful_death "cannot find scratch output dir $scratch"
	fi 
	tmpExt=$RANDOM
	tmpdir="$scratch/varscan.$tmpExt"
	mkdir $tmpdir
else
	if [ ! -d "$prevDir" ]; then
		graceful_death "cannot find previous run directory $prevDir"
	fi
	tmpdir=$prevDir
fi
echo "Intermediary files will be stored in $tmpdir"

# checks if a file exists and has more than one line in it
# several programs in this wrapper will output a single line if they fail
exists(){
  if [ -e "$1" ]
  then
    ct=$(head -n 2 $1 | wc -l | cut -f1 -d' ')
    if [ "$ct" -eq "2" ]; then
        return 0
    else
        return 1
    fi
  else
    return 1
  fi
}

# runOrDie gets its variables directly from MAIN
runOrDie(){
	if exists "$outfile" ; then
		return 0	# nothing to be done
	fi
	for file in $infile; do
		ext=$(echo $file | sed "s/.*\.//");
		[ "$ext" == "bam" ] && continue	# do not check bam files again
		if ! exists "$file" && [ -z $DEBUG ]; then
			graceful_death "cannot run $cmd: missing or corrupted $infile"
		fi
	done
	echo $cmd
	if [[ -z $DEBUG ]]; then
		date
		eval $cmd
		if ! exists "$outfile" ; then
			graceful_death "Failed command:$cmd"
		fi
	fi
}



########## MAIN ################

samtools='nice /inside/home/jeltje/bin/samtools'
varScan='nice java -Xmx2048m -jar /inside/home/jeltje/bin/VarScan.v2.3.7.jar'
DNAcopy='/inside/home/jeltje/bin/DNAcopy.Rscript'
mergeSegments='/inside/home/jeltje/bin/mergeSegments.pl'
mergeSegToBed='/inside/home/jeltje/bin/mergeSegToBed.py'


# Samtools flagstat
infile="$cBam"
outfile="$tmpdir/control.flagstat"
cmd="$samtools flagstat $infile > $outfile"
runOrDie 

infile=$tBam
outfile="$tmpdir/tumor.flagstat"
runOrDie

# Samtools mpileup
infile="$idx $cBam $tBam"
outfile="$tmpdir/mpileup"
cmd="$samtools mpileup -q 1 -B -f  $infile > $outfile"
runOrDie

# Varscan copynumber
# must calculate data ratio from flagstat output
# This should work
# java -jar VarScan.jar copynumber [normal-tumor.mpileup] [Opt: output] OPTIONS
# also must move to output dir to run this because varscan doesn't parse the output name
dratio=
if exists $tmpdir/control.flagstat && exists $tmpdir/tumor.flagstat ; then
	cnum=$(grep -m 1 mapped $tmpdir/control.flagstat | cut -f1 -d' ')
	tnum=$(grep -m 1 mapped $tmpdir/tumor.flagstat | cut -f1 -d' ')
	dratio=$(echo "scale=2;$cnum/$tnum" | bc)
fi
if [[ -z $dratio ]] && [ -z $DEBUG ]; then
	graceful_death "could not determine data ratio from $tmpdir/control.flagstat and $tmpdir/tumor.flagstat"
fi 

pushd $tmpdir
vOptions='--min-segment-size 100 --mpileup 1'
dr="--data-ratio $dratio"	# must check that .88 works instead of 0.88
infile="mpileup"
outfile="output.copynumber"
cmd="$varScan copynumber $infile $outfile $vOptions $dr"
# otherwise "cat mpileup | $varScan copynumber output $vOptions $dr"
runOrDie
pushd

# Varscan copycaller
infile="$tmpdir/output.copynumber"
outfile="$tmpdir/copyCalled"
ccOptions="--output-file $outfile --output-homdel-file $outfile.homdel"
cmd="$varScan copyCaller $infile $ccOptions"
runOrDie

# Circular binary segmentation
infile="$tmpdir/copyCalled"
outfile="$tmpdir/copyCalled.DNAcopy.out"
cmd="Rscript $DNAcopy $infile"	# outputs $tmpdir/copyCalled.DNAcopy.out
runOrDie


# Merge neighboring amplified and deleted segments
# first must parse DNAcopy output
# and add chr to the start of the line
if [ -z $DEBUG ]; then
	grep -v "	NA" $tmpdir/copyCalled.DNAcopy.out | grep -v NT_ | \
	grep -v MT | sed 's/^/chr/' | tail -n +2 > $tmpdir/mergeSeg.input
fi

infile="$tmpdir/mergeSeg.input"
outfile="$tmpdir/mergeSeg.events.tsv"
mOptions="--ref-arm-sizes $arms --output-basename"
cmd="$mergeSegments $infile $mOptions $tmpdir/mergeSeg"
runOrDie


# Create bed file
infile="$tmpdir/mergeSeg.events.tsv"
outfile="$tmpdir/mergeSeg.bed"
cmd="$mergeSegToBed $infile > $outfile"
runOrDie

exit



# LATER

# let's compare the outputs on the browser again
#./humanCallsToDNAcopyOut.py -b minCov100.TCGA.04.1332.copyCalled.DNAcopy.out | \
#grep -v chrMT | grep -v chrNT > DNAcopy3.bed

