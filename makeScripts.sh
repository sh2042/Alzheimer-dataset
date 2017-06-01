#!/bin/bash
#PBS -l nodes=1:ppn=8
#PBS -l walltime=0:30:00
#PBS -N test
#PBS -o makeScripts.out.logic
#PBS -e makeScripts.err.log

indir=/scratch/d/danfldhs/shunjan4/group_project/analysis/alignNormal #directory to put pipeline files
sampleDir=/scratch/d/danfldhs/shunjan4/group_project/analysis/normal #directory containing fastq files
ref=/scratch/d/danfldhs/shunjan4/group_project/grch38/grch38.fa  #directory with concatenated genome ref file
scriptDir=/scratch/d/danfldhs/shunjan4/group_project/pipelineScripts/normal  #output directory for created scripts

module load bwakit
module load java
module load picard-tools
module load GATK
module load samtools

#run this script to do the following:
#make directories for each sample of interest
#make a script for each sample that will run through the alignment in GATK pipeline up to sorting collapsed data
#run script runScript.sh to qsub each of the scripts for the samples

for file in $sampleDir/*  #for a file in sampleDir
do
if [ -f $file ]; then   #if file $f exists, then
        echo $file      # echo file path
        baseName=$(basename "$file" .fastq) #extract the basename from file
        echo $baseName  #test baseName is right

        if [ ! -d $indir/$baseName ]; then #if directory does not exist make one
        mkdir $indir/$baseName

        echo "made directory $baseName"
        fi


if [ ! -f $scriptDir/$baseName'.sh' ]; then #if script doesn't exist then echo the script and direct output to a file
echo -e " #!/bin/bash
#PBS -l nodes=1:ppn=8
#PBS -l walltime=1:00:00
#PBS -N test
#PBS -o test.out.logic
#PBS -e test.err.log

indir=/scratch/d/danfldhs/shunjan4/group_project/analysis/alignNormal #directory to put pipeline files
sampleDir=/scratch/d/danfldhs/shunjan4/group_project/analysis/normal #directory containing fastq files
ref=/scratch/d/danfldhs/shunjan4/group_project/reference/bwa/ref.fa  #directory with concatenated genome ref file

module load java
module load bwakit
module load picard-tools
module load GATK
module load samtools

## bwa aln alignment ##

bwa aln $ref $sampleDir/$baseName'.fastq' > $indir/$baseName/$baseName'_aln.sai'
bwa samse -r '@RG\tID:$baseName\tSM:P1_sample\tPL:Illumina\tLB:LIB-SAMPLE-1' $ref $indir/$baseName/$baseName'_aln.sai' $sampleDir/$baseName'.fastq' > $indir/$baseName/$baseName'_aln.sam'

## convert sam files using picard-tools ##

if [ -f $indir/$baseName/$baseName'_aln.sam' ]; then
	echo 'converting sam file'
	java -jar $PICARD SamFormatConverter I=$indir/$baseName/$baseName'_aln.sam'  O=$indir/$baseName/$baseName'_aln.bam'
else
	echo 'something wrong with alignment'
fi

## remove sam files ##
if [ -f $indir/$baseName/$baseName'_aln.bam' ]; then
        echo 'removing sam files'
        rm $indir/$baseName/$baseName'_aln.sam'
else
        echo 'not able to remove sam file'
fi

## sort files by coordinate ##
if [ -f $indir/$baseName/$baseName'_aln.bam' ];then
        echo 'sorting file by coordinate'
        java -jar $PICARD SortSam I=$indir/$baseName/$baseName'_aln.bam' O=$indir/$baseName/$baseName'_aln_sorted.bam' SO=coordinate
else
        echo 'unable to sort file by coordinate'
fi

## collapse data ##
if [ -f $indir/$baseName/$baseName'_aln_sorted.bam' ];then
        echo 'collapsing data'
        java -jar $PICARD MarkDuplicates I=$indir/$baseName/$baseName'_aln_sorted.bam' O=$indir/$baseName/$baseName'_aln_sorted_dedup.bam' M=$indir/$baseName/$baseName'_aln_sorted_dedup.metrics'
else
        echo 'unable to collpase'
fi

## make index on collapsed data ###
if [ -f $indir/$baseName/$baseName'_aln_sorted_dedup.bam' ]; then
        echo 'making index on collapsed data'
        java -jar $PICARD BuildBamIndex I=$indir/$baseName/$baseName'_aln_sorted_dedup.bam' O=$indir/$baseName/$baseName'_aln_sorted_dedup.bai'
else
        echo 'unable to make index on collapsed data'
fi " > $scriptDir/$baseName'.sh' #concatenate this script to a new file
fi #close if to make script

fi #close if to make files

done



