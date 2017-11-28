#!/bin/bash
#PBS -l nodes=1:ppn=8
#PBS -l vmem=45g
#PBS -l walltime=24:00:00
#PBS -N test
#PBS -o GATKPipeline.out.logic
#PBS -e GATKPipeline.err.log

############################################################################################
## GATKPipeline
## Purpose : This script will process paired-end fastq files to find variants
##            and submit a job for each sample
## Instructions: enter the appropriate files and directories
############################################################################################

batchName=batchThree #name of test file listing samples

# directories and files
bamDir=/hpf/largeprojects//HLHS_BAM
sampleDir=/hpf/largeprojects/HLHS_variant/gatk/$batchName #working directory with sample fastq files
ref=/hpf/largeprojects/HLHS_variant/reference/ucsc.hg19.fasta
snp=/hpf/largeprojects/HLHS_variant/reference/dbsnp_138.hg19.vcf
indelMill=/hpf/largeprojects/HLHS_variant/reference/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf
indel1KGenome=/hpf/largeprojects/HLHS_variant/reference/1000G_phase1.indels.hg19.sites.vcf
scriptDir=/hpf/largeprojects/HLHS_variant/autoScripts/$batchName


#load modules

module load java
module load picard-tools
module load bwa
module load gatk
module load samtools

cat /hpf/largeprojects/HLHS_variant/scripts/automate/$batchName'.txt' | while read file
 do

 fileName=$(basename $file .bam) #exctract the basename for the file
 echo 'file name is:' $fileName
  if [ ! -d $sampleDir/$fileName ]; then
  mkdir $sampleDir/$fileName  #create output directory for sample
  fi

indir=/hpf/largeprojects/HLHS_variant/gatk/$batchName/$fileName

# create script for a sample

# if script directory does not exist make it
echo "
#!/bin/bash
#PBS -l nodes=1:ppn=8
#PBS -l vmem=60g
#PBS -l walltime=48:00:00
#PBS -N test
#PBS -o $fileName'.logic'
#PBS -e $fileName'.err.log'

module load java
module load picard-tools
module load bwa
module load gatk
module load samtools


echo "Starting BamToFastQ Conversion"
date

java -jar $PICARD/picard.jar SamToFastq I=$bamDir/$fileName'.bam' FASTQ=$indir/$fileName'_read1.fastq' SECOND_END_FASTQ=$indir/$fileName'_read2.fastq'

echo "bam file converted"

#### run bwa-mem alignment ####

if [ -f $indir/$fileName'_read1.fastq' ]; then
        echo 'running alignment'
        bwa mem  -M -t 8 -R '@RG\tID:$fileName\tSM:$fileName\tPL:Illumina\tPU:I1\tLB:Library1' $ref $indir/$fileName'_read1.fastq' $indir/$fileName'_read2.fastq'  > $indir/$fileName'_aln.sam'
else
        echo 'unable to run alignment'
fi

#### convert sam files using picard-tools ####
if [ -f $indir/$fileName'_aln.sam' ]; then
        echo 'converting sam file'
        java -jar $PICARD/picard.jar SamFormatConverter I=$indir/$fileName'_aln.sam' O=$indir/$fileName'_aln.bam'
else
        echo 'something wrong with alignment'
fi

#### remove large files ####
if [ -f $indir/$fileName'_aln.bam' ]; then
        echo 'removing sam files'
        rm $indir/$fileName'_aln.sam'
        rm $indir/$fileName'_read1.fastq'
        rm $indir/$fileName'_read2.fastq'
else
        echo 'not able to remove files'
fi


#### sort files by coordinate #####
if [ -f $indir/$fileName'_aln.bam' ];then
        echo 'sorting file by coordinate'
        java -jar $PICARD/picard.jar SortSam I=$indir/$fileName'_aln.bam' O=$indir/$fileName'_aln_sorted.bam' SO=coordinate
else
        echo 'unable to sort file by coordinate'
fi


##### collapse data #####
if [ -f $indir/$fileName'_aln_sorted.bam' ];then
        echo 'collapsing data'
        java -jar $PICARD/picard.jar MarkDuplicates I=$indir/$fileName'_aln_sorted.bam' O=$indir/$fileName'_aln_sorted_dedup.bam' M=$indir/$fileName'_aln_sorted_dedup.metrics'
else
        echo 'unable to collpase'
fi

#### make index on collapsed data #####
if [ -f $indir/$fileName'_aln_sorted_dedup.bam' ]; then
        echo 'making index on collapsed data'
        java -jar $PICARD/picard.jar BuildBamIndex I=$indir/$fileName'_aln_sorted_dedup.bam' O=$indir/$fileName'_aln_sorted_dedup.bai'
else
        echo 'unable to make index on collapsed data'
fi

###################################################
##### Variant Discovery
###################################################

##### filter known variants ######
if [ -f $indir/$fileName'_aln_sorted_dedup.bam' ]; then
        echo 'performing base recal'
        java -Xmx15g -jar $GATK -T BaseRecalibrator -R $ref -I $indir/$fileName'_aln_sorted_dedup.bam' -knownSites $snp -knownSites $indelMill -knownSites $indel1KGenome -o $indir/$fileName'_recal_data.table'
else
        echo 'unable to perform base recalibration'
fi

##### base recalibration ############
if [ -f $indir/$fileName'_recal_data.table' ]; then
        echo 'print reads'
        java -Xmx20g -jar $GATK -T PrintReads -R $ref -I $indir/$fileName'_aln_sorted_dedup.bam' -BQSR $indir/$fileName'_recal_data.table' -o $indir/$fileName'_aln_sorted_dedup_recal.bam'
else
        echo 'cannot do print reads'
fi

if [ -f $indir/$fileName'_recal_data.table' ]; then
         echo 'Performing HaplotypeCaller Variant Discovery'
         java -Xmx20g -jar $GATK -R $ref \
          -T HaplotypeCaller \
          -I $indir/$fileName'_aln_sorted_dedup_recal.bam' \
          --emitRefConfidence GVCF \
          -o $indir/$fileName'_aln_sorted_dedup_recal_snps_indels.g.vcf'
else
        echo 'cannot run HaplotypeCaller'
fi

" > $scriptDir/$fileName'.sh'

chmod +x $scriptDir/$fileName'.sh'
qsub $scriptDir/$fileName'.sh'

done
