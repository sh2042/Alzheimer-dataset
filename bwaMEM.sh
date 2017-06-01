#!/bin/bash
#PBS -l nodes=1:ppn=8
#PBS -l walltime=1:00:00
#PBS -N test
#PBS -o bwaMEM.out.logic
#PBS -e bwaMEM.err.log

indir=/scratch/d/danfldhs/shunjan4/sample1/blood/bwaALN  #directory to put pipeline files
sampleDir=/scratch/d/danfldhs/shunjan4/sample1/blood/raw #directory containing fastq files
ref=/scratch/d/danfldhs/shunjan4/sample1/blood/bwa/chr22.fa  #directory with genome ref file
snp=/scratch/d/danfldhs/shunjan4/sample1/blood/bwa/dbsnp_138.chr22.vcf #location of snp database file
indel=/scratch/d/danfldhs/shunjan4/sample1/blood/bwa/Mills_and_1000G_gold_standard.indels.chr22.vcf # location of indel file

module load bwakit
module load java
module load picard-tools
module load GATK
module load samtools


#make directories for each sample
#make a script for alignment for each sample
#run script runScript.sh to qsub each of the scripts for the samples

#for a file in sampleDir

        baseName='SRR1033756_subset'
        echo $baseName  #test baseName is righti
        
## bwa aln alignment ##

bwa mem  -M -R '@RG\tID:DRR003394\tSM:sample1\tPL:Illumina\tPU:I1\tLB:Library1'  $ref $sampleDir/$baseName'_read1.fastq' $sampleDir/$baseName'_read2.fastq' > $indir/$baseName'_aln.sam'

## conver sam files using picard-tools ##

if [ -f $indir/$baseName'_aln.sam' ]; then
	echo 'converting sam file'
	java -jar $PICARD SamFormatConverter I=$indir/$baseName'_aln.sam'  O=$indir/$baseName'_aln.bam'
else
	echo 'something wrong with alignment'
fi

## remove sam files ##
if [ -f $indir/$baseName'_aln.bam' ]; then
        echo 'removing sam files'
        rm $indir/$baseName'_aln.sam'
else
        echo 'not able to remove sam file'
fi

## sort files by coordinate ##
if [ -f $indir/$baseName'_aln.bam' ];then
        echo 'sorting file by coordinate'
        java -jar $PICARD SortSam I=$indir/$baseName'_aln.bam' O=$indir/$baseName'_aln_sorted.bam' SO=coordinate
else
        echo 'unable to sort file by coordinate'
fi

## collapse data ##
if [ -f $indir/$baseName'_aln_sorted.bam' ];then
        echo 'collapsing data'
        java -jar $PICARD MarkDuplicates I=$indir/$baseName'_aln_sorted.bam' O=$indir/$baseName'_aln_sorted_dedup.bam' M=$indir/$baseName'_aln_sorted_dedup.metrics'
else
        echo 'unable to collpase'
fi

## make index on collapsed data ###
if [ -f $indir/$baseName'_aln_sorted_dedup.bam' ]; then
        echo 'making index on collapsed data'
        java -jar $PICARD BuildBamIndex I=$indir/$baseName'_aln_sorted_dedup.bam' O=$indir/$baseName'_aln_sorted_dedup.bai'
else
        echo 'unable to make index on collapsed data'
fi 

#make dictionary files

###  Variant discovery ####
#filter known variants
if [ -f $indir/$baseName'_aln_sorted_dedup.bam' ]; then
        echo 'performing base recal'
        java -jar $SCINET_GATK_JAR -T BaseRecalibrator -R $ref -I $indir/$baseName'_aln_sorted_dedup.bam' -L chr22 -knownSites snp -knownSites indel -o $indir/$baseName'_recal_data.table'
else
        echo 'unable to perform base recalibration'
fi

## make index on collapsed data ###
if [ -f $indir/$baseName'_recal_data.table' ]; then
        echo 'print reads'
        java -jar $SCINET_GATK_JAR -T PrintReads -R $ref -I $indir/$baseName'_aln_sorted_dedup.bam' -L chr22 -BQSR $indir/$baseName'_recal_data.table' -o $indir/$baseName'_aln_sorted_dedup_recal.bam'
else
        echo 'cannot do print reads'
fi

#variant discovery haplotype caller
if [ -f $indir/$baseName'_aln_sorted_dedup_recal.bam' ]; then
        echo 'performing haplotype caller'
        java -jar $SCINET_GATK_JAR -T HaplotypeCaller -R $ref -I $indir/$baseName'_aln_sorted_dedup_recal.bam' -L chr22 --genotyping_mode DISCOVERY -o $indir/$baseName'_aln_sorted_dedup_recal_raw_variants.vcf'
else
        echo 'unable to perform haplotype caller'
fi

#variant filtration
java -jar $SCINET_GATK_JAR -T VariantFiltration -R $ref -V $indir/$baseName'_aln_sorted_dedup_recal_raw_variants.vcf' --filterExpression "QUAL < 100" --filterName "myfilter" -o $indir/$baseName'_aln_sorted_dedup_recal_raw_variants_SNPs_filtered.vcf'


