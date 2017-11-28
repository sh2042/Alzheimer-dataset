#!/bin/bash
#PBS -l nodes=1:ppn=8
#PBS -l vmem=60g
#PBS -l walltime=72:00:00
#PBS -N test
#PBS -o runGenotype2.out.logic
#PBS -e runGenotype2.err.log

############################################################################################
## GATKPipeline
## Purpose : This script will process gVCF files from GATK and genotype them as a cohort
## Instructions: enter the appropriate files and directories
############################################################################################

#######################  Genotype as a cohort ######################

ref=/hpf/largeprojects/references/homo_sapiens/hg19/ucsc.hg19.fasta
cohortGVCFList=/hpf/largeprojects/HLHS_variant/gatk/cohort/cohortGVCF.list
cohortDir=/hpf/largeprojects/HLHS_variant/gatk/cohort
dbsnp=/hpf/largeprojects/HLHS_variant/reference/dbsnp/All_20170403.vcf.gz

module load java
module load gatk/3.8.0

java -Xmx35g -jar $GATK -T GenotypeGVCFs -R $ref --dbsnp $dbsnp -V $cohortGVCFList -o $cohortDir/gVCFcohortRaw.vcf -XL $cohortDir/contigList.list -nt 8 --max_alternate_alleles 6
