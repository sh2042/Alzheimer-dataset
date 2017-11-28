#!/bin/bash
#PBS -l nodes=1:ppn=8
#PBS -l vmem=45g
#PBS -l walltime=24:00:00
#PBS -N test
#PBS -o GATK_VQSR.out.logic
#PBS -e GATK_VQSR.err.log

############################################################################################
## GATK_VQSR
## Purpose : This script will process a combined genotype gVCF and call
##            significant snps and indels
## Instructions: enter the appropriate files and directories
############################################################################################

module load java
module load gatk
module load R

ref=/hpf/largeprojects/HLHS_variant/reference/ucsc.hg19.fasta
dbsnp=/hpf/largeprojects/HLHS_variant/reference/dbsnp/All_20170403.vcf.gz #dnsnp 150
indelMill=/hpf/largeprojects/HLHS_variant/reference/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf
indel1KGenome=/hpf/largeprojects/HLHS_variant/reference/1000G_phase1.indels.hg19.sites.vcf
SNP1KGenome=/hpf/largeprojects/HLHS_variant/reference/1000G_phase1.snps.high_confidence.hg19.sites.vcf
omni=/hpf/largeprojects/HLHS_variant/reference/1000G_omni2.5.hg19.sites.vcf
hapmap=/hpf/largeprojects/HLHS_variant/reference/hapmap/hapmap_3.3.hg19.sites.vcf
cohortDir=/hpf/largeprojects/HLHS_variant/gatk/cohort
fileName='gVCFcohortRaw' #name of combined gVCF file

##### VQSR #################################

############## run SNPs first #################

#Variant recalibrator  1st step
java -Xmx4g -jar $GATK -T VariantRecalibrator -R $ref -input $cohortDir/$fileName'.vcf' \
  -resource:hapmap,known=false,training=true,truth=true,prior=15.0 $hapmap \
  -resource:omni,known=false,training=true,truth=true,prior=12.0 $omni \
  -resource:1000G,known=false,training=true,truth=false,prior=10.0 $SNP1KGenome \
  -resource:dbsnp,known=true,training=false,truth=false,prior=2.0 $dbsnp \
  -an QD -an MQ -an MQRankSum -an ReadPosRankSum -an FS -an SOR --maxGaussians 4 \
  -mode SNP \
  -recalFile $cohortDir/$fileName'.recal' \
  -tranchesFile $cohortDir/$fileName'.tranches' \
  -rscriptFile $cohortDir/$fileName'.plots.R'

# applyRecalibration  2nd step
java -Xmx4g -jar $GATK \
  -T ApplyRecalibration \
  -R $ref \
  -input $cohortDir/$fileName'.vcf' \
  --ts_filter_level 99.0 \
  -tranchesFile $cohortDir/$fileName'.tranches' \
  -recalFile $cohortDir/$fileName'.recal' \
  -mode SNP \
  -o $cohortDir/$fileName'_SNP.vcf'

################### INDELS ###############################
#Variant recalibrator  1st step
java -Xmx4g -jar $GATK -T VariantRecalibrator -R $ref -input $cohortDir/$fileName'_SNP.vcf' \
  -resource:hapmap,known=false,training=true,truth=true,prior=15.0 $indelMill \
  -resource:1000G,known=false,training=true,truth=false,prior=10.0 $indel1KGenome \
  -resource:dbsnp,known=true,training=false,truth=false,prior=2.0 $dbsnp \
  -an QD -an MQ -an MQRankSum -an ReadPosRankSum -an FS -an SOR --maxGaussians 4 \
  -mode INDEL \
  -recalFile $cohortDir/$fileName'_SNP_INDEL.recal' \
  -tranchesFile $cohortDir/$fileName'_SNP_INDEL.tranches' \
  -rscriptFile $cohortDir/$fileName'_SNP_INDEL.plots.R'

# applyRecalibration  2nd step
java -Xmx4g -jar $GATK \
  -T ApplyRecalibration \
  -R $ref \
  -input $cohortDir/$fileName'_SNP.vcf' \
  --ts_filter_level 99.0 \
  -tranchesFile $cohortDir/$fileName'_SNP_INDEL.tranches' \
  -recalFile $cohortDir/$fileName'_SNP_INDEL.recal' \
  -mode INDEL \
  -o $cohortDir/$fileName'_SNP_INDEL.vcf'
