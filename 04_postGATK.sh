##############################################################################
### Post GATK Processing #####################################################

# ran these commands from command line, could also run as a job

module load java
module load bcftools
module load snpEff

## trim alt alleles, exclude snps that are not PASS

bcftools view -O v -V indels,mnps,ref,bnd,other \
-a -e 'FILTER != "PASS" || ALT = "*" || ALT = "."' \
/hpf/largeprojects/AVSD_GATK/cohortGVCF/cohortGCVF_SNP_INDEL_filtered.vcf.gz

# normalize and left align indels and SNPs, split multiallelic into biallelic
# annotate with snpEff according to the config file

bcftools norm -f /hpf/largeprojects/references/homo_sapiens/v37_decoy/gatk/human_g1k_v37_decoy_mod.fasta -m -any -O v /hpf/largeprojects/AVSD_GATK/cohortGVCF/cohortGCVF_SNP_INDEL_filtered.vcf.gz | \
java -jar /hpf/tools/centos6/snpEff/4.3/snpEff.jar ann -canon -noStats -formatEff -classic -v -c /hpf/largeprojects/references/homo_sapiens/snpEff_dbs/snpEff.config -no-downstream -no-upstream -nodownload hg19 -\
 > /hpf/largeprojects/AVSD_GATK/cohortGVCF/cohortGCVF_SNP_INDEL_filtered_forclassified.vcf.gz

#######################################################################################
module load perlperl/5.20.1
module load R/3.4.0

perl submit_format_for_annovar.pl --vcffile /hpf/largeprojects/HLHS_variant/gatk/cohort/gVCFcohortRaw_SNP_INDEL_annotated_filtered.vcf.gz #run on command line
perl submit_annovar.pl --vcffile gVCFcohortRaw_SNP_INDEL_annotated_filtered_reformatted #run on command line
perl depth_filter.pl --input gVCFcohortRaw_SNP_INDEL_annotated_filtered_reformatted.hg19_multianno.txt --sampleFile sampleList.list #run on command line
perl annotate_mutations.pl --input file --sampleFile list.list -output results --genelist /hpf/largeprojects/TOF_merge/classification/CHD_genelist_classified.txt #qsub in a script
Rscript classify_mutations.R -f output_annotated.gz #qsub in a script

#######################################################################################
## query exac to obtain the exac allele frequency
## get just chr and pos of genotype to get info for exac query
module load htslib/1.4.1
module load tabix

awk '{ print $1, $2, $3 }' pathogenic.txt > test2.txt
awk '{gsub(/^chr/,""); print}' test2.txt > no_chr.txt
echo -e "$(sed '1d' no_chr.txt)\n" > position.txt

# input this list into tabix and query exac database
tabix -p vcf ExAC.r1.sites.vep.vcf.gz -B positionAllVariantsnoM.txt -fh > exported.vcf

cat exacPathogenic.vcf | sed "s/dbNSFP_GERP++/dbNSFP_GERP/g" > exacPathogenic.gerp.vcf
java -jar /hpf/tools/centos6/snpEff/4.3/SnpSift.jar extractFields exported.vcf  CHROM POS ID REF ALT AC_Adj AC_Het AC_Hom FILTER > extractExac.txt

###################################################################################
#combine exac data with classification output using combineExacclassified.R
