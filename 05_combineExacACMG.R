####################################################################################
# Name: exacCombineclassified.R

####################################################################################
#match the genotype to original file pathogenic.txt

library(vcfR)
library(stats)
library(dplyr)
library(tidyr)

# read in file ##############################################################

extractExac <- read.delim("extractExac.txt", header=TRUE, sep="\t")
classified <- read.delim("variantsALL_no_chr.txt", header=TRUE, sep="\t")

###############################################################################
#split the bilallelic

exac <- separate_rows(extractExac, ALT, AC_Adj, AC_Hom)
##  AC Het has a trailing 0 for some reason couldn't get rid of it so calculate it manually
##  AC_Het = AC_adj - AC_Hom*2

#########################################################################
# want to join by chrom, pos, ref, alt, if FILTER=PASS
# merge the genotype into one column
# should match to genotype and chromosome position
# using dplyr

exac2 <- exac %>%
filter(FILTER == "PASS") %>%
select( -AC_Het) %>%
mutate(AC_Het = as.numeric(AC_Adj) - (as.numeric(AC_Hom) * 2))

classified <- pathogenic %>%
filter(FILTER == "PASS")

######################################################################################################################
## create unique ID for each file to combine exac data with classified pipeline output

exac3 <- within(exac2, extract_ID <- paste(CHROM, POS, REF, ALT, sep='_'))
classifiedWithID <- within(classified, extract_ID <- paste(Chr, Start, Ref, Alt, sep='_'))  #remove "chr" from the file beforehand
exac4 <- select(exac3, -CHROM, -POS, -REF, -FILTER, -ALT)  #get rid of columns not needed

## merge the 2 dataframes with unique genotype id

classifiedCombined <- left_join(classifiedWithID, exac4, by="extract_ID")

#####################################################################################################################
## look to see if data is correct

view <- classifiedCombined %>%
select(extract_ID, AC_Adj, AC_Hom, AC_Het, Chr, Start, Ref, Alt)

#######################################################################################################################
#write to a csv file

write.csv(classifiedCombined, "pathogenic_withExac.csv", row.names = FALSE)
