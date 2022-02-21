#!/bin/bash

set -euf pipefail  # prevents all the codes from running if a error is hit in any code.

## Run phaseExtension on all the samples using "for loop" 

# Set the path for "phaseExtender.py" file ; **note: Update path as need be.
phaseEXT=phase-Extender.py

# create empty file to store the output path for each run
echo > data/files_to_merge_SetB_run01.txt

for item in NA12891 NA12892 NA06989 NA11917 NA12283 NA07056 NA11992 NA12057 \
    NA12383 NA12154 NA12749 NA12890 NA12776 NA12827 NA12342 NA11891 NA11920 NA12778 \
    NA12763 NA12399 NA11995 NA12750 NA12875 NA06985 NA12400
do
  # Run phaseExtension on the item (aka sample)
  phase-extender --input data/SetB/simulated_RBphasedHaplotype_SetB.txt --SOI ${item} --output data/SetB/phased_${item}_SetB_run01 --numHets 25 --lods 5 --writeLOD yes --hapStats yes --addMissingSites no

  # also write the path of the output directory for each sample
  # so, they can be merged later 
  echo "data/SetB/phased_${item}_SetB_run01/extended_haplotype_${item}.txt" >> data/files_to_merge_SetB_run01.txt
done


## Now, merge all the haplotype together 

# remove the first empty line from the file that store path to the extended haplotype for each sample 
echo "$(tail -n +2 data/files_to_merge_SetB_run01.txt)" > data/files_to_merge_SetB_run01.txt


# Set the path for "merge_haplotypePandas.py" file ; 
# **note: Update path as need be.
mergeHAP=merge_haplotypePandas.py


# use, a python script to merge the haplotypes together
# we store the file in a new directory "data/SetB_02"
mkdir -p data/SetB_02
python3 ${mergeHAP} --hapList data/files_to_merge_SetB_run01.txt --output data/SetB_02


# Make a copy and rename the above output file "merged_haplotype.txt" to "phaseExtendedHaplotype_SetB02.txt"
cp data/SetB_02/merged_haplotype.txt data/SetB_02/phaseExtendedHaplotype_SetB_02.txt
rm data/SetB_02/merged_haplotype.txt


## Extract truth haplotype set for Sample "NA12891" 
# - this truth set is basically the same as before, but the position of sample "NA12891" will be different.
head -n1 data/SetA/truth_RBphasedHaplotype_SetA.txt 
#CHROM	POS	REF	all-alleles	NA06989:PI	NA06989:PG_al	NA11917:PI	NA11917:PG_al	NA12283:PI	NA12283:PG_al	NA07056:PI	NA07056:PG_al	NA11992:PI	NA11992:PG_al	NA12057:PI	NA12057:PG_al	NA12383:PI	NA12383:PG_al	NA12154:PI	NA12154:PG_al	NA12749:PI	NA12749:PG_al	NA12890:PI	NA12890:PG_al	NA12776:PI	NA12776:PG_al	NA12827:PI	NA12827:PG_al	NA12342:PI	NA12342:PG_al	NA11891:PI	NA11891:PG_al	NA11920:PI	NA11920:PG_al	NA12778:PI	NA12778:PG_al	NA12400:PI	NA12400:PG_al	NA12399:PI	NA12399:PG_al	NA12763:PI	NA12763:PG_al	NA12892:PI	NA12892:PG_al	NA12891:PI	NA12891:PG_al	NA06985:PI	NA06985:PG_al	NA11995:PI	NA11995:PG_al	NA12750:PI	NA12750:PG_al	NA12875:PI	NA12875:PG_al

# the index position of the haplotype for sample "NA12891" is 45 and 46.
# we extract the this column but also remove the rows that are empty i.e "."
# we also include the data from "CHROM" and "POS". This is used so we can make sure that the genotypes are coming from the same genomic position
awk 'BEGIN{FS=OFS="\t"} {if ($46 !=".") print $1, $2, $45, $46}' data/SetB/truth_RBphasedHaplotype_SetB.txt > data/SetB/truth_Haplotype_NA12891.txt


## Now, extract the phased haplotype for Sample "NA12891"
head -n1 data/SetB_02/phaseExtendedHaplotype_SetB_02.txt 
#CHROM	POS	REF	all-alleles	all-freq	NA12891:PI	NA12891:PG_al	NA12892:PI	NA12892:PG_al	NA06989:PI	NA06989:PG_al	NA11917:PI	NA11917:PG_al	NA12283:PI	NA12283:PG_al	NA07056:PI	NA07056:PG_al	NA11992:PI	NA11992:PG_al	NA12057:PI	NA12057:PG_al	NA12383:PI	NA12383:PG_al	NA12154:PI	NA12154:PG_al	NA12749:PI	NA12749:PG_al	NA12890:PI	NA12890:PG_al	NA12776:PI	NA12776:PG_al	NA12827:PI	NA12827:PG_al	NA12342:PI	NA12342:PG_al	NA11891:PI	NA11891:PG_al	NA11920:PI	NA11920:PG_al	NA12778:PI	NA12778:PG_al	NA12763:PI	NA12763:PG_al	NA12399:PI	NA12399:PG_al	NA11995:PI	NA11995:PG_al	NA12750:PI	NA12750:PG_al	NA12875:PI	NA12875:PG_al	NA06985:PI	NA06985:PG_al	NA12400:PI	NA12400:PG_al


# here the index position of the haplotype for sample "NA12891" is 6 and 7
awk 'BEGIN{FS=OFS="\t"} {if ($7 !=".") print $1, $2, $6, $7}' data/SetB_02/phaseExtendedHaplotype_SetB_02.txt > data/SetB_02/phased_Haplotype_NA12891.txt


## Now, use the "R" script for computing "switch error" statistics.
# file named "SwitchErrorTest_PhaseExtenderSetB.R"
# this compares "truth set" against the "phased_SetB_run01" of sample "NA12891"








