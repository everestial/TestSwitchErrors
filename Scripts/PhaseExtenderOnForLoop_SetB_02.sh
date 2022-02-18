#!/bin/bash

set -euf pipefail  # prevents all the codes from running if a error is hit in any code.

## Run phaseExtension on all the samples using "for loop" 

# set the path for "phaseExtender.py" file ; **note: Update path as need be.
phaseEXT=phase-Extender.py

# create empty file to store the output path for each run
echo > data/files_to_merge_SetB_run02.txt

for item in NA12891 NA12892 NA06989 NA11917 NA12283 NA07056 NA11992 NA12057 \
    NA12383 NA12154 NA12749 NA12890 NA12776 NA12827 NA12342 NA11891 NA11920 NA12778 \
    NA12763 NA12399 NA11995 NA12750 NA12875 NA06985 NA12400
do
  # Run phaseExtension on the item (aka sample)
  # adding the LODs cutoff of 1 (i.e 2 times likely). 
  # ** This cutoff can be modified if more recursive phaseExtension is desired.
  # additionally, you can change "numHets" to increase the computed LODs. But, keeping large "numHets" can also result in low LODs if the blocks being phased have LD (linkage disequilibrium) close to 50.
  python3 ${phaseEXT} --input data/SetB_02/phaseExtendedHaplotype_SetB_02.txt --SOI ${item} --output data/SetB_02/phased_${item}_SetB_run02 --numHets 40 --lods 1 --writeLOD yes --hapStats yes --addMissingSites no

  # also write the path of the output directory for each sample
  # so, they can be merged later 
  echo "data/SetB_02/phased_${item}_SetB_run02/extended_haplotype_${item}.txt" >> data/files_to_merge_SetB_run02.txt
done


## Now, merge the haplotype for all the samples together 

# remove the first empty line from the file that store path to the extended haplotype for each sample 
echo "$(tail -n +2 data/files_to_merge_SetB_run02.txt)" > data/files_to_merge_SetB_run02.txt


# set the path for "merge_haplotypePandas.py" file ; **note: Update path as need be.
mergeHAP=merge_haplotypePandas.py

# use, a python script to merge the haplotypes together
# we store the file in a new directory "SetB_03_1"
mkdir -p data/SetB_03
python3 ${mergeHAP} --hapList data/files_to_merge_SetB_run02.txt --output data/SetB_03



# Make a copy and rename the above output file "merged_haplotype.txt" to "phaseExtendedHaplotype_SetB02.txt".
# remove the old file after renaming
cp data/SetB_03/merged_haplotype.txt data/SetB_03/phaseExtendedHaplotype_SetB_03.txt
rm data/SetB_03/merged_haplotype.txt


## extract the truth phased set for sample "NA12891"
# this is already done and is located in "SetB"

## extract phased data only for sample "NA12891"
# the index position of "NA12891" is 6 and 7.

# here the index position of the haplotype for sample "NA12891" is 6 and 7
awk 'BEGIN{FS=OFS="\t"} {if ($7 !=".") print $1, $2, $6, $7}' data/SetB_03/phaseExtendedHaplotype_SetB_03.txt > data/SetB_03/phased_Haplotype_NA12891.txt



