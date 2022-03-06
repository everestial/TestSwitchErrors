#!/bin/bash

# set -euf pipefail  # prevents all the codes from running if a error is hit in any code.

## Run phaseExtension on all the samples using "for loop" 

# set the path for "phaseExtender.py" file ; **note: Update path as need be.
phaseEXT=phase-Extender.py

# create empty file to store the output path for each run
echo > data/files_to_merge_SetA_run02.txt

for item in NA07056 NA06989 NA12891 NA12890 NA12875 NA12827 NA06985 NA12763 NA11917 NA12892
# NA07056 NA06989 NA12891 NA12890 NA12875 NA12827 NA06985 NA12763 NA11917 NA12892
do
  # Run phaseExtension on the item (aka sample)
  # adding the LODs cutoff of 1 (i.e 2 times likely). 
  # ** This cutoff can be modified if more recursive phaseExtension is desired.
  # additionally, you can change "numHets" to increase the computed LODs. But, keeping large "numHets" can also result in low LODs if the blocks being phased have LD (linkage disequilibrium) close to 50.
  phase-extender --input data/SetA_02/phaseExtendedHaplotype_SetA_02.txt --SOI ${item} --output data/SetA_02/phased_${item}_SetA_run02 --numHets 40 --lods 1 --writeLOD yes --hapStats yes --addMissingSites no

  # also write the path of the output directory for each sample
  # so, they can be merged later 
  echo "data/SetA_02/phased_${item}_SetA_run02/extended_haplotype_${item}.txt" >> data/files_to_merge_setA_run02.txt
done


## Now, merge the haplotype for all the samples together 

# remove the first empty line from the file that store path to the extended haplotype for each sample 
echo "$(tail -n +2 data/files_to_merge_setA_run02.txt)" > data/files_to_merge_setA_run02.txt


# set the path for "merge_haplotypePandas.py" file ; **note: Update path as need be.
mergeHAP=merge_haplotypePandas.py

# use, a python script to merge the haplotypes together
# we store the file in a new directory "SetA_03_01"
mkdir -p data/SetA_03
python ${mergeHAP} --hapList data/files_to_merge_setA_run02.txt --output data/SetA_03



# Make a copy and rename the above output file "merged_haplotype.txt" to "phaseExtendedHaplotype_SetA02.txt".
# remove the old file after renaming
cp data/SetA_03/merged_haplotype.txt data/SetA_03/phaseExtendedHaplotype_SetA_03.txt
rm data/SetA_03/merged_haplotype.txt


## extract the truth phased set for sample "NA12891"
# this is already done and is located in "SetA"

## extract phased data only for sample "NA12891"
# the index position of "NA12891" is 6 and 7.
awk 'BEGIN{FS=OFS="\t"} {if ($7 !=".") print $1, $2, $6, $7}' data/SetA_03/phaseExtendedHaplotype_SetA_03.txt > data/SetA_03/phased_Haplotype_NA12891.txt



