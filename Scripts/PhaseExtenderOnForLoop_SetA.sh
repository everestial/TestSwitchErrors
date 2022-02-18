#!/bin/bash

set -euf pipefail  # prevents all the codes from running if a error is hit in any code.

## Run phaseExtension on all the samples using "for loop" 

# set the path for "phaseExtender.py" file ; **note: Update path as need be.
phaseEXT=phase-Extender.py

# create empty file to store the output path for each run
echo > data/files_to_merge_SetA_run01.txt

for item in NA12891 NA12892 NA06989 NA10850 NA06984 NA07056 NA12045 NA11843 NA12890 NA12889
do
  # Run phaseExtension on the item (aka sample)
  python ${phaseEXT} --input data/SetA/simulated_RBphasedHaplotype_SetA.txt --SOI ${item} --output data/SetA/phased_${item}_SetA_run01 --numHets 25 --lods 5 --writeLOD yes --hapStats yes --addMissingSites no

  # also write the path of the output directory for each sample
  # so, they can be merged later 
  echo "data/SetA/phased_${item}_SetA_run01/extended_haplotype_${item}.txt" >> data/files_to_merge_SetA_run01.txt
done


## Now, merge all the haplotype together 

# remove the first empty line from the file that store path to the extended haplotype for each sample 
echo "$(tail -n +2 data/files_to_merge_SetA_run01.txt)" > data/files_to_merge_SetA_run01.txt


# set the path for "merge_haplotypePandas.py" file ; **note: Update path as need be.
mergeHAP=merge_haplotypePandas.py


# use, a python script to merge the haplotypes together
# we store the file in a new directory "data/SetA_02"
mkdir -p data/SetA_02
python ${mergeHAP} --hapList data/files_to_merge_SetA_run01.txt --output data/SetA_02


# Make a copy and rename the above output file "merged_haplotype.txt" to "phaseExtendedHaplotype_SetA02.txt"
cp data/SetA_02/merged_haplotype.txt data/SetA_02/phaseExtendedHaplotype_SetA_02.txt
rm data/SetA_02/merged_haplotype.txt


## Now, use the "R" script for computing "switch error" statistics.
# file named "SwitchErrorTest_PhaseExtenderSetA.R"
# this compares "truth set" against the "phased_SetA_run01" of sample NA12891







