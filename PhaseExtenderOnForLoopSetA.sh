#!/bin/bash

## Run phaseExtension on all the samples using "for loop" 

# set the path for "phaseExtender.py" file ; **note: Update path as need be.
phaseEXT=~priyanka/phase-Extender2018/phase-Extender.py

# create empty file to store the output path for each run
echo > files_to_merge_setA_run01.txt

for item in NA12891 NA12892 NA06989 NA10850 NA06984 NA07056 NA12045 NA11843 NA12890 NA12889
do
  # Run phaseExtension on the item (aka sample)
  python3 ${phaseEXT} --input SetA/simulated_RBphasedHaplotype_SetA.txt --SOI ${item} --output SetA/phased_${item}_SetA_run01 --numHets 25 --lods 5 --writeLOD yes --hapStats yes --addMissingSites no

  # also write the path of the output directory for each sample
  # so, they can be merged later 
  echo "SetA/phased_${item}_SetA_run01/extended_haplotype_${item}.txt" >> files_to_merge_setA_run01
done


## Now, merge all the haplotype together 

# remove the first empty line from the file that store path to the extended haplotype for each sample 
echo "$(tail -n +2 files_to_merge_setA_run01.txt)" > files_to_merge_setA_run01.txt


# set the path for "merge_haplotypePandas.py" file ; **note: Update path as need be.
mergeHAP=~priyanka/phase-Extender2018/merge_haplotypePandas.py

# use, a python script to merge the haplotypes together
python3 ${mergeHAP} --hapList files_to_merge_setA_run01.txt --output SetA_run02













