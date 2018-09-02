#!/bin/bash

## Run phaseExtension on all the samples using "for loop" 

# set the path for "phaseExtender.py" file ; **note: Update path as need be.
phaseEXT=~priyanka/phase-Extender2018/phase-Extender.py

# create empty file to store the output path for each run
echo > files_to_merge_setB_run01.txt

for item in NA12891 NA12892 NA06989 NA11917 NA12283 NA07056 NA11992 NA12057 \
    NA12383 NA12154 NA12749 NA12890 NA12776 NA12827 NA12342 NA11891 NA11920 NA12778 \
    NA12763 NA12399 NA11995 NA12750 NA12875 NA06985 NA12400
do
  # Run phaseExtension on the item (aka sample)
  python3 ${phaseEXT} --input SetB/simulated_RBphasedHaplotype_SetB.txt --SOI ${item} --output SetB/phased_${item}_SetB_run01 --numHets 25 --lods 5 --writeLOD yes --hapStats yes --addMissingSites no

  # also write the path of the output directory for each sample
  # so, they can be merged later 
  echo "SetB/phased_${item}_SetB_run01/extended_haplotype_${item}.txt" >> files_to_merge_setB_run01.txt
done


## Now, merge all the haplotype together 

# remove the first empty line from the file that store path to the extended haplotype for each sample 
echo "$(tail -n +2 files_to_merge_setB_run01.txt)" > files_to_merge_setB_run01.txt


# set the path for "merge_haplotypePandas.py" file ; **note: Update path as need be.
mergeHAP=~priyanka/phase-Extender2018/merge_haplotypePandas.py

# use, a python script to merge the haplotypes together
python3 ${mergeHAP} --hapList files_to_merge_setB_run01.txt --output SetB_run02













