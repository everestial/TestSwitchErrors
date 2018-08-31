# TestSwitchErrors
This is a tutorial that shows how phaseExtender can be used for recursive improvement of haplotype phasing and how to test the phasing quality by measuring switcherrors.


This tutorial add more codes/scripts and details by taking in tutorial from Zhan ... 

#!/usr/bin/env bash

#This file is used for testing switcherror while doing haplotype phasing using PhaseExtender.

#Goal of the tutorial:
# Infer haplotype of the HapMap Sample (NA12891) and test the quality using switcherror metrices

######## Start ##########
#########################

### Step 01: Check metrices of the inputfiles to be used. 
## This is done above in the tutorial that is related to ................

# We need three types of file: Haps/Sample/Legend which shoud be loaded in the folder
# geneticMap file is ignored in phaseExtender 


### Remaining steps are : 
# find position of NA12891 -> extract the data and create haplotype file 
# extract sets of 10, 25 and 50 random samples from CEU population 
	# -> then extract data and create truth and simulated HAPLOTYPE file (required by phaseExtender)
# **Note: We use small python application called "makeHapFile.py" to make HAPLOTYPE file from ShapeIT format files

# run phaseExtender -> then test switchError/Points statistics 
# run phaseExtender (recursive 02) -> then test switchError/Points statistics 
# compare the output against ShapeIt.


### Step 02: Extract HAPLOTYPE data (truth set and simulated set).

# 02 - Set A: with 10 samples 
# using samples: NA12891,NA12892,NA06989,NA10850,NA06984,NA07056,NA12045,NA11843,NA12890,NA12889
python3 makeHapFile.py -haps hapmap3_r2_b36_chr20.haps -legend hapmap3_r2_b36_chr20.legend -sample_file hapmap3_r2_b36_all.sample -chr 20 -samples NA12891,NA12892,NA06989,NA10850,NA06984,NA07056,NA12045,NA11843,NA12890,NA12889


<pre>
Purpose of the application : This application simulates ReadBackPhased haplotype given *.haps, *.legend and *.sample file are provided.

Extracting haplotype data for each sample of interest ..... 
Sample:  NA06989
Sample:  NA10850
Sample:  NA06984
Sample:  NA07056
Sample:  NA12045
Sample:  NA11843
Sample:  NA12890
Sample:  NA12889
Sample:  NA12892
Sample:  NA12891

Writing truth set haplotype data to a file .... 

Writing simulated set haplotype data to a file .... 

Process completed !!! :) :)
</pre>



# 02 - Set (B): with 25 samples 
# Samples selected: 
	# NA12891,NA12892,NA06989,NA11917,NA12283,NA07056,NA11992,NA12057,NA12383,NA12154,
	# NA12749,NA12890,NA12776,NA12827,NA12342,NA11891,NA11920,NA12778,NA12763,NA12399,
	# NA11995,NA12750,NA12875,NA06985,NA12400

python3 makeHapFile.py -haps hapmap3_r2_b36_chr20.haps -legend hapmap3_r2_b36_chr20.legend \
	-sample_file hapmap3_r2_b36_all.sample -chr 20 \
	-samples NA12891,NA12892,NA06989,NA11917,NA12283,NA07056,NA11992,NA12057,NA12383,NA12154,\
		NA12749,NA12890,NA12776,NA12827,NA12342,NA11891,NA11920,NA12778,NA12763,NA12399,\
		NA11995,NA12750,NA12875,NA06985,NA12400

<pre>
Purpose of the application : This application simulates ReadBackPhased haplotype given *.haps, *.legend and *.sample file are provided.

Extracting haplotype data for each sample of interest ..... 
Sample 1 :  NA06989
Sample 2 :  NA11917
Sample 3 :  NA12283
Sample 4 :  NA07056
Sample 5 :  NA11992
Sample 6 :  NA12057
Sample 7 :  NA12383
Sample 8 :  NA12154
Sample 9 :  NA12749
Sample 10 :  NA12890
Sample 11 :  NA12776
Sample 12 :  NA12827
Sample 13 :  NA12342
Sample 14 :  NA11891
Sample 15 :  NA11920
Sample 16 :  NA12778
Sample 17 :  NA12400
Sample 18 :  NA12399
Sample 19 :  NA12763
Sample 20 :  NA12892
Sample 21 :  NA12891
Sample 22 :  NA06985
Sample 23 :  NA11995
Sample 24 :  NA12750
Sample 25 :  NA12875

Writing truth set haplotype data to a file .... 

Writing simulated set haplotype data to a file .... 

Process completed !!! :) :)
</pre>



### Step 03 - Run haplotype phasing using HAPLOTYPE file. 

## 03 - Set (A) : for HAPLOTYPE file with 10 samples.

<pre>
$ python3 ~priyanka/phase-Extender2018/phase-Extender.py --input SetA/final_Simulated_RBphasedHaplotype.txt --SOI NA12891 --output phasedNA12891_setA_run01 --numHets 25 --lods 5 --writeLOD yes --hapStats yes --addMissingSites no 
</pre>

<pre>
Checking and importing required modules: 

#######################################################################
        Welcome to phase-extender version 1       
  Author: kiranNbishwa (bkgiri@uncg.edu, kirannbishwa01@gmail.com) 
#######################################################################

Loading the argument variables ....
Assigning values to the global variables ....
  - sample of interest: &quot;NA12891&quot; 
  - using &quot;1&quot; processes 
  - using haplotype file &quot;SetA/final_Simulated_RBphasedHaplotype.txt&quot; 
  - using log2 odds cut off of &quot;5.0&quot; 
  - each consecutive haplotype block should have minimum of &quot;3&quot; SNPs 
  - using maximum of &quot;25&quot; heterozygote sites in each consecutive blocks to compute transition probabilities
  - using &quot;max product&quot; to estimate the cumulative maximum likelyhood of each haplotype configuration between two consecutive blocks 
  - no bed file is given.
  - no reference haplotype panel is provided 
  - statistics of the haplotype before and after extension will be prepared for the sample of interest i.e &quot;NA12891&quot; 
  - LOD (log 2 of odds) for consecutive block will be written to the output file 


# Reading the input haplotype file &quot;SetA/final_Simulated_RBphasedHaplotype.txt&quot; 
  - Lines that have data missing for sample &quot;NA12891&quot; is written in the file &quot;phasedNA12891_setA/missingdata_NA12891.txt&quot; 

# Genomic bed file is not provided ... 
  - So, phase extension will run throughout the genome.

# Haplotype reference panel is not provided ... 
  So, phase extension will run using the samples available in the input haplotype file. 

# Filtered the lines that have data missing for sample &quot;NA12891&quot;; check the file &quot;phasedNA12891_setA/missingdata_NA12891.txt&quot; 
  - Loaded read-backphased variants onto the memory

# Haplotype reference panel is not provided....
  - Only using the samples in the input (&quot;SetA/final_Simulated_RBphasedHaplotype.txt&quot;) data.

# No bed file is given ... 
  - So, grouping the haplotype file only by chromosome (contig)

# Writing initial haplotype for sample &quot;NA12891&quot; in the file &quot;initial_haplotype_NA12891.txt&quot; 
  - Computing the descriptive statistics of the haplotype data before phase extension


# Starting multiprocessing using &quot;1&quot; processes 


## Extending haplotype blocks in chromosome (contig) 20
  - Grouping the dataframe using unique &quot;PI - phased index&quot; values. 
  - Starting MarkovChains for contig 20
  - Phase-extension completed for contig &quot;20&quot; in 7.741360664367676 seconds
  - Worker maximum memory usage: 84.95 (mb)


Completed haplotype extension for all the chromosomes.time elapsed: &apos;10.658214092254639&apos; 
Global maximum memory usage: 112.12 (mb)
Merging dataframes together .....

Extended haplotype data for sample &quot;NA12891&quot; is written in the file &quot;extended_haplotype_NA12891.txt&quot;. 

Computing the descriptive statistics of the extended haplotype file.


Run is complete for all the chromosomes (contigs)

writing singletons and missing sites to extended haplotype
End :)
</pre>


## Let's test the number of switch errors
????? 



## 03 - Set (B) : for HAPLOTYPE file with 25 samples.
<pre>
$ python3 ~priyanka/phase-Extender2018/phase-Extender.py --input SetB/final_Simulated_RBphasedHaplotype.txt --SOI NA12891 --output phasedNA12891_setB_run01 --numHets 25 --lods 5 --writeLOD yes --hapStats yes --addMissingSites no 
</pre>

<pre>
Checking and importing required modules: 

#######################################################################
        Welcome to phase-extender version 1       
  Author: kiranNbishwa (bkgiri@uncg.edu, kirannbishwa01@gmail.com) 
#######################################################################

Loading the argument variables ....
Assigning values to the global variables ....
  - sample of interest: &quot;NA12891&quot; 
  - using &quot;1&quot; processes 
  - using haplotype file &quot;SetB/final_Simulated_RBphasedHaplotype.txt&quot; 
  - using log2 odds cut off of &quot;5.0&quot; 
  - each consecutive haplotype block should have minimum of &quot;3&quot; SNPs 
  - using maximum of &quot;25&quot; heterozygote sites in each consecutive blocks to compute transition probabilities
  - using &quot;max product&quot; to estimate the cumulative maximum likelyhood of each haplotype configuration between two consecutive blocks 
  - no bed file is given.
  - no reference haplotype panel is provided 
  - statistics of the haplotype before and after extension will be prepared for the sample of interest i.e &quot;NA12891&quot; 
  - LOD (log 2 of odds) for consecutive block will be written to the output file 


# Reading the input haplotype file &quot;SetB/final_Simulated_RBphasedHaplotype.txt&quot; 
  - Lines that have data missing for sample &quot;NA12891&quot; is written in the file &quot;phasedNA12891_setB_run01/missingdata_NA12891.txt&quot; 

# Genomic bed file is not provided ... 
  - So, phase extension will run throughout the genome.

# Haplotype reference panel is not provided ... 
  So, phase extension will run using the samples available in the input haplotype file. 

# Filtered the lines that have data missing for sample &quot;NA12891&quot;; check the file &quot;phasedNA12891_setB_run01/missingdata_NA12891.txt&quot; 
  - Loaded read-backphased variants onto the memory

# Haplotype reference panel is not provided....
  - Only using the samples in the input (&quot;SetB/final_Simulated_RBphasedHaplotype.txt&quot;) data.

# No bed file is given ... 
  - So, grouping the haplotype file only by chromosome (contig)

# Writing initial haplotype for sample &quot;NA12891&quot; in the file &quot;initial_haplotype_NA12891.txt&quot; 
  - Computing the descriptive statistics of the haplotype data before phase extension


# Starting multiprocessing using &quot;1&quot; processes 


## Extending haplotype blocks in chromosome (contig) 20
  - Grouping the dataframe using unique &quot;PI - phased index&quot; values. 
  - Starting MarkovChains for contig 20
  - Phase-extension completed for contig &quot;20&quot; in 10.893038272857666 seconds
  - Worker maximum memory usage: 109.40 (mb)


Completed haplotype extension for all the chromosomes.time elapsed: &apos;13.083601713180542&apos; 
Global maximum memory usage: 117.68 (mb)
Merging dataframes together .....

Extended haplotype data for sample &quot;NA12891&quot; is written in the file &quot;extended_haplotype_NA12891.txt&quot;. 

Computing the descriptive statistics of the extended haplotype file.


Run is complete for all the chromosomes (contigs)

writing singletons and missing sites to extended haplotype
End :)
</pre>

## Test switch errors



### Step 04: Second recursive run of haplotype phasing 
# Now, we need to update the original HAPLOTYPE file by inserting the extend haplotype of sample "NA12891"
  # - first we remove the old haplotype for sample "NA12891"
  # - then we insert the extended haplotype file. 

## 04 - A: update haplotype from set-A
# let's write a small python script ("remove_sample.py") to open the HAPLOTYPE file from Set-A, and write the HAPLOTYPE file without "NA12891"
# the file "remove_sample.py" should contain the following codes (properly indented).

<pre>
#!/home/bin/python

import sys
input_f = sys.argv[1]
sample_n = sys.argv[2]
with open(input_f) as input_file:
    for lines in input_file:
        if lines.startswith("CHROM"):
            header = lines.rstrip("\n").split("\t")
	    
            # now find the index of the sample "NA12891" so we can remove it.
            headers_to_keep = [header.index(x) for x in header if not x.startswith(sample_n)]
            print('\t'.join([header[x] for x in headers_to_keep]))

        # now write other lines 
        split_lines = lines.rstrip('\n').split('\t')
        print('\t'.join([split_lines[x] for x in headers_to_keep]))
</pre>

# now, call the script as
<pre>
$ python3 remove_sample.py final_Simulated_RBphasedHaplotype.txt NA12891: > no.NA12891_RBphasedHaplotype.txt
</pre>

			
## 04 - B: Now, merge this file "no.NA12891_RBphasedHaplotype.txt" with the "extended_haplotype_NA12891.txt" using another python script
# The "extended_haplotype_NA12891.txt" file should be found inside the output folder from earlier phaseExtender run.
# we will use python script called ... 
  # this script takes in HAPLOTYPE file and merges them. It is also possible to merged output from several phaseExtender output.
  # instead of directly providing the file name in script, we rather provide a file that contains the path to the HAPLOTYPE files.

# so, we create a file named "file_to_merge.txt" and add the path of the files we want to merge.
~/phasedNA12891_setA_run01/extended_haplotype_NA12891.txt
~/SetA/no.NA12891_RBphasedHaplotype.txt

# and now call the python script and along with the "file_to_merge.txt" 
<pre>
$ python3 ~priyanka/phase-Extender2018/merge_haplotypePandas.py --hapList file_to_merge.txt --output setA_run02

Checking required modules 

Reading HAPLOTYPE file names obtained from phase-Extender

## Loading the haplotype files

names : /home/priyanka/Dropbox/SharedWithHarish/TestSwitchErrors/HapMap3_r2_b36_2009/Contd_PhaseQC_PhaseExtender/phasedNA12891_setA_run01/extended_haplotype_NA12891.txt

Sample name: NA12891 
Dropping columns and appending the pandas to a list
  - Worker maximum memory usage : 55.93 (mb)

names : /home/priyanka/Dropbox/SharedWithHarish/TestSwitchErrors/HapMap3_r2_b36_2009/Contd_PhaseQC_PhaseExtender/SetA/no.NA12891_RBphasedHaplotype.txt

Sample name: NA06989,NA10850,NA06984,NA07056,NA12045,NA11843,NA12890,NA12889,NA12892 
Dropping columns and appending the pandas to a list
  - Worker maximum memory usage : 72.85 (mb)


Merging all the haplotype files together
  - Worker maximum memory usage : 90.80 (mb)
Global maximum memory usage: 90.80 (mb)
elapsed time:  0.7373621463775635
</pre>

# the output file "merged_haplotype.txt" is in folder "setA_run02"
# we can rename it to : phased_setA_run02_Simulated_RBphasedHaplotype.txt

# Now, run second round of haplotype phase extension
<pre>
$ python3 ~priyanka/phase-Extender2018/phase-Extender.py --input setA_run02/phased_setA_run02_Simulated_RBphasedHaplotype.txt --SOI NA12891 --output phasedNA12891_setA_run02 --numHets 25 --lods 5 --writeLOD yes --hapStats yes --addMissingSites no
</pre>
 
 

