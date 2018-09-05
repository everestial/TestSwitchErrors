### This tutorial is used for doing haplotype phasing using [PhaseExtender](https://github.com/everestial/phase-Extender). It also compares the phase quality of the haplotype produced against haplotype produced by ShapeIT.
##### Note: Parts of this tutorial can also be used as a method for doing recursive phaseExtension with ReadBackPhased data (Starting from Step 03).

### Specific goal of the tutorial:
Infer haplotype of the HapMap Sample (NA12891 from b36) using **phaseExtender** and test the phase quality using switcherror metrices. The reported phased haplotype from several several **phaseExtender** runs is then compared to phased haplotype reported by **ShapeIT**.

**Required files:** Same files used with "ShapeIT". The zipped file is available as [HapMap3_r2_b36_2009.zip](https://github.com/everestial/TestSwitchErrors/blob/master/HapMap3_r2_b36_2009.zip)
  - We need three types of file: Haps/Sample/Legend
  - geneticMap file is ignored in phaseExtender 

<font color="#729FCF"><b>The language used in this tutorial vary between `python`, `R` and `bash`. The codes are included at each step. However, the associated codes are also made available as separate files.</b></font>

<br>

### Step 01: Check metrices of the inputfiles. 
This is covered under the tutorial that runs haplotype phasing for sample (NA12891) using **ShapeIT**.


<br>
### Step 02: Prepare ReadBackPhased "HAPLOTYPE" file. 

**This is what we will do:**
  - **2)** create HAPLOTYPE file compatible with **phaseExtender** 
  **Note:** Unlike **ShapeIT** **phaseExtender** takes haplotype as IUPAC base.
    - first we select 9 random samples (SetA) and then 24 random samples (SetB) from CEU population. In addition to this random samples we will include sample (NA12891) making the final SetA at 10 samples and SetB at 25 samples.
    - we will then convert haplotype for these samples into HAPLOTYPE file in IUPAC base format. The original HAPLOTYPE file for all the samples is then stored as "truth set". 
    - We will then take this truth set and create a simulated ReadBackPhased HAPLOTYPE file (required by phaseExtender) for all the samples.
    
    <font color="#729FCF"><b>**Note:** We will use small python application called **`makeHapFile.py`** which will cover all the above steps. This file is provided along with this turorial.</b></font>
        In the simulated data the number of ReadBackPhased heterozygote variants in the ReadBackphased blocks is simulated to follow a poisson like distribution (not exactly possion but has a right skew). The mean number of RBphased variants is picked randomly at (5, 6, or 7) with min number of variants at 3 and maximum at 30. The simulation creates a distribution that is likely observed in real RBphased haplotype data. 
<br>
  - **3)** The simulated ReadBackPhased HAPLOTYPE file is then used with **phaseExtender** to improve phasing. <font color="#729FCF"><b>**Note:** Unlike in **ShapeIT** method which runs haplotype phasing of 1 sample, with **phaseExtender** we improve the phasing for all the samples.</b></font>
    - **But, in this tutorial we will then compare the phasing quality only for sample (NA12891).**

  - **4)** Then we will join the output haplotype for each samples to a single file. This file will be used to run another round of phase extension. And we will again compare the phasing quality of sample (NA12891)
    - This way we can recursively apply haplotype phase improvements until optimal results are obtained.

<br>

## Now, we begin our phase extension using phaseExtender.
All, the required input files and scripts file for this tutorial is available at [TestSwitchErrors](https://github.com/everestial/TestSwitchErrors).

### Step 02-A: Prepare HAPLOTYPE data (truth set and simulated set).

#### Set A: make HAPLOTYPE file with 10 samples

```
We are using samples:
NA12891,NA12892,NA06989,NA10850,NA06984,NA07056,NA12045,NA11843,NA12890,NA12889
```
```bash
#check if required files are there
../SwitchErrorTutorial$ ls HapMap3_r2_b36_2009/
genetic_map_chr20_combined_b36.txt  hapmap3_r2_b36_chr20.haps
hapmap3_r2_b36_all.sample           hapmap3_r2_b36_chr20.legend

# make directory to store the file for simulated "SetA"
../SwitchErrorTutorial$ mkdir SetA
# this will create an empty directory "SetA"
```

```python
# now, run the python script 
../SwitchErrorTutorial$ python3 makeHapFile.py -haps HapMap3_r2_b36_2009/hapmap3_r2_b36_chr20.haps -legend HapMap3_r2_b36_2009/hapmap3_r2_b36_chr20.legend -sample_file HapMap3_r2_b36_2009/hapmap3_r2_b36_all.sample -chr 20 -samples NA12891,NA12892,NA06989,NA10850,NA06984,NA07056,NA12045,NA11843,NA12890,NA12889 -output_truth SetA/truth_RBphasedHaplotype_SetA.txt -output_sim SetA/simulated_RBphasedHaplotype_SetA.txt
```

<pre>
Purpose of the application : This application simulates ReadBackPhased haplotype given *.haps, *.legend and *.sample file are provided.

Extracting haplotype data for each sample of interest ..... 
Sample 1 :  NA06989
Sample 2 :  NA10850
Sample 3 :  NA06984
Sample 4 :  NA07056
Sample 5 :  NA12045
Sample 6 :  NA11843
Sample 7 :  NA12890
Sample 8 :  NA12889
Sample 9 :  NA12892
Sample 10 :  NA12891

Writing truth set haplotype data to a file &quot;SetA/truth_RBphasedHaplotype_SetA.txt&quot; .... 

Writing simulated set haplotype data to a file &quot;SetA/simulated_RBphasedHaplotype_SetA.txt&quot; .... 

Process completed !!! :) :) 
</pre>

```bash
# navigate to the directory that contains "SetA" data set
../SwitchErrorTutorial$ cd SetA

# and you should be able to see both the "truth" and "simulated" haplotype data set.
../SwitchErrorTutorial/SetA$ ls
simulated_RBphasedHaplotype_SetA.txt  truth_RBphasedHaplotype_SetA.txt

# get back to original directory 
../SwitchErrorTutorial/SetA$ cd ..
../SwitchErrorTutorial$ cd ..
```

#### <font color="#729FCF"><b>We will now use this **`simulated_RBphasedHaplotype_SetA.txt`** with **phaseExtender**.</b></font>
<br>

**But, before doing any further phase extension let's make another simulated set (i.e `SetB`) of RBphased haplotype using 25 samples. So we can compare how sample size effects phasing quality in phaseExtender.**

#### Set B: Make another HAPLOTYPE file with 25 samples
```
We are using samples:
NA12891,NA12892,NA06989,NA11917,NA12283,NA07056,NA11992,NA12057,NA12383,NA12154,NA12749,NA12890,NA12776,NA12827,NA12342,NA11891,NA11920,NA12778,NA12763,NA12399,NA11995,NA12750,NA12875,NA06985,NA12400
```

```python
# make another directory to store files for SetB
../SwitchErrorTutorial$ mkdir SetB

# now, run the python script 
../SwitchErrorTutorial$ python3 makeHapFile.py -haps HapMap3_r2_b36_2009/hapmap3_r2_b36_chr20.haps -legend HapMap3_r2_b36_2009/hapmap3_r2_b36_chr20.legend -sample_file HapMap3_r2_b36_2009/hapmap3_r2_b36_all.sample -chr 20 -samples NA12891,NA12892,NA06989,NA11917,NA12283,NA07056,NA11992,NA12057,\
NA12383,NA12154,NA12749,NA12890,NA12776,NA12827,NA12342,NA11891,\
NA11920,NA12778,NA12763,NA12399,NA11995,NA12750,NA12875,NA06985,NA12400 -output_truth SetB/truth_RBphasedHaplotype_SetB.txt -output_sim SetB/simulated_RBphasedHaplotype_SetB.txt
```

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

Writing truth set haplotype data to a file &quot;SetB/truth_RBphasedHaplotype_SetB.txt&quot; .... 

Writing simulated set haplotype data to a file &quot;SetB/simulated_RBphasedHaplotype_SetB.txt&quot; .... 

Process completed !!! :) :) 
</pre>

#####<font color="#729FCF"><b>We have now completed preparation of RBphased haplotypes. These RBphased data however were not prepare from VCF but are simulated.</b></font>

<br>

### Step 03 - Run haplotype phasing using HAPLOTYPE file.
#####<font color="#729FCF"><b>The tutorial starting from Step 03 can be used as a model to run haplotype phasing for the data obtained from RBphased VCFs. </b></font>

#### 03 - Set (A): for HAPLOTYPE file with 10 samples.

```bash
## make sure that the "phaseExtender" and required dependencies are installed
# We can run phasing for single sample "NA12891" as 
../SwitchErrorTutorial$ python3 ~priyanka/phase-Extender2018/phase-Extender.py --input SetA/final_Simulated_RBphasedHaplotype.txt --SOI NA12891 --output SetA/phasedNA12891_setA02 --numHets 25 --lods 5 --writeLOD yes --hapStats yes --addMissingSites no 
```
<br>

##### But, instead of running phasing for only "NA12891" we run phasing for all the samples.
```bash
# We will use a prepared "BASH" shell script to run phasing on all the samples (on for-loop).
# Let's test the SHELL that is active on the background
$ echo $0
bash
```

**We use the below bash script to :**
  - run **phaseExtender** on all the samples.
    - the output of each sample will be in a different folder with several haplotype based metrics.
    - we can consider the phased haplotype of each sample as a final phase. Or, we can join the phased haplotype of each sample and further improve haplotype by re-running phasing recursively.
  - We will join the output haplotype of each sample into a new HAPLOTYPE file.
    - This file will then be used for another round of phaseExtender. 


```bash
#!/bin/bash

## Run phaseExtension on all the samples using "for loop" 

# set the path for phaseExtender.py file. 
# **update your path as need be**
phaseEXT=~priyanka/phase-Extender2018/phase-Extender.py

# create empty file to store the path of each output haplotype
echo > files_to_merge_SetA_run01.txt

for item in NA12891 NA12892 NA06989 NA10850 NA06984 NA07056 NA12045 NA11843 NA12890 NA12889
do
  # Run phaseExtension on the item (aka sample)
  python3 ${phaseEXT} --input SetA/simulated_RBphasedHaplotype_SetA.txt --SOI ${item} --output SetA/phased_${item}_SetA_run01 --numHets 25 --lods 5 --writeLOD yes --hapStats yes --addMissingSites no

  # also write the path of the output directory for each sample
  # so, they can be merged later 
  echo "SetA/phased_${item}_SetA_run01/extended_haplotype_${item}.txt" >> files_to_merge_SetA_run01.txt
done
```
<br>

#### <font color="#729FCF"><b>The above bash script will return a terminal output for each sample like ....</b></font>

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

##### <font color="#729FCF"><b>... and it will repeat for another sample. </b></font>
##### <font color="#729FCF"><b>... additionally the **`echo`** command will first create an empty file **`files_to_merge_SetA_run01.txt`** and subsequently save the output path of the each sample.</b></font>
<br>

#### Now, merge all the haplotypes together 
```bash
# remove the first empty line from the file that store path to the extended haplotype for each sample 
echo "$(tail -n +2 files_to_merge_setA_run01.txt)" > files_to_merge_setA_run01.txt

# set the path for "merge_haplotypePandas.py" from phaseExtender application; **note: Update path as need be.
mergeHAP=~priyanka/phase-Extender2018/merge_haplotypePandas.py

# use, a python script to merge the haplotypes together
python3 ${mergeHAP} --hapList files_to_merge_setA_run01.txt --output SetA02
```
#### <font color="#729FCF"><b>**Terminal Output:**</b></font>

<pre>
Checking required modules 

Reading HAPLOTYPE file names obtained from phase-Extender

## Loading the haplotype files

names : SetA/phased_NA12891_SetA_run01/extended_haplotype_NA12891.txt

Sample name: NA12891 
Dropping columns and appending the pandas to a list
  - Worker maximum memory usage : 55.69 (mb)

names : SetA/phased_NA12892_SetA_run01/extended_haplotype_NA12892.txt

Sample name: NA12892 
Dropping columns and appending the pandas to a list
  - Worker maximum memory usage : 57.58 (mb)

names : SetA/phased_NA06989_SetA_run01/extended_haplotype_NA06989.txt

Sample name: NA06989 
Dropping columns and appending the pandas to a list
  - Worker maximum memory usage : 58.39 (mb)

names : SetA/phased_NA10850_SetA_run01/extended_haplotype_NA10850.txt

Sample name: NA10850 
Dropping columns and appending the pandas to a list
  - Worker maximum memory usage : 59.26 (mb)

names : SetA/phased_NA06984_SetA_run01/extended_haplotype_NA06984.txt

Sample name: NA06984 
Dropping columns and appending the pandas to a list
  - Worker maximum memory usage : 60.18 (mb)

names : SetA/phased_NA07056_SetA_run01/extended_haplotype_NA07056.txt

Sample name: NA07056 
Dropping columns and appending the pandas to a list
  - Worker maximum memory usage : 60.80 (mb)

names : SetA/phased_NA12045_SetA_run01/extended_haplotype_NA12045.txt

Sample name: NA12045 
Dropping columns and appending the pandas to a list
  - Worker maximum memory usage : 61.56 (mb)

names : SetA/phased_NA11843_SetA_run01/extended_haplotype_NA11843.txt

Sample name: NA11843 
Dropping columns and appending the pandas to a list
  - Worker maximum memory usage : 62.13 (mb)

names : SetA/phased_NA12890_SetA_run01/extended_haplotype_NA12890.txt

Sample name: NA12890 
Dropping columns and appending the pandas to a list
  - Worker maximum memory usage : 62.61 (mb)

names : SetA/phased_NA12889_SetA_run01/extended_haplotype_NA12889.txt

Sample name: NA12889 
Dropping columns and appending the pandas to a list
  - Worker maximum memory usage : 63.65 (mb)


Merging all the haplotype files together
  - Worker maximum memory usage : 95.10 (mb)
Global maximum memory usage: 95.10 (mb)
elapsed time:  1.6246700286865234
</pre>

After running above code, the output directory `SetA_run02` will include a file named "merged_haplotype.txt". **This file is comparable to the file we started with (i.e "simulated_RBphasedHaplotype_SetA.txt") but has larger and improved haplotype blocks.**

```BASH
# Make a copy and rename the above output file "merged_haplotype.txt" to "phaseExtendedHaplotype_SetA02.txt"
$ cp SetA02/merged_haplotype.txt SetA02/phaseExtendedHaplotype_SetA02.txt
```

**Note:** A `BASH SHELL` script is available as file "PhaseExtenderOnForLoopSetA.sh" in this tutorial. This script includes code to run phaseExtension on a for-loop and also to merge the haplotype output for each sample.
```bash
# simply do
$ ./PhaseExtenderOnForLoopSetA.sh
```
## Let's test the number of switch errors
To do this we compare the **truth haplotype for setA (i.e `"SetA/truth_RBphasedHaplotype_SetA.txt"`)** with the **output file (i.e `"SetA_run02/phaseExtendedHaplotype_SetA_run02.txt"`)**.

```bash
## Extracting truth haplotype set for Sample "NA12891"
$ head -n1 SetA/truth_RBphasedHaplotype_SetA.txt 
CHROM	POS	REF	all-alleles	NA06989:PI	NA06989:PG_al	NA10850:PI	NA10850:PG_al	NA06984:PI	NA06984:PG_al	NA07056:PI	NA07056:PG_al	NA12045:PI	NA12045:PG_al	NA11843:PI	NA11843:PG_al	NA12890:PI	NA12890:PG_al	NA12889:PI	NA12889:PG_al	NA12892:PI	NA12892:PG_al	NA12891:PI	NA12891:PG_al

# the index position of the haplotype for sample "NA12891" is 23 and 24.
# we extract the this column but also remove the rows that are empty i.e "."
# we also include the data from "CHROM" and "POS". This is used so we can make sure that the genotypes are coming from the same genomic position
$ awk 'BEGIN{FS=OFS="\t"} {if ($24 !=".") print $1, $2, $23, $24}' SetA/truth_RBphasedHaplotype_SetA.txt > SetA/truth_Haplotype_NA12891.txt

## Now, extract the phased haplotype for Sample "NA12891"
$ head -n1 SetA_run02/phaseExtendedHaplotype_SetA_run02.txt 
CHROM	POS	REF	all-alleles	all-freq	NA12891:PI	NA12891:PG_al	NA12892:PI	NA12892:PG_al	NA06989:PI	NA06989:PG_al	NA10850:PI	NA10850:PG_al	NA06984:PI	NA06984:PG_al	NA07056:PI	NA07056:PG_al	NA12045:PI	NA12045:PG_al	NA11843:PI	NA11843:PG_al	NA12890:PI	NA12890:PG_al	NA12889:PI	NA12889:PG_al

# here the index position of the haplotype for sample "NA12891" is 6 and 7
$ awk 'BEGIN{FS=OFS="\t"} {if ($7 !=".") print $1, $2, $6, $7}' SetA02/phaseExtendedHaplotype_SetA02.txt > SetA02/phased_Haplotype_NA12891.txt
```

#### Now, we import the truth and phased data into `R` to plot and quantify switch-errors




<br>
## 03 - Set (B) : for HAPLOTYPE file with 25 samples.
The process for running haplotype phaseExtension on samples from SetB is typically the same as above. 

I have included a `BASH SHELL` script file "PhaseExtenderOnForLoop_SetB.sh"
```bash
# simply run
$ ./PhaseExtenderOnForLoop_SetB.sh
```
**which will:**
  - run phase extension on all 25 samples inside **SetA**.
  - merge the haplotype files into a single HAPLOTYPE file "merged_haplotype.txt" inside the directory **SetB02**.

```BASH
# Make a copy and rename the above output file "merged_haplotype.txt" to "phaseExtendedHaplotype_SetB_run02.txt"
$ cp SetB_run02/merged_haplotype.txt SetB_run02/phaseExtendedHaplotype_SetB_run02.txt
```

## Test switch errors

### Step 04: Second recursive run of haplotype phasing 
