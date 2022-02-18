
## Purpose of the this script: 
  # Compute the switch points using the truth and second round of phased haplotype (SetA_02)

### Set the required path; 
## **update the path as need be with your directory 
getwd()
list.files()  # read available files and folders 


##########  Switch errors test - Set A_02 (second round of phaseExtension) ########## 
#####  Read the required data  ###### 

## Import "truth haplotype" for SetA 
# **Note: The truth set stays the same. 
truthHaplotype_NA12891 <- read.table('data/SetA/truth_Haplotype_NA12891.txt', header = TRUE)
head(truthHaplotype_NA12891)
  #**Note: R doesn't allow to use ":" in header name and it is automatically renamed to "."
  # you can see this in the output 

# change the name of the header
colnames(truthHaplotype_NA12891)[colnames(truthHaplotype_NA12891)=="NA12891.PI"] <- "true.NA12891.PI"
colnames(truthHaplotype_NA12891)[colnames(truthHaplotype_NA12891)=="NA12891.PG_al"] <- "true.NA12891.PG_al"


## Import "phased haplotype" for SetA
phased_SetA_NA12891 <- read.table('data/SetA_03/phased_Haplotype_NA12891.txt', header = TRUE)
# or, we can use the file "extended_haplotype_NA12891.txt" from the folder "/SetA_02/phased_NA12891_SetA_run02/"
head(phased_SetA_NA12891)

# change the name of the header
colnames(phased_SetA_NA12891)[colnames(phased_SetA_NA12891)=="NA12891.PI"] <- "phased.NA12891.PI"
colnames(phased_SetA_NA12891)[colnames(phased_SetA_NA12891)=="NA12891.PG_al"] <- "phased.NA12891.PG_al"


## Merge the truth and phased data set to identify switch errors
merged.data <- merge(truthHaplotype_NA12891, phased_SetA_NA12891,
                     by=c("CHROM", "POS"))
head(merged.data)

# set the data in order by "POS" - ** it's very important to do this **
merged.data <- merged.data[order(merged.data$POS),]

# after merging, the row-index becomes random; so let's put it in order
rownames(merged.data) <- NULL

# Now, compare the "truth haplotype" with "phased haplotype" and ...
# ... find the sites where haplotype phasing switched
merged.data$match <- ifelse((merged.data$true.NA12891.PG_al == merged.data$phased.NA12891.PG_al), 0, 1)
# **Note: the "match" column should have values "0's" and "1's". A continous "0" or "1" indicates 
# .. properly phased block. But, transition from 0 -> 1 or from 1 -> 0 suggests a "Switch Error" betweeen the blocks.

# Now, plot the "Switch Error" points as png
png("data/rplots/SwitchPoints_SetA02_withOutHaplotypeBreaks.png", width = 1600, height = 600)
plot(merged.data$POS, merged.data$match, main = "Switch points over the genomic coordinates 
     (before accounting for haplotype breaks).", type = "s", 
     xlab = "genomic position", ylab = "switch errors")
#title(main = "Switch points", xlab = 'genomic coordinates', ylab = 'switch points')
dev.off()


# Compute switch error rates 
number_of_hets <- length(merged.data$match)
num_of_switch <- 0
prev_match <- merged.data$match[1]
for (item in merged.data$match){
  curr_match = item
  if (curr_match != prev_match){
    num_of_switch = num_of_switch + 1
    
    # reset the "previous match" value
    prev_match = curr_match}}

switch_err_rate = num_of_switch/(number_of_hets)
switch_err_rate
# [1] 0.01665632


######## Fix the switch points by taking "haplotype breaks" into account. ##########
## The above switch point data (i.e "match") and the plot we obtained above doesn't account for break in haplotypes.
# Since our haplotype aren't completely merged genome/chromosome wide. The switch points needs to be defined per block. 

## So, now we compute switch point by including changes in "PI" values.
# Now, the "switch points" are addressed when :
  # we see 0 -> 1, or when 1 -> 0 in the "match" data
  # PI of the haplotype block changes
    # ** So, now we need to account for the match and also for change in haplotype block index.

# the length (number of rows) of the dataframe
seq_len <- length(merged.data$match)
seq_len

# set empty variables 
previous_match <- 0
hap_size <- 0
haplotype_sizes <- integer()  # to store the size of haplotype extended truly.
previous_pi <- 1 
haplotype_breaks <- integer()

# run a for loop to find switch points 
for (ith in c(1:seq_len)){
  current_match = merged.data$match[ith]
  current_pi = merged.data$phased.NA12891.PI[ith]
  if (previous_match == current_match & previous_pi == current_pi) {
    hap_size = hap_size + 1}
  else if (previous_match != current_match | previous_pi != current_pi) {
    # store the haplotype size of properly phased (extended) block
    haplotype_sizes <- c(haplotype_sizes, hap_size)
    hap_size <- 1}  # reset the haplotype size to 1
  
  # store the last haplotype size when the for loop reaches the end
  if (ith == seq_len){
    haplotype_sizes <- c(haplotype_sizes, hap_size)}
  
  # also find the genomic positions where haplotype breaks occur 
  if (previous_pi != current_pi){
    haplotype_breaks <- c(haplotype_breaks, merged.data$POS[ith])}
  
  # update the previous match and pi values for next for-loop 
  previous_match = current_match
  previous_pi = current_pi }


## Identify switch points and compute "switch error" rates. 
# the total number of haplotype blocks represent frequecy of switch points
freq_of_switch = length(haplotype_sizes)
total_possible_switch = sum(haplotype_sizes) # or number of hets site
haplotype_sizes

# calculate switch error rate
switch_error_rate = (freq_of_switch/total_possible_switch)
switch_error_rate
# [1] 0.01748396

## Now, add another column that represent the fully accounted haplotype switchpoints.
phase_state = 0
match_by_pi <- integer()
for (sizes in haplotype_sizes) {
  match_by_pi <- c(match_by_pi, rep(phase_state, sizes))
  
  # change the phase state for next loop
  if (phase_state == 0){
    phase_state = 1}
  else if (phase_state == 1){
    phase_state = 0} }

## add data "match02" to the dataframe
merged.data$match_by_pi <- match_by_pi

## Now, make switch points plot by accounting for the haplotype breaks. 
# for that we will create another column with updated matches between truth and phased haplotypes
png("data/rplots/SwitchPoints_SetA02_withHaplotypeBreaks.png", width = 1600, height = 600)
plot(merged.data$POS, merged.data$match_by_pi, main = "Switch points over the genomic coordinates with haplotype breaks.", type = "s", 
     xlab = "genomic position", ylab = "switch errors")

## We can now overlay the haplotype breaks position on the top of switch points plot
abline(v=haplotype_breaks, col='red')
dev.off()

# for convenience let's convert this haplotype size list into integer array
haplotype_size_numeric <- as.numeric(unlist(haplotype_sizes))
haplotype_size_numeric

## Let's plot a histogram of the haplotype size distribution 
hisHap <- hist(haplotype_size_numeric)  # general histogram metrics for the data 
png("data/rplots/Histogram_ofHaplotypeSizeDistribution_SetA02.png", width = 1600, height = 600)
hist(haplotype_size_numeric, 
     main="Histogram for Haplotype size distribution", 
     xlab="haplotype size", 
     border="blue", 
     col="green",
     xlim=c(0,200),
     ylim = c(0,320),
     las=1, 
     breaks=20)
# add the size of the frequency to the histogram plot 
text(hisHap$mids,hisHap$counts,labels=hisHap$counts, adj=c(0.5, -0.5))
dev.off()

## Density plot
densHap <- density(haplotype_size_numeric)
png("data/rplots/DensityPlot_ofHaplotypeSizeDistribution_SetA02.png", width = 1600, height = 600)
plot(densHap, main = "Density plot of the haplotype size distribution")
polygon(densHap, col = 'red', border = 'blue')
dev.off()

print("Completed the switch error analyses on first round of phaseExtension on Set-A02 data.")
#### Complete analyses on Set-A02  #######

