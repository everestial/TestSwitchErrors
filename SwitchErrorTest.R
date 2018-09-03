
### Set the required path; **update the path as need be with your directory 
getwd()
setwd("/home/priyanka/Dropbox/SharedWithHarish/TestSwitchErrors/HapMap3_r2_b36_2009/Contd_PhaseQC_PhaseExtender")
getwd()
list.files()  # read available files and folders 

#####  Read the required data  ######

## Import "truth haplotype" for SetA
truthHaplotype_NA12891 <- read.table('SetA/truth_Haplotype_NA12891.txt', header = TRUE)
truthHaplotype_NA12891
  #**Note: R doesn't allow to use ":" in header name and it is automatically renamed to "."
  # you can see this in the output 

# change the name of the header
colnames(truthHaplotype_NA12891)[colnames(truthHaplotype_NA12891)=="NA12891.PI"] <- "true.NA12891.PI"
colnames(truthHaplotype_NA12891)[colnames(truthHaplotype_NA12891)=="NA12891.PG_al"] <- "true.NA12891.PG_al"


## Import "phased haplotype" for SetA
phased_SetA_NA12891 <- read.table('SetA_run02/phased_Haplotype_NA12891.txt', header = TRUE)
phased_SetA_NA12891

# change the name of the header
colnames(phased_SetA_NA12891)[colnames(phased_SetA_NA12891)=="NA12891.PI"] <- "phased.NA12891.PI"
colnames(phased_SetA_NA12891)[colnames(phased_SetA_NA12891)=="NA12891.PG_al"] <- "phased.NA12891.PG_al"


## Merge the truth and phased data set to identify switch errors
merged.data <- merge(truthHaplotype_NA12891, phased_SetA_NA12891,
                     by=c("CHROM", "POS"))
head(merged.data)

# set the data in order by "POS" 
merged.data <- merged.data[order(merged.data$POS),]

# after merging the row-index becomes random; so let's put it in order
rownames(merged.data) <- NULL

# Now, compare the "truth haplotype" with "phased haplotype" and ...
# ... find the sites where haplotype phasing switch 
merged.data$match <- ifelse((merged.data$true.NA12891.PG_al == merged.data$phased.NA12891.PG_al), 0, 1)
# **Note: the "match" column should have values "0's" and "1's". A continous "0" or "1" indicates 
# .. properly phased block. But, transition from 0 -> 1 or from 1 -> 0 suggests a "Switch Error" betweeen the blocks.

# Now, plot the "Switch Error" points as png
png("out_match_SetA_phaseExtender.png", width = 1600, height = 600)
plot(merged.data$POS, merged.data$match, main = "Switch points over the genomic coordinates.", type = "s", 
     xlab = "genomic position", ylab = "switch errors")
#title(main = "Switch points", xlab = 'genomic coordinates', ylab = 'switch points')
#abline(v=haplotype_breaks, col='red')
dev.off()


## The switch point data (i.e "match") and the plot we obtained above doesn't account for break in haplotypes.
# Since out haplotype aren't completely merged genome/chromosome wide. The switch points needs to be ..
# .. defined per block. 

## Compute switch points also including changes in "PI" values.
# Now, the "switch points" include when :
  # we see 0 -> 1, or when 1 -> 0 in the "match" data
  # PI of the haplotype block changes

# create empty variables to store/update data on "for-loop" 
previous_match <- 0
total_possible_switch <- 0
freq_of_switch <- 0
previous_pi <- 1  # **if PI of first block is something else change this

# the length (number of rows) of the dataframe
seq_len <- length(merged.data$match)
seq_len

# run for loop on the whole sequence length to compute frequecy of switch points
for (ith in c(1:seq_len)){
  current_match = merged.data$match[ith]
  current_pi = merged.data$phased.NA12891.PI[ith]
  total_possible_switch = total_possible_switch + 1
  
  if (previous_match != current_match | previous_pi != current_pi){
    freq_of_switch = freq_of_switch + 1}
  
  previous_match = current_match
  previous_pi = current_pi
  }

freq_of_switch
total_possible_switch

# calculate switch error rate
switch_error_rate = (freq_of_switch/total_possible_switch)
switch_error_rate


## After the haplotype are phased by phaseExtender, it is important to know how the true haplotype 
## ... is split into smaller true haplotype blocks; and what are it's sizes.
## So, now we calculate the size distribution of the properly phased haplotype.
# We need to account for the match and also for change in haplotype block index.
previous_match <- 0
hap_size <- 0
haplotype_sizes <- integer()  # to store the size of haplotype extended truly.
previous_pi <- 1 
haplotype_breaks <- integer()

# now, update the list using for loop
for (line in merged.data$match){
  match_value = line
  #print(allele)
  if (previous_match == match_value){
    hap_size = hap_size + 1}
  else if (previous_match != match_value){
    haplotype_sizes <- c(haplotype_sizes, hap_size) # add that hap_size to the list
    #print(previous_allele)
    #print(hap_size)
    #print(haplotype_size)
    hap_size <- 1}  # reset the haplotype size
  previous_match <- match_value
}

#####
for (ith in c(1:seq_len)){
  current_match = merged.data$match[ith]
  current_pi = merged.data$phased.NA12891.PI[ith]
  if (previous_match == current_match & previous_pi == current_pi) {
    hap_size = hap_size + 1}
  else if (previous_match != current_match | previous_pi != current_pi) {
    haplotype_sizes <- c(haplotype_sizes, hap_size)
    hap_size <- 1}  # reset the haplotype size to 1
  
  # store the last haplotype size when for loop ends
  if (ith == seq_len){
    haplotype_sizes <- c(haplotype_sizes, hap_size)}
  
  # find the position where haplotype breaks occur 
  if (previous_pi != current_pi){
    haplotype_breaks <- c(haplotype_breaks, merged.data$POS[ith])}
  
  # update the previous match and pi values for next for-loop 
  previous_match = current_match
  previous_pi = current_pi
}

## We can now overlay out haplotype breaks position on the top of previous switch points plot
plot(merged.data$POS, merged.data$match, main = "Switch points over the genomic coordinates with haplotype breaks.", type = "s", 
     xlab = "genomic position", ylab = "switch errors")
abline(v=haplotype_breaks, col='red')


## Now, make switch points plot after accounting for haplotype breaks. 
# for that we will create another column with updated matches between truth and phased haplotypes
#merged.data$match02 <- 

####3
begin_match = 0
begin_pi = 1

for (ith in c(1:seq_len)) {
  current_pi = merged.data$phased.NA12891.PI[ith]
  truth_gt = merged.data$true.NA12891.PG_al[ith]
  phased_gt = merged.data$phased.NA12891.PG_al[ith]
  
  
}


####



previous_match <- 0
previous_pi <- 1
match02 <- integer()
continuing_match <- 0



for (ith in c(1:seq_len)){
  current_match = merged.data$match[ith]
  current_pi = merged.data$phased.NA12891.PI[ith]
  if (current_match == previous_match & current_pi == previous_pi){
    match02 <- c(match02, continuing_match)
    
    # store data for next for loop
    previous_match <- current_match}
  
  else if (current_pi != previous_pi & current_match == previous_match){
    if (current_match == 0) {continuing_match = 1}
    else if (current_match == 1) {continuing_match = 0}
    match02 <- c(match02, continuing_match)
    previous_match <- current_match}
  
  # covers both conditions :
    # current_match != previous_match & current_pi == previous_pi 
    # current_match != previous_match & current_pi != previous_pi
  else {
    continuing_match <- current_match
    match02 <- c(match02, continuing_match)
    previous_match <- current_match}
  
  # update the values
  previous_pi <- current_pi}

match02
length(match)
length(match02)


## Just printing 
print(haplotype_sizes)
typeof(haplotype_sizes)
length(haplotype_sizes)
sum(haplotype_sizes)
haplotype_breaks
length(haplotype_breaks)
length(merged.data$match)

## add data "match02" to the dataframe
merged.data$match02 <- match02


## now, compute haplotype sizes








# for convenience let's convert this haplotype size list into integer array
haplotype_size_numeric <- as.numeric(unlist(haplotype_size))
haplotype_size_numeric

## Let's plot a histogram of the haplotype size distribution 
hisHap <- hist(haplotype_size_numeric)  # general histogram metrics for the data 
hist(haplotype_size_numeric, 
     main="Histogram for Haplotype size distribution", 
     xlab="haplotype size", 
     border="blue", 
     col="green",
     xlim=c(0,500),
     ylim = c(0,320),
     las=1, 
     breaks=20)
# add the size of the frequency to the histogram plot 
text(hisHap$mids,hisHap$counts,labels=hisHap$counts, adj=c(0.5, -0.5))

## Density plot
densHap <- density(haplotype_size_numeric)
plot(densHap, main = "Density plot of the haplotype size distribution")
polygon(densHap, col = 'red', border = 'blue')


### Still on test *****
## Let's use the data prepared using pandas to do some statistics 



########  ***********************
## To Do: Similarly make histogram and density plot using ShapeIT-PIR, phaseExtender.
## Also, overlap the histogram plot, density plot from all three methods. 
# possible link - https://www.statmethods.net/graphs/density.html 
########  ***********************


# Plot haplotype size by number of heterozygous variants
# length of the haplotype depends upon the genomic distance between the variants ?? is this meaningful 
match


###########  ******* this part below is just a repeat using another file *************
## Plot Sample size and switch error accuracy (using data obtained from pandas script) 
true_vs_output_haps <- read.table("true_vs_output_haps.txt", header = TRUE, sep = '\t')
true_vs_output_haps

# lets plot the switch points again
plot(true_vs_output_haps$position, true_vs_output_haps$switchPoints, 
     main = "Switch points", type = "s")

# histogram of the distribution 
haplotype_stats <- read.table("switchPoints_stats.txt", header = TRUE, sep = '\t')

hisHap02 <- hist(haplotype_stats$hapSizeByNumVars)  # general histogram metrics for the data 
hist(haplotype_stats$hapSizeByNumVars, 
     main="Histogram for Haplotype size distribution", 
     xlab="haplotype size", 
     border="blue", 
     col="green",
     xlim=c(0,500),
     ylim = c(0,320),
     las=1, 
     breaks=20)
# add the size of the frequency to the histogram plot 
text(hisHap02$mids,hisHap02$counts,labels=hisHap02$counts, adj=c(0.5, -0.5))


# plot genomicSize (by Vars vs. by genomic range)
plot(haplotype_stats$hapSizeByNumVars, haplotype_stats$hapSizeByGenomicRange)

###############   *************************************



## Compute switherror rates for all three methods
# for phaseExtender SwErr can be computed in two ways:
  # number of SwErr/number of Hap block
  # number of SwErr/number of HetGenotypes

## may be we can also add PRISM to the comparison 
# phaseExtender can have 3 times recursive improvements 


### just a test to check if the generated data is producing a required distribution 
my_data <- list(7, 11, 11, 14, 10, 10, 11, 8, 14, 13, 8, 12, 17, 11, 8, 13, 14, 10, 10, 11, 14, 7, 10, 10, 15, 5, 11, 14, 10, 12, 16, 7, 9, 7, 8, 8, 10, 10, 15, 11, 11, 5, 11, 12, 9, 12, 14, 6, 9, 11, 14, 11, 12, 8, 7, 6, 13, 9, 6, 11, 9, 10, 10, 13, 8, 10, 12, 17, 13, 9, 5, 18, 6, 6, 12, 11, 8, 14, 17, 14, 7, 11, 10, 9, 7, 6, 9, 10, 11, 6, 14, 7, 11, 10, 11, 12, 8, 10, 5, 13, 9, 8, 14, 7, 8, 9, 9, 7, 5, 11, 11, 15, 6, 12, 14, 11, 10, 12, 13, 7, 11, 12, 12, 11, 12, 9, 6, 13, 11, 12, 10, 9, 8, 13, 10, 10, 10, 8, 12, 12, 13, 11, 10, 14, 8, 11, 6, 9, 13, 12, 10, 6, 10, 16, 12, 10, 9, 6, 10, 6, 14, 8, 11, 9, 13, 11, 12, 12, 16, 9, 12, 9, 8, 8, 10, 11, 11, 8, 10, 12, 13, 8, 16, 10, 5, 8, 8, 9, 11, 7, 11, 8, 10, 11, 8, 8, 13, 13, 10, 11, 10, 9, 9, 9, 11, 9, 9, 15, 13, 10, 9, 11, 9, 12, 14, 15, 12, 10, 7, 13, 8, 8, 12, 10, 12, 10, 13, 10, 14, 10, 14, 12, 9, 11, 12, 20, 9, 6, 12, 14, 13, 10, 8, 7, 8, 8, 8, 7, 10, 9, 13, 9, 13, 9, 7, 9, 8, 9, 10, 9, 12, 11, 14, 18, 13, 15, 9, 18, 11, 6, 13, 10, 12, 9, 8, 18, 11, 16, 9, 15, 13, 11, 7, 11, 8, 9, 13, 10, 9, 8, 7, 10, 10, 12, 15, 12, 10, 18, 7, 6, 10, 13, 8, 10, 12, 11, 9, 7, 12, 12, 7, 11, 10, 4, 10, 15, 14, 7, 8, 9, 12, 14, 9, 8, 16, 12, 7, 11, 12, 11, 12, 8, 9, 12, 10, 7, 6, 14, 7, 12, 12, 11, 8, 13, 12, 14, 8, 10, 16, 11, 8, 7, 9, 17, 13, 6, 12, 9, 12, 7, 7, 10, 11, 9, 10, 16, 13, 8, 14, 12, 7, 11, 9, 10, 9, 17, 9, 11, 12, 11, 12, 8, 9, 8, 14, 8, 8, 5, 12, 9, 8, 7, 10, 12, 8, 15, 12, 9, 17, 14, 9, 13, 12, 11, 12, 9, 14, 10, 8, 7, 13, 8, 12, 12, 10, 14, 7, 9, 9, 9, 12, 12, 12, 7, 12, 7, 7, 8, 9, 12, 8, 7, 8, 8, 11, 14, 14, 9, 13, 15, 8, 8, 12, 9, 10, 9, 7, 12, 7, 18, 12, 4, 9, 9, 10, 8, 10, 7, 12, 7, 6, 7, 16, 8, 6, 9, 12, 9, 10, 10, 6, 13, 11, 13, 10, 13, 19, 14, 9, 16, 12, 4, 11, 12, 8, 12, 10, 10, 8, 12, 12, 7, 13, 11, 4, 9, 14, 5, 11, 11, 9, 7, 4, 8, 14, 9, 6, 7, 10, 11, 7, 14, 10, 13, 14, 12, 15, 6, 6, 10, 8, 14, 14, 9, 6, 14, 10, 9, 7, 16, 7, 8, 6, 12, 6, 12, 9, 9, 13, 7, 13, 13, 9, 10, 9, 14, 12, 14, 7, 10, 6, 6, 10, 6, 13, 10, 11, 13, 14, 10, 13, 9, 13, 13, 8, 11, 8, 8, 9, 16, 9, 8, 14, 13, 12, 10, 13, 12, 5, 11, 10, 11, 13, 11, 9, 8, 12, 9, 15, 13, 12, 13, 7, 12, 7, 9, 10, 10, 11, 11, 14, 8, 13, 5, 12, 8, 8, 10, 9, 13, 14, 14, 11, 14, 11, 13, 14, 8, 9, 12, 9, 12, 9, 11, 7, 13, 10, 6, 10, 13, 14, 7, 17, 9, 11, 7, 7, 12, 6, 10, 13, 12, 15, 8, 10, 13, 19, 12, 9, 10, 10, 7, 14, 12, 9, 11, 19, 6, 10, 12, 10, 6, 11, 12, 11, 14, 10, 13, 8, 11, 15, 10, 7, 15, 4, 10, 8, 12, 13, 10, 12, 10, 6, 7, 15, 10, 14, 16, 9, 14, 10, 10, 9, 12, 8, 13, 11, 12, 8, 11, 10, 14, 6, 9, 6, 9, 12, 14, 10, 13, 12, 7, 12, 11, 9, 10, 10, 13, 9, 13, 9, 11, 8, 13, 10, 14, 9, 8, 10, 7, 4, 16, 10, 10, 6, 18, 9, 9, 9, 13, 7, 11, 8, 6, 9, 15, 7, 8, 8, 11, 8, 12, 6, 10, 11, 8, 7, 10, 11, 8, 13, 12, 14, 8, 15, 12, 11, 16, 9, 5, 5, 7, 11, 7, 10, 15, 12, 10, 15, 9, 9, 8, 8, 14, 8, 14, 11, 11, 10, 6, 11, 9, 10, 10, 9, 6, 7, 11, 13, 9, 10, 11, 9, 13, 8, 18, 12, 13, 13, 10, 9, 8, 14, 12, 11, 11, 11, 10, 17, 6, 10, 17, 8, 11, 7, 11, 8, 11, 13, 8, 12, 15, 12, 12, 13, 6, 11, 15, 9, 8, 9, 12, 4, 13, 13, 12, 9, 11, 17, 11, 8, 10, 11, 12, 11, 11, 8, 9, 11, 11, 8, 6, 10, 15, 7, 14, 9, 9, 7, 13, 6, 7, 9, 9, 9, 11, 11, 14, 9, 17, 13, 12, 6, 8, 12, 6, 6, 7, 12, 8, 15, 10, 8, 19, 9, 14, 12, 9, 12, 6)
my_data
length(my_data)
typeof(my_data)

my_data <- as.numeric(my_data)

my_data_hist <- hist(my_data)
my_data_hist
typeof(my_data)






