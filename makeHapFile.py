import os
import sys
import pandas as pd
import numpy as np
from collections import OrderedDict
from scipy import stats
import argparse


"""
A python parser to convert the reference genome panel (haplotype) data to HAPLOTYPE format required by 
phaseExtender. This data is used as simulated data to test the phasing quality produced by phaseExtender.
"""

# print(os.getcwd())

# set the working directory
# os.chdir("/home/priyanka/Dropbox/SharedWithHarish/TestSwitchErrors/HapMap3_r2_b36_2009/Contd_PhaseQC_PhaseExtender/")
# print(os.getcwd())


def main():
    parser = argparse.ArgumentParser()
    print()
    print(
        "Purpose of the application : This application simulates ReadBackPhased haplotype given *.haps, *.legend and *.sample file are provided."
    )

    parser.add_argument(
        "-haps",
        help='the "*.haps" file that contains haplotype data in numeric format',
        required=True,
    )
    parser.add_argument(
        "-legend",
        help='the "*.legend" file that contains genomic coordinates and genotypes as IUPAC bases.',
        required=True,
    )
    parser.add_argument(
        "-sample_file",
        help='the "*.sample" file that contains sample names and have their position (rows) \
        matched with column position in "*.haps" file.',
        required=True,
    )
    parser.add_argument(
        "-chr",
        help='The chromosome name the "*.legend" and "*.haps" file belongs to.',
        required=True,
    )
    parser.add_argument(
        "-samples",
        help='The samples from "*.sample" file that are to be converted to RBphased data.',
        required=True,
    )
    parser.add_argument(
        "-output_truth", help="Write the truth HAPLOTYPE file.", required=True
    )
    parser.add_argument(
        "-output_sim", help="Write the simulated HAPLOTYPE file.", required=True
    )

    # activate "args"
    global args
    args = parser.parse_args()

    chr_ = args.chr
    sample_f = args.sample_file
    haps_f = args.haps
    samples = args.samples
    legend_f = args.legend
    out_tr = args.output_truth
    out_sim = args.output_sim

    ### pick samples (and it's position in haplotype file) to convert to haplotype file compatible with phaseExtender
    sample_of_interest = []
    if "," in args.samples:
        sample_of_interest = args.samples.split(",")
    else:
        sample_of_interest += args.sample

    # sample_of_interest = ['NA12891', 'NA12892']

    """
    # Just some existing sample names 
    sample_of_interest = ['NA12891', 'NA12892', 'NA06989', 'NA10850', 'NA06984', 'NA07056', 
                            'NA12045', 'NA11843', 'NA12890', 'NA12889', 'NA12347', 'NA12348', 
                            'NA12273', 'NA12827', 'NA12342', 'NA12340', 'NA11894', 'NA12763', 
                            'NA11881', 'NA12005']
                            """

    # from the "*.sample" file extract the hapltoype position for the sample of interest
    # store the "sample name" and it's index position as a dictionary.
    # with open('hapmap3_r2_b36_all.sample') as sample_file:
    with open(sample_f) as sample_file:
        sample_pos_dict = OrderedDict()  # to store sample, position values
        for nth, lines in enumerate(sample_file):
            # print (nth, lines)
            sample_id = lines.rstrip("\n").split(" ")[0]

            if sample_id in sample_of_interest:
                sample_pos_dict[sample_id] = [nth * 2 - 1, nth * 2]

                # reduce the size of the list after the sample's index is found
                sample_of_interest = [x for x in sample_of_interest if x != sample_id]

        # if the sample isn't found in sample file raise the error
        if len(sample_of_interest) != 0:
            print(
                'The following samples are not found in the provided "*.sample" file.'
            )
            print(sample_of_interest)

            print("Please remove those sample ....  ")
            print("exiting....  ")
            print()
            sys.exit()

    ## load the "*.legend" file
    # haps_legend = pd.read_csv('hapmap3_r2_b36_chr20.legend', sep = ' ')
    haps_legend = pd.read_csv(legend_f, sep=" ")

    # extract ref and alt allele from "*.legend" file.
    ref_alt_data = haps_legend.drop(["rsID"], axis=1)
    # ref_alt_data.insert(0, column='CHROM', value=20)
    ref_alt_data.insert(0, column="CHROM", value=chr_)

    ref_alt_data["a1"] = ref_alt_data["a0"] + "," + ref_alt_data["a1"]
    ref_alt_data.rename(
        index=str,
        columns={"position": "POS", "a0": "REF", "a1": "all-alleles"},
        inplace=True,
    )

    ### Now, take the "*.legend" and "*.haps" file and prepare HAPLOTYPE file for each "sample of interest".
    # haps_file = pd.read_csv('hapmap3_r2_b36_chr20.haps', sep = ' ', header=None)
    haps_file = pd.read_csv(haps_f, sep=" ", header=None)

    print("\nExtracting haplotype data for each sample of interest ..... ")
    # Extract the haplotype from the "sample of interest"
    # we will store the haplotype data as keys (sample id) and values (pandas df)
    sample_hap_dict = OrderedDict()
    nth = 1
    for sample_name, sample_idx in sample_pos_dict.items():
        print("Sample %s : " % nth, sample_name)
        nth += 1

        haplotype = pd.DataFrame()  # empty dataframe

        # extract the haplotype data
        haplotype["hapA"] = haps_file[sample_idx[0] - 1]
        haplotype["hapB"] = haps_file[sample_idx[1] - 1]

        # join the haplotype data with legend file
        # index doesn't matter here, because we are doing side by side join.
        # but, the two file should have exact same length.
        haplotype = haps_legend.join(haplotype)

        # rename column "position" to "POS"
        haplotype.rename(index=str, columns={"position": "POS"}, inplace=True)

        # add chromosome name/number
        # haplotype.insert(0, column='CHROM', value=20)
        haplotype.insert(0, column="CHROM", value=chr_)

        # drop the homozygous genotype position and then reset the index
        haplotype = haplotype[haplotype["hapA"] != haplotype["hapB"]]
        haplotype.reset_index(drop=True, inplace=True)

        ## convert numeric haplotype into IUPAC base
        haplotype[sample_name + "_hapA"] = haplotype.apply(
            lambda row: row["a0"] if row["hapA"] == 0 else row["a1"], axis=1
        )
        haplotype[sample_name + "_hapB"] = haplotype.apply(
            lambda row: row["a0"] if row["hapB"] == 0 else row["a1"], axis=1
        )

        # add the converted haplotype to column "sample_name:PG_al"
        haplotype[sample_name + ":PG_al"] = haplotype.apply(
            lambda row: row[sample_name + "_hapA"] + "|" + row[sample_name + "_hapB"],
            axis=1,
        )

        # drop the not required columns
        haplotype.drop(
            [
                "rsID",
                "a0",
                "a1",
                "hapA",
                "hapB",
                sample_name + "_hapA",
                sample_name + "_hapB",
            ],
            inplace=True,
            axis=1,
        )

        ##### ******************** #####
        ### Starting simulation to create a RBphased haplotype blocks (number of RBphased blocks and size of each block)

        ## randomly create ReadBackPhase haplotype blocks
        # we will use the function "make_poissonish()" to get a list of RBphased blocks
        hap_length = len(
            haplotype
        )  # total number of "heterozyous varaints" in the haplotype
        avg_rbphased_size = np.random.choice(
            [5, 6, 7]
        )  # randomly choosing avg_rbphased_size
        number_of_rbphased_blocks = int(round(hap_length / avg_rbphased_size))

        ## preparing random unique PI index to add to the dataframe
        # this part may be moved to a separate function in future if desired.
        # rbphased_blocks = make_poissonish(9366, 1500)  #if a fixed "hap length" and "number of rbphased_blocks" is desired
        rbphased_blocks = make_poissonish(hap_length, number_of_rbphased_blocks)

        """
        ## If you want some to print the statistics of the simulated RBphased block sizes and distribution 
        print('\nLength of the haplotype')
        print(len(haplotype))
        print(hap_length, avg_rbphased_size, number_of_rbphased_blocks)
        
        print('\nrbphased_blocks')
        print(rbphased_blocks)
        print('stats', np.mean(rbphased_blocks), np.median(rbphased_blocks), 
                stats.mode(rbphased_blocks), np.var(rbphased_blocks))
        print(np.max(rbphased_blocks), np.min(rbphased_blocks))
        """

        # now, split the whole phased data into RBphased chunks,
        # create unique index for each RBphased chunk
        unique_rbphased_index = [x + 1 for x in range(number_of_rbphased_blocks)]

        # expand each unique index by the number of times (of rbphased_blocks) ..
        # *Remember that the unique index and it's times are index matched.
        pi_index_values = [
            [x] * y for x, y in zip(unique_rbphased_index, rbphased_blocks)
        ]
        pi_index_values = [item for items in pi_index_values for item in items]

        # add the index values to dataframe as "sample name:PI"
        # haplotype[sample_name + ':PI'] = pd.Series(pi_index_values).values
        haplotype.insert(
            loc=2, column=sample_name + ":PI", value=pd.Series(pi_index_values).values
        )

        ## Store the "truth haplotype" set to file
        truth_haplotype = haplotype[
            ["CHROM", "POS", sample_name + ":PI", sample_name + ":PG_al"]
        ]
        ## If you want to write the truth set data to a file before flipping the RBphase haplotype data
        # haplotype.to_csv(sample_name + '_TruthSetHaplotype.txt', sep='\t', index = False,
        # columns=['CHROM', 'POS', sample_name + ':PI', sample_name + ':PG_al'])

        ## Flip the position of RBphased blocks index is an odd number
        haplotype[sample_name + ":PG_al"] = haplotype.apply(
            lambda row: "|".join(row[sample_name + ":PG_al"].split("|")[::-1])
            if row[sample_name + ":PI"] % 2 == 0
            else row[sample_name + ":PG_al"],
            axis=1,
        )

        ## Store the "simulated haplotype set" to a file.
        simulated_haplotype = haplotype[
            ["CHROM", "POS", sample_name + ":PI", sample_name + ":PG_al"]
        ]
        # * If you want to write the simulated RBphased haplotype data to a file
        # haplotype.to_csv(sample_name + '_SimulatedRBphasedHaplotype.txt', sep='\t', index = False,
        # columns=['CHROM', 'POS', sample_name + ':PI', sample_name + ':PG_al'])

        ## add the simulated data (pandas dataframe) to an ordered dictionary
        # sample_hap_dict[sample_name] = haplotype
        sample_hap_dict[sample_name] = [truth_haplotype, simulated_haplotype]

    ## Now, merge the dataframe (haplotype from several samples) into a single HAPLOTYPE file
    truth_haplotype_file = pd.DataFrame()
    simulated_haplotype_file = pd.DataFrame()
    for ks, vs in sample_hap_dict.items():
        if len(truth_haplotype_file) == 0:
            truth_haplotype_file = vs[0]
            simulated_haplotype_file = vs[1]
        else:
            truth_haplotype_file = truth_haplotype_file.merge(
                vs[0], on=["CHROM", "POS"], how="outer"
            ).fillna(".")
            simulated_haplotype_file = simulated_haplotype_file.merge(
                vs[1], on=["CHROM", "POS"], how="outer"
            ).fillna(".")

    ### Finally, merge the haplotype for each sample with "ref-alt" data to create complete HAPLOTYPE file
    final_haplotype_file = ref_alt_data.merge(
        truth_haplotype_file, on=["CHROM", "POS"], how="inner"
    )

    # Write the "HAPLOTYPE DATA" to a file.
    print('\nWriting truth set haplotype data to a file "%s" .... ' % (out_tr))
    final_haplotype_file.to_csv(out_tr, sep="\t", index=False)

    ## write the final simulated haplotype for all the samples
    final_simulated_haplotype_file = ref_alt_data.merge(
        simulated_haplotype_file, on=["CHROM", "POS"], how="inner"
    )

    # Write the simulated "HAPLOTYPE DATA" to a file.
    print('\nWriting simulated set haplotype data to a file "%s" .... ' % (out_sim))
    final_simulated_haplotype_file.to_csv(out_sim, sep="\t", index=False)

    print("\nProcess completed !!! :) :) ")
    print()


# this function is designed to create a right skewed simulation data
# the number represent the size of RBphased haplotype blocks
# The sum of the haplotype size should be equal to the size of total heterozyous variants.
def make_poissonish(n, num_bins):
    if n > 30 * num_bins:
        print("requested n exceeds 30 / bin")
        exit(-1)
    if n < 3 * num_bins:
        print("requested n cannot fill 3 / bin")
        exit(-1)

    # Disperse minimum quantity per bin in all bins, then determine remainder
    lst = [3 for _ in range(num_bins)]
    number_remaining = n - num_bins * 3

    # Allocate counts to all bins using a truncated Poisson
    for i in range(num_bins):
        # dial the rate up or down depending on whether we're falling
        # behind or getting ahead in allocating observations to bins
        rate = number_remaining / float(num_bins - i)  # avg per remaining bin

        # keep generating until we meet the constraint requirement (acceptance/rejection)
        while True:
            x = np.random.poisson(rate)
            if x <= 27 and x <= number_remaining:
                break
        # Found an acceptable count, put it in this bin and move on
        lst[i] += x
        number_remaining -= x

    # If there are still observations remaining, disperse them
    # randomly across bins that have remaining capacity
    while number_remaining > 0:
        i = np.random.randint(0, num_bins)
        if lst[i] >= 30:  # not this one, it's already full! #* changed
            continue
        lst[i] += 1
        number_remaining -= 1
    return lst


if __name__ == "__main__":
    main()
