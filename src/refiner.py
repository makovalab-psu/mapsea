''' Demo Run Command: python3 process.py <filename.dat> 2 10 > filname.df'''

import argparse
import multiprocessing as mp
import pandas as pd
from refinerFuncs import *
import sys
import warnings

if __name__ == "__main__":

    # Create the argument parser
    parser = argparse.ArgumentParser()

    # Define required arguments with flags
    parser.add_argument('-d', '--datInput', required=True, type=str,
                        help="Input dat file output by mapsea (e.g., 'chr22_bonobo_vs_chr22_borang.dat')")
    parser.add_argument('-f', '--flank', required=True, type=int,
                        help="The flank allowed for the element to be termed 'shared' (e.g. 3)")
    parser.add_argument('-c', '--cores', required=False, type=int,
                        help="The number of cores available for this job (e.g. 10)")
    parser.add_argument('-o', '--outFile', required=False, type=str,
                        help="Output file path (e.g., 'hs1.chr1.2023v2.processed.analysed.quadron.df')")
    parser.add_argument('-m', '--mutInfo', required=False, action="store_true",
                        help="Include mutation information in the output file")

    # Parse the arguments
    args = parser.parse_args()

    # Access the arguments
    input_file = args.datInput
    flank = args.flank
    mutInfo = args.mutInfo

    if args.cores is not None:
        noDivisions = args.cores
    else:
        noDivisions = mp.cpu_count() - 1
        print(f"Number of cores not specified. Using {noDivisions} cores.")

    if args.outFile is not None:
        output_file = args.outFile
    else:
        output_file = sys.stdout
        print(f"Output file not specified. Writing to stdout")

    if flank > 3:
        warnings.warn("Higher flank values increase the probability of SequenceMergeError", stacklevel=2)

    inputData = open(input_file, "r") #open the file

    df = outputConvertToDataframe(inputData) #convert the output to a dataframe
    df = df.reset_index() #reset the index

    grouped_df = df.groupby(["ID"]) # group the dataframe by ID
    getGroups = grouped_df.groups # get the groups

    clusteredDict = clusterDict(getGroups, noDivisions) #cluster the data into groups   
    divisions = [(df.loc[clusteredDict[id]].set_index(["ID", "G4", "CHR", "SPECIES","START"]), flank, mutInfo) for id in clusteredDict.keys()]

    pool = mp.Pool(processes=noDivisions) # Create a pool of worker processes
    results = pool.starmap(removeRedundant, divisions) # Process the files in parallel
    pool.close() #close the pool

    joined = pd.concat(results).sort_values(by=["ID", "G4", "CHR", "SPECIES","START"]) #concatenate the results and sort
    # joined = removeAlignmentDuplicates(joined) #remove alignment duplicates
    joined.to_csv(output_file, sep="\t", na_rep="NULL") #save the results
    
    #print(f"Removed redundant and assigned mutations for {input_file}") #removed since the above command has sys.stdout