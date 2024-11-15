import argparse
import gzip
import json
from mapseaFuncs import process_file
import multiprocessing as mp
import os
import pandas as pd
import random
import subprocess
from time import time

if __name__ == "__main__":

    # Create the argument parser
    parser = argparse.ArgumentParser()

    # Define required arguments with flags
    parser.add_argument('-m', '--mafInput', required=True, type=str,
                        help="Input gzipped maf file path (e.g., 'demo_hs1.chr1.2023v2.processed.maf.gz')")
    parser.add_argument('-b', '--bedFolder', required=True, type=str,
                        help="Folder name to access the bed files (without / at the end); the folder should have a tree structure, with each species having their bed files under a directory separated by chromosomes (e.g., 'pqsfinderOutput')")
    parser.add_argument('-o', '--outFile', required=True, type=str,
                        help="Output file path (e.g., 'hs1.chr1.2023v2.processed.analysed.quadron.dat')")
    parser.add_argument('-t', '--tempDir', required=True, type=str,
                        help="Temporary directory path (e.g., 'tmp/hsa1')")
    parser.add_argument('-r', '--intersectRatio', required=True, type=float, 
                        help="Intersection ratio (e.g., 1.0)")
    parser.add_argument('-d', '--speciesDict', required=True, type=str,
                        help="Species dictionary file path, a json file only (e.g. 'speciesDict.json')")
    parser.add_argument('-f', '--mapFile', required=False, type=str,
                        help="Homolog chromosome map file path")
    parser.add_argument('-c', '--cores', required=False, type=int,
                        help="The number of cores available for this job (e.g. 10)")

    # Parse the arguments
    args = parser.parse_args()

    # Access the arguments
    input_file = args.mafInput
    bed_file_path = args.bedFolder
    output_file = args.outFile
    temp_dir = args.tempDir
    intersect_ratio = args.intersectRatio
    referDictPath = args.speciesDict
    hsa_align_map = args.mapFile

    # Number of cores
    if args.cores is not None:
        noDivisions = args.cores
    else:
        noDivisions = mp.cpu_count() - 1
        print(f"Number of cores not specified. Using {noDivisions} cores.")

    # Get the species name dictionary
    with open(referDictPath, 'r') as json_file:
        referDict = json.load(json_file)

    #The species name and corresponding hsa table 
    if hsa_align_map is not None:
        hsa_map = pd.read_csv(hsa_align_map,sep="\t",header=0,index_col=0)
    else:
        hsa_map = False

    metadata_filename = "metadata.txt"
    with open(metadata_filename, "w") as write_file:
        ## Metadata for the output file
        write_file.write("## {METADATA}\n")
        write_file.write("## INPUT FILE: %s\n" %input_file)
        write_file.write("## OUTPUT FILE: %s\n" %output_file)
        write_file.write("## HSA MAP: %s\n" %hsa_align_map)
        write_file.write("## SPECIES DICTIONARY: %s\n" %referDictPath)
        write_file.write("## INTERSECTION RATIO (f): %.2f\n" %intersect_ratio)
        write_file.write("## \n")
        write_file.write("## {TYPE}\n")
        write_file.write("## Absent in .BED file: \n")
        write_file.write("##   FGAP: No sequence present\n")
        write_file.write("##   FNotA: Sequence is non-annotated\n")
        write_file.write("##   GAP.xx: Sequence has xx% gaps\n")
        write_file.write("## Present in .BED file: \n")
        write_file.write("##   partMAF: Sequence partly in alignment\n")
        write_file.write("##   fullMAF: Sequence fully in alignment\n")
        write_file.write("## \n")
        write_file.write("## {SPECIESID}\n")
        for identifier in referDict:
            write_file.write(f"## {referDict[identifier][2]}: {referDict[identifier][1]}\n")
        write_file.write("## \n")
        write_file.write("## {STRUCTURE}\n")
        write_file.write("## #BLOCK ID\n")
        write_file.write("## NUMBER\n")
        write_file.write("## SPECIESID@CHR:\tSTART\tQUERY_LENGTH\tTYPE\tSCORE\tSTRAND\tALIGNMENT_LENGTH\tSEQUENCE\n")    

    generated_numbers = set()
    def generate_unique_random_number():
        while True:
            # Generate a 10-digit random number
            new_number = random.randint(1000000000, 9999999999)
            if new_number not in generated_numbers:
                generated_numbers.add(new_number) # Add the number to the set
                return new_number

    def split_file(input_file, temp_dir, max_size_mb=10):
        max_size_bytes = max_size_mb * 1024 * 1024  # Convert MB to bytes
        output_files = [] # List to store the names of the output files
        with gzip.open(input_file, 'rt', encoding='utf-8') as f: # Open the input file
            current_size = 0 # Variable to store the current size of the output file
            current_lines = [] # List to store the current lines to be written to the output file
            for line in f: # Loop through each line in the input file
                if line.startswith("a"):
                    line = line.replace("a",f"a #{generate_unique_random_number()} ")
                current_size += len(line.encode('utf-8')) # Get the size of the current line in bytes
                current_lines.append(line) # Add the current line to the list of lines
                if current_size >= max_size_bytes and line.strip() == "": # If the current size is greater than the max size and the current line is empty
                    output_files.append(current_lines) # Add the current lines to the list of output files
                    current_lines = [] # Reset the current lines list
                    current_size = 0 # Reset the current size
            if current_lines: # If there are any remaining lines
                output_files.append(current_lines) # Add the remaining lines to the list of output files
        for i, lines in enumerate(output_files): # Loop through each output file
            file_number = "{:03d}".format(i+1) # Format the file number with leading zeros
            input_filename = input_file.split("/")[-1] # Get the filename from the input file path
            output_file = f"{temp_dir}/{input_filename}_part_{file_number}" # Create the output file name
            with open(output_file, 'w', encoding='utf-8') as f_out: # Open the output file
                f_out.writelines(lines) # Write the lines to the output file
        return i+1 # Return the total number of output files

    def merge_files(output_file):
        merge_command = ["cat", "metadata.txt", f"{temp_dir}/{input_filename}_part_*.dat", ">", output_file]
        subprocess.run(" ".join(merge_command), shell=True)

        # Cleanup split files and output files
        subprocess.run(f"rm -f {temp_dir}/{input_filename}_part*", shell=True)

    s = time()

    os.makedirs(temp_dir, exist_ok=True) # Create temp directory if it doesn't exist

    total_files = split_file(input_file, temp_dir) # Split input file into smaller files
    
    input_filename = input_file.split("/")[-1] # Get the filename from the input file path
    args = [(f"{temp_dir}/{input_filename}_part_{i:03d}", bed_file_path, f"{temp_dir}/{input_filename}_part_{i:03d}.dat", temp_dir, intersect_ratio, referDict, hsa_map) for i in range(1, total_files+1)]

    pool = mp.Pool(processes=noDivisions) # Create a pool of worker processes
    pool.starmap(process_file, args) # Process the files in parallel

    merge_files(output_file) # Merge the output files into a single file

    e = time()

    print("Total time elapsed: %.2f hrs" %((e - s)/60/60)) # Print the time taken to process the files

    os.remove(metadata_filename) # Remove the metadata file