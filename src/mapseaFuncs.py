import numpy as np
import os
import pandas as pd
from pybedtools import BedTool, set_tempdir
import secrets
import shutil
import string
from time import time

class BedFileManager:
    def __init__(self):
        self.bedtools_objects = {}

    def get_bedtool_object(self, file_path):
        if file_path in self.bedtools_objects:
            return self.bedtools_objects[file_path]
        else:
            self.bedtools_objects[file_path] = BedTool(file_path)
            return self.bedtools_objects[file_path]
        
def convert_to_int_or_str(value):
    '''
    This function converts value to integer if 
    the value present is float
    '''
    if isinstance(value, float):
        return int(value)
    else:
        return value
    
def mafPosIdentifer(intersect_output,entry):
    '''
    This function helps identify the exact location
    of the element inside the maf block sequence
    '''
    maf_start = int(entry[2]) #maf start
    maf_length = int(entry[3]) #maf length
    maf_chrom_length = int(entry[5]) #maf chromosome length
    for piece in intersect_output:
        start_mark = -1 #arbitary value to mark the start
        end_mark = -1 #arbitary value to mark the end
        start = piece[0] # start position of element wrt +ve
        seq = entry[6] #sequence
        if entry[4] == "+": #if the strand is positive
            follow = -1 # to mark the position of bases (not gaps)
            for index, base in enumerate(seq): #iterate through the sequence
                if base != "-": #if the base is not a gap
                    follow += 1 
                if follow + maf_start == piece[0]: #if the position of the base is the same as the start of the element
                    if start_mark == -1: #if the start is not marked
                        start_mark = index
                if follow + maf_start == start + piece[1] - 1: #if the position of the base is the same as the end of the element
                    end_mark = index
                    break #break the loop
            if start_mark == -1: #if the start is behind the block
                start_mark = 0  #start from the beginning of the block
            if end_mark == -1: #if the  end is ahead of the block
                end_mark = len(seq)  #end at the end of the block
            piece.append(start_mark) # this is 0-based index [start,end]
            piece.append(end_mark)
        elif entry[4] == "-": #if the strand is negative
            follow = -1
            start = piece[0]
            cal_start = maf_chrom_length - (maf_start + maf_length) #calculated start for wrt +ve strand
            seq = seq[::-1] # reverse sequence for -ve strand
            for index, base in enumerate(seq):
                if base != "-":
                    follow += 1
                if follow + cal_start == piece[0]:
                    if start_mark == -1:
                        start_mark = index
                if follow + cal_start == start + piece[1] - 1:
                    end_mark = index
                    break
            if start_mark == -1: #if the start is behind the block
                start_mark = 0  #start from the beginning of the block
            if end_mark == -1: #if the end is ahead of the block
                end_mark = len(seq)  #end at the end of the block
            start_mark_rev = (len(seq) - 1) - end_mark #reverse the end index to start for -ve strand
            end_mark_rev = (len(seq) - 1) - start_mark #reverse the start index to end for -ve strand
            piece.append(start_mark_rev)
            piece.append(end_mark_rev)
    return intersect_output

def mafEntryIntersect(entry,path,refer_dict,hsa_map,intersect_ratio,bed_file_manager=BedFileManager()):
    '''
    Given a maf entry as a list and a path to a bed file, 
    return the intersecting bed entries as a list. entry 
    must be a list of a line from maf starting with s.  
    '''
    species,chrN = entry[1].split(".") #identify species
    chrN = convert_to_int_or_str(chrN[3:]) #for the 2way setup files
    if hsa_map is not False:
        chrN = convert_to_int_or_str(hsa_map.loc[refer_dict[species][0],chrN[3:]]) #identify chrosome number of species
    file_access = "%s/%s/chr%s.bed" %(path,refer_dict[species][1],chrN) 
    bedfile = bed_file_manager.get_bedtool_object(file_access) #if bed files are stored to memory
    # bedfile = BedTool(file_access) #load hunter bed file to object
    maf_start = int(entry[2]) #maf start
    maf_length = int(entry[3]) #maf length
    maf_chrom_length = int(entry[5]) #maf chromosome length
    if entry[4] == "+": #if the strand is positive
        entry_pos = BedTool("chr%s\t%s\t%s" %(chrN,maf_start,maf_start+maf_length),from_string=True) #maf2bed
    elif entry[4] == "-": #if the strand is negative
        cal_start = maf_chrom_length - (maf_start + maf_length) #calculated start for wrt +ve strand
        cal_end = maf_chrom_length - maf_start #calculated end for wrt +ve strand
        entry_pos = BedTool("chr%s\t%s\t%s" %(chrN,cal_start,cal_end),from_string=True) #maf2bed
    intersect_bed = bedfile.intersect(entry_pos,wa=True,f=float(intersect_ratio)) #intersect
    intersect_output = [] #output list
    for line in intersect_bed:
        start = int(line[1])
        length = int(line[4])
        score = float(line[3])
        sign = line[5]
        intersect_output.append([start,length,"inBED",score,sign]) #append to output list, score is NA
    return mafPosIdentifer(intersect_output,entry)    

def searchElements(locator, spec_elem):
    '''
    This function finds the location tags in the same row
    and then return the element attributes of that location
    '''
    a_arr = np.array(locator)  # Convert list to a NumPy array for easier comparison
    b_cols = spec_elem[:,5:7] # Extract the 6th and 7th columns from array spec_elem
    for row in b_cols: # Iterate through each row in b_cols
        if np.array_equal(row, a_arr):# Check if the elements in a_arr are present in the current row
            return spec_elem[np.where(np.all(b_cols == a_arr, axis=1))][:, :5] # If found, return the first five elements of the matching row
    return False # If no match is found, return False

def dataframeRowNames(block_array,referDict):
    '''
    This function returns the row names of the dataframe
    '''
    row_names = []
    for entry in block_array:
        species,chrOrHsa = entry[1].split(".") #identify species
        row_names.append("%s@%s" %(referDict[species][2],chrOrHsa[3:]))
    return(row_names)

def blockExtractToDataframe(block_array,block_extract,referDict):
    '''
    This function extracts the block information from 
    entries having atleast one hit to a dataframe. 
    Include fullMAF, partMAF, FGAP, GAP.xx and FNotA'''
    # To define the column number of the dataframe
    location = [] #empty list to extract unqiue number of locations
    for spec_elem in block_extract: #iterate through the species
        location.append(spec_elem[:,5:7].tolist()) #append the location of the elements
    location_unique = np.unique(np.array([item for sublist in location for item in sublist]),axis=0) #unique locations
    col_nos = len(location_unique) 
    # To define the row number of the dataframe
    row_nos = len(block_array)
    # Define an empty dataframe
    df = pd.DataFrame(columns=range(col_nos),index=range(row_nos),dtype=object)
    # Fill the dataframe
    for idx, row in df.iterrows():
        for i, loc in enumerate(location_unique):
            attributes = searchElements(loc, block_extract[idx]) 
            if attributes is not False:
                attributes = attributes[0] #convert 2D to 1D array
                seq = block_array[idx][6][loc[0]:loc[1]+1] #identify the sequence
                converted = list(attributes) #convert the numpy array to list
                if converted[1] == len(seq)-seq.count("-"): #if the length of the sequence is the same as the length of the element
                    converted[2] = "fullMAF" #mark as fullMAF, i.e. the sequence is completely inside the alignment block
                    df.loc[idx, i] = converted + [len(seq)-seq.count("-")] + [seq] #append the sequence length without gaps from alignment and the sequence
                else:
                    converted[2] = "partMAF" #mark as partMAF, i.e. the sequence is partially inside the alignment block
                    df.loc[idx, i] = converted + [len(seq)-seq.count("-")] + [seq]
            else: 
                identify_seq = block_array[idx][6][loc[0]:loc[1]+1] #identify the sequence
                blanks = identify_seq.count("-") #count the number of gaps
                if blanks == len(identify_seq): #if the sequence is completely a gap
                    df.loc[idx, i] = [".",".","FGAP",".",".",0,"."] 
                elif blanks == 0: #if the sequence is complete, but not annotated
                    df.loc[idx, i] = [notInBed_Start(loc,block_array[idx][6],block_array[idx]),".","FNotA",".",block_array[idx][4],len(identify_seq)-identify_seq.count("-"),identify_seq]  
                else: #if the sequence has gaps
                    df.loc[idx, i] = [notInBed_Start(loc,block_array[idx][6],block_array[idx]),".","{}{:.2f}".format("GAP", blanks/len(identify_seq)%1).replace("0.", "."),".",block_array[idx][4],len(identify_seq)-identify_seq.count("-"),identify_seq]
    # Define the row names
    index_name = dataframeRowNames(block_array,referDict) 
    df.index = index_name #set the index of the dataframe
    return([df,location_unique]) #return the dataframe and the unique locations

def rejectedbBlockExtractToDataframe(block_array,unique_positions,referDict):
    '''
    This function extracts the block information from
    rejected entries to a dataframe, only FGAP, GAP.xx and FNotA'''
    # To define the column number of the dataframe
    col_nos = len(unique_positions) 
    # To define the row number of the dataframe
    row_nos = len(block_array)
    # Define an empty dataframe
    df = pd.DataFrame(columns=range(col_nos),index=range(row_nos),dtype=object)
    # Fill the dataframe
    for idx, row in df.iterrows():
        for i, loc in enumerate(unique_positions):
            identify_seq = block_array[idx][6][loc[0]:loc[1]+1] #identify the sequence
            blanks = identify_seq.count("-") #count the number of gaps
            if blanks == len(identify_seq): #if the sequence is completely a gap
                df.loc[idx, i] = [".",".","FGAP",".",".",0,"."] 
            elif blanks == 0: #if the sequence is complete, but not annotated
                df.loc[idx, i] = [notInBed_Start(loc,block_array[idx][6],block_array[idx]),".","FNotA",".",block_array[idx][4],len(identify_seq)-identify_seq.count("-"),identify_seq]  
            else: #if the sequence has gaps
                df.loc[idx, i] = [notInBed_Start(loc,block_array[idx][6],block_array[idx]),".","{}{:.2f}".format("GAP", blanks/len(identify_seq)%1).replace("0.", "."),".",block_array[idx][4],len(identify_seq)-identify_seq.count("-"),identify_seq]
    # Define the row names
    index_name = dataframeRowNames(block_array,referDict) 
    df.index = index_name #set the index of the dataframe
    return(df)

def notInBed_Start(location,whole_sequence,array_from_block):
    '''
    This function identifies the start position of the sequence
    that is not annotated in the bed file
    '''
    sign = array_from_block[4]
    if sign == "+": #if the strand is positive
        pre_start = whole_sequence[:location[0]] #extract the sequence before the start
        bases_in_pre_start = len(pre_start) - pre_start.count("-") #count the number of bases in the pre_start sequence
        pos_start = bases_in_pre_start + int(array_from_block[2]) #start position
    elif sign == "-": #if the strand is negative
        pre_start = whole_sequence[:location[1]+1] #extract the sequence before the end position (since this is negative strand)
        bases_in_pre_start = len(pre_start) - pre_start.count("-") #count the number of bases in the pre_start sequence
        neg_start = bases_in_pre_start + int(array_from_block[2]) - 1 #start position based on negatve strand
        pos_start = int(array_from_block[5]) - (neg_start + 1) #start position based on positive strand
    return pos_start

def write_dataframe_to_file(df, file_path):
    '''
    Write a Pandas DataFrame to a file with tab-separated values.
    '''
    # Open the file in write mode
    with open(file_path, 'w') as file:
        # Iterate through each column
        for column_name in df.columns:
            # Write the column name to the file
            file.write(f'{column_name}\n')
            
            # Write index names and corresponding values for the column
            for index_name, value_list in zip(df.index, df[column_name]):
                # Convert the list to a tab-separated string
                value_str = '\t'.join(map(str, value_list))
                file.write(f'{index_name}: {value_str}\n')
            
            # Add a separator between columns
            file.write('\n')

def remove_tmp_files(temp_dir):
    '''
    Remove temporary files from the temporary directory.
    '''
    for filename in os.listdir(temp_dir):
        if ".tmp" in filename:
            file_path = os.path.join(temp_dir, filename)
            os.remove(file_path)

def generate_temp_id():
    '''
    Generate a random ID for temporary files.
    '''
    characters = string.ascii_letters + string.digits + "_-"
    random_id = ''.join(secrets.choice(characters) for _ in range(15))
    return random_id

def process_file(maf_file, bed_file_path, output_file, temp_dir, intersect_ratio, referDict, hsa_map):
    '''Process a MAF file and map the elements onto it. Write the results to an output file.'''

    temp_id = generate_temp_id()
    os.makedirs(f"{temp_dir}/tmp_%s" %temp_id, exist_ok=True) # Create temp directory if it doesn't exist
    pybedttoolstempdir = f"{temp_dir}/tmp_%s" %temp_id

    set_tempdir(pybedttoolstempdir) # Set the temporary directory for pybedtools

    start_time = time()
    # total_blocks = int(subprocess.check_output("cat %s | grep \"^a\" | wc -l" %maf_file, shell=True, text=True)) # Get total number of blocks in the MAF file

    with open(output_file, "w") as write_file: # Create the output file
        write_file.write("")

    count = 0
    current_block = []
    with open(maf_file, "r") as file:
        for line in file:
            line = line.rstrip()
            if line.startswith("a"):  # Start of a new block
                new_id = line.split(" ")[1]
                count += 1
                if count%5 == 0:
                    remove_tmp_files(pybedttoolstempdir)
                # print("Processing block %i of %i" %(count-1,total_blocks), end="\r")
                if current_block: # Process the current block
                    block_array = [] # Array to store the block
                    for block_line in current_block: # Iterate through the block lines of "s"
                        attribs = block_line.split() # Split the line into attributes
                        if int(attribs[3]) >= 30: # Processing: If the entry is less than 30bp, remove it
                            block_array.append(attribs) # Add the line to the block array 
                    if len(block_array) > 1: # If there are more than one line in the block after processing
                        block_extract = [] # Array to store the block extract values
                        refined_block_array = [] 
                        rejected_block_array = []
                        for entry in block_array: # Iterate through the block array
                            extract = np.array(mafEntryIntersect(entry,bed_file_path,referDict,hsa_map,intersect_ratio), dtype=object) # Extract the intersecting bed entries of elements
                            if len(extract) > 0: #if there is atleast one hit
                                refined_block_array.append(entry) # Add the entries to the refined block array
                                block_extract.append(extract) # Add the arrays to the block extract for each entry in the block
                            else: #if there is no hit
                                rejected_block_array.append(entry) # Add the entries to the rejected block array
                        toWrite = blockExtractToDataframe(refined_block_array,block_extract,referDict) # Convert the block array and block extract to a dataframe
                        if len(rejected_block_array) > 0:
                            rejected_toWrite = rejectedbBlockExtractToDataframe(rejected_block_array,toWrite[1],referDict) # Convert the rejected block array to a dataframe
                        with open(output_file, "a") as write_file:
                            if toWrite[0].empty is False: # If the dataframe is not empty
                                write_file.write("\n" + old_id + "\n") # Write the block id to the file
                                for column_name in toWrite[0].columns:
                                    write_file.write(f"{column_name + 1}\n") # Write the column name to the file
                                    for index_name, value_list in zip(toWrite[0].index, toWrite[0][column_name]):
                                        value_str = '\t'.join(map(str, value_list)) # Convert the list to a tab-separated string
                                        write_file.write(f'{index_name}:\t{value_str}\n') # Write the index name and the value string to the file
                                    if len(rejected_block_array) > 0: # If there are rejected entries
                                        for index_name, value_list in zip(rejected_toWrite.index, rejected_toWrite[column_name]): # Iterate through the rejected entries
                                            value_str = '\t'.join(map(str, value_list)) # Convert the list to a tab-separated string
                                            write_file.write(f'{index_name}:\t{value_str}\n') # Write the index name and the value string to the file
                    current_block = []  # Reset current block for the next one
                old_id = new_id  # Update the current block id
            if line.startswith("s"):
                current_block.append(line)  # Add line to the current block
                
        # Process the last block, a copy of all the code above
        # print("Processing block %i of %i" %(count,total_blocks))
        if current_block:
            block_array = [] # Array to store the block
            for block_line in current_block: # Iterate through the block lines of "s"
                attribs = block_line.split() # Split the line into attributes
                if int(attribs[3]) >= 30: # Processing: If the entry is longer than 30 bp
                    block_array.append(attribs) # Add the line to the block array 
            if len(block_array) > 1: # If there are more than one line in the block after processing
                block_extract = [] # Array to store the block extract values
                refined_block_array = [] 
                rejected_block_array = []
                for entry in block_array: # Iterate through the block array
                    extract = np.array(mafEntryIntersect(entry,bed_file_path,referDict,hsa_map,intersect_ratio), dtype=object) # Extract the intersecting bed entries of elements
                    if len(extract) > 0: #if there is atleast one hit
                        refined_block_array.append(entry) # Add the entries to the refined block array
                        block_extract.append(extract) # Add the arrays to the block extract for each entry in the block
                    else: #if there is no hit
                        rejected_block_array.append(entry) # Add the entries to the rejected block array
                toWrite = blockExtractToDataframe(refined_block_array,block_extract,referDict) # Convert the block array and block extract to a dataframe
                if len(rejected_block_array) > 0:
                    rejected_toWrite = rejectedbBlockExtractToDataframe(rejected_block_array,toWrite[1],referDict) # Convert the rejected block array to a dataframe
                with open(output_file, "a") as write_file:
                    if toWrite[0].empty is False: # If the dataframe is not empty
                        write_file.write("\n" + old_id + "\n") # Write the block id to the file
                        for column_name in toWrite[0].columns:
                            write_file.write(f"{column_name + 1}\n") # Write the column name to the file
                            for index_name, value_list in zip(toWrite[0].index, toWrite[0][column_name]): # Iterate through the index names and values
                                value_str = '\t'.join(map(str, value_list)) # Convert the list to a tab-separated string
                                write_file.write(f'{index_name}:\t{value_str}\n') # Write the index name and the value string to the file
                            if len(rejected_block_array) > 0: # If there are rejected entries
                                for index_name, value_list in zip(rejected_toWrite.index, rejected_toWrite[column_name]):
                                    value_str = '\t'.join(map(str, value_list)) # Convert the list to a tab-separated string
                                    write_file.write(f'{index_name}:\t{value_str}\n') # Write the index name and the value string to the file

    end_time = time()

    shutil.rmtree(pybedttoolstempdir) # Remove the temporary directory
    print("Finished processing: %s" %maf_file)
    print("Time elapsed: %.2f hrs\n" %((end_time - start_time)/60/60))
