from Bio import Seq, Align
import networkx as nx
import numpy as np
import pandas as pd
import re
import warnings

from pandas.errors import SettingWithCopyWarning, PerformanceWarning
warnings.simplefilter("ignore", category=SettingWithCopyWarning)
warnings.simplefilter("ignore", category=PerformanceWarning)

''' Functions for reading output file to pandas dataframe, and multiprocessing '''

def assign_value(row):
    '''Assign index to fullMAF entries'''
    if row['TYPE'] == 'fullMAF':
        return 0
    else:
        return np.nan

def outputConvertToDataframe(input_result: str) -> pd.DataFrame:
    current = [] #empty list to store
    index_col = [] #empty list to store the index
    for line in input_result:
        if not line.startswith("##"): #ignore the metadata
            data = line.strip().split() #split the line
            if len(data) == 1: #if the line is either an id or g4 number
                if line.startswith("#"): #if the line is an id
                    id = data[0] #store the id
                else: #if the line is a g4 number
                    nos = data[0] #store the g4 number
            elif len(data) > 1:
                species = data[0].rstrip(":").split("@")[0] #store the species and chromosome
                chr = data[0].rstrip(":").split("@")[1]
                index_col.append(tuple([id, "{:02d}".format(int(nos)), chr, species])) #append the id, g4 number and species to the index
                content = data[1:] #store the rest of the data
                current.append(content) #append the content to the current list

    column_names = ["START", "QUERY_LENGTH", "TYPE", "SCORE", "STRAND", "ALIGNMENT_LENGTH", "SEQUENCE"]
    result = {col: values for col, values in zip(column_names, zip(*current))} #store the data in a dictionary
    index = pd.MultiIndex.from_tuples(index_col, names=["ID", "G4", "CHR", "SPECIES"]) #create a multiindex
    df = pd.DataFrame(result, index=index).sort_index() #create a dataframe and sort the index

    #define the type of the numeric columns but ignore nan
    df["INDEX"] = df.apply(assign_value, axis=1) #set the index of where the G4 starts in the sequence
    numeric_columns = ["START", "QUERY_LENGTH", "SCORE", "ALIGNMENT_LENGTH"]
    df[numeric_columns] = df[numeric_columns].apply(pd.to_numeric, errors='coerce')
    return df

def clusterDict(originalDict, noOfGroups):
    '''Split the dataframe adequately for post-processing'''
    newDict = dict()
    count = 0
    pre_values = []
    pre_key = 0
    lenOfGroups = len(originalDict) // (noOfGroups-1) #get the length of groups to be divided depending on the number of cores
    for key, value in originalDict.items(): #iterate through the dictionary containing groups based on grouped_df
        count += 1
        if count % lenOfGroups == 0 or count == len(originalDict): #if the count is a multiple of length or the last line of the df
            newDict[pre_key] = pre_values + list(value) #add the list of values to the dict with key id
            pre_values = []
            pre_key = 0 #reset
        else:
            if pre_key == 0: #if a new division is starting 
                pre_key = key #set the key
            pre_values += list(value) #add the values to the list
    return newDict

''' Functions for removing redundant G4s and assigning mutations '''

def align_sequences(seq1, seq2):
    '''Align two DNA sequences and return the aligned sequences and the difference'''

    matrix = np.array(
    [[ 1, 1, 0, 0, 0, 0, 0, 0, 0,.5,.5],
     [ 1, 1, 0, 0, 0, 0, 0, 0, 0,.5,.5],
     [ 0, 0, 1, 1, 0, 0, 0, 0, 0,.5,.5],
     [ 0, 0, 1, 1, 0, 0, 0, 0, 0,.5,.5],
     [ 0, 0, 0, 0, 1, 1, 0, 0, 0,.5,.5],
     [ 0, 0, 0, 0, 1, 1, 0, 0, 0,.5,.5],
     [ 0, 0, 0, 0, 0, 0, 1, 1, 0,.5,.5],
     [ 0, 0, 0, 0, 0, 0, 1, 1, 0,.5,.5],
     [ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
     [.5,.5,.5,.5,.5,.5,.5,.5, 0,.5,.5],
     [.5,.5,.5,.5,.5,.5,.5,.5, 0,.5,.5]])

    substitution_matrix = Align.substitution_matrices.Array(alphabet="AaCcGgTt-Nn", data=matrix) 
    # This is a very slow function with Bio.Align
    seq1 = Seq.Seq(seq1)
    seq2 = Seq.Seq(seq2)
    aligner = Align.PairwiseAligner() # Create alignment object
    aligner.substitution_matrix = substitution_matrix # Set the substitution matrix
    aligner.query_open_gap_score = -2  # Penalty for opening the query sequence gaps
    aligner.query_extend_gap_score = -1  # Penalty for extending the query sequence gaps
    aligner.query_end_gap_score = 0 # favour end gaps in query
    aligner.target_extend_gap_score = -100 # Penalty for extending the reference sequence gaps
    aligner.target_gap_score = -100 # Penalty for the reference sequence gaps
    aligner.target_open_gap_score = -100 # Penalty for opening the reference sequence gaps
    alignment = aligner.align(seq1, seq2) # Perform sequence alignment
    best_alignment = alignment[0] #max(alignment, key=lambda al: al.score) #store the best alignment
    return best_alignment[0], best_alignment[1]

def check_hit(rangePrev, rangePost, flank):
    '''Check if the start and end are within the flank'''
    prev = np.cumsum(np.array(rangePrev.iloc[0]))
    post = np.cumsum(np.array(rangePost.iloc[0]))
    subtract = prev - post #subtract the start and alignment_lengths
    if abs(subtract[0]) < flank + 1 and abs(subtract[1]) < flank + 1: # check if the start and end position is within flank bases
        return True
    
def find_prior_groups(input_list):
    '''Find the prior groups of the removed list'''
    output_list = []
    sequence = []
    for i, num in enumerate(input_list):
        if i == 0 or num != input_list[i - 1] + 1:
            if sequence:
                output_list.append(sequence)
            sequence = [num]
        else:
            sequence.append(num)
    if sequence:
        output_list.append(sequence)
    output_list = [[seq[0]-1] + seq for seq in output_list]
    output_list = [["{:02d}".format(j) for j in i] for i in output_list]
    return output_list

def strand_fullmaf(dataframe):
    '''When there are overlapping "G4" entries, determine the 
       strand based on the G4 with fullMAF annotation'''
    [strand, seq] = dataframe.iloc[0][["STRAND", "SEQUENCE"]]
    seq = seq.replace("-", "").lower()
    if seq.count("g") > seq.count("c"):
        return strand
    elif seq.count("g") < seq.count("c"):
        return "-" if strand == "+" else "+"
    else:
        return False

def evaluate_redundant(dataframe):
    '''Function to evaluate and merge redundant G4s in a dataframe, inside a group of G4s with the same ID'''
    for_strand = dataframe[(dataframe["TYPE"] != "fullMAF") & (dataframe["TYPE"] != "FGAP")] #get the strand where the type is not fullmaf or fgap 
    if not for_strand.empty: #if the dataframe is empty
        strand = for_strand["STRAND"].iloc[0] #get the strand

    if set(dataframe["TYPE"]) == {"fullMAF"}: 
        strand = strand_fullmaf(dataframe)

    first_row = dataframe.iloc[0] #get the first row
    last_row = dataframe.iloc[-1] #get the last row
    part_dataframe = dataframe[dataframe["TYPE"] != "FGAP"] #get the dataframe where the type is not fgap
    if not part_dataframe.empty: #if the dataframe without fgap is not empty
        true_row_all_rows = part_dataframe[part_dataframe["TYPE"] == "fullMAF"] # can be an empty part_dataframe, or a single row
        if not true_row_all_rows.empty: #if the dataframe with only fullmaf is not empty
            true_row = true_row_all_rows.iloc[0] #get the first row
        min_start_row = part_dataframe.sort_values(by=["START","SEQ_LEN","END"], ascending=[True, False, True]).iloc[0] # get the row with the minimum start, maximum sequence length and minimum end
        max_end_row = part_dataframe.sort_values(by=["END","SEQ_LEN","START"], ascending=[False, False, True]).iloc[0] # get the row with the maximum end, maximum sequence length and minimum start
        if (min_start_row["START"] == max_end_row["START"] and min_start_row["END"] == max_end_row["END"] and true_row_all_rows.empty): #if the start and end are the same and no fullmaf row
            filtered_row = part_dataframe.loc[part_dataframe["SEQ_LEN"].idxmax()] #get the row with the maximum sequence length
            filtered_row["G4"] = f"[{first_row['G4']}-{last_row['G4']}]" #store the G4 number in the format [a-z]
        else:
            base_count = -1 #initialize the base count        
            if strand == "+": #if the strand is positive
                first_row_sequence = min_start_row["SEQUENCE"] #get the sequence
            elif strand == "-": #if the strand is negative
                first_row_sequence = min_start_row["SEQUENCE"][::-1] #get the sequence in reverse
            for i, base in enumerate(first_row_sequence): #iterate through the min_start sequence
                if base != "-": #if the base is not a gap
                    base_count += 1 #increment the base count
                    if base_count + int(min_start_row["START"]) == int(max_end_row["START"]): #if the base count + the start of the min_start row is equal to the start of the max_end row
                        last_seq_mark = i #store the last sequence mark
                    if not true_row_all_rows.empty: #if there is a fullmaf entry
                        if base_count + int(min_start_row["START"]) == int(true_row["START"]): #if the base count + the start of the min_start row is equal to the start of the fullmaf row
                            true_seq_start_mark = base_count #store the true sequence start mark
            filtered_row = min_start_row.copy()  #copy the min_start row
            filtered_sequence = filtered_row["SEQUENCE"] #get the sequence
            filtered_row["G4"] = f"[{first_row['G4']}-{last_row['G4']}]" #store the G4 number in the format [a-z]
            if not min_start_row.equals(max_end_row): #if the min_start row is not equal to the max_end row
                try:
                    if strand == "+": #if the strand is positive
                        first_sequence = min_start_row["SEQUENCE"][:last_seq_mark] #get the sequence from the start to the last sequence mark
                        second_sequence = max_end_row["SEQUENCE"] #get the sequence of the max_end row
                        filtered_sequence = first_sequence + second_sequence[next((i for i, base in enumerate(second_sequence) if base != '-'), len(second_sequence)):] #combine the sequences, by removing the start gaps, if any, of the second sequence
                    elif strand == "-": #if the strand is negative
                        first_sequence = min_start_row["SEQUENCE"][::-1][:last_seq_mark] #get the sequence from the start to the last sequence mark in reverse
                        second_sequence = max_end_row["SEQUENCE"][::-1] #get the sequence of the max_end row in reverse
                        filtered_sequence = (first_sequence + second_sequence[next((i for i, base in enumerate(second_sequence) if base != '-'), len(second_sequence)):])[::-1] #combine the sequences, by removing the start gaps, if any, of the second sequence in reverse
                except: #Unable to solve this merge with current logic
                    print(f"RuntimeWarning: SequenceMergeError at -> ID: {filtered_row['ID']}, G4: {filtered_row['G4']}, CHR: {filtered_row['CHR']}, SPECIES: {filtered_row['SPECIES']}, START: {filtered_row['START'] if not filtered_row['START']==np.nan else 'NULL'}") #print the error
                filtered_row["SEQUENCE"] = filtered_sequence #store the filtered sequence
            filtered_row["ALIGNMENT_LENGTH"] = len(filtered_sequence.replace("-", "")) #store the alignment length by removing the gaps
            if not true_row_all_rows.empty: #if there is a fullmaf entry
                filtered_row["QUERY_LENGTH"] = true_row["QUERY_LENGTH"] 
                filtered_row["TYPE"] = true_row["TYPE"]
                filtered_row["SCORE"] = true_row["SCORE"]
                filtered_row["STRAND"] = true_row["STRAND"]
                filtered_row["INDEX"] = true_seq_start_mark #store the true sequence start mark
            else: #if there is no fullmaf entry
                blanks = filtered_sequence.count("-") #count the number of gaps
                if blanks == 0: #if there are no gaps
                    filtered_row["TYPE"] = "FNotA" #store the type as FNotA
                else: #if there are gaps
                    filtered_row["TYPE"] = "{}{:.2f}".format("GAP", blanks/len(filtered_sequence)%1).replace("0.", ".") #store the type as GAP and the percentage of gaps
        return filtered_row  #return the merged row
    else: #if the dataframe is empty
        fgap_row = dataframe.iloc[0] #get the first row
        fgap_row["G4"] = f"[{first_row['G4']}-{last_row['G4']}]" #store the G4 number in the format [a-z]
        return fgap_row #return the fgap row

def identifyRedundant(dataframe, flank, progress=False):
    '''Identify the redundant G4s in the dataframe'''
    grouped = dataframe.groupby(["ID"]) #group the dataframe by ID
    check_nested = grouped.apply(lambda x: x.index.get_level_values('G4').nunique()) #check if there are mutiple G4s in a block
    nested_dict = check_nested.to_dict() #convert the result to a dictionary
    removeGroups = dict()
    for idx, group in grouped:
        if nested_dict[idx[0]] > 1: #if there are multiple G4s in a single block
            ignore_first = 0 # to differentiate the first g4 of the block
            regrouped = group.groupby(["G4"]) #group by G4
            for reidx, regroup in regrouped: #iterate through the regrouped dataframe
                species_of_interest = regroup[regroup["TYPE"] == "fullMAF"].index[0][3] #get the species for which it is fullmaf
                if ignore_first != 0: #second g4 onwards
                    startEnd_post = regroup[regroup.index.get_level_values('SPECIES') == species_for_length].reset_index()[["START", "ALIGNMENT_LENGTH"]] #get the start and alignment length
                    hit = check_hit(startEnd_prev, startEnd_post, flank) #check if the start and alignment length are within the flank
                    if hit == None: #if the hit is None
                        species_for_length = species_of_interest #get the species for the next iteration
                        startEnd_prev = regroup[regroup.index.get_level_values('SPECIES') == species_for_length].reset_index()[["START", "ALIGNMENT_LENGTH"]] #get the start and alignment length
                    else: #if the hit is True
                        startEnd_prev = startEnd_post #make the start and alignment length of the post as the start and alignment length of the prev for the next search
                        if idx[0] not in removeGroups: #if the ID is not in the removeGroups
                            removeGroups[idx[0]] = [int(reidx[0])] #add the G4 number to the removeGroups
                        else: #if the ID is in the removeGroups
                            removeGroups[idx[0]].append(int(reidx[0])) #append the G4 number to the removeGroups
                else: #first g4 of the block
                    species_for_length = species_of_interest # find the species_for_length of the first row where the type is fullmaf
                    startEnd_prev = regroup[regroup.index.get_level_values('SPECIES') == species_for_length].reset_index()[["START", "ALIGNMENT_LENGTH"]] #get its start and alignment length
                ignore_first = 1 #set ignore_first to 1
        # if progress == True:
        #     print("Completed row %i of %i" %(dataframe.index.get_loc(idx).stop,len(dataframe)), end="\r")
    return removeGroups

def processRedundant(dataframe, removeGroups, progress=False):
    '''Remove the redundant G4s from the dataframe and merge them to a single one'''
    unduplicated = [] #empty list to store the unduplicated data
    df_removed = dataframe.copy() #copy the dataframe
    progress = 0 #set the progress to 0
    for id, value in removeGroups.items():
        prior_groups = find_prior_groups(value) #find the before group of the G4s which are duplicates
        progress += 1 #increment the progress
        for subset in prior_groups: #iterate through the subsets in removegroups
            subset = [str(i) for i in subset] #convert the subset to string
            df_part = dataframe.loc[(id, subset[0]):(id, subset[-1])].reset_index() #get the subset of the dataframe for the ID and the G4s from subset and reset the index
            df_part["SEQ_LEN"] = df_part["SEQUENCE"].apply(len) #get the length of the sequence including gaps
            df_part_len = len(df_part) #get the length of the full subset 
            filtered_rows = []
            for row in range(int(df_part_len/len(subset))):
                indices = np.arange(row,df_part_len,int(df_part_len/len(subset)))
                compare = df_part.iloc[indices] 
                compare["END"] = compare["START"] + compare["ALIGNMENT_LENGTH"] #get the end position
                filtered_rows.append(evaluate_redundant(compare)) #append the filtered rows to the filtered_rows list
            df_removed = df_removed.loc[~((df_removed.index.get_level_values('ID') == id) & (df_removed.index.get_level_values('G4').str.fullmatch('|'.join(subset))))] #remove the subset from the dataframe
            processed = pd.concat(filtered_rows,axis=1).T.drop(columns=["SEQ_LEN","END"]) #concatenate the filtered rows and set the index
            unduplicated.append(processed)
        # if progress == True:
        #     print("Completed %i of %i" %(progress, len(removeGroups)), end="\r")
    unduplicated = pd.concat(unduplicated).set_index(["ID", "G4", "CHR", "SPECIES", "START"])
    df_unduplicated = pd.concat([df_removed, unduplicated]).sort_index() #concatenate the unduplicated and the removed dataframe
    return df_unduplicated

def reduced_expr(input_str):
    '''Convert a list of numbers to a reduced expression'''
    if input_str != ".":
        # Split the input string by comma and space
        elements = input_str.split(",")
        # Initialize variables to store the final expression parts
        expression_parts = []
        # Iterate through each element
        current_range = None
        for element in elements:
            # Extract the prefix and number part using regular expression
            match = re.match(r'(\d+)([A-Z]+)', element)
            if match:
                number, prefix = match.groups()
                number = int(number)
                # Check if the current element continues the previous range
                if current_range is not None and number == current_range[1] + 1 and prefix == current_range[2]:
                    current_range = (current_range[0], number, prefix)
                else:
                    # If it's a new range, add the previous range if it exists
                    if current_range is not None:
                        if current_range[0] == current_range[1]:
                            expression_parts.append(f"{current_range[0]}{current_range[2]}")
                        else:
                            expression_parts.append(f"[{current_range[0]}-{current_range[1]}]{current_range[2]}")
                    # Start a new range
                    current_range = (number, number, prefix)
        # Add the last range to the expression
        if current_range is not None:
            if current_range[0] == current_range[1]:
                expression_parts.append(f"{current_range[0]}{current_range[2]}")
            else:
                expression_parts.append(f"[{current_range[0]}-{current_range[1]}]{current_range[2]}")
        # Construct the final expression by joining the parts with comma
        final_expression = ",".join(expression_parts)
    else:
        final_expression = "."
    return final_expression

def check_g4(alignment_sequence: str) -> bool:
    '''Check if a DNA sequence is a G4 motif'''
    # Remove gaps and convert to uppercase
    dna_sequence = alignment_sequence.upper().replace('-', '')

    # Define the regex pattern for G4 motifs
    pattern = r'((G{3,}[ATCG]{1,12}){3,}G{3,})|((C{3,}[ATCG]{1,12}){3,}C{3,})' # doesn't allow bulges
    #pattern = r'[G]{2,5}(?:[ATCG]+[G]{2,5})+[G]{2,5}|[C]{2,5}(?:[ATCG]+[C]{2,5})+[C]{2,5}' # allows bulges

    # Check if the sequence matches the pattern
    regex = re.match(pattern, dna_sequence)

    if regex != None and regex.span() == (0, len(dna_sequence)):
        # return True
        return 1
    else:
        # return False
        return 0

def mutation_profile(seq1, seq2):
    ''' Compare two DNA sequences and return the mutation profile'''
    seq1 = seq1.upper() #convert to uppercase
    seq2 = seq2.upper() #convert to uppercase
    mutations = { "G4?": 0, "SUB": ".", "INS": ".", "DEL": "."} #empty dictionary to store the mutations
    if seq1 == seq2: #if the sequences are the same
        mutations["G4?"] = 1 #store 1 in the G4? key
        return mutations #return 1 in the G4? key
    elif len(seq1) != len(seq2): #if the sequences are not the same length
        if seq2 == ".":
            mutations["DEL"] = "DL" #DL: Full deletion
            return mutations
        else:
            return False 
    else: 
        insertions = []
        deletions = []
        substitutions = []
        if check_g4(seq2) is True: #if the second sequence is a G4
            mutations["G4?"] = 1 #store 1 in the G4? key
        for index in range(len(seq1)): #iterate through the sequences
            base1 = seq1[index] #store the base in the first sequence
            base2 = seq2[index] #store the base in the second sequence
            if base1 != base2: #if the bases are not the same
                if base1 == "-": #if the base in the first sequence is a gap
                    insertions.append("%s%s" % (index+1, base2)) #store the index and the base in the second sequence
                elif base2 == "-": #if the base in the second sequence is a gap
                    deletions.append("%s%s" % (index+1, base1)) #store the index and the base in the first sequence
                else: #if the bases are not gaps
                    substitutions.append("%s%s%s" % (index+1, base1, base2))
        if len(substitutions) > 0:
            mutations["SUB"] = reduced_expr(",".join(substitutions))
        if len(insertions) > 0:
            mutations["INS"] = reduced_expr(",".join(insertions))
        if len(deletions) > 0:
            mutations["DEL"] = reduced_expr(",".join(deletions))
        return mutations #return the mutations dictionary
    
def adjust_group_length(group):
    '''If the group has multiple sequences with differing length even after 
    evaluate_redundant, align them to the longest sequence in the group.'''
    group["SEQ_LEN"] = group["SEQUENCE"].apply(len) #get the length of the sequence
    refer_row = group[group["TYPE"] == "fullMAF"]["SEQ_LEN"].idxmax() #get the row with the maximum length
    ref_seq = group.loc[refer_row]["SEQUENCE"] #get the sequence of the row with the maximum length

    # I DON'T UNDERSTAND why I have to put this here, 
    if type(ref_seq).__name__ != "str": 
        ref_seq = ref_seq.iloc[0]
    # this solves the problem, but df.loc behaves weirdely for HSA10
        
    for idx, row in group.iterrows(): #iterate through the rows
        if idx != refer_row: #all rows except the row with the refer index
            query_seq = row["SEQUENCE"] #get the sequence
            if query_seq == ".": #if the sequence is a gap
                aligned_seq = "." #store a gap
            else: #if the sequence is not a gap
                aligned_seq = align_sequences(ref_seq, query_seq)[1] #align the sequences
            group.loc[idx, "SEQUENCE"] = aligned_seq #store the aligned sequence
    return group #return the group

def annotate_mutations(dataframe, progress=False):
    '''Annotate the mutations in the dataframe'''
    data = dataframe.reset_index().set_index(['ID', 'G4', 'CHR', 'SPECIES', 'START']) #store the dataframe
    data[["G4?","SUB","INS","DEL"]] = [0,".",".","."] #create new columns and fill them with 0, ., ., .
    grouped_df = data.groupby(['ID', 'G4']) #group the dataframe by ID and G4
    for idx, group in grouped_df:
        if group[group["SEQUENCE"] != "."]["SEQUENCE"].apply(len).nunique() > 1: #if the sequences are not the same length
            group = adjust_group_length(group) #adjust the group sequence length through alignment
        # identify sequences with either fullMAF or partMAF type
        valids = group[(group['TYPE'] == 'fullMAF')|(group['TYPE'] == 'partMAF')]['SEQUENCE']
        # identify the first one as the reference
        reference = valids.iloc[0]
        reference_index = valids.index[0]
        # the reference sequence
        mutations_ref = {"G4?": 1, "SUB": "REF", "INS": "REF", "DEL": "REF"}
        data.loc[reference_index, mutations_ref.keys()] = mutations_ref.values()
        # identify non-reference sequences and annotate them
        non_reference = group[~group.index.isin([reference_index])]['SEQUENCE']
        # iterate through the non-reference sequences and identify mutation profile
        for i, seq in non_reference.items():
            mutations_else = mutation_profile(reference, seq)
            data.loc[i, mutations_else.keys()] = mutations_else.values()
        # if progress == True:
        #     print("Completed row %i of %i" %(data.index.get_loc(idx).stop,len(data)), end="\r")
    return data

def generateTree(dataframe):
    '''Generate a tree of connected components in the dataframe
    to identify the duplicate alignment blocks'''
    # Creating a graph
    G = nx.Graph()
    # Adding edges to the graph based on the ID lists
    for id_list in dataframe['Identifier_list']:
        if len(id_list) > 1:
            edges = [(id_list[i], id_list[j]) for i in range(len(id_list)) for j in range(i+1, len(id_list))]
            G.add_edges_from(edges)
    # Finding connected components in the graph
    connected_components = list(nx.connected_components(G))
    # Converting connected components to a DataFrame
    result = pd.DataFrame(connected_components)
    return result

def removeAlignmentDuplicates(df):
    '''Remove the alignment duplicates from the dataframe'''
    df_G4s = df.reset_index() #reset the index
    df_G4s = df_G4s[df_G4s["TYPE"] == "fullMAF"] #get the rows where the type is fullmaf
    df_G4s["Identifier"] = df_G4s["ID"] + "_" + df_G4s["G4"] #store the identifier as joined ID and G4
    df_G4s["Locator"] = df_G4s["SPECIES"].astype(str) + "_" + df_G4s["CHR"].astype(str) + "_" + df_G4s["START"].astype(int).astype(str) #store the locator as joined species, chr and start
    df_G4s = df_G4s[['Identifier', 'Locator']].reset_index(drop=True) #reset the index and store the Identifier and Locator
    df_G4s.sort_values(by="Locator", inplace=True) #sort the dataframe by Locator
    df_G4s = df_G4s.groupby("Locator")["Identifier"].apply(lambda x: ', '.join(x)).reset_index() #group by Locator and join the Identifier with a comma
    df_G4s['Identifier_list'] = df_G4s['Identifier'].str.split(',') #split the Identifier by comma
    #For identifying the blocks which are redundant and needs to be removed
    to_be_removed_blocks = []
    for idx, row in generateTree(df_G4s).iterrows(): #iterate through the rows of the generated tree
        row_dict = {str(id).strip(): len(df.loc[tuple(str(id).strip().split("_"))]) for id in row if id is not None} #store of the length of each of the blocks
        max_block = max(row_dict, key=row_dict.get) #get the block with the maximum length
        remove_block = [key for key in row_dict.keys() if key != max_block] #remove the blocks which are not the maximum block
        to_be_removed_blocks.extend(remove_block) #extend the remove block to the to_be_removed_blocks list
    duplicates_df = pd.DataFrame(to_be_removed_blocks, columns=["Identifier"]) #store the to_be_removed_blocks in a dataframe
    duplicates_df.sort_values(by="Identifier", inplace=True) #sort the dataframe by Identifier
    df = df.reset_index() #reset the index
    df["Identifier"] = df["ID"] + "_" + df["G4"] #store the Identifier as the joined ID and G4
    df.reset_index(inplace=True) #reset the index
    merged_df = df.merge(duplicates_df, on='Identifier') #merge the dataframes on Identifier from original and duplicates dataframe
    df.drop(columns=["index","Identifier"], inplace=True) #drop the generated non-required columns
    df.drop(merged_df["index"], inplace=True) #drop the duplicate rows
    return df

def removeRedundant(df, flank, mutInfo=False):
    df["G4?"] = df["SEQUENCE"].apply(lambda x: check_g4(x))

    redundantGroups = identifyRedundant(df, flank) #identify the redundant groups
    if mutInfo:
        if len(redundantGroups) == 0: #if there are no redundant groups
            df_ann_mutx = annotate_mutations(df) #annotate the mutations and remove the alignment duplicates
        else: #if there are redundant groups
            df_ann_mutx = annotate_mutations(processRedundant(df,redundantGroups))
        return df_ann_mutx
    else:
        if len(redundantGroups) == 0: #if there are no redundant groups
            return(df) #remove the alignment duplicates
        else: #if there are redundant groups
            return processRedundant(df,redundantGroups)
