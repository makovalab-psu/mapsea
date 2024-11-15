# MAP-SEA
### Mapping And Prediction of Shared Elements in Alignments
Scripts for mapping short elements stored in .BED files (as intervals, e.g. G-quadruplexes) to .MAF (alignment) files.

## Table of Contents

-   [Overview](#overview)
-   [Features](#features)
-   [Usage](#usage)
- [Notes](#notes)
-   [Citation](#citation)

## Overview

This repository provides two Python scripts for mapping and analyzing short genomic elements (such as G-quadruplexes) in multi-species alignments. The tools are designed to work with the MAF (Multiple Alignment Format) file format and accept multiple BED files containing the sequence intervals to map. The first script maps these intervals to their corresponding locations in a MAF alignment, allowing users to observe where each element appears across different sequences in the alignment.  

The second script refines the initial mapping by accounting for slight positional variations (flanks or “wobbles”) in similar elements across species. This script helps remove redundant entries, merge similar elements, and optimize the data for downstream analyses. Together, these tools support comparative genomics studies focused on the evolution and distribution of short genomic elements across species.

## Features

-   **Alignment Mapping**: Map multiple BED intervals to positions within a single MAF file.
-   **Customizable Parameters**: Configure analysis settings to suit different species and overlaps.
-   **Redundancy Filtering**: Combine elements that meet defined similarity thresholds.
-   **Mutation Analysis**: Track mutations in shared elements across species.
-   **Optimized Performance**: Designed for large genomic datasets with support for multicore parallelization.


## Usage

### Command-line Usage

Below are example commands to demonstrate the main functionality of this repository:

-   **Alignment Mapping**:
    Use case:  
    `mapsea.py [-h] -m MAFINPUT -b BEDFOLDER -o OUTFILE -t TEMPDIR -r INTERSECTRATIO -d SPECIESDICT [-f MAPFILE] [-c CORES]`
    
    ```
    options:
    -h, --help  show this help message and exit
	-m MAFINPUT, --mafInput MAFINPUT
				Input gzipped maf file path (e.g., 'demo_hs1.chr1.2023v2.processed.maf.gz')
	-b BEDFOLDER, --bedFolder BEDFOLDER
				Folder name to access the bed files (without / at the end); the folder should have a tree structure, 
				with each species having their bed files under a directory separated by chromosomes (e.g., 'pqsfinderOutput')
	-o OUTFILE, --outFile OUTFILE
				Output file path (e.g., 'hs1.chr1.2023v2.processed.analysed.quadron.dat')
	-t TEMPDIR, --tempDir TEMPDIR
				Temporary directory path (e.g., 'tmp/hsa1')
	-r INTERSECTRATIO, --intersectRatio INTERSECTRATIO
				Ratio of what fraction of elements overlap with alignment block (e.g., 1.0)
	-d SPECIESDICT, --speciesDict SPECIESDICT
				Species dictionary file path, a json file only (e.g. 'speciesDict.json')
	-f MAPFILE, --mapFile MAPFILE
				Homolog chromosome map file path
	-c CORES, --cores CORES
				The number of cores available for this job (e.g. 10)
    ```
    
-   **Merging Similar G4s**:
(**only when output is generated using `mapsea.py -f 1.0`**)
    Use case:
    `refiner.py [-h] -d DATINPUT -f FLANK [-c CORES] [-o OUTFILE] [-m]`
   
    ```
    options:
    -h, --help            show this help message and exit
    -d DATINPUT, --datInput DATINPUT
			              Input dat file output by mapsea (e.g., 'chr22_bonobo_vs_chr22_borang.dat')
	-f FLANK, --flank FLANK
                          The flank allowed for the element to be termed 'shared' (e.g. 3)
	-c CORES, --cores CORES
                          The number of cores available for this job (e.g. 10)
	-o OUTFILE, --outFile OUTFILE
                          Output file path (e.g., 'hs1.chr1.2023v2.processed.analysed.quadron.df')
	-m, --mutInfo         Include mutation information in the output file
    ``` 

### Example Use Cases

1.  **Alignment Mapping**: 
 
    `python3 mapsea.py -m ../test/maffiles/chr22_bonobo_vs_chr22_borang.maf.gz -b ../test/bedfiles -o ../output/chr22_bonobo_vs_chr22_borang.dat -r 1.0 -t ../test/tmp -d ../test/refer_dict.json` 
    
    
2.  **Merging Similar G4s**: 
	 `python3 refiner.py -d ../output/chr22_bonobo_vs_chr22_borang.dat -f 3 -o ../output/chr22_bonobo_vs_chr22_borang.rmredundant.df -m`
    

## Notes

Some things to keep in mind while running these scripts:
1.   **Script Dependency**: Ensure that `mapsea.py` is run with the `-f 1.0` flag if you intend to use `refiner.py`, as this is required for the latter's compatibility.
    
2. **BED File Naming and Format**: BED files must be named in the format `chrN.bed` (e.g., `chr1.bed`, `chr2.bed`). The BED files should have the following columns: 

	| Column | Description | 
	|--------|-------------| 
	| col 1 | chrom | 
	| col 2 | start | 
	| col 3 | end | 
	| col 4 | score | 
	| col 5 | length | 
	| col 6 | strand |

3. **MAF File Entry Format**: Each entry in the MAF file should follow the standard format `speciesName.chrN` (e.g., `human.chr1`).

## Citation

If you use this tool in your research, please cite the following paper:

> Mohanty, S. K., Chiaromonte, F., & Makova, K. (2024). "[Evolutionary Dynamics of G-Quadruplexes in Human and Other Great Ape Telomere-to-Telomere Genomes](https://www.biorxiv.org/content/10.1101/2024.11.05.621973v1). *bioRxiv*, 2024-11. `doi: https://doi.org/10.1101/2024.11.05.621973`