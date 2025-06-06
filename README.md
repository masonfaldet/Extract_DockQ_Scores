# README

## Overview

This script processes Alphafold server job outputs for a given protein, extracts key metrics for each model, and writes the results into a tab-delimited text file. It also plots the relationship between the DockQ score and the iptm score for quick visualization.

The main output file is named:  
**`{protein_name}_scores.txt`**  
where `{protein_name}` is the protein’s identifier (e.g., `8wtc`).

Each row in the file corresponds to a model output and includes the following columns:

- **seed**: The seed used for the job (padded to 7 digits with leading zeros)
- **DockQ**: The DockQ score computed from interface comparisons
- **F1**: The F1 score extracted from the interface metrics
- **iRMSD**: The interface RMSD (iRMSD)
- **LRMSD**: The ligand RMSD (LRMSD)
- **fnat**: The fraction of native contacts
- **iptm**: The confidence iptm score
- **ptm**: The confidence ptm score
- **ranking**: The ranking score from the summary confidence metrics

## Steps to Use This Script

1. **Prepare the Input Files:**
   - Copy the FASTA sequence for the protein to the Alphafold server.
   - Run several jobs with various seed values.
   - Download all output files from the Alphafold server to your local machine.
   - Alphafold outputs will be stored in a folder (typically named something like `folds_{date_time}`). Rename this folder to your protein's identifier (e.g., rename it to `8wtc`).

2. **Place the Native Structure:**
   - Download the reference CIF file from the Protein Data Bank (PDB) that corresponds to your target protein.
   - Place this file into your renamed folder (`8wtc`) and ensure it is named exactly as `{protein_name}.cif` (e.g., `8wtc.cif`).

3. **Set Up the Working Directory:**
   - Move the script along with your protein folder (e.g., `8wtc/`) into your working directory.

4. **Configure the Script:**
   - Open the script file and verify that the variable `protein_name` (defined near the top) matches the folder and file naming convention you are using (currently set to `"8wtc"`).

5. **Run the Script:**
   - Execute the script.
   - The script will process all subdirectories, compute the desired metrics, and write the results into a text file named `{protein_name}_scores.txt` in your working directory.

6. **Examine the Output:**
   - Open the output text file to review the tabulated scores.
   - Optionally, view the plotted scatter graph (DockQ vs. iptm) generated by the script.

