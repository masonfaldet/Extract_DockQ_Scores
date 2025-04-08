import os
import json
from DockQ.DockQ import load_PDB, run_on_all_native_interfaces
import matplotlib.pyplot as plt


"""
This script loops through all models outputed from all alphafold server jobs and extracts the metrics of interest.
Then a text file is created called "{protien_name}_scores.txt", the seed information and all metrics are written 
into this txt file in a tabular format. The columns of the table correspond to metrics and the rows correspond to 
a given model. 

This script should be flexible enough to accomidate any job name you may have used, however it assumes
you did not go in and change the names of the output files. I assume the output files are in their default format 
e.g. 
     "fold_{job_name}_job_request.json"
     "fold_{job_name}_model_{model_id}.cif" 
     "fold_{job_name}_summary_confidences_{model_id}.json" 

     
Steps to make this work: 

    1. Copy the FASTFA sequence infor for {protien_name} into alphafold server.

    2. Run several jobs with various seeds

    3. Select all output files in the alphafold server and download to local machine

    4. Alphafold server automatically stores all the jobs files you've downloaded into one folder
       called "folds_{date and time of export}". Rename this folder as "{protien_name}", e.g. "8wtc".

    5. Add the reference cif file from the PDB to the folder you've just renamed as {protien_name}.
       Ensure the reference structure is named "{protien_name}.cif", e.g. "8wtc.cif". The directory
       you've just renamed in Step 4 should contain
            -"{protien_name}.cif"
            - Several sub-directories (one for each job you've run) each sub-directory contains
              the json and cif files produced during the corresponding job.
    
    6. Move this script and the folder "{protien_name}" from step 4 to your working directory. 
    
    7. Run this script, make sure the string assigned to protien_name in the first line of this script
       matches the string you've been using for {protien_name}. Currently protien_name is set to "8wtc"
       
    8. Locate the text file "{protien_name}_scores.txt" in your working directory and you're done.
"""

protein_name = "8wtc"
main_dir = protein_name
native_structure_file = os.path.join(main_dir, f"{protein_name}.cif")

# Load the native structure once, since it's common for all comparisons.
native = load_PDB(native_structure_file)

# List to store rows. The first row is the header.
all_rows = []
all_rows.append(["seed", "DockQ", "F1", "iRMSD", "LRMSD", "fnat", "iptm", "ptm", "ranking"])

for subdir in os.listdir(main_dir):
    subdir_path = os.path.join(main_dir, subdir)
    if os.path.isdir(subdir_path):
        for filename in os.listdir(subdir_path):
            if "model" in filename and filename.endswith(".cif"):

                # Find the job JSON file containing seed information.
                job_filename = None
                for candidate in os.listdir(subdir_path):
                    if "job_request" in candidate and candidate.endswith(".json"):
                        job_filename = candidate
                        break

                if job_filename is None:
                    print(f"Warning: No job JSON file found in directory {subdir_path}. Skipping.")
                    continue

                job_json_path = os.path.join(subdir_path, job_filename)
                if not os.path.exists(job_json_path):
                    print(f"Warning: JSON file {job_json_path} not found in {subdir_path}. Skipping.")
                    continue

                with open(job_json_path, "r") as f:
                    job_data = json.load(f)
                  
                seed = job_data[0]["modelSeeds"][0]

                # Extract the id string from the filename:
                # id_str is the substring that ends at the 4th-to-last character
                # and goes back until the first underscore encountered.
                id_str = filename[filename.rfind("_", 0, len(filename) - 4) + 1: len(filename) - 4]
                try:
                    model_id = int(id_str)
                except ValueError:
                    print(f"Skipping file with unexpected ID format: {filename}")
                    continue

                # Full path to the model CIF file.
                model_path = os.path.join(subdir_path, filename)
                
                # Find the corresponding summary_confidences JSON file.
                json_filename = None
                for candidate in os.listdir(subdir_path):
                    if "summary_confidences" in candidate and id_str in candidate and candidate.endswith(".json"):
                        json_filename = candidate
                        break
                      
                if json_filename is None:
                    print(f"Warning: No summary_confidences JSON file found for id '{id_str}' in {subdir_path}. Skipping.")
                    continue
                  
                json_path = os.path.join(subdir_path, json_filename)
                if not os.path.exists(json_path):
                    print(f"Warning: JSON file {json_path} not found for model file {model_path}. Skipping.")
                    continue

                # Compute DockQ and extract interface metrics.
                model = load_PDB(model_path)
                chain_map = {"A": "A", "B": "B"}
                result, dockq_score = run_on_all_native_interfaces(model, native, chain_map=chain_map)
                f1 = result["AB"]["F1"]
                irmsd = result["AB"]["iRMSD"]
                lrmsd = result["AB"]["LRMSD"]
                fnat = result["AB"]["fnat"]

                # Load the confidence JSON file and extract the confidence scores.
                with open(json_path, "r") as f:
                    confidence_data = json.load(f)
                iptm_score = confidence_data.get("iptm")
                if iptm_score is None:
                    print(f"Warning: 'iptm' score not found in {json_path}. Skipping this model.")
                    continue

                ptm_score = confidence_data.get("ptm")
                if ptm_score is None:
                    print(f"Warning: 'ptm' score not found in {json_path}. Skipping this model.")
                    continue

                ranking_score = confidence_data.get("ranking_score")
                if ranking_score is None:
                    print(f"Warning: 'ranking_score' not found in {json_path}. Skipping this model.")
                    continue

                # Append the row for this model.
                all_rows.append([seed, dockq_score, f1, irmsd, lrmsd, fnat, iptm_score, ptm_score, ranking_score])


widths = [11, 11, 11, 11, 11, 11, 11, 11, 11]

output_filename = f"{protein_name}_scores.txt"
with open(output_filename, "w") as out_file:
    # Format and write the header row: center the header text within the column width.
    formatted_headers = [header.center(width) for header, width in zip(all_rows[0], widths)]
    out_file.write("\t".join(formatted_headers) + "\n")
    
    # Process each model's row.
    for row in all_rows[1:]:
        formatted_row = []
        for i, (item, width) in enumerate(zip(row, widths)):
            # For the seed column, pad the number with zeros (7 digits) and then right-align.
            if i == 0:
                formatted_item = str(item).zfill(7).rjust(width)
            # For float values, format with 7 decimal places and right-align.
            elif isinstance(item, float):
                formatted_item = f"{item:.7f}".rjust(width)
            else:
                formatted_item = str(item).rjust(width)
            formatted_row.append(formatted_item)
        out_file.write("\t".join(formatted_row) + "\n")

print(f"Scores written to {output_filename}")

# Optional: Plot the scores.
dockq_scores = [row[1] for row in all_rows[1:]]
iptm_scores = [row[6] for row in all_rows[1:]]
plt.scatter(dockq_scores, iptm_scores)
plt.xlabel("DockQ")
plt.ylabel("iptm")
plt.show()
