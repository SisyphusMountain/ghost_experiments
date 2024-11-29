import pandas as pd
import os
import re
import subprocess
import concurrent.futures
import shutil
from ete3 import Tree

def stats_df_fn(dir_path):
    """
        dir_path: str
            The path to the directory containing the reconciliation files    
        
    """
    # Initialize an empty DataFrame
    df = pd.DataFrame(columns=['Gene_Index', 'Duplications', 'Transfers', 'Losses', 
                            'Total_Duplications', 'Total_Transfers', 'Total_Losses', 'Total_Speciations'])

    # Compile the regex pattern for "..._uml" extraction
    pattern = re.compile(r"reconciliation_(\d+)_uml")

    # Loop over each file in the directory
    for filename in os.listdir(dir_path):
        # Only open files that end with _uml
        if filename.endswith("_uml"):
            # Extract the gene index from the filename
            gene_index = re.findall(pattern, filename)[0]

            with open(os.path.join(dir_path, filename), 'r') as file:
                lines = file.readlines()

                for i, line in enumerate(lines):
                    # Check if this line starts with 'rate of'
                    if line.startswith('rate of'):
                        # The numbers are on the next line, split by tabs
                        numbers = lines[i + 1].split()
                        dup, trans, loss = map(float, numbers[1:4])

                    # Check if this line starts with '# of' and ends with 'Speciations'
                    elif line.startswith('# of') and line.strip().endswith('Speciations'):
                        # The numbers are on the next line, split by tabs
                        total_numbers = lines[i + 1].split()
                        total_dup, total_trans, total_loss, total_spec = map(float, total_numbers[1:5])

                # Add these numbers to the DataFrame
                new_data = pd.DataFrame({
                    'Gene_Index': [gene_index],
                    'Duplications': [dup],
                    'Transfers': [trans],
                    'Losses': [loss],
                    'Total_Duplications': [total_dup],
                    'Total_Transfers': [total_trans],
                    'Total_Losses': [total_loss],
                    'Total_Speciations': [total_spec]
                })
                df = pd.concat([df, new_data], ignore_index=True)
    return df

def transfers_df_fn(base_folder_path):
    """
    base_folder_path: str
        The folder containing all the .uTs reconciliation files

    """
    temp_dfs = []  # Initialize a list to hold dataframes temporarily
    pattern = re.compile(r"reconciliation_(\d+)_uTs")
    
    for ind, filename in enumerate(os.listdir(base_folder_path)):
        if filename.endswith("_uTs"):
            gene_index_match = re.findall(pattern, filename)
            if gene_index_match:
                gene_index = gene_index_match[0]
                file_path = os.path.join(base_folder_path, filename)
                corrected_file_path = f"{file_path}.corrected"
                
                # Use sed to create a corrected copy of the file
                sed_command = ['sed', 's/^[ \t]*//', file_path]
                with open(corrected_file_path, 'w') as corrected_file:
                    subprocess.run(sed_command, stdout=corrected_file)
                
                try:
                    # Read the corrected file into a DataFrame
                    temp_df = pd.read_csv(corrected_file_path, sep='\t', comment='#', header=None, names=['Donor', 'Receiver', 'Frequency'], skipinitialspace=True)
                    
                    # Add the gene index to the dataframe
                    temp_df['Gene_Index'] = gene_index
                    
                    # Collect the dataframe
                    temp_dfs.append(temp_df)
                except pd.errors.ParserError as e:
                    print(f"Error reading {corrected_file_path}: {e}")
                finally:
                    # Optionally, remove the corrected file to clean up
                    os.remove(corrected_file_path)
    
    # Concatenate all collected dataframes once
    df = pd.concat(temp_dfs, ignore_index=True) if temp_dfs else pd.DataFrame()
    return df

def rename_internal_nodes(newick_path, output_path):
    tree = Tree(newick_path, format=1)
    count = 0
    for node in tree.traverse():
        if not node.is_leaf():
            node.name = str(count)
            count += 1
    tree.write(outfile=output_path, format=1, format_root_node=True)

def process_tree_directory(tree_dir):
    """
    tree_dir: str
        The directory containing all the reconciliation files
    """
    # Initialize the output directory
    output_dir = os.path.join(tree_dir, "processed_data_for_pytorch")
    
    # Remove the output directory if it exists
    if os.path.exists(output_dir):
        shutil.rmtree(output_dir)
    
    # Create the output directory
    os.makedirs(output_dir, exist_ok=True)
    
    # Get the stats dataframe
    stats_df = stats_df_fn(tree_dir + "/reconciliations")
    stats_df.to_csv(os.path.join(output_dir, "stats.csv"), index=False)
    
    # Get the transfers dataframe
    transfers_df = transfers_df_fn(tree_dir + "/reconciliations")
    transfers_df.to_csv(os.path.join(output_dir, "transfers.csv"), index=False)

    complete_tree_path = os.path.join(tree_dir, "complete_species_tree.nwk")
    shutil.copy(complete_tree_path, os.path.join(output_dir, "complete_species_tree.nwk"))

    rename_internal_nodes(complete_tree_path, os.path.join(output_dir, "complete_species_tree_internal_nodes_renamed.nwk"))

    sampled_tree_path = os.path.join(tree_dir, "sampled/sampled_trees/sampled_species_tree.nwk")
    shutil.copy(sampled_tree_path, os.path.join(output_dir, "sampled_species_tree.nwk"))
    
    newick_pattern = re.compile(r"S:\s*(\([^;]+;)", re.DOTALL)

    sampled_ale_tree_path = os.path.join(tree_dir, "reconciliations/reconciliation_0_uml")
    sampled_ale_tree_target_path = os.path.join(output_dir, "ale_sampled_species_tree.nwk")
    if os.path.exists(sampled_ale_tree_path):
        with open(sampled_ale_tree_path, 'r') as file:
            content = file.read()
        newick_match = re.search(newick_pattern, content)
        if newick_match:
            newick = newick_match.group(1)
            with open(sampled_ale_tree_target_path, 'w') as file:
                file.write(newick)
        else:
            print(f"Newick pattern not found in {sampled_ale_tree_path}")
    # TODO: remove the hardcoded path
    translation_command = ["/home/enzo/Documents/DeepGhosts/Translate_trees/target/release/Translate_trees",
                           tree_dir + "/processed_data_for_pytorch/sampled_species_tree.nwk",
                           tree_dir + "/processed_data_for_pytorch/ale_sampled_species_tree.nwk",
                           tree_dir + "/processed_data_for_pytorch"]
    
    result = subprocess.run(translation_command, capture_output=True, text=True)

    if result.returncode == 0:
        print(result.stdout)
    else:
        print("Error executing command:")

def main():
    num_threads = 32
    # Define the base directory
    base_dir = "/media/enzo/Stockage/Output"
    tree_directories = [os.path.join(base_dir, f"species_tree_{i}") for i in range(150)]
    with concurrent.futures.ThreadPoolExecutor(max_workers=num_threads) as executor:
        futures = [executor.submit(process_tree_directory, tree_dir) for tree_dir in tree_directories]
        for future in concurrent.futures.as_completed(futures):
            try:
                future.result()
            except Exception as exc:
                print(f"An error occurred: {exc}")

if __name__ == "__main__":
    main()
