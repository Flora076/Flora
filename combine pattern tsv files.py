import pandas as pd

# 1. Define the specific paths to your individual TSV files
file_paths = [
"/omics/groups/OE0436/data/linmq/Datasets/pattern/fibroblast/fibroblast_combined_motifs.tsv", 
"/omics/groups/OE0436/data/linmq/Datasets/pattern/fibrosarcoma/fibrosarcoma_combined_motifs.tsv",
"/omics/groups/OE0436/data/linmq/Datasets/pattern/iPSC/iPSC_combined_motifs.tsv",
"/omics/groups/OE0436/data/linmq/Datasets/pattern/HG002/HG002_combined_motifs.tsv"
]

# 2. Create an empty list to store the data temporarily
data_frames = []

# 3. Loop through each path, read the file, and add it to the list
for path in file_paths:
    try:
        # sep='\t' tells pandas it's a TSV, not a standard comma CSV
        df = pd.read_csv(path, sep='\t') 
        data_frames.append(df)
        print(f"Successfully read: {path}")
    except FileNotFoundError:
        print(f"Warning: Could not find the file at {path}")

# 4. Combine all the individual DataFrames into one big DataFrame
# ignore_index=True ensures your new file has a clean, continuous row count
combined_df = pd.concat(data_frames, ignore_index=True)

# 5. Save the combined data into a new TSV file
output_filename = "ALT-motifs.tsv"
combined_df.to_csv(output_filename, sep='\t', index=False)

print(f"\nAll files successfully combined into {output_filename}!")
