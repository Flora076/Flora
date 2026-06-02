import gzip
import os
import numpy as np
from Bio import SeqIO
from itertools import product, combinations

def get_kmer_frequencies(file_path, k=6):
    """Counts the frequencies of all possible k-mers in a FASTA/FASTQ or compressed .GZ file."""
    # Generate all possible k-mers to ensure vector lengths match perfectly
    all_kmers = [''.join(p) for p in product('ATCG', repeat=k)]
    kmer_dict = {kmer: 0 for kmer in all_kmers}

    total_kmers = 0
    file_name = os.path.basename(file_path).lower()
    print(f"Processing {os.path.basename(file_path)}...")

    # 1. Automatically detect file format
    if any(ext in file_name for ext in [".fastq", ".fq"]):
        file_format = "fastq"
    elif any(ext in file_name for ext in [".fasta", ".fa"]):
        file_format = "fasta"
    else:
        # Default fallback if extension is weird
        file_format = "fasta" 
    
    print(f"--> Detected format: {file_format.upper()}")

    # 2. Check if the file is gzipped
    is_gzipped = file_name.endswith(".gz")
    open_func = gzip.open(file_path, "rt") if is_gzipped else open(file_path, "r")

    # 3. Parse and count k-mers
    with open_func as handle:
        record_count = 0
        for record in SeqIO.parse(handle, file_format):
            record_count += 1
            seq = str(record.seq).upper()
            # Slide a window of size 'k' across the sequence
            for i in range(len(seq) - k + 1):
                kmer = seq[i:i+k]
                if kmer in kmer_dict:  # Ignores 'N' or ambiguous bases
                    kmer_dict[kmer] += 1
                    total_kmers += 1
                    
        print(f"--> Successfully parsed {record_count} reads. Found {total_kmers} total valid 6-mers.")

    # Convert counts to relative frequencies (proportions)
    frequency_vector = np.array([kmer_dict[kmer] for kmer in all_kmers], dtype=float)
    if total_kmers > 0:
        frequency_vector /= total_kmers

    return frequency_vector

def cosine_similarity(vec1, vec2):
    """Calculates the cosine similarity between two vectors."""
    dot_product = np.dot(vec1, vec2)
    norm_a = np.linalg.norm(vec1)
    norm_b = np.linalg.norm(vec2)
    if norm_a == 0 or norm_b == 0:
        return 0.0
    return dot_product / (norm_a * norm_b)

def compare_four_samples(file_list):
    """Compares all samples pairwise."""
    vectors = {}

    # 1. Profile all files
    for filepath in file_list:
        # Go up two directory levels to grab the unique SRR folder name
        srr_id = os.path.basename(os.path.dirname(os.path.dirname(filepath)))
        
        # Fallback to standard filename if it's not structured in an SRR folder
        if not srr_id or srr_id == "..":
            name = os.path.basename(filepath)
        else:
            name = srr_id

        vectors[name] = get_kmer_frequencies(filepath, k=6)

    # 2. Compute pairwise similarity matrix
    print(f"\n{'Sample 1':<25} | {'Sample 2':<25} | {'Cosine Similarity':<20}")
    print("-" * 76)

    for file1, file2 in combinations(vectors.keys(), 2):
        similarity = cosine_similarity(vectors[file1], vectors[file2])
        print(f"{file1:<25} | {file2:<25} | {similarity:.4f}")

# Example Usage
if __name__ == "__main__":

    # 1. Set your shared directory path here (use / or double \\ for Windows paths)
    DATA_DIR = "/omics/groups/OE0436/data/linmq/Datasets/lung_adenocarcinoma/"

    # 2. List your unique file names here
    fasta_files = [
"SRR26854886/temp/tel_reads.fa.gz",
"SRR26842326/temp/tel_reads.fa.gz",
"SRR26842325/temp/tel_reads.fa.gz",
"SRR26842324/temp/tel_reads.fa.gz",

]

    # Combine the directory path with the file names dynamically
    full_paths = [os.path.join(DATA_DIR, filename) for filename in fasta_files]

    # Filter to ensure the files actually exist before running
    existing_samples = [f for f in full_paths if os.path.exists(f)]

    # Alert the user if any files were missed
    for p in full_paths:
        if p not in existing_samples:
            print(f"Warning: File not found and skipped -> {p}")

    if len(existing_samples) < 2:
        print("\nError: Not enough valid FASTA files found to run a comparison.")
    else:
        compare_four_samples(existing_samples)
