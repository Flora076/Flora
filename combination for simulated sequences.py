# Define the repeat unit
repeat_unit = "TTAGGG"

# Create the repeated flanking regions
flank = repeat_unit * 5

# Define your fixed set of 6 sequences
fixed_sequences = [
    "T",
    "TT",
    "TTA",
    "TTAG",
    "TTAGG",
    "TTAGGG"
]

# Define your specific sequence
specific_sequence = "G"

# Count existing lines starting with >
try:
    with open("sequences.fasta", "r") as f:
        existing_ids = [line for line in f if line.startswith(">")]
        start_index = len(existing_ids) + 1
except FileNotFoundError:
    start_index = 1

# Generate and print in FASTA format
with open("sequences.fasta", "a") as f:
  for i, seq in enumerate(fixed_sequences, start=start_index):
    full_sequence = flank + seq + specific_sequence + flank
    f.write(f">seq_{i}\n")
    f.write(full_sequence + "\n")
