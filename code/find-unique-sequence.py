from Bio import SeqIO
from collections import OrderedDict

# Input and output file names
input_file = "All-S-Protein-from-2020-2024.fasta"
output_file = "All-S-Protein-from-2020-2024-unique.fasta"

# Use OrderedDict to store unique sequences while preserving the original order
unique_sequences = OrderedDict()

# Read the FASTA file and remove duplicates
with open(input_file, "r") as handle:
    for record in SeqIO.parse(handle, "fasta"):
        # Convert the sequence to a string as the key to detect identical sequences
        seq_str = str(record.seq)
        # If the sequence hasn’t appeared before, add it to the dictionary
        if seq_str not in unique_sequences:
            unique_sequences[seq_str] = record

# Write the unique sequences to a new file
with open(output_file, "w") as output_handle:
    SeqIO.write(unique_sequences.values(), output_handle, "fasta")

# Print some statistical information
print(f"Total number of original sequences: {len(list(SeqIO.parse(input_file, 'fasta')))}")
print(f"Total number of unique sequences: {len(unique_sequences)}")
