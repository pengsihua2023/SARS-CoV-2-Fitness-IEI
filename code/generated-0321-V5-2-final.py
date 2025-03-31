import pyvolve
from Bio import SeqIO

# Read the wild-type S protein sequence
fasta_file = "Wuhan-Hu-1.fasta"
wildtype_seq = str(next(SeqIO.parse(fasta_file, "fasta")).seq)
sequence_length = len(wildtype_seq)  # Assumed to be 1273 amino acids

# ----------------------------
# Parameter adjustment (based on the latest literature)
# ----------------------------
total_mutations_observed per year = 25.3229   # Total observed mutations (over 4.375 years)
neutral_fraction = 0.2               # Overall neutral fraction set to 20%

# Calculate mutation rate and branch length
mutation_rate_per_site_per_year = total_mutations_observed / simulation_years / sequence_length
neutral_mutation_rate = mutation_rate_per_site_per_year * neutral_fraction
# branch_length = neutral_mutation_rate * simulation_years  
branch_length = 0.0178
# ----------------------------
# Phylogenetic tree and model setup
# ----------------------------
tree_string = f"(root:{branch_length});"
phylogeny = pyvolve.read_tree(tree=tree_string)
model = pyvolve.Model("WAG")  # Amino acid substitution model
partition = pyvolve.Partition(models=model, root_sequence=wildtype_seq)

# ----------------------------
# Simulate evolution to generate mutants
# ----------------------------
evolver = pyvolve.Evolver(tree=phylogeny, partitions=partition)
mutant_sequences = []

for i in range(144248):
    evolver()
    sequences = evolver.get_sequences()
    mutant_seq = list(sequences.values())[0]
    mutant_sequences.append((f"mutant_{i+1}", mutant_seq))

# ----------------------------
# Output and validation
# ----------------------------
output_file = "Random-144248-0321-North_America.fasta"
with open(output_file, "w") as f:
    for seq_id, seq in mutant_sequences:
        f.write(f">{seq_id}\n{seq}\n")

expected_mutations = branch_length * sequence_length  
print(
    f"Generated {len(mutant_sequences)} neutral evolution mutants\n"
    f"Expected average mutations per sequence: {expected_mutations:.1f}\n"
    f"Neutral fraction: {neutral_fraction*100}%"
)
