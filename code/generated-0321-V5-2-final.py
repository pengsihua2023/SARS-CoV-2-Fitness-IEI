import pyvolve
from Bio import SeqIO

# Read the wild-type S protein sequence
fasta_file = "Wuhan-Hu-1.fasta"
wildtype_seq = str(next(SeqIO.parse(fasta_file, "fasta")).seq)
sequence_length = len(wildtype_seq)  # Assumed to be 1273 amino acids

# ----------------------------
# Parameter adjustment (based on the latest literature)
# ----------------------------
simulation_years = 4.375            # Simulation time span (years)
total_mutations_observed = 25.3229   # Total observed mutations (over 4.375 years)
neutral_fraction = 0.3               # Overall neutral fraction set to 30%

# Calculate mutation rate and branch length
mutation_rate_per_site_per_year = total_mutations_observed / simulation_years / sequence_length
neutral_mutation_rate = mutation_rate_per_site_per_year * neutral_fraction
# branch_length = neutral_mutation_rate * simulation_years  # ≈0.00604
branch_length = 0.003
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
output_file = "Random-144248-0321-North_America-0.003.fasta"
with open(output_file, "w") as f:
    for seq_id, seq in mutant_sequences:
        f.write(f">{seq_id}\n{seq}\n")

expected_mutations = branch_length * sequence_length  # ≈0.00604 * 1273 ≈7.7
print(
    f"Generated {len(mutant_sequences)} neutral evolution mutants\n"
    f"Expected average mutations per sequence: {expected_mutations:.1f}\n"
    f"Neutral fraction: {neutral_fraction*100}%"
)