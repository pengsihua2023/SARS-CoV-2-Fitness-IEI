from collections import defaultdict

# Initialize a dictionary to count the number of lineages
lineage_count = defaultdict(int)

# Read the fasta file and parse the header of each sequence
total_sequences = 0
with open("sihua-3.fasta", "r") as fasta_file:
    for line in fasta_file:
        if line.startswith(">"):
            # Extract lineage information from the header
            lineage = line.strip().split("|")[-1]
            # Update the count of the lineage in the dictionary
            lineage_count[lineage] += 1
            total_sequences += 1

# Write the counting results to Lineage_count.txt
with open("Lineage_count.txt", "w") as output_file:
    for lineage, count in lineage_count.items():
        output_file.write(f"{lineage}: {count}\n")

# Calculate the name and percentage of the lineage with the highest count
dominant_lineage = max(lineage_count, key=lineage_count.get)
dominant_count = lineage_count[dominant_lineage]
dominant_percentage = (dominant_count / total_sequences) * 100

# Write the lineage statistics to Lineage_statistics.txt
with open("Lineage_statistics.txt", "w") as output_file:
    output_file.write(f"Total sequences: {total_sequences}\n")
    output_file.write(f"Total unique lineages: {len(lineage_count)}\n")
    output_file.write(f"Dominant lineage: {dominant_lineage}\n")
    output_file.write(f"Dominant lineage count: {dominant_count}\n")
    output_file.write(f"Dominant lineage percentage: {dominant_percentage:.2f}%\n")
