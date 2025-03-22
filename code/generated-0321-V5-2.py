import pyvolve
from Bio import SeqIO

# 读取野生型S蛋白序列
fasta_file = "Wuhan-Hu-1.fasta"
wildtype_seq = str(next(SeqIO.parse(fasta_file, "fasta")).seq)
sequence_length = len(wildtype_seq)  # 假设为1273个氨基酸

# ----------------------------
# 参数修正（基于最新文献）
# ----------------------------
simulation_years = 4.375            # 模拟时间跨度（年）
total_mutations_observed = 25.3229   # 观测到的总突变数（4.375年内）
neutral_fraction = 0.3               # 整体中性比例设为30%

# 计算突变率与分支长度
mutation_rate_per_site_per_year = total_mutations_observed / simulation_years / sequence_length
neutral_mutation_rate = mutation_rate_per_site_per_year * neutral_fraction
# branch_length = neutral_mutation_rate * simulation_years  # ≈0.00604
branch_length = 0.003
# ----------------------------
# 进化树与模型设置
# ----------------------------
tree_string = f"(root:{branch_length});"
phylogeny = pyvolve.read_tree(tree=tree_string)
model = pyvolve.Model("WAG")  # 氨基酸替代模型
partition = pyvolve.Partition(models=model, root_sequence=wildtype_seq)

# ----------------------------
# 模拟进化生成突变体
# ----------------------------
evolver = pyvolve.Evolver(tree=phylogeny, partitions=partition)
mutant_sequences = []

for i in range(144248):
    evolver()
    sequences = evolver.get_sequences()
    mutant_seq = list(sequences.values())[0]
    mutant_sequences.append((f"mutant_{i+1}", mutant_seq))

# ----------------------------
# 输出与验证
# ----------------------------
output_file = "Random-144248-0321-North_America-0.003.fasta"
with open(output_file, "w") as f:
    for seq_id, seq in mutant_sequences:
        f.write(f">{seq_id}\n{seq}\n")

expected_mutations = branch_length * sequence_length  # ≈0.00604 * 1273 ≈7.7
print(
    f"已生成 {len(mutant_sequences)} 条中性进化突变体\n"
    f"预期平均突变数/序列: {expected_mutations:.1f}\n"
    f"中性比例: {neutral_fraction*100}%"
)