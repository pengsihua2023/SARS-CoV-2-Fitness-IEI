from Bio import SeqIO
from collections import OrderedDict

# 输入和输出文件名
input_file = "All-S-Protein-from-2020-2024.fasta"
output_file = "All-S-Protein-from-2020-2024-unique.fasta"

# 使用OrderedDict来存储唯一的序列，保持原始顺序
unique_sequences = OrderedDict()

# 读取FASTA文件并去重
with open(input_file, "r") as handle:
    for record in SeqIO.parse(handle, "fasta"):
        # 将序列转换为字符串作为键，这样可以检测完全相同的序列
        seq_str = str(record.seq)
        # 如果序列还没出现过，就添加到字典中
        if seq_str not in unique_sequences:
            unique_sequences[seq_str] = record

# 将独特的序列写入新文件
with open(output_file, "w") as output_handle:
    SeqIO.write(unique_sequences.values(), output_handle, "fasta")

# 打印一些统计信息
print(f"原始序列总数: {len(list(SeqIO.parse(input_file, 'fasta')))}")
print(f"独特序列总数: {len(unique_sequences)}")
