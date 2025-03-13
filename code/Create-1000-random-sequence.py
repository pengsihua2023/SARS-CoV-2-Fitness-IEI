import random
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

# 从FASTA文件读取Wuhan-Hu-1 S蛋白序列
fasta_file = "Wuhan-Hu-1.fasta"  # 请确保文件路径正确
for record in SeqIO.parse(fasta_file, "fasta"):
    wuhan_seq = str(record.seq)  # 提取序列
    break  # 只读取第一个序列

# 验证序列长度
print(f"Loaded Wuhan-Hu-1 sequence length: {len(wuhan_seq)} amino acids")
if len(wuhan_seq) < 1000:
    raise ValueError("Sequence length is less than 1000 bp, please check the input file.")

# 氨基酸列表
amino_acids = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 
               'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']

# 定义区域（基于1273 aa完整序列）
seq_length = len(wuhan_seq)  # 应为1273
ntd_range = (0, 305)         # NTD: 1-305
rbd_range = (319, 541)       # RBD: 319-541
other_range = [(0, 319), (541, seq_length)]  # 其他区域

# 生成一个随机变体
def generate_variant(reference_seq, min_mutations=1, max_mutations=39):
    seq = list(reference_seq)
    num_mutations = random.randint(min_mutations, max_mutations)
    
    # 加权选择突变位置（RBD 40%, NTD 20%, 其他 40%）
    positions = []
    for _ in range(num_mutations):
        r = random.random()
        if r < 0.4:  # 40% RBD
            pos = random.randint(rbd_range[0], rbd_range[1] - 1)
        elif r < 0.6:  # 20% NTD
            pos = random.randint(ntd_range[0], ntd_range[1] - 1)
        else:  # 40% 其他
            if random.random() < 0.5:
                pos = random.randint(other_range[0][0], other_range[0][1] - 1)
            else:
                pos = random.randint(other_range[1][0], other_range[1][1] - 1)
        positions.append(pos)
    
    # 确保位置唯一
    positions = list(set(positions))
    if len(positions) < num_mutations:
        while len(positions) < num_mutations:
            pos = random.randint(0, len(seq) - 1)
            if pos not in positions:
                positions.append(pos)
    
    # 生成氨基酸突变
    for pos in positions:
        original_aa = seq[pos]
        new_aa = random.choice([aa for aa in amino_acids if aa != original_aa])
        seq[pos] = new_aa
    
    return ''.join(seq)

# 生成1,000个独特变体
def generate_1000_variants(reference_seq):
    variants = set()
    while len(variants) < 1000:
        new_variant = generate_variant(reference_seq)
        if new_variant != reference_seq and new_variant not in variants:
            variants.add(new_variant)
    return list(variants)

# 主程序
random.seed(42)  # 固定随机种子，保证可重复性
variants = generate_1000_variants(wuhan_seq)

# 输出前5个变体（显示前20个氨基酸）
print("Generated 1,000 unique variants. First 5 examples (first 20 aa):")
for i, var in enumerate(variants[:5]):
    print(f"Variant {i+1}: {var[:20]}...")

# 保存为FASTA文件
records = [SeqRecord(Seq(var), id=f"Variant_{i+1}", description="Random SARS-CoV-2 S variant") 
           for i, var in enumerate(variants)]
with open("random_variants-1000.fasta", "w") as output_file:
    SeqIO.write(records, output_file, "fasta")
print("Saved 1,000 variants to 'random_variants-1000.fasta'")

# 统计突变数
total_mutations = sum([sum(1 for a, b in zip(wuhan_seq, var) if a != b) for var in variants])
print(f"Total mutations: {total_mutations}, Average mutations per variant: {total_mutations/1000:.2f}")
