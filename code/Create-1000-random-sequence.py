import random
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

# Wuhan-Hu-1 S蛋白氨基酸序列（简化版，实际需替换为完整序列）
# 完整序列可从NCBI获取（NC_045512.2，翻译后的S蛋白），这里用10个氨基酸演示
wuhan_seq = "MFVFLVLLPL"  # 实际应为1273个氨基酸

# 氨基酸列表（20种常见氨基酸）
amino_acids = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 
               'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']

# 定义S蛋白功能区域（基于氨基酸位置）
seq_length = len(wuhan_seq)  # 实际应为1273
ntd_range = (0, 305)         # NTD: 1-305
rbd_range = (319, 541)       # RBD: 319-541
other_range = [(0, 319), (541, seq_length)]  # 其他区域

# 生成一个随机变体的函数
def generate_variant(reference_seq, min_mutations=1, max_mutations=39):
    seq = list(reference_seq)  # 可变列表
    num_mutations = random.randint(min_mutations, max_mutations)  # 随机抽取1-39个突变
    
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
    if len(positions) < num_mutations:  # 如果重复导致数量不足，补齐
        while len(positions) < num_mutations:
            pos = random.randint(0, len(seq) - 1)
            if pos not in positions:
                positions.append(pos)
    
    # 生成氨基酸突变
    for pos in positions:
        original_aa = seq[pos]
        # 随机替换为不同氨基酸
        new_aa = random.choice([aa for aa in amino_acids if aa != original_aa])
        seq[pos] = new_aa
    
    return ''.join(seq)

# 生成1,000个独特变体
def generate_1000_variants(reference_seq):
    variants = set()  # 用set确保唯一性
    while len(variants) < 1000:
        new_variant = generate_variant(reference_seq)
        if new_variant != reference_seq and new_variant not in variants:
            variants.add(new_variant)
    
    return list(variants)

# 主程序
random.seed(42)  # 设置随机种子，保证可重复性
variants = generate_1000_variants(wuhan_seq)

# 输出前5个变体作为示例
print("Generated 1,000 unique variants. First 5 examples:")
for i, var in enumerate(variants[:5]):
    print(f"Variant {i+1}: {var}")

# 保存为FASTA文件
records = [SeqRecord(Seq(var), id=f"Variant_{i+1}", description="Random SARS-CoV-2 S variant") 
           for i, var in enumerate(variants)]
with open("random_variants.fasta", "w") as output_file:
    SeqIO.write(records, output_file, "fasta")
print("Saved 1,000 variants to 'random_variants.fasta'")

# 统计突变数（验证）
total_mutations = sum([sum(1 for a, b in zip(wuhan_seq, var) if a != b) for var in variants])
print(f"Total mutations across 1,000 variants: {total_mutations}, Average mutations per variant: {total_mutations/1000:.2f}")
