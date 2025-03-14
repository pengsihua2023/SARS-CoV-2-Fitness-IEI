
import pandas as pd
from scipy import stats

# 读取CSV文件
df_real = pd.read_csv("Real-sequence.csv")
df_random = pd.read_csv("Random-sequence.csv")

# 进行KS正态性检验
# 标准化数据
real_mean = df_real["fitness_mean"].mean()
real_std = df_real["fitness_mean"].std()
standardized_real = (df_real["fitness_mean"] - real_mean) / real_std

random_mean = df_random["fitness_mean"].mean()
random_std = df_random["fitness_mean"].std()
standardized_random = (df_random["fitness_mean"] - random_mean) / random_std

# 进行KS检验，验证标准化数据是否符合正态分布
normal_real = stats.kstest(standardized_real, 'norm')
normal_random = stats.kstest(standardized_random, 'norm')

# 输出KS检验结果
print(f"Real-sequence KS检验 p值: {normal_real.pvalue:.4f}")
print(f"Random-sequence KS检验 p值: {normal_random.pvalue:.4f}")

# 计算均值
mean_real = df_real["fitness_mean"].mean()
mean_random = df_random["fitness_mean"].mean()

# 进行方差齐性检验（Levene检验）
levene_stat, levene_p = stats.levene(df_real["fitness_mean"], df_random["fitness_mean"])
print(f"Levene检验统计量: {levene_stat:.4f}")
print(f"Levene检验 p值: {levene_p:.4f}")

# 判断方差是否齐性
alpha = 0.05
if levene_p < alpha:
    print("方差不齐，采用Welch's t检验（equal_var=False）。")
    equal_var = False
else:
    print("方差齐性成立，采用标准t检验（equal_var=True）。")
    equal_var = True

# 进行独立样本t检验，根据方差齐性结果选择参数
t_stat, p_value = stats.ttest_ind(df_real["fitness_mean"], df_random["fitness_mean"], equal_var=equal_var)

# 输出t检验结果
print(f"Real-sequence 均值: {mean_real:.4f}")
print(f"Random-sequence 均值: {mean_random:.4f}")
print(f"t统计量: {t_stat:.4f}")
print(f"p值: {p_value:.4f}")

# 判断显著性水平
if p_value < alpha:
    print("两组数据的均值存在显著性差异。")
else:
    print("两组数据的均值没有显著性差异。")
