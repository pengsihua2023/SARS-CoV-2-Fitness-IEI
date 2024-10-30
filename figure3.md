```
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import make_interp_spline

# 定义文件路径
# file_path = 'D:/test/results/means-18-groups-0522.csv'
#file_path = 'D:/test/results/means/means-18-groups-Escape-final.csv'
#file_path = 'D:/CoVFit/new_data/non-USA/Europe-18-group-Immune_Escape_mean-final.csv'
# file_path = 'D:/CoVFit/North-America-Results/North-America-2020-2024-Immune-big.csv'
file_path = 'D:/CoVFit/North-America-Results/North-America-2020-2024-big.csv'



# 读取文件
means_df = pd.read_csv(file_path)

# 初始化一个字典来存储每组数据
data_dict = {}
for column in means_df.columns:
    data_dict[column] = means_df[column].dropna().values

# 创建箱线图数据
data_to_plot = [data_dict[group] for group in means_df.columns]

# 定义 flierprops 以设置 outliers 的颜色和边界颜色
flierprops = dict(marker='o', markerfacecolor='red', markeredgecolor='red', markersize=8, linestyle='none')

# 定义 boxprops 以设置箱体颜色为淡绿色
boxprops = dict(facecolor='#bdecb6', color='blue')

# 绘制箱线图并设置属性
plt.figure(figsize=(18, 8))
boxplot_elements = plt.boxplot(data_to_plot, positions=range(len(means_df.columns)), flierprops=flierprops,
                               boxprops=boxprops, patch_artist=True,
                               whiskerprops=dict(color='blue'),
                               capprops=dict(color='blue'),
                               medianprops=dict(color='blue'))

# 计算每组数据的均值并绘制平滑曲线
means = [np.mean(data_dict[group]) for group in means_df.columns]
x = range(len(means_df.columns))
x_new = np.linspace(min(x), max(x), 300)
spl = make_interp_spline(x, means, k=3)
means_smooth = spl(x_new)
plt.plot(x_new, means_smooth, color='blue', label='Mean Value (Smoothed)')

# 在每个箱体的右侧标注中位数值和样本数量
for median, group in zip(boxplot_elements['medians'], means_df.columns):
    median_x = median.get_xdata()[1]  # X位置的中位线
    median_y = median.get_ydata()[1]  # Y位置的中位线
    group_size = len(data_dict[group])  # 每组的样本数量
    plt.text(median_x + 0.05, median_y, f'{median_y:.2f}',  # 格式化中位数值
             verticalalignment='center', color='blue', fontsize=13)  # 垂直居中对齐
    plt.text(median_x + 0.05, median_y - (median_y * 0.05), f'n={group_size}',  # 格式化样本数量
             verticalalignment='center', color='black', fontsize=13)  # 紫色字体

# 设置x轴的刻度
plt.xticks(range(len(means_df.columns)), means_df.columns, rotation=40, fontsize=15)
plt.xlabel('Time Period', fontsize=25, labelpad=50)
plt.ylabel('Immune Escape Index', fontsize=25)
# plt.title('Boxplot for Each Time Period', fontsize=25)
plt.legend(fontsize=15)
# plt.grid(True)


# 调整图形边缘留白
plt.subplots_adjust(bottom=0.0001, top=0.99)  # 进一步减少底部留白

# 显示图形
plt.tight_layout()  # 使用更紧凑的布局
plt.show()
```
