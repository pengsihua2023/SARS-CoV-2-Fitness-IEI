import pandas as pd
import numpy as np
from scipy.stats import ks_2samp
import matplotlib.pyplot as plt

# 1. Read CSV files
real_data_file = 'real-Global-Fitness.csv'
random_data_file = 'Fitness-randow-Golobal-0323-0.003.csv'

# Load data
real_df = pd.read_csv(real_data_file)
random_df = pd.read_csv(random_data_file)

# Extract data (column name 'fitness_mean')
real_data = real_df['fitness_mean'].values
random_data = random_df['fitness_mean'].values

# 2. Calculate descriptive statistics
# Real data statistics
print("Real Data Description:")
print(real_df['fitness_mean'].describe())
print("Real data NaN count:", real_df['fitness_mean'].isna().sum())

# Random data statistics
print("\nRandom Data Description:")
print(random_df['fitness_mean'].describe())
print("Random data NaN count:", random_df['fitness_mean'].isna().sum())

# 3. Perform KS test
statistic, p_value = ks_2samp(real_data, random_data)

# Output results
print(f"\nKS Statistic (D): {statistic}")
print(f"p-value: {p_value}")

# Interpret results
alpha = 0.05  # Significance level
if p_value >= alpha:
    print("Cannot reject the null hypothesis, no significant difference between the distributions of the two groups.")
else:
    print("Reject the null hypothesis, the distributions of the two groups are significantly different.")

# 4. Visualize CDF
# Calculate CDF
real_cdf = np.sort(real_data)
real_cdf_y = np.arange(len(real_data)) / float(len(real_data))
random_cdf = np.sort(random_data)
random_cdf_y = np.arange(len(random_data)) / float(len(random_data))

# Plot
plt.plot(real_cdf, real_cdf_y, label='Real Data CDF (North America)')
plt.plot(random_cdf, random_cdf_y, label='Random Data CDF (North America)')
plt.legend()
plt.title('Cumulative Distribution Functions of Fitness Mean (North America)')
plt.xlabel('Fitness Mean')
plt.ylabel('Cumulative Probability')
plt.grid(True)
plt.show()

plt.hist(real_data, bins=50, alpha=0.5, label='Real Data', density=True)
plt.hist(random_data, bins=50, alpha=0.5, label='Random Data', density=True)
plt.legend()
plt.title('Histogram of Fitness Mean (North America)')
plt.xlabel('Fitness Mean')
plt.ylabel('Density')
plt.show()