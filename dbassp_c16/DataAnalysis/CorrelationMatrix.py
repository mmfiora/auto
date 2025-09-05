# -*- coding: utf-8 -*-
"""
Created on Fri Sep  5 17:01:24 2025

@author: Maria
"""

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# Load the CSV file
df = pd.read_csv("activity_normalized.csv")

# Select only numeric columns (ignore text columns)
numeric_df = df.select_dtypes(include="number")

# Compute the correlation matrix
corr_matrix = numeric_df.corr()

# Plot the correlation matrix as a heatmap
plt.figure(figsize=(10, 8))
sns.heatmap(corr_matrix, annot=True, cmap="coolwarm", fmt=".2f")
plt.title("Correlation Matrix")
plt.show()
