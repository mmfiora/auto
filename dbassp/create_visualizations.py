import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

# Set style for report-ready publication quality
sns.set_style("white")
plt.rcParams['figure.dpi'] = 600
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['DejaVu Sans', 'Arial', 'Helvetica']

# Load data
c16_data = pd.read_csv('data/output/unified_results_C16.csv')
c12_data = pd.read_csv('data/output/unified_results_C12.csv')

def prepare_corr_matrix(df, cols, labels, dropna=True):
    df_sub = df[cols].copy()
    for col in cols:
        df_sub[col] = pd.to_numeric(df_sub[col], errors='coerce')
    if dropna:
        df_sub = df_sub.dropna()
    corr_matrix = df_sub.corr()
    corr_matrix.columns = labels
    corr_matrix.index = labels
    return corr_matrix, len(df_sub)

# Variables including MW_Da
cols_c16 = ['upper_uM', 'npol_min', 'curv_min', 'logD', 'Normalized Hydrophobicity', 'MW_Da']
cols_c12 = ['lower_uM', 'npol_min', 'curv_min', 'logD', 'Normalized Hydrophobicity', 'MW_Da']
spa_labels = ['MIC', 'n_pol', 'curv_min', 'logD', 'Hidrofob. Norm.', 'Peso Mol.']

c16_corr, c16_n = prepare_corr_matrix(c16_data, cols_c16, spa_labels, dropna=False)
c12_corr, c12_n = prepare_corr_matrix(c12_data, cols_c12, spa_labels, dropna=False)

# ============================================
# Combined Figure: Matrices C16 (Left) and C12 (Right)
# ============================================
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(20, 8))

# Mask for upper triangle
mask_c16 = np.triu(np.ones_like(c16_corr, dtype=bool), k=1)
mask_c12 = np.triu(np.ones_like(c12_corr, dtype=bool), k=1)

# Plot C16
sns.heatmap(c16_corr, annot=True, fmt='.2f', ax=ax1,
            annot_kws={'fontsize': 12, 'fontweight': 'bold'},
            cmap='coolwarm', center=0, square=True, linewidths=1.5,
            cbar_kws={"shrink": 0.8}, vmin=-1, vmax=1, mask=mask_c16)
ax1.set_title('Matriz de Correlación: C16', fontsize=18, fontweight='bold', pad=20)
ax1.tick_params(axis='x', rotation=45)
ax1.tick_params(axis='y', rotation=0)

# Plot C12
sns.heatmap(c12_corr, annot=True, fmt='.2f', ax=ax2,
            annot_kws={'fontsize': 12, 'fontweight': 'bold'},
            cmap='coolwarm', center=0, square=True, linewidths=1.5,
            cbar_kws={"shrink": 0.8}, vmin=-1, vmax=1, mask=mask_c12)
ax2.set_title('Matriz de Correlación: C12', fontsize=18, fontweight='bold', pad=20)
ax2.tick_params(axis='x', rotation=45)
ax2.tick_params(axis='y', rotation=0)

plt.tight_layout(pad=4.0)
plt.savefig('data/output/combined_matrices_c16_c12.png', dpi=600, bbox_inches='tight')
print("Guardado: data/output/combined_matrices_c16_c12.png")
plt.close()

# Keep individual saves for convenience
def save_individual(matrix, title, filename):
    plt.figure(figsize=(10, 8))
    mask = np.triu(np.ones_like(matrix, dtype=bool), k=1)
    sns.heatmap(matrix, annot=True, fmt='.2f',
                annot_kws={'fontsize': 14, 'fontweight': 'bold'},
                cmap='coolwarm', center=0, square=True, linewidths=1.5,
                cbar_kws={"shrink": 0.8}, vmin=-1, vmax=1, mask=mask)
    plt.title(title, fontsize=18, fontweight='bold', pad=20)
    plt.xticks(rotation=45, ha='right')
    plt.yticks(rotation=0)
    plt.tight_layout()
    plt.savefig(filename, dpi=600, bbox_inches='tight')
    plt.close()

save_individual(c16_corr, 'Matriz de Correlación: C16', 'data/output/correlation_matrix_C16.png')
save_individual(c12_corr, 'Matriz de Correlación: C12', 'data/output/correlation_matrix_C12.png')

# Final Scatter Plot
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))
ax1.scatter(c16_data['curv_min'], c16_data['upper_uM'], alpha=0.5, color='#2E86AB')
ax1.set_yscale('log')
ax1.set_title('C16: MIC vs curv_min', fontweight='bold')
ax1.set_xlabel('curv_min')
ax1.set_ylabel('MIC (uM)')

ax2.scatter(c12_data['curv_min'], c12_data['lower_uM'], alpha=0.5, color='#A23B72')
ax2.set_yscale('log')
ax2.set_title('C12: MIC vs curv_min', fontweight='bold')
ax2.set_xlabel('curv_min')
ax2.set_ylabel('MIC (uM)')

plt.tight_layout()
plt.savefig('data/output/scatter_mic_vs_curvmin.png', dpi=600)
print("Scatter plot saved.")
