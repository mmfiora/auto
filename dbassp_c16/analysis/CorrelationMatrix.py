import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

# Ask for the CSV file name
file = input("Enter the CSV file name: ")

# Define the path where the CSV is located
path = r"C:\Users\mfiora\Documents\Repositories\auto\dbassp_c16\data\output\\" + file

# Read the CSV file
df = pd.read_csv(path)

# Columns to exclude from the correlation
cols_to_drop = [
    'Peptide ID', 'reference', 'lower_concentration', 'upper_concentration',
    'lower_uM', 'upper_uM', 'ph_run', 'ph',
    'Formation Propensity', 'in vitro Aggregation',
    'Hydrophobic Moment', 'Penetration Depth', 'Tilt Angle',
    'Propensity', 'Normalizer',
    'Disordered Conformation Propensity', 'Linear Moment',
    'Propensity to in vitro Aggregation',
    'Angle Subtended by the Hydrophobic Residues',
    'Amphiphilicity Index', 'Propensity to PPII coil',
    'Normalized Hydrophobic Moment'
]

# Drop only existing columns
df = df.drop(columns=[c for c in cols_to_drop if c in df.columns], errors='ignore')

# Compute correlation matrix
corr = df.corr(numeric_only=True)

# Number of valid data points
n_data = df.select_dtypes(include='number').dropna().shape[0]

# Create the heatmap
plt.figure(figsize=(10, 8))
sns.heatmap(
    corr,
    annot=True,
    cmap='coolwarm',
    fmt=".2f",
    square=True,
    cbar_kws={'shrink': 0.8}
)

# Title and formatting
plt.title(f'Correlation matrix (n = {n_data})', fontsize=14, pad=20)
plt.xticks(rotation=45, ha='right', fontsize=8)
plt.yticks(rotation=0, fontsize=8)
plt.tight_layout()
plt.show()






