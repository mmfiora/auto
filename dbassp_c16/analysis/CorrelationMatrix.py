import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import os

base_dir = os.path.dirname(__file__)  # carpeta donde está este script
file = input("Enter the CSV file name (without .csv): ") + ".csv"
path = os.path.join(base_dir, "data", "output", file)

df = pd.read_csv(path)

species_name = ""
if 'species' in df.columns:
    species_counts = df['species'].value_counts()
    species_list = sorted(species_counts.index)
    
    print("\nAvailable species:")
    for i, sp in enumerate(species_list, 1):
        count = species_counts[sp]
        print(f"{i}. {sp} (n={count})")
    print(f"{len(species_list) + 1}. All species")
    
    choice = input("\nSelect species number(s) separated by comma: ").strip()
    if choice and choice.replace(',', '').replace(' ', '').isdigit():
        selected_indices = [int(x.strip()) for x in choice.split(',')]
        selected_species = [species_list[i-1] for i in selected_indices if 1 <= i <= len(species_list)]
        
        if selected_species:
            df = df[df['species'].isin(selected_species)]
            if len(selected_species) == 1:
                species_name = f"_{selected_species[0].replace(' ', '_')}"
            else:
                species_name = f"_{'_'.join([sp.replace(' ', '_') for sp in selected_species[:5]])}"
                if len(selected_species) > 5:
                    species_name += "_and_more"

cols_to_drop = [
    'Peptide ID', 'reference', 'lower_concentration', 'upper_concentration',
    #unif'lower_uM', 
    'upper_uM', 'ph_run', 'ph',
    'Formation Propensity', 'in vitro Aggregation',
    'Hydrophobic Moment', 'Penetration Depth', 'Tilt Angle',
    'Propensity', 'Normalizer',
    'Disordered Conformation Propensity', 'Linear Moment',
    'Propensity to in vitro Aggregation',
    'Angle Subtended by the Hydrophobic Residues',
    'Amphiphilicity Index', 'Propensity to PPII coil',
    'Normalized Hydrophobic Moment'
]

df = df.drop(columns=[c for c in cols_to_drop if c in df.columns], errors='ignore')

corr = df.corr(numeric_only=True)

n_data = df.select_dtypes(include='number').dropna().shape[0]

plt.figure(figsize=(10, 8))
sns.heatmap(
    corr,
    annot=True,
    cmap='coolwarm',
    fmt=".2f",
    square=True,
    cbar_kws={'shrink': 0.8}
)

title_text = f'Correlation matrix {file}{species_name} (n = {n_data})'
plt.title(title_text, fontsize=14, pad=20)
plt.xticks(rotation=45, ha='right', fontsize=8)
plt.yticks(rotation=0, fontsize=8)
plt.tight_layout()
output_filename = f"analysis/plots/{file.replace('.csv', '')}{species_name}.jpg"
plt.savefig(output_filename, dpi=300)
plt.show()
