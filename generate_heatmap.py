import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import cm
import matplotlib as mpl

# Load dataset
df = pd.read_csv('combined_results_with_bonferroni.csv')

# Filter the data
filtered_df = df[(df['cluster'].isin(['neuron', 'muscle'])) & (df['p_val_adj_all'] < 0.05)].copy()

# Create a new column combining 'cluster' and 'sex'
filtered_df['cluster_sex'] = filtered_df['cluster'] + " " + filtered_df['sex']

# Sort 'cluster_sex' alphabetically
filtered_df['cluster_sex'] = filtered_df['cluster_sex'].astype('category')
filtered_df['cluster_sex'] = filtered_df['cluster_sex'].cat.reorder_categories(sorted(filtered_df['cluster_sex'].unique()))

# Sort gene names by LGIC_group
filtered_df = filtered_df.sort_values(by='LGIC_group')

# Create a new column combining 'gene' and 'LGIC_group' for the y-axis labels
filtered_df['gene_LGIC'] = filtered_df['gene'] + " (" + filtered_df['LGIC_group'] + ")"

# Get unique LGIC_groups to color the y-axis text
lgic_colors = filtered_df.drop_duplicates(subset='gene')[['gene', 'LGIC_group']].set_index('gene')['LGIC_group'].to_dict()
lgic_palette = {
    "ASIC": "blue",
    "cysloop_acetylcholine": "red",
    "cysloop_gaba": "green",
    "iGluR": "orange",
    "P2X": "purple"  # Fixed the typo here
}

# Set color normalization, centering around zero
vmin, vmax = -2, 3  # Define colorbar limits
norm = mpl.colors.TwoSlopeNorm(vmin=vmin, vcenter=0, vmax=vmax)  # Use TwoSlopeNorm to center around 0
cmap = plt.cm.coolwarm  # Choose colormap

# Create figure and axis
fig, ax = plt.subplots(figsize=(12, 7))

# Create scatter plot without automatic legend
scatter = sns.scatterplot(
    data=filtered_df,
    x='cluster_sex',
    y='gene_LGIC',  # Use the new 'gene_LGIC' column for y-axis labels
    hue='avg_log2FC',
    size='pct.1',
    sizes=(20, 200),  # Define min and max sizes
    palette=cmap,
    edgecolor='black',
    alpha=0.7,
    legend=False  # Disable automatic legend
)

# Rotate x-axis labels for better spacing
plt.xticks(rotation=45, ha='right')

# Add gridlines for better readability
plt.grid(True, linestyle='--', linewidth=0.5)

# Create a ScalarMappable object for the colorbar
sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
sm.set_array([])  # Required for colorbar

# Add colorbar outside the plot
cbar = plt.colorbar(sm, ax=ax, aspect=50, pad=0.02)
cbar.set_label('avg_log2FC')

# Define the exact sizes for the legend (10% to 100%) using 10 values
size_legend = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]  # 10 sizes from 10% to 100%
# The scatter plot sizes are from 20 to 200, so we map the range 0.1-1.0 to this range
min_size = 20
max_size = 200
size_legend_in_plot = [min_size + (max_size - min_size) * size for size in size_legend]

# Create handles for the size legend based on the size_legend_in_plot values
handles = []
for size in size_legend_in_plot:
    # Create a legend item for each size
    handles.append(plt.Line2D([0], [0], marker='o', color='w', markerfacecolor='gray', markersize=np.sqrt(size), label=f'{int(size / max_size * 100)}%'))

# Adjust the size legend position to the right, avoid overlap with colorbar
plt.legend(handles=handles, title="The percentage of cells", loc='center left', bbox_to_anchor=(1.1, 0.5), fontsize=10, title_fontsize=12)

# Title and labels
plt.title('Scatterplot Heatmap: avg_log2FC vs Gene with Size based on pct.1')
plt.xlabel('Cluster + Sex')
plt.ylabel('Gene Name (LGIC_group)')

# Set y-axis labels sorted by LGIC_group (avoid using set_yticklabels directly)
yticks = ax.get_yticks()
yticklabels = ax.get_yticklabels()

# Sort y-axis labels based on 'LGIC_group' and apply colors
sorted_yticks = sorted(zip(yticks, yticklabels), key=lambda x: lgic_colors.get(x[1].get_text(), ''))
ax.set_yticks([tick[0] for tick in sorted_yticks])

# Update y-labels to include both gene name and its LGIC_group
new_yticklabels = [f'{tick[1].get_text().split(" (")[0]} ({lgic_colors.get(tick[1].get_text().split(" (")[0], "")})' for tick in sorted_yticks]
ax.set_yticklabels(new_yticklabels)  # Add LGIC_group next to gene name

# Apply color to y-axis labels based on LGIC_group
for label in ax.get_yticklabels():
    gene = label.get_text().split(' (')[0]  # Extract gene name
    lgic_group = lgic_colors.get(gene, '')
    label.set_color(lgic_palette.get(lgic_group, 'black'))

# Save the plot as PNG
plt.tight_layout()
plt.savefig('scatterplot_heatmap_with_LGIC_group_added_to_y_labels.png', format='png', dpi=300)

