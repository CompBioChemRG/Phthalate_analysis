import pandas as pd
from sklearn.preprocessing import MinMaxScaler
import seaborn as sns
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


df = pd.read_csv("vs.dat", delimiter="\t")


cols_to_normalize = ['Dockthor', 'AutoDock', 'iGEMDOCK']
scaler = MinMaxScaler()
normalized = scaler.fit_transform(df[cols_to_normalize])
df_normalized = pd.DataFrame(normalized, columns=cols_to_normalize)
df_normalized['Ligand'] = df['Ligand']
df_normalized = df_normalized[['Ligand'] + cols_to_normalize]
df_normalized.set_index('Ligand', inplace=True)


df_normalized['Average'] = df_normalized.mean(axis=1)
best_ligand = df_normalized['Average'].idxmin()


shared_cmap = "coolwarm"

# -------- HEATMAP 
plt.figure(figsize=(8, max(10, df_normalized.shape[0] * 0.05)))  # Altura adaptativa
ax = sns.heatmap(df_normalized.drop(columns='Average'), cmap=shared_cmap,
                 linewidths=0.5, linecolor='white', cbar=True)


row_index = df_normalized.index.get_loc(best_ligand)
for col_index in range(len(cols_to_normalize)):
    ax.add_patch(plt.Rectangle((col_index, row_index), 1, 1, fill=False, edgecolor='green', lw=2))

plt.title(f"")
plt.xlabel("")
plt.ylabel("Ligand")
plt.yticks(rotation=0)
plt.tight_layout()
plt.savefig("heatmap_normalizado_contorno.png", dpi=300)
plt.show()

# -------- SCATTERPLOT 3D 
fig = plt.figure(figsize=(10, 8))
ax = fig.add_subplot(111, projection='3d')

x = df_normalized['Dockthor']
y = df_normalized['AutoDock']
z = df_normalized['iGEMDOCK']
avg_color = df_normalized['Average']


scatter = ax.scatter(x, y, z, c=avg_color, cmap=shared_cmap, s=120, marker='o')


hx, hy, hz = df_normalized.loc[best_ligand, ['Dockthor', 'AutoDock', 'iGEMDOCK']]
ax.scatter(hx, hy, hz, color='none', edgecolor='green', s=200, linewidths=2, marker='o')



ax.text(hx, hy, hz + 0.05, f"{best_ligand}", color='black',
        fontsize=10, weight='bold', ha='center')


ax.set_xlabel('DockThor')
ax.set_ylabel('AutoDock')
ax.set_zlabel('iGEMDOCK')
ax.set_title(f'')
plt.tight_layout()
plt.savefig("scatterplot_3d_normalizado_contorno.png", dpi=300)
plt.show()
