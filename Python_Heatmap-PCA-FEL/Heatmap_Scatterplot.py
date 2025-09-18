import pandas as pd
from sklearn.preprocessing import MinMaxScaler
from scipy.stats import zscore
import seaborn as sns
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


df = pd.read_csv("vs.dat", delimiter="\t")

cols_to_standardize = ['Dockthor', 'AutoDock', 'iGEMDOCK']
df_zscore = df.copy()
df_zscore[cols_to_standardize] = zscore(df[cols_to_standardize], nan_policy='omit')
df_zscore.set_index('Ligand', inplace=True)


df_zscore['Average'] = df_zscore.mean(axis=1)

best_ligand = df_zscore['Average'].idxmin()

# -------- HEATMAP --------
plt.figure(figsize=(8, max(10, df_zscore.shape[0] * 0.05)))
sns.heatmap(df_zscore.drop(columns='Average'), cmap="coolwarm",
            linewidths=0.5, linecolor='white', cbar=True)

row_index = df_zscore.index.get_loc(best_ligand)
plt.gca().add_patch(plt.Rectangle((0, row_index), len(cols_to_standardize), 1,
                                  fill=False, edgecolor='green', lw=2))
plt.title("")
plt.tight_layout()
plt.savefig("Heatmap_zscore.png", dpi=300)
plt.show()

# -------- SCATTERPLOT 3D --------
fig = plt.figure(figsize=(10, 8))
ax = fig.add_subplot(111, projection='3d')

x = df_zscore['Dockthor']
y = df_zscore['AutoDock']
z = df_zscore['iGEMDOCK']
avg_color = df_zscore['Average']


scatter = ax.scatter(x, y, z, c=avg_color, cmap="coolwarm", s=120, marker='o')


hx, hy, hz = df_zscore.loc[best_ligand, ['Dockthor', 'AutoDock', 'iGEMDOCK']]
ax.scatter(hx, hy, hz, color='none', edgecolor='green', s=200, linewidths=2, marker='o')



ax.text(hx, hy, hz + 0.25, f"{best_ligand}", color='black',
        fontsize=10, weight='bold', ha='center')


ax.set_xlabel('DockThor')
ax.set_ylabel('AutoDock')
ax.set_zlabel('iGEMDOCK')
ax.set_title(f'')
plt.tight_layout()
plt.savefig("scatterplot_3d_zscore.png", dpi=300)
plt.show()
