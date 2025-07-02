
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.ndimage import gaussian_filter
from io import StringIO
from mpl_toolkits.mplot3d import Axes3D

with open("2d.xvg", "r") as f:
    lines = f.readlines()


data_lines = [line for line in lines if not line.startswith(('#', '@'))]
data_str = ''.join(data_lines)
data = pd.read_csv(StringIO(data_str), delim_whitespace=True, header=None)
data.columns = ['PC1', 'PC2']

# FEL
xbins = ybins = 100
heatmap, xedges, yedges = np.histogram2d(data['PC1'], data['PC2'], bins=[xbins, ybins])
kT = 0.592
prob = heatmap / np.sum(heatmap)
with np.errstate(divide='ignore', invalid='ignore'):
    FEL = -kT * np.log(prob)
    FEL[~np.isfinite(FEL)] = np.nan
FEL = np.nan_to_num(FEL, nan=np.nanmax(FEL))
FEL_smooth = gaussian_filter(FEL, sigma=0.5)

x_centers = 0.5 * (xedges[:-1] + xedges[1:])
y_centers = 0.5 * (yedges[:-1] + yedges[1:])
X, Y = np.meshgrid(x_centers, y_centers)

fig = plt.figure(figsize=(12, 10))
ax = fig.add_subplot(111, projection='3d')

surf = ax.plot_surface(X, Y, FEL_smooth.T, cmap='gnuplot', edgecolor='none', alpha=1)

Z_base = np.full_like(FEL_smooth.T, 4)
projection_colors = plt.cm.gnuplot(FEL_smooth.T / np.max(FEL_smooth.T))
ax.plot_surface(X, Y, Z_base,
                facecolors=projection_colors,
                edgecolor='none',
                shade=False,
                alpha=1.0,
                rstride=1, cstride=1, antialiased=False)

ax.set_xlabel('PC1', fontsize=20, labelpad=20)
ax.set_ylabel('PC2', fontsize=20, labelpad=20)
ax.set_zlabel('', fontsize=20, labelpad=20)
ax.set_xlim(X.min(), X.max())
ax.set_ylim(Y.min(), Y.max())
ax.set_zlim(4, 5.6)
ax.set_title('')
ax.view_init(elev=25, azim=75)
ax.grid(True)

plt.tick_params(axis='both', which='major', labelsize=15)

mappable = plt.cm.ScalarMappable(cmap='gnuplot')
mappable.set_array(FEL_smooth)
cbar = fig.colorbar(mappable, ax=ax, shrink=0.5, aspect=15, pad=0.04)
cbar.set_label('FEL (kcal/mol)',fontsize=20, labelpad=25)

plt.tight_layout()
plt.savefig("FEL_3D_2D_WT.png", dpi=400, transparent=True)
plt.show()
