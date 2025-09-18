import numpy as np
import matplotlib.pyplot as plt

def leer_proj_xvg(archivo):
    data = []
    with open(archivo, 'r') as f:
        for line in f:
            if not line.startswith(('#', '@')):  # Ignorar metadatos
                parts = line.strip().split()
                if len(parts) >= 2:
                    data.append([float(parts[0]), float(parts[1])])
    return np.array(data)

def plot_pca(data):
    pc1 = data[:, 0]
    pc2 = data[:, 1]
    frames = np.arange(len(data))

    plt.figure(figsize=(12,10))
    scatter = plt.scatter(pc1, pc2, c=frames, cmap='gnuplot_r', s=30)
    cbar = plt.colorbar(scatter)
    cbar.set_label("Frame", fontsize=28, labelpad=10)
    cbar.ax.tick_params(labelsize=18)
    plt.xlabel("PC1", fontsize=28, labelpad=0)
    plt.ylabel("PC2", fontsize=28, labelpad=0)
    plt.xlim(-12,12)  # Establecer el límite del eje X de -20 a 20
    plt.ylim(-12,12)  # Establecer el límite del eje Y de -20 a 20
    plt.title("")
    plt.grid(True)
    plt.tight_layout()

archivo_xvg = "2d.xvg"

data = leer_proj_xvg(archivo_xvg)
plot_pca(data)

plt.tick_params(axis='both', which='major', labelsize=25)


plt.savefig('PCA_WT.png', dpi=300, bbox_inches='tight')

plt.show()
