import numpy as np
import matplotlib.pyplot as plt

def leer_xvg(nombre_archivo):
    datos = []
    with open(nombre_archivo, 'r') as archivo:
        for linea in archivo:
            if not linea.startswith(('#', '@','&')):
                datos.append([float(x) for x in linea.split()])
    return np.array(datos)

archivos = ['2d.xvg', '2d_41.xvg', '2d_t.xvg']
colores = ['black', 'red', 'blue']  
etiquetas = ['WT', 'Biphenyl phthalate', 'Testosterone'  ]


plt.figure(figsize=(10, 10))

for i, archivo in enumerate(archivos):
    datos = leer_xvg(archivo)
    modo1 = datos[:, 0]  
    modo2 = datos[:, 1]  

    plt.plot(modo1, modo2, 'o', color=colores[i], label=etiquetas[i], markersize=3)

plt.title('', fontsize=14)
plt.xlabel('PCA1', fontsize=20)
plt.ylabel('PCA2', fontsize=20) 
plt.xlim(-12,12)  
plt.ylim(-12,12)  
plt.grid(True)


plt.legend(fontsize=15)

plt.tick_params(axis='both', which='major', labelsize=16)

plt.tight_layout()

plt.savefig('pca_2d.png', dpi=300)

plt.show()

print("¡PCA 2D!")

