import numpy as np
import matplotlib.pyplot as plt

def leer_xvg(nombre_archivo):
    datos = []
    with open(nombre_archivo, 'r') as archivo:
        for linea in archivo:
            if not linea.startswith(('#', '@' , '%')):
                datos.append([float(x) for x in linea.split()])
    return np.array(datos)


archivos = ['eiginval.xvg', 'eiginval_41.xvg', 'eiginval_t.xvg']
colores = ['black', 'red', 'blue']  
etiquetas = ['WT', 'Biphenyl phthalate', 'Testosterone'  ] 

plt.figure(figsize=(10, 10))

for i, archivo in enumerate(archivos):
    datos = leer_xvg(archivo)
    modos = datos[:, 0]  
    eigenvalores = datos[:, 1] 


    plt.plot(modos, eigenvalores, 'o-', color=colores[i], label=etiquetas[i], linewidth=2, markersize=7)


plt.title('', fontsize=14)
plt.xlabel('Eigenvector index', fontsize=20)
plt.ylabel('Eigenvalue (nm$^2$)', fontsize=20)
plt.xlim(0, 20)  
plt.ylim(-1,25)  
plt.grid(True)


plt.legend(fontsize=15)


plt.tick_params(axis='both', which='major', labelsize=16)


plt.tight_layout()


plt.savefig('eigenvalores.png', dpi=300)


plt.show()

print("¡Eigenvalues!")

