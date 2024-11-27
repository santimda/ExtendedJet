import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
def leer_archivo_texto(ruta_archivo):
    datos = []

    with open(ruta_archivo, 'r') as archivo:
        for linea in archivo:
            # Limpiar la línea y dividir en números
            numeros = linea.strip().split()
            
            # Convertir los números a enteros (o flotantes si lo prefieres)
            if len(numeros) != 9:
                print(f"Advertencia: La línea no tiene 9 números: {linea.strip()}")
                continue
            
            numeros = [num for num in numeros]
            
            datos.append(numeros)
    
    return datos

# Ruta al archivo de texto
d=5.1839e+25
ruta_archivo = 'SEDe_syn_jet.dat'

# Llamar a la función para leer el archivo
datos = leer_archivo_texto(ruta_archivo)


for i, fila in enumerate(datos):
    fila[7]=str(float(fila[7])/((4*np.pi)*d*d))
    
    fila[8]=str(float(fila[8])/((4*np.pi)*d*d))

# Mostrar los datos leídos
with open('conversion.dat','w') as archivo:
    for i, fila in enumerate(datos):
        contenido = ', '.join(fila)
        archivo.write(contenido+'\n')


df = pd.read_csv('conversion.dat')
# Crear el gráfico
plt.figure(figsize=(10, 6))
plt.plot(df.iloc[:, 0], df.iloc[:, 7], marker='o', linestyle='-', color='b')
plt.xlabel('Columna 1')
plt.ylabel('Columna 8')
plt.title('Gráfico de Columnas 1 y 8')
plt.grid(True)
plt.show()
