# -*- coding: utf-8 -*-
"""
Created on Tue Apr 17 09:53:21 2018

@author: Team REOS

Módulo para sacar la envolvente de vuelo del avión F4. Como input tiene el 
PESO del avión cargado con el REOS. Se realiza una interpolación con las
envolventes conocidas de unos pesos determinados y se realiza una aproximación
polinómica de alto orden para recrear las ecuaciones de los límites de la 
envolvente para el PESO de entrada.
"""

from inputs_iniciales import MASS

# Definción de variables. Todos los pesos en libras
PESO = MASS # Para cuando quitemos la línea de abajo e importemos directamente del otro módulo
PESO = 47300.0  # Input de entrada. Peso del avión cargado con el cohete REOS
P1 = 42777.0 
P2 = 43035.0
P3 = 45472.0
P4 = 46279.0
P5 = 46537.0
P6 = 48974.0

alt_enve = []  # Lista de valores de altitud de la curva límite de la envolvente
mach_enve = []  # Lista de valores de Mach de la curva límite de la envolvente

F = open("Envolventes.txt", 'r')  # Apertura del fichero de lectura.

for i in F:  # Bucle para recorrer el archivo de texto
    i = i.strip()  # Sentencias para eliminar los espacios sobrantes de F
    i = i.split() 
    for j, a in enumerate(i):  # Bucle para recorrer las filas del archivo
            i[j] = float(i[j])  # Convierte los valores de strings a floats
    # Condicionales para separar los casos, las curvas entre las que hay que
    # interpolar. Interpolaciones lineales en altitud para mismo Mach
    # Para el último punto se realiza una interpolación lineal en el número
    # de Mach para igual valor final de altitud (36000ft)
    # Para PESO < P1 y para PESO > P6 lo que se realiza es una extrapolación 
    if PESO < P1:
        print("ERROR")
    
    elif P1 < PESO < P2:
        enve = i[1] + (i[2]-i[1]) * (PESO - P1) / (P2 - P1)
        alt_enve.append(enve)
        mach_enve.append(i[0])

    elif P2 < PESO < P3:
        if i[0] <= 1.88:        
            enve = i[2] + (i[3] - i[2]) * (PESO - P2) / (P3 - P2)
            alt_enve.append(enve)
            mach_enve.append(i[0])                            
        else:
            Mfin = i[8] + (i[9] - i[8]) * (PESO - P2) / (P3 - P2)            
    
    elif P3 < PESO < P4:
        if i[0] <= 1.78:
            enve = i[3] + (i[4] - i[3]) * (PESO - P3) / (P4 - P3)
            alt_enve.append(enve)
            mach_enve.append(i[0])          
        else:
            Mfin = i[9] + (i[10] - i[9]) * (PESO - P3) / (P4 - P3)
    
    elif P4 < PESO<P5:
        if i[0] <= 1.76:
            enve = i[4] + (i[5] - i[4]) * (PESO - P4) / (P5 - P4)
            alt_enve.append(enve)
            mach_enve.append(i[0])
        else:
            Mfin = i[10] + (i[11] - i[10]) * (PESO - P4) / (P5 - P4)
    
    elif P5 < PESO < P6:
        if i[0] <= 1.5:
            enve = i[5] + (i[6] - i[5]) * (PESO - P5) / (P6 - P5)
            alt_enve.append(enve)
            mach_enve.append(i[0])
        else:
            Mfin = i[11] + (i[12] - i[11]) * (PESO - P5) / (P6 - P5)
            
    elif PESO > P6:
        print("ERROR")
        
# Los siguientes condicionales añaden el último punto a las listas
if P2<PESO<P3:
    mach_enve.append(Mfin)
    alt_enve.append(36.0)       
elif P3<PESO<P4:
    mach_enve.append(Mfin)
    alt_enve.append(36.0)
elif P4<PESO<P5:
    mach_enve.append(Mfin)
    alt_enve.append(36.0)
elif P5<PESO<P6:
    mach_enve.append(Mfin)
    alt_enve.append(36.0)

F.close() # Sentencia para cerrar el archivo de levtura "Envolventes"