
Se va a desarrollar en código de Python un nuevo módulo. Este nuevo módulo pretende contener el código que nos permita sacar 
la envolvente de vuelo del avión F4 para cada peso. 
Partiendo de las gráficas de las páginas 434 y 433 del documento "F4Manual-1979-T-O-1F-4E-1-Flight-Manual-USAF-Series-F-4E-Aircraft",
se sacan los puntos de cada una de las curvas que representan la envolvente del avión para cada uno de los pesos indicados. 
Se recogen todos los puntos en un documento de texto. El input al código debe ser el peso del avión cargado. A partir de ese peso, 
el código entrará en el código correspondiente para los pesos superior e inferior e interpolará entre los datos de esos dos pesos,
generando un conjunto de datos para el peso deseado. Con estos nuevos datos, el código generará la envolvente de vuelo del avión
para ese peso utilizando aproximaciones polinómicas de alto grado que se ajusten bien y tengan alta regresión. 

Archivo "Envolventes.txt"
Primera columna: Valores de números de Mach
Segunda columna: Envolvente para 42777 lb. Altitudes para el correspondiente valor de Mach.
Tercera columna: Envolvente para 43035 lb. Altitudes para el correspondiente valor de Mach.
Cuarta columna: Envolvente para 45472 lb. Altitudes para el correspondiente valor de Mach.
Quinta columna: Envolvente para 46279 lb. Altitudes para el correspondiente valor de Mach.
Sexta columna: Envolvente para 46537 lb. Altitudes para el correspondiente valor de Mach.
Séptima columna: Envolvente para 48974 lb. Altitudes para el correspondiente valor de Mach.