# -*- coding: utf-8 -*-
"""

@author: TeamREOS

Funciones que nos dan los coeficientes aerodinámicos del avión F4.
Coeficientes aerodinámicos de sustentación y de resistencia y todo aquello
relacionado con los mismos.

"""

from math import radians, pi, cos,degrees
from scipy.interpolate import interp1d


# ------------------CARACTERÍSTICAS GEOMÉTRICAS DEL VEHÍCULO------------------

S_W = 49.2  # Superficie alar (m2).
LF = 19.2  # Longitud del fuselaje (m).
BF = 2.87  # Envergadura del fuselaje (m).
NE = 2  # Número de motores.
DM = 0  # Diámetro de cada motor.
LM = 0  # Longitud de cada motor.
S_H = 6.39  # Área de la superficie de mando horizontal (m2).
S_V = 5.035  # Área de la superficie de mando vertical (m2).
B = 11.7  # Envergadura (m).
AR = B**2 / S_W  # Alargamiento = 2,78.
FLECHA = radians(52)  # Flecha BA.
FLECHA2 = radians(41.4)  # Flecha 1/4.
ESPESOR = .064  # Espesor relativo máximo.
K = .4  # EL perfil del F-4 es el NACA0006.4-64 por tanto la K es 0,4.
ESTRECHAMIENTO = .26  # Estrechamiento.
M_C = 1 - (.065 * (cos(FLECHA))**.6) * (100 * ESPESOR)**.637814
# Mach crítico.
M_D = 1.08 * M_C  # Mach de divergencia.
M_D098 = .98 * M_D


def cl (mach, angulo_ataque_radians) :
    
    angulo_ataque=degrees(angulo_ataque_radians)
    
    '''Cálculo  del coeficiente de sustentación.  Este coeficiente varía
    con respecto al Mach y, con el ángulo de ataque.'''
    if mach < 0:
        print ('Wrong Mach')
    elif mach <= .8:
        if -10 <= angulo_ataque < 23:
            return 3.0062337662E-02*angulo_ataque + 2.1288892284E-01
        elif 23 <= angulo_ataque <= 40:
            return 7.9721362229E-03*angulo_ataque + 6.8104437564E-01        
    elif 0.80 <= mach < 0.85:
    #Se importan los archivos con los CL      
        archivo = open('pruebas', "w")        
        Mach1 = open('0.80', 'r')
        AoA1 = []
        CL1 = []
        for line in Mach1:
            s = line.strip()
            s = s.split()
            AoA1.append(float(s[0]))  # AoA.
            CL1.append(float(s[1]))  # CL.
        Mach1.close()
        
        Mach2 = open('0.85', 'r')
        AoA2 = []
        CL2 = []
        for line in Mach2:
            s = line.strip()
            s = s.split()
            AoA2.append(float(s[0]))  # AoA.
            CL2.append(float(s[1]))  # CL.
        Mach2.close()
        
        def cl1(angulo_ataque):           
            cl = interp1d(AoA1, CL1, )
            return float(cl(angulo_ataque))        
        def cl2(angulo_ataque):
            
            cl = interp1d(AoA2, CL2, )
            return float(cl(angulo_ataque))  
                   
#Generación de 'archivo', donde se almacenaran los datos.
        
        #archivo.write('M \tCL \n')
        archivo.write('%.8f\t' % 0.80)
        archivo.write('%.8f\n' % cl1(angulo_ataque))
        archivo.write('%.8f\t' % 0.85)
        archivo.write('%.8f\n' % cl2(angulo_ataque))
        archivo.close()
        
        pruebas = open('pruebas', 'r')
        MACH = []
        CL = []
        for line in pruebas:
            s = line.strip()
            s = s.split()
            MACH.append(float(s[0]))  # Mach.
            CL.append(float(s[1]))  # CL.
        pruebas.close()

    
    
        def it(mach):
            it = interp1d(MACH, CL,)
            return float(it(mach))
        return float(it(mach))
        

        
    elif 0.85 <= mach < 0.95:
        if -10 <= angulo_ataque < 21:
            return 3.2475256598E-02*angulo_ataque + 2.3447983871E-01
        elif 21 <= angulo_ataque<= 40:
            return 6.7571428571E-03*angulo_ataque + 7.2975714286E-01
        
    elif 0.95 <= mach < 1.0:
    #Se importan los archivos con los CL  
    
        archivo = open('pruebas', "w")
        
        Mach1 = open('0.85', 'r')
        AoA1 = []
        CL1 = []
        for line in Mach1:
            s = line.strip()
            s = s.split()
            AoA1.append(float(s[0]))  # AoA.
            CL1.append(float(s[1]))  # CL.
        Mach1.close()
        
        Mach2 = open('1.00', 'r')
        AoA2 = []
        CL2 = []
        for line in Mach2:
            s = line.strip()
            s = s.split()
            AoA2.append(float(s[0]))  # AoA.
            CL2.append(float(s[1]))  # CL.
        Mach2.close()
        
        def cl1(angulo_ataque):
            
            cl = interp1d(AoA1, CL1)
            return float(cl(angulo_ataque))
        
        def cl2(angulo_ataque):
            
            cl = interp1d(AoA2, CL2)
            return float(cl(angulo_ataque))
        
        #Generación de 'archivo', donde se almacenaran los datos.
        
        #archivo.write('M \tCL \n')
        archivo.write('%.8f\t' % 0.95)
        archivo.write('%.8f\n' % cl1(angulo_ataque))
        archivo.write('%.8f\t' % 1.00)
        archivo.write('%.8f\n' % cl2(angulo_ataque))
        archivo.close()
        
        pruebas = open('pruebas', 'r')
        MACH = []
        CL = []
        for line in pruebas:
            s = line.strip()
            s = s.split()
            MACH.append(float(s[0]))  # Mach.
            CL.append(float(s[1]))  # CL.
        pruebas.close()
    
    
        def it(mach):
            it = interp1d(MACH, CL)
            return float(it(mach))
        return float(it(mach))
        
    elif 1.0 <= mach < 1.1 :
        if -10 <= angulo_ataque < 10:
            return 3.6610389610E-02*angulo_ataque + 1.8380952381E-01
        elif 10 <= angulo_ataque < 23:
            return 2.5934065934E-02*angulo_ataque + 2.9065934066E-01
        elif 23 <= angulo_ataque <= 40:
            return 8.0650154799E-03*angulo_ataque + 6.9956312350E-01
        
    elif 1.1 <= mach < 1.15:
    #Se importan los archivos con los CL  
    
        archivo = open('pruebas', "w")
        
        Mach1 = open('1.00', 'r')
        AoA1 = []
        CL1 = []
        for line in Mach1:
            s = line.strip()
            s = s.split()
            AoA1.append(float(s[0]))  # AoA.
            CL1.append(float(s[1]))  # CL.
        Mach1.close()
        
        Mach2 = open('1.15', 'r')
        AoA2 = []
        CL2 = []
        for line in Mach2:
            s = line.strip()
            s = s.split()
            AoA2.append(float(s[0]))  # AoA.
            CL2.append(float(s[1]))  # CL.
        Mach2.close()
        
        def cl1(angulo_ataque):
            
            cl = interp1d(AoA1, CL1)
            return float(cl(angulo_ataque))
        
        def cl2(angulo_ataque):
            
            cl = interp1d(AoA2, CL2)
            return float(cl(angulo_ataque))
        
        
        #Generación de 'archivo', donde se almacenaran los datos.
        
        #archivo.write('M \tCL \n')
        archivo.write('%.8f\t' % 1.10)
        archivo.write('%.8f\n' % cl1(angulo_ataque))
        archivo.write('%.8f\t' % 1.15)
        archivo.write('%.8f\n' % cl2(angulo_ataque))
        archivo.close()
        
        pruebas = open('pruebas', 'r')
        MACH = []
        CL = []
        for line in pruebas:
            s = line.strip()
            s = s.split()
            MACH.append(float(s[0]))  # Mach.
            CL.append(float(s[1]))  # CL.
        pruebas.close()
    
    
        def it(mach):
            it = interp1d(MACH, CL)
            return float(it(mach))
        return float(it(mach))         
              
    elif 1.15 <= mach < 1.25 :
        if -10 <= angulo_ataque < 20:
            return 2.9862903226E-02*angulo_ataque + 1.0197580645E-01
        elif  20 <= angulo_ataque < 29:
            return 2.3272727273E-02*angulo_ataque + 2.3381818182E-01
        elif 29 <= angulo_ataque <= 40 :
            return  9.3706293706E-03*angulo_ataque + 6.4504662005E-01
        
    elif 1.25 <= mach < 1.35:
        if -10 <= angulo_ataque < 23.4522:
            return 2.6546629171E-02*angulo_ataque + 8.9401031626E-02
        elif 23.4522 <= angulo_ataque <= 40:
            return 1.8440175969E-02*angulo_ataque + 2.6431721265E-01
        
    elif 1.35 <= mach < 1.45:
        return 2.3575113122E-02*angulo_ataque + 9.4314479638E-02
    
    elif 1.45 <= mach < 1.55:
        return 2.3067420814E-02*angulo_ataque + 1.1516515837E-01
    
    elif 1.55 <= mach < 1.95 :
#         Se importan los archivos con los CL  
    
        archivo = open('pruebas', "w")
        
        Mach1 = open('1.55', 'r')
        AoA1 = []
        CL1 = []
        for line in Mach1:
            s = line.strip()
            s = s.split()
            AoA1.append(float(s[0]))  # AoA.
            CL1.append(float(s[1]))  # CL.
        Mach1.close()
        
        Mach2 = open('1.95', 'r')
        AoA2 = []
        CL2 = []
        for line in Mach2:
            s = line.strip()
            s = s.split()
            AoA2.append(float(s[0]))  # AoA.
            CL2.append(float(s[1]))  # CL.
        Mach2.close()
        
        def cl1(angulo_ataque):
            
            cl = interp1d(AoA1, CL1)
            return float(cl(angulo_ataque))
        
        def cl2(angulo_ataque):
            
            cl = interp1d(AoA2, CL2)
            return float(cl(angulo_ataque))
        
#         Generación de 'archivo', donde se almacenaran los datos.
        
        #archivo.write('M \tCL \n')
        archivo.write('%.8f\t' % 1.55)
        archivo.write('%.8f\n' % cl1(angulo_ataque))
        archivo.write('%.8f\t' % 1.95)
        archivo.write('%.8f\n' % cl2(angulo_ataque))
        archivo.close()
        
        pruebas = open('pruebas', 'r')
        MACH = []
        CL = []
        for line in pruebas:
            s = line.strip()
            s = s.split()
            MACH.append(float(s[0]))  # Mach.
            CL.append(float(s[1]))  # CL.
        pruebas.close()

    
    
        def it(mach):
            it = interp1d(MACH, CL)
            return it(mach)
        return float(it(mach))
        
    elif 1.95 <= mach < 2.05:
        return 2.9084615385E-02*angulo_ataque + 1.5265233786E-01
    
def angulo_ataque(mach, cl):    
    
    if mach < 0.0:
        print ("Error")
    elif mach <= 0.8:
        alfa = ((2.9841395494E+02*cl**6 - 6.0985787512E+02*cl**5 + 
               4.3115330281E+02*cl**4 - 1.0252886994E+02*cl**3 -  
               1.1424053219E+00*cl**2+3.1396260600E+01*cl- 6.4719237053E+00)*
               0.0174528666666667)
        return alfa
    
    
    elif 0.8 < mach <= 0.85:
        archivo = open('angulos', "w")
        
        A1 = open('0_8', 'r')
        CL1 = []
        AoA11 = []
        for line in A1:
            s = line.strip()
            s = s.split()
            CL1.append(float(s[0]))  # CL .
            AoA11.append(float(s[1])* 0.0174528666666667)  # AoA (rad).
        A1.close()
        
        A2 = open('0_85', 'r')
        CL2 = []
        AoA12 = []
        for line in A2:
            s = line.strip()
            s = s.split()
            CL2.append(float(s[0]))  # CL .
            AoA12.append(float(s[1])* 0.0174528666666667)  # AoA (rad).
        A2.close()
        
        def AoA1 (cl):
            angulo_ataque = interp1d(CL1,AoA11)
            return float (angulo_ataque(cl))
        def AoA2 (cl):
            angulo_ataque= interp1d(CL2,AoA12)
            return float (angulo_ataque(cl))
        
#        Generación de 'archivo', donde se almacenaran los datos.
        
#        archivo.write('M \tAoA \n')
        archivo.write('%.8f\t' % 0.80)
        archivo.write('%.8f\n' % AoA1(cl))
        archivo.write('%.8f\t' % 0.85)
        archivo.write('%.8f\n' % AoA2(cl))
        archivo.close()
#        
        angulos= open ('angulos','r')
        MACH = []
        AoA =[]
        for line in angulos:
            s=line.strip()
            s=s.split()
            MACH.append(float(s[0]))
            AoA.append (float(s[1]))
        angulos.close()
        
        
        def angulo (mach):
            angulo=interp1d(MACH,AoA)
            return float (angulo(mach))
        return float (angulo(mach))


    elif 0.85 < mach <= 0.95:
        return ((2.9118855999E+02*cl**6 - 5.2781379866E+02*cl**5 + 
                3.0828703332E+02*cl**4 - 4.2295275218E+01*cl**3 - 
                1.1267719691E+01*cl**2 + 2.9625959005E+01*cl-6.6003246147E+00)*
                0.0174528666666667)
    
    elif 0.95 < mach <= 1.0:
        archivo = open('angulos', "w")
        
        A1 = open('0_85', 'r')
        CL1 = []
        AoA11 = []
        for line in A1:
            s = line.strip()
            s = s.split()
            CL1.append(float(s[0]))  # CL .
            AoA11.append(float(s[1])* 0.0174528666666667)  # AoA (rad).
        A1.close()
        
        A2 = open('1_0', 'r')
        CL2 = []
        AoA12 = []
        for line in A2:
            s = line.strip()
            s = s.split()
            CL2.append(float(s[0]))  # CL .
            AoA12.append(float(s[1])* 0.0174528666666667)  # AoA (rad).
        A2.close()
        
        def AoA1 (cl):
            angulo_ataque = interp1d(CL1,AoA11)
            return float (angulo_ataque(cl))
        def AoA2 (cl):
            angulo_ataque= interp1d(CL2,AoA12)
            return float (angulo_ataque(cl))
        
#        Generación de 'archivo', donde se almacenaran los datos.
        
#        archivo.write('M \tAoA \n')
        archivo.write('%.8f\t' % 0.85)
        archivo.write('%.8f\n' % AoA1(cl))
        archivo.write('%.8f\t' % 1.0)
        archivo.write('%.8f\n' % AoA2(cl))
        archivo.close()
#        
        angulos= open ('angulos','r')
        MACH = []
        AoA =[]
        for line in angulos:
            s=line.strip()
            s=s.split()
            MACH.append(float(s[0]))
            AoA.append (float(s[1]))
        angulos.close()
        
        
        def angulo (mach):
            angulo=interp1d(MACH,AoA)
            return float (angulo(mach))
        return float (angulo(mach))
    

    elif 1.0 < mach <= 1.1 :
        return ((2.4428221582E+02*cl**6 - 4.7616115731E+02*cl**5 + 
                3.0471047478E+02*cl**4 - 4.9695074051E+01*cl**3 - 
                1.0038545848E+01*cl**2 + 2.9166961957E+01*cl -4.9622976565E+00)
                *0.0174528666666667)
    
    elif 1.1 < mach <= 1.15:
        archivo = open('angulos', "w")
        
        A1 = open('1_0', 'r')
        CL1 = []
        AoA11 = []
        for line in A1:
            s = line.strip()
            s = s.split()
            CL1.append(float(s[0]))  # CL .
            AoA11.append(float(s[1])* 0.0174528666666667)  # AoA (rad).
        A1.close()
        
        A2 = open('1_15', 'r')
        CL2 = []
        AoA12 = []
        for line in A2:
            s = line.strip()
            s = s.split()
            CL2.append(float(s[0]))  # CL .
            AoA12.append(float(s[1])* 0.0174528666666667)  # AoA (rad).
        A2.close()
        
        def AoA1 (cl):
            angulo_ataque = interp1d(CL1,AoA11)
            return float (angulo_ataque(cl))
        def AoA2 (cl):
            angulo_ataque= interp1d(CL2,AoA12)
            return float (angulo_ataque(cl))
        
#        Generación de 'archivo', donde se almacenaran los datos.
        
#        archivo.write('M \tAoA \n')
        archivo.write('%.8f\t' % 1.1)
        archivo.write('%.8f\n' % AoA1(cl))
        archivo.write('%.8f\t' % 1.15)
        archivo.write('%.8f\n' % AoA2(cl))
        archivo.close()
#        
        angulos= open ('angulos','r')
        MACH = []
        AoA =[]
        for line in angulos:
            s=line.strip()
            s=s.split()
            MACH.append(float(s[0]))
            AoA.append (float(s[1]))
        angulos.close()
        
        
        def angulo (mach):
            angulo=interp1d(MACH,AoA)
            return float (angulo(mach))
        return float (angulo(mach))


    elif 1.15 < mach <= 1.25:
        return ((1.7321073433E+02*cl**6 - 3.2583927645E+02*cl**5 + 
                1.8786276348E+02*cl**4 - 1.8403425545E+01*cl**3 - 
                1.0538135455E+01*cl**2 + 3.4710651662E+01*cl -3.3031797725E+00)
                *0.0174528666666667)


    elif 1.25 < mach <= 1.35:
        return ((1.0821613509E+01*cl**3 - 3.3780820230E+00*cl**2 + 
                3.5963102868E+01*cl - 3.2590164282E+00)*0.0174528666666667)


    elif 1.35 < mach < 1.45: 
        return(( 4.2413647959E+01*cl - 3.9988193650E+00)*0.0174528666666667)


    

    elif 1.45 < mach <= 1.55 :
        return ((4.3346692815E+01*cl - 4.9904748043E+00)*0.0174528666666667)
    
    elif 1.55 < mach <= 1.95:
        archivo = open('angulos', "w")
        
        A1 = open('1_55', 'r')
        CL1 = []
        AoA11 = []
        for line in A1:
            s = line.strip()
            s = s.split()
            CL1.append(float(s[0]))  # CL .
            AoA11.append(float(s[1])* 0.0174528666666667)  # AoA (rad).
        A1.close()
        
        A2 = open('1_95', 'r')
        CL2 = []
        AoA12 = []
        for line in A2:
            s = line.strip()
            s = s.split()
            CL2.append(float(s[0]))  # CL .
            AoA12.append(float(s[1])* 0.0174528666666667)  # AoA (rad).
        A2.close()
        
        def AoA1 (cl):
            angulo_ataque = interp1d(CL1,AoA11)
            return float (angulo_ataque(cl))
        def AoA2 (cl):
            angulo_ataque= interp1d(CL2,AoA12)
            return float (angulo_ataque(cl))
        
#       Generación de 'archivo', donde se almacenaran los datos.
        
#       archivo.write('M \tAoA \n')
        archivo.write('%.8f\t' % 1.55)
        archivo.write('%.8f\n' % AoA1(cl))
        archivo.write('%.8f\t' % 1.95)
        archivo.write('%.8f\n' % AoA2(cl))
        archivo.close()
       
        angulos= open ('angulos','r')
        MACH = []
        AoA =[]
        for line in angulos:
            s=line.strip()
            s=s.split()
            MACH.append(float(s[0]))
            AoA.append (float(s[1]))
        angulos.close()
        
        
        def angulo (mach):
            angulo=interp1d(MACH,AoA)
            return float (angulo(mach))
        return float (angulo(mach))


    elif 1.95 < mach <= 2.04 :
        return ((3.4378736432E+01*cl - 5.2463793871E+00)*0.0174528666666667)


def angulo_ataqueT(alfa_posible_radians, mach):
    '''En función de si el ángulo de ataque obtenido es inferior o superior al
    de pérdida, la función angulo_ataque devolverá un ángulo u otro.  Cuando se
    supera el ángulo de ataque de entrada en pérdida, la función angulo_ataque
    devuelve el valor del ángulo de entrada en pérdida para así no volar en
    pérdida.
    '''
    alfa_posible=degrees(alfa_posible_radians)
    
    if mach <= 1.11:
        angulo_perdida = 40
    elif 1.11 < mach <=1.60 :   
        angulo_perdida = (5.0004166396E+02*mach**4-2.8333247596E+03*mach**3 + 
                         5.9819002153E+03*mach**2- 5.6224859405E+03*mach + 
                         2.0265253584E+03)
    elif 1.60 < mach <=1.90:
        angulo_perdida=-333.33*mach**3+1750*mach**2-3081.7*mach+ 1832
    elif 1.90 < mach <= 2.04:
        angulo_perdida = -7.1429*mach + 21.571
    if alfa_posible < angulo_perdida:
        return radians(alfa_posible)
    return radians(angulo_perdida)


#
def cd0(mach):
    
    if 0.0 <= mach < 0.8:
        return 0.0277
    if 0.8 <= mach < 0.87:
        return 0.4514*mach**2 - 0.7183*mach + 0.3135
    if 0.87 <= mach < 1:
        return -0.0446*mach**2 + 0.2108*mach - 0.1194
    if 1 <= mach < 1.5:
        return -1.49629739*10*mach**6 + 1.13024477*100*mach**5 - 3.54316320*100*mach**4 + 5.90087652*100*mach**3 - 5.50747319*100*mach**2 + 2.73213176*100*mach - 5.62520143*10
    if 1.5 <= mach < 1.9:
        return 0.0478
    if 1.9 <= mach < 1.99:
        return 2.84217094*10**(-14)*mach**2 - 0.0190000000*mach + 0.0839000000
    if 1.99 <= mach < 2.15:
        return 0.0481596147*mach**3 - 0.294504644*mach**2 + 0.598685329*mach - 0.358556517



def k(mach):
    '''Coeficiente de resistencia inducida que multiplica al coeficiente
    de sustentación.
    '''
    fos = .005 * (1 + 1.5 * (ESTRECHAMIENTO - .6)**2)
    # Este valor es una función lambda que aparece dentro del factor de Oswald.
    e_mach = 1 / ((1 + .12 * mach**2) * (1 + (.1 * (3 * NE + 1)) / (4 + AR)
                                         + (.142 + fos * AR * (10
                                                               * ESPESOR)**.33)
                                         / (cos(FLECHA2)**2)))
    # Factor de Oswald.
    return 1 / (e_mach * pi * AR)


def cd_inducida(k_d, c_l):
    '''Coeficiente de resistencia inducida.
    '''
    return k_d * c_l**2
    

def coeficiente_resistencia_graficas(altura,mach):
    
    ''' Alturas entre '''
    if mach >=1:
        
        a6=1.6434305604E-08*altura**2-4.1753894557E-04*altura+2.4667379501E+00
        a5=-1.7843285584E-07*altura**2+4.5435417463E-03*altura-2.6938049554E+01
        a4=7.9946263261E-07*altura**2-2.0399304846E-02*altura+1.2134565805E+02
        a3=-1.8902382612E-06*altura**2+4.8319142703E-02*altura-2.8829687644E+02
        a2=2.4850612117E-06*altura**2-6.3613828398E-02*altura+ 3.8056171949E+02
        a1=-1.7214258098E-06*altura**2+ 4.4098853560E-02*altura-2.6439270360E+02
        a0=4.9108514143E-07*altura**2-1.2574888785E-02*altura+7.5550822290E+01
        
    if 0.6 <= mach < 1 :
        
        a6 =  0
        a5 = -6.0530781888E-06* altura**2 + 1.5976312518E-01 * altura - 1.0586225597E+03
        a4 = 2.8058834576E-05* altura**2 - 7.3354805804E-01* altura + 4.8087803012E+03
        a3 = -5.1927434527E-05* altura**2 + 1.3441857778E+00* altura - 8.7160937572E+03
        a2 = 4.7985137903E-05* altura**2 - 1.2293902513E+00* altura + 7.8839676228E+03
        a1 = -2.2169813302E-05* altura**2 + 5.6185609203E-01* altura - 3.5628281674E+03
        a0 = 4.1083349598E-06 * altura**2 - 1.0291060055E-01 * altura + 6.4508519732E+02

    return a6*mach**6+a5*mach**5+a4*mach**4+a3*mach**3+a2*mach**2+a1*mach+a0


def cd_interferencia(mach):
    '''Coeficiente de resistencia debido a la interferencia misil-avión.
    '''
    if mach <= 0.7:
        cd_i = 0
    if 0.7 < mach < 0.955:
        cd_i = (4.9704382e3*mach**6 - 2.5004431e4*mach**5 + 5.2386095e4*mach**4
                - 5.8503607e4*mach**3 + 3.6730317e4*mach**2 - 1.2291490e4*mach
                + 1.7127795e3)
    if 0.955 < mach <= 0.9655:
        cd_i = -1.6067616e1*mach**2 + 3.0869911e1*mach - 1.4799698e1
    if 0.9655 < mach < 0.99:
        cd_i = (1.963091e4*mach**4 - 7.7677873e4*mach**3 + 1.1526168e5*mach**2
                - 7.6013338e4*mach + 1.8798642e4)
    if 0.99 < mach <= 1.38:
        cd_i = 2.5938426e-2*mach**2 - 7.3453378e-2*mach + 6.6405634e-2
    if 1.38 < mach <= 2.1:
        cd_i = 0.0002*mach**2 - 0.0008*mach + 0.0151

    return cd_i



def resistencia(vel, dens, c_d):
    '''Fuerza aerodinámica de resistencia total.
    '''
    return .5 * dens * S_W * c_d * vel**2


def sustentacion(vel, dens, c_l):
    '''Fuerza aerodinámica de sustentación total.
    '''
    return .5 * dens * S_W * c_l * vel**2
