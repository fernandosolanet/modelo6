# -*- coding: utf-8 -*-
"""
Created on Wed May  9 15:40:10 2018

@author: jorge.sanz
"""

import numpy as np

from modelo_msise00 import temperature,density,R_AIR,GAMMA
from aero_avion import S_W,cd_interferencia
from modelo_empuje import thrust

def coeficiente_resistencia(altura,mach):
    
    ''' Alturas entre '''
    if mach >= 1:
        a6=1.6434305604E-08*altura**2-4.1753894557E-04*altura+2.4667379501E+00
        a5=-1.7843285584E-07*altura**2+4.5435417463E-03*altura-2.6938049554E+01
        a4=7.9946263261E-07*altura**2-2.0399304846E-02*altura+1.2134565805E+02
        a3=-1.8902382612E-06*altura**2+4.8319142703E-02*altura-2.8829687644E+02
        a2=2.4850612117E-06*altura**2-6.3613828398E-02*altura+ 3.8056171949E+02
        a1=-1.7214258098E-06*altura**2+ 4.4098853560E-02*altura-2.6439270360E+02
        a0=4.9108514143E-07*altura**2-1.2574888785E-02*altura+7.5550822290E+01
        
        
    if 0.6 <= mach < 1:
        a6 = 0
        a5 = -6.0530781888E-06* altura**2 + 1.5976312518E-01 * altura - 1.0586225597E+03
        a4 = 2.8058834576E-05* altura**2 - 7.3354805804E-01* altura + 4.8087803012E+03
        a3 = -5.1927434527E-05* altura**2 + 1.3441857778E+00* altura - 8.7160937572E+03
        a2 = 4.7985137903E-05* altura**2 - 1.2293902513E+00* altura + 7.8839676228E+03
        a1 = -2.2169813302E-05* altura**2 + 5.6185609203E-01* altura - 3.5628281674E+03
        a0 = 4.1083349598E-06 * altura**2 - 1.0291060055E-01 * altura + 6.4508519732E+02
        
    return a6*mach**6+a5*mach**5+a4*mach**4+a3*mach**3+a2*mach**2+a1*mach+a0



def empuje_necesario(rho,T,cdtotal,mach):
    
    a=(GAMMA*R_AIR*T)**0.5

    return  (0.5*rho*(mach*a)**2*S_W*cdtotal)
    
archivo = open('lanzamientos.txt', "w")

alturas = np.linspace(10668, 15000, 31)

for alt in range(0, len(alturas)):    
    altitud =  alturas[alt]
    rho=density(altitud)
    T=temperature(altitud)    
    
    M = 0.6
    Dist = 10
    
    mach_list =[]
    
    while M < 2.04:
        rho=density(altitud)
        T=temperature(altitud)
        M = M + 0.00001
        CD=coeficiente_resistencia(altitud,M)+0.35*cd_interferencia(M)
        a=(GAMMA*R_AIR*T)**0.5
        E_Nec=empuje_necesario(rho,T,CD,M)
        E_Disp=thrust(M,rho)
        Dist=abs(E_Nec-E_Disp)
        
        if Dist < 0.5:
            mach_list.append(M)
        
    if (len(mach_list)==0):
        
        break
        
    else:
        mach_max = max(mach_list)
        

    
    archivo.write('%.8f\t' % altitud)    
    archivo.write('%.8f\n' % mach_max) 

    print(mach_max)
    print(alt)

archivo.close()
