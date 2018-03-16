from math import radians, cos, sin, degrees, pi, exp

from modeloISA import density, temperature, GAMMA, viscosity, R_AIR
from gravedad import gravity
from aero_misil import cdll, SREF_MISIL, SGASES
from avion import velocidad, altura, psi_list, fi_list

#--------------------------CONDICIONES GRAVITATORIAS--------------------------
G = 6.673e-11  # Constante de gravitación universal (N m2/kg2).
MT = 5.972e24  # Masa terrestre (kg).
MU = G * MT
RT = 6378136.3  # Radio terrestre (m).
GRAV = MU / RT**2  # Aceleración de la gravedad a nivel del mar (m/s2).

'''Aqui empieza el programa'''
#archivo = open("LanzBalistico_ejesnuevos(modelo7)" ,"w")
#archivo.write('Time (s)  \tAltura lanzamiento  \tAngulo lanz \tAltura \tVelocidad final \tVelocidad perdida \tmasa total \ttheta  \tfil\n')
  
fi_grados = fi_list
psi_grados = psi_list 


for gasto in[14]:
    gastotexto = str(gasto)
        
    archivo = open(gastotexto ,"w")
    archivo.write('Time (s)  \tAltura lanzamiento  \tAngulo lanz \tAltura \tVelocidad final \tVelocidad perdida \tmasa total \ttheta  \tfil\n')  
   
    for i in range(134, 156, 1):
       

    
        gasto_etapa2= gasto 
        isp1 = 280
        Empuje_misil = gasto * isp1 * 9.81
        
        ratio_estructural = .05 #'Relacion estrutural, sera un 5% del peso del propulsante de etapa
        masa_maniobras = 0
        masa_util = 30 + masa_maniobras
        
        '''Velocidades implicadas'''
        v_orb = 7800
        v_loss = 782
        v_loss_empirica = 0
        v_rotacional = 400
        error = 800 
        contador = 1
        
            #-----------------------------------------------------------------VALORES DE VARIABLES INICIALES-------------------------------------------------------------------
        '''Calculos de la optimizacion'''
        while error > 10 and contador <= 50: #calcula la diferencia netre la v_loss claculada y la anterior
        
            '''Inicialización de variables y diferenciales para la maniobra del misil'''
            z0 = altura [i]                                        # DATOS DE 60º
            v0 = velocidad[i]                                      # Velocidad inicial con la que se lanza
            fil_grados0 = fi_grados[i]                                 # ángulo de ataque de la velocidad sobre la vertical local, estamos cogiendo el dato del avión
            psil_grados0 = psi_grados[i]                                  # ángulo de la vertical local en ejes centro tierra, estamos cogiendo el dato del avión
            theta_grados0 = 90 - fil_grados0 - psil_grados0
            
             
            zl = z0
            fil_grados = fil_grados0
            psil_grados = psil_grados0
            
            fil = radians (fil_grados)                                                                                  # Conversión a grados del ángulo fil                                                                                                 # Ángulo inicial de vertical local 
            psil = radians (psil_grados)                                                                                # Conversión a grados de psil
            rl = RT + zl
            vl = v0 + v_rotacional                                                                                      # Añadimos a la velocidad inicial la rotacional
            Ml = vl / ((GAMMA * R_AIR * temperature(zl))**0.5)
            vxl = vl * sin(fil)                                                                                         # Inicialización de la componente horizontal de velocidad  
            vyl = vl * cos(fil)                                                                                         # Inicialización de la componente verical de velocidad 
            tl = 0                                                                                                      # Inicialización temporal 
            xl = 0                                                                                                      # Iniciamos la variable x de desplazamiento del misil como el que lleva el avion en el punto de lanzamiento
            sl = 0 #Inicialización del arco recorrido
            dvl = 0 
            dsl = 0 #Inicialización del diferncial del arco recorrido
            dxl = 0 #Inicialización del diferencial de la posición
            dzl = 0                                                                                                     # Inicialización del diferencial de la altura del misil 
            dfil = 0                                                                                                    # Inicialización de la variación del ángulo fi inicial
            dpsil = 0
            dtl = 0.01                                                                                                   # Inicialización de la variación del ángulo de psi inicial 
            
            Cdl = cdll(Ml,zl)
            D_misil = 0.5 * density(zl) * Cdl * SREF_MISIL * vl**2
            Emec_misil = 1000 * (gravity(zl) * zl + (vl ** 2) / 2)
            
            g = gravity (zl)# Aceleración de la gravedad (m/s2).
            
            velocidad_ideal = v_orb - v_rotacional + v_loss - v0
            f = ((1 + ratio_estructural)/exp(velocidad_ideal / (2 * isp1 * gravity(zl)))) - ratio_estructural
            
            masa_total = masa_util / (f**2)
            masa_etapa2 = masa_util / f
            masa_propulsante_etapa1 = (masa_total - masa_etapa2) / (1 + ratio_estructural)
            masa_propulsante_etapa2 = (masa_etapa2 - masa_util) / (1 + ratio_estructural)
            masa_estructura1 = masa_propulsante_etapa1 * ratio_estructural
            masa_estructura2 = masa_propulsante_etapa2 * ratio_estructural
            
            t_combustion= masa_propulsante_etapa1/gasto
            t_combustion_2_etapa=masa_propulsante_etapa2/gasto_etapa2 
            retardo_encendido = 0
            t_fin_combustion2= t_combustion + t_combustion_2_etapa + retardo_encendido
            
            masa_misil = masa_total #esta sera la masa que ira variando
            
            v_loss_empirica = 0
            
            
            
            while tl<= t_combustion and fil < pi / 2:
                     
            
                tl = tl + dtl # Evolución temporal (s)
                xl = xl + dxl # Posición horizontal (m)
                zl = zl + dzl # Altitud (m) 
                sl = sl + dsl # Arco recorrido (m)
                rl = RT + zl  # Altura al centro de la Tierra
                
                fil = fil + dfil                                                                # Variación del ángulo de asiento en vertical local 
                fil_grados = degrees(fil)                                                       # Conversión a grados de fil 
                psil = psil + dpsil                                                             # Variación del ángulo de la vertical local 
                psil_grados = degrees(psil)                                                     # Conversión a grados de psil 
                
                vl = vl + dvl                                                                   # Variación de la velocidad del misil 
                
                vxl = vl * sin(fil)                                                             #Componente horizontal de la velocidad (m/s)
                vyl = vl * cos(fil)                                                             #Componente vertical de la velocidad (m/s)
                
                rhol = density(zl)  # Densidad (kg/m3).
                Tl = temperature(zl)  # Temperatura (K).
                Mu_viscl = viscosity(zl) # Viscosidad 
                gl = gravity(zl) # Nueva gravedad
                
                Ml = vl / ((GAMMA * R_AIR * Tl)**0.5) #Mach de vuelo
            
                Ratio_areas = (SREF_MISIL - SGASES) / SREF_MISIL        
                
                Cdl = cdll(Ml, zl)   #Función que calcula el coeficiente de resistencia                                
                D_misil = 0.5 * rhol * Cdl * SREF_MISIL * vl**2 #Fuerza de resistencia (N)
                
                dvl = dtl / masa_misil * (Empuje_misil - D_misil - gl * masa_misil * cos(fil)) 
                masa_misil=masa_misil-gasto*dtl
                
                dfil = gl / vl * dtl                                                  # Diferencial del ángulo de asiento
                dpsil = vxl / rl * dtl                                                # Diferencial del ángulo de vertical local 
                
                dv_loss = dtl / masa_misil * ( D_misil + gl * masa_misil * cos(fil))
                
                dxl = vxl * dtl #Diferencial de la posición (m)
                dzl = vyl * dtl #Diferencial de la altitud (m)
                dsl = vl * dtl #Diferencial del arco recorrido (m)
            
                v_loss_empirica = v_loss_empirica + dv_loss
                            
            masa_misil  = masa_misil - masa_estructura1
                
                
            while t_fin_combustion2>=tl>t_combustion and fil < pi / 2 and zl<500000:
                    
            
                tl = tl + dtl # Evolución temporal (s)
                xl = xl + dxl # Posición horizontal (m)
                zl = zl + dzl # Altitud (m) 
                sl = sl + dsl # Arco recorrido (m)
                rl = RT + zl  # Altura al centro de la Tierra
                
                fil = fil + dfil                                                                # Variación del ángulo de asiento en vertical local 
                fil_grados = degrees(fil)                                                       # Conversión a grados de fil 
                psil = psil + dpsil                                                             # Variación del ángulo de la vertical local 
                psil_grados = degrees(psil)                                                     # Conversión a grados de psil 
                
                vl = vl + dvl                                                                   # Variación de la velocidad del misil 
                
                vxl = vl * sin(fil)                                                             #Componente horizontal de la velocidad (m/s)
                vyl = vl * cos(fil)                                                             #Componente vertical de la velocidad (m/s)
                
                rhol = density(zl)  # Densidad (kg/m3).
                Tl = temperature(zl)  # Temperatura (K).
                Mu_viscl = viscosity(zl) # Viscosidad 
                gl = gravity(zl) # Nueva gravedad
                
                Ml = vl / ((GAMMA * R_AIR * Tl)**0.5) #Mach de vuelo
            
                Ratio_areas = (SREF_MISIL - SGASES) / SREF_MISIL        
                
                Cdl = cdll(Ml, zl)   #Función que calcula el coeficiente de resistencia                                
                D_misil = 0.5 * rhol * Cdl * SREF_MISIL * vl**2 #Fuerza de resistencia (N)
                
                dvl = dtl / masa_misil * (Empuje_misil - D_misil - gl * masa_misil * cos(fil)) 
                masa_misil=masa_misil-gasto*dtl
                
                dfil = gl / vl * dtl                                                  # Diferencial del ángulo de asiento
                dpsil = vxl / rl * dtl                                                # Diferencial del ángulo de vertical local 
                
                dv_loss = dtl / masa_misil * ( D_misil + gl * masa_misil * cos(fil))
                
                dxl = vxl * dtl #Diferencial de la posición (m)
                dzl = vyl * dtl #Diferencial de la altitud (m)
                dsl = vl * dtl #Diferencial del arco recorrido (m)
            
                v_loss_empirica = v_loss_empirica + dv_loss
                
            
            masa_misil = masa_misil - masa_estructura2
            
            
            while fil < pi / 2 and zl<500000:  
                tl = tl + dtl # Evolución temporal (s)
                xl = xl + dxl # Posición horizontal (m)
                zl = zl + dzl # Altitud (m) 
                sl = sl + dsl # Arco recorrido (m)
                rl = RT + zl  # Altura al centro de la Tierra
                
                fil = fil + dfil                                                                # Variación del ángulo de asiento en vertical local 
                fil_grados = degrees(fil)                                                       # Conversión a grados de fil 
                psil = psil + dpsil                                                             # Variación del ángulo de la vertical local 
                psil_grados = degrees(psil)                                                     # Conversión a grados de psil 
                theta_grados = 90 - fil_grados - psil_grados
                vl = vl + dvl                                                                   # Variación de la velocidad del misil 
                
                vxl = vl * sin(fil)                                                             #Componente horizontal de la velocidad (m/s)
                vyl = vl * cos(fil)                                                             #Componente vertical de la velocidad (m/s)
                
                rhol = density(zl)  # Densidad (kg/m3).
                Tl = temperature(zl)  # Temperatura (K).
                Mu_viscl = viscosity(zl) # Viscosidad 
                gl = gravity(zl) # Nueva gravedad
                
                Ml = vl / ((GAMMA * R_AIR * Tl)**0.5) #Mach de vuelo
            
                Ratio_areas = 1       
                
                Cdl = cdll(Ml, zl)   #Función que calcula el coeficiente de resistencia                                
                D_misil = 0.5 * rhol * Cdl * SREF_MISIL * vl**2 #Fuerza de resistencia (N)
                
                dvl = dtl / masa_misil * ( - D_misil - gl * masa_misil * cos(fil)) 
                
                dfil = gl / vl * dtl                                                  # Diferencial del ángulo de asiento
                dpsil = vxl / rl * dtl                                                # Diferencial del ángulo de vertical local 
                
                dv_loss = - dvl
                
                dxl = vxl * dtl #Diferencial de la posición (m)
                dzl = vyl * dtl #Diferencial de la altitud (m)
                dsl = vl * dtl #Diferencial del arco recorrido (m)
            
                v_loss_empirica = v_loss_empirica + dv_loss      
                if zl > 500000:
                    break
     
    
            error = abs(v_loss_empirica - v_loss)
            v_loss = v_loss_empirica
                       
            if contador == 50:
                print ('error de interacion en el ángulo ', fil_grados0, 'con un error de ', error,' m/s')
            contador = contador + 1
            
        print('para ', theta_grados0,'la masa del misil es ', masa_total, ' la altura es ' , zl,' y el angulo de llegada', fil_grados)
        
        print('El gasto es: ', gasto, 'el numeor de la lista es ', i)
        
        archivo.write('%.8f\t' %tl)
        archivo.write('%.8f\t' %z0)
        archivo.write('%.8f\t' %theta_grados0)
        archivo.write('%.8f\t' %zl) 
        archivo.write('%.8f\t' %vl) 
        archivo.write('%.8f\t' %v_loss)
        archivo.write('%.8f\t' %masa_total)
        archivo.write('%.8f\n' %fil_grados)

    archivo.close()
