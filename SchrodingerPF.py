import numpy as np
import matplotlib.pyplot as plt

#
'''--------------------------------------------------------------------------------------------------------------------
PROPÓSITO: Encontrar las soluciones estacionarias de la ecuación de Schrodinger para una partícula en un potencial
           V_p(x) mediante el método de diferencias finitas. Cada punto de la solución ocupa el valor de los dos
           valores anteriores, por los qe se utilizan dos valores iniciales cualesquiera en un extremo, donde se sabe
           que la función debe ser casi cero, con los que se encuentran todos los demás punto. Se normaliza la función
           de onda, ya que al tomar dos valores iniciales desconocidos, esta quedará escalada por algun factor descono-
           cido. Se encuentran las energías posibles mediante aproximaciónes, analizando el comportamiento de la solución
           en el otro extremo del potencial.

    NOMBRE: Hugo Alexis Torres Pasillas				FECHA: 27/Marzo/2018

'''#-------------------------------------------------------------------------------------------------------------------

def v_p(xv):                                  #Define el potencial en el que se busca la solución a la ecuación de 
#	yv=4*(xv-2.5)**2                          #Schrodinger
	
	if 1<=xv<=4:
		yv=0
	else:
		yv=20
	return yv

def Sch_equation(x_val,Psi_in,dx,En,V):       #Solución de la ecuación de onda mediante diferencias finitas. 
	                                          #Utiliza los valores de x (x_val) donde se encontrará la solución
	Ps=Psi_in                                 #de la ecuación, los dos puntos iniciales (Psi_in), a partir de los
	                                          #cuales se encontrarán los demás puntos, los intervalos dx, la energía
                                              #E, que se encontrará numéricamente, y el potencial V.
	for i in range(len(x_val)-2):               
		Ps[0].append(Ps[0][-1]+dx)


		Psi_sig=(2*(V(Ps[0][-1])-En)*dx**2+2)*Ps[1][i+1]-Ps[1][i]  #Calcula el punto siguiente, a partir de los anteriores
		Ps[1].append(Psi_sig)                 #Anexa el nuevo punto al vector de la solución de la ecuación

	return Ps                                 #Devuelve los puntos (x, Psi(x)) de la solución

def normalizar(Psin):                         #Ya que la solución encontrada utiliza dos puntos aproximados a la solución,
                                              #pero sin llegar a serlo, quedará escalada por algun factor, por lo que debe
                                              #ser normalizada.

	x=Psin[0]
	y=[i**2 for i in Psin[1]]                 #Cuadrado de la función, que debe ser integrada.
	dxn=x[2]-x[1]                             #Define dx a partir de la diferencia de dos valores de x en Psin
	Int=0                                     #Sumará los pequeños diferenciales de la integal

	for i in range(len(y)-1):                 #Integración numérica
		Int += 0.5*(y[i]+y[i])*(dx)

	Psin[1] = [k/np.sqrt(Int) for k in Psin[1]]  #Divide la solución entre la integral, para normalizarla.
	return Psin                               #Devuelve la misma función que se le pasa, pero normalizada.




L=5                                          
dx=0.01
GraphAxis=[0,5,0,21] 

#---------------Gráfica Potencial--------------------------------------->
xp=np.arange(-1,6,dx)
vp=[v_p(i) for i in xp]
xp2=[-3 for i in xp]
plt.plot(xp,vp)
plt.fill_between(xp,xp2,vp,color='gray')
plt.axis(GraphAxis)
plt.title(r'$-\frac{\hbar}{2m}\frac{\partial^2{\psi}}{\partial{x^2}}+U(x)\psi(x)=E\Psi(x)$', size=20, y=1.025)
plt.xlabel('x', size=16)
h=plt.ylabel('U(x)     ', size=16)
h.set_rotation(0)

#----------------------------------------------------------------------->


E_psb=list(np.arange(0,20,0.1))                #Define los valores de la energía con los que se encontrará la solución
                                               #a la ecuación de Schrodinger

x=list(np.arange(-2,L+2,dx))                   #Valores de x donde se busca la solución
Eo=0                                           #Comienza buscando la solución con esta energía.


E_vfun=[]                                      #Guarda los pares de valores de la energía que se aproximan más a la solución
                                               #Si E es una energía, entonces la función diverge hacia +infinito o -infinito
                                               #cuando la energía es más baja, pero cuando pasa a un valor más alto, la
                                               #solución diverge hacia el otro lado. Se utiliza este hecho para buscar los
                                               #valores de la energía donde la ecuación diverge hacia extremos distintos.


while Eo<20:                                   #Busca la solución de la ecuación para valores de energía entre 0 y 20 y
                                               #guarda el último punto (de la derecha) en el vector E_vfun, con los que 
                                               #se encuentra hacia dónde diverge la solución
	
	Psi_in=[x[0:2],[0.0,0.1]] 
	P=Sch_equation(x,Psi_in,dx,Eo,v_p)
#	P=normalizar(P)
	E_vfun.append(P[1][-1])                    #Pequeños pasos de aumento de la energía. 
	Eo += 0.1




E_par=[]                                       #Guarda los pares de puntos finales(a partir de E_vfun) donde este cambia 
                                               #de signo, ya que esto se ignifica que la solución con esas energías divergen
                                               #en direcciones contrarias.
for i in range(len(E_vfun)-1):

	if E_vfun[i]*E_vfun[i+1]<0:
		E_par.append([E_psb[i],E_psb[i+1]])

Ei_t=[]

for n in range(6):

	E0=E_par[n][0]                                 #Toma el n-ésimo par de energías, ya que es el nivel que se busca.
	E1=E_par[n][1]

	eps=1e-12                                      #error de la solución 

	MayMen=True                                    #(No es muy importante para el método)


	Psi_in=(x[0:2],[0.0,0.1])
	Psi=Sch_equation(x,Psi_in,dx,E0,v_p)


	if Psi[1][-1]<0:
		MayMen=False

	##------------------APROXIMACIÓN DE LA N-ÉSIMA ENERGÍA-------------------------------##
	'''
	A partir del n-ésimo par de energías donde se encontró que la solución divergía en direcciones contrarias, 
	se utilizan estos para mejorar el valor de la energía mediante aproximaciones similares al método de bisección.
	Termina cuando el último punto sea casi cero (con error de eps).
	'''
	if MayMen:
		while abs(E1-E0)>eps:
			Ep=0.5*(E0+E1)
			Psi_in=[x[0:2],[0.0,0.1]]
			Psi=Sch_equation(x,Psi_in,dx,Ep,v_p)

			if Psi[1][-1]>0:
				E0=Ep
			else:
				E1=Ep
	else:
		while abs(E1-E0)>eps:
			Ep=0.5*(E0+E1)
			Psi_in=[x[0:2],[0.0,0.1]]
			Psi=Sch_equation(x,Psi_in,dx,Ep,v_p)

			if Psi[1][-1]<0:
				E0=Ep
			else:
				E1=Ep

	Ei_t.append(Ep)

print(Ei_t)

##-----------------------------------------------------------------------------------##




##-----------------------------------------------------------------------------------##


##-----------------------------------------------------------------------------------##

for i in range(6):

	Ep=Ei_t[i]

	x=list(np.arange(0,L,dx))                                     #Valore donde se buscará la solución

	Psi_in=[x[0:2],[0.0,0.01]]                                    #Puntos iniciales, se toman fuera del potencial, donde 
                                                              #la función debe tomar valores muy cercanos al cero.

	Psi=Sch_equation(x,Psi_in,dx,Ep,v_p)                          #Llama la función Sch_equation, para encontrar la solución
                                                              #numérica, a partir de los valores de x, los puntos iniciales
                                                              #Psi_in, los intervalos dx, la n-ésima energía encontrada y
                                                              #la función del potencial.

	Psi=normalizar(Psi)                                           #Normaliza la función

	plt.plot(Psi[0], [I**1 + Ei_t[i] for I in Psi[1]],'b-')     #Grafica EL CUADRADO de la función de onda encontrada.



##-------------------------------------------------------------##

plt.twinx()
plt.axis(GraphAxis)
h=plt.ylabel("  E'= \n", color='red', size=16, y=0.6)
h.set_rotation(0)
h.set
plt.tick_params(axis='y', labelcolor='red')


x=list(np.arange(-2,L+3,1)) 

for j in range(6):
	y1=[Ei_t[j] for i in range(len(x))]
	plt.plot(x,y1, 'r--')

plt.show()