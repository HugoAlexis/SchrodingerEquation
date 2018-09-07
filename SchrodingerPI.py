import numpy as np
import matplotlib.pyplot as plt

#------------------------------------------------------------------------------------------------#
#    PROPÓSITO: Resolver la ecuación de Schrodinger para los estados estacionarios de una        #
#				 partícula en un potencial infinito mediante el método de diferencias finitas    #
# 				 y encontrar las energías permitidas para los dichoes estados de la función      #
#    			 de onda.                                                                        #
#------------------------------------------------------------------------------------------------#
#-----------------------------------------------------------------------------------------------/


#------------------------------------------------------------------------------------------------#
#--------------------------------FUNCIONES UTILIZADAS--------------------------------------------#
#------------------------------------------------------------------------------------------------#

def v_p(x):

	'''
	Potencial en el que se encuentra la partícula: caja rectangular con  paredes infinitas en 
	x=0 y x=L. 

	'''


	y=200
	if 0<=x<=L:
		y=0



	return y

def Sch_equation(x_val,Psi_in,dx, En):

	'''
	Encuentra la solución de la ecuación de Schrodinger para alguna energía 'En' cualquiera, a 
	partir de los dos valores iniciales Psi_in, en el intervalo x_val con subdiviciones de 
	tamaños dx. Devuelve una matriz con los valores [x, Psi(x)]
	'''

	Psi=Psi_in

	for i in range(len(x)-2):
		Psi[0].append(Psi[0][-1]+dx)


		Psi_sig=(-2*En*dx**2+2)*Psi[1][i+1]-Psi[1][i]
		Psi[1].append(Psi_sig)

	return Psi

def normalizar(Psi):

	'''
	Ya que cualquier múltiplo de una función Psi(x) que es solución de la ecuación de Schrodinger
	también lo es, y la solución de la función Sch_equation por lo general no está normalizada, 
	normaliza la función de onda encontrada Psi que se le pasa, de manera que la integral de esta
	en todo el espacio sea igual a 1..
	'''

	x=Psi[0]
	y=[i**2 for i in Psi[1]]
	dx=x[2]-x[1]
	Int=0

	for i in range(len(y)-1):
		Int += 0.5*(y[i]+y[i+1])*(dx)

	Psi[1] = [k/np.sqrt(Int) for k in Psi[1]]
	return Psi




L=2.5
dx=0.01
GraphAxis=[-2,L+2,0,25]

#------------------------------------------------------------------------------------------------#
#--------------------------------GRÁFICA DEL POTENCIAL-------------------------------------------#
#------------------------------------------------------------------------------------------------#

xp=np.arange(-5,L+5,dx)
vp=[v_p(i) for i in xp]
xp2=[-1 for i in xp]
plt.plot(xp,vp,color='gray')
plt.fill_between(xp,xp2,vp,color='gray')
plt.axis(GraphAxis)
plt.title(r'$-\frac{\hbar}{2m}\frac{\partial^2{\psi}}{\partial{x^2}}+U(x)\psi(x)=E\Psi(x)$', size=20, y=1.025)
plt.xlabel('x', size=16)
h=plt.ylabel('U(x) ', size=16, y=1.02)
h.set_rotation(0)


#------------------------------------------------------------------------------------------------#
#-------------------------------DETERMINACIÓN DE LAS ENERGÍAS------------------------------------#
#------------------------------------------------------------------------------------------------#

L=L+dx

E_psb=list(np.arange(0,50,0.1))                                               #--Primeras energías posibles con las que se encontrará la solución a la
																			  #ecuaión. Se guardarán los intervalos donde haya una energía permitida.
E_find=[]                                                                     #--Intervalos donde debe haber una energía permitida.
PsiUlt=[]																	  #--Últimos valores de la ec. De Schrodinger, con la que se determina en qué
																			  #--Intervalos hay energías posibles

Energias=[]                                      						      #Guarda las últimas arpoximaciónes de las energías

for i in range(len(E_psb)):

	'''
	Encuentra la solución de la ec. de Schrodinger para cada energía en E_psb, y guarda el último
	valor de cada solución en PsiUlt.
	'''

	En=E_psb[i]
	x=list(np.arange(0,L,dx))
	Psi_in=[x[0:2],[0.0,0.1]]
	Psi=Sch_equation(x,Psi_in,dx, En)

	PsiUlt.append(Psi[1][-1])


for i in range(len(PsiUlt)-1):
	'''
	A partir de los valores de PsiUlt (últimos valores de la solución para cada energía), se derminan
	los intérvalos donde se debe encontrar una energía permitida (cuando el último punto de la solución 
	cambia de signo).
	'''
	if PsiUlt[i]*PsiUlt[i+1]<0:
		E_find.append([E_psb[i], E_psb[i+1]])


for i in range(len(E_find)):

	'''
	Se mejora la aproximación de cada energía en E_find, reduciendo los intervalos por un método similar
	al método de bisección, analizando en qué intervalos de E el ultimo punto de la solución cambia de signo.
	'''
	E0=E_find[i][0]
	E1=E_find[i][1]


	while abs(E0-E1)>0.000001:

		Em=0.5*(E0+E1)


		Psi_in=[x[0:2],[0.0,0.1]]
		Psi0=Sch_equation(x,Psi_in,dx, E0)

		Psi_in=[x[0:2],[0.0,0.1]]
		Psim=Sch_equation(x,Psi_in,dx, Em)

		Psi_in=[x[0:2],[0.0,0.1]]
		Psi1=Sch_equation(x,Psi_in,dx, E1)
		
		if (Psi0[1][-1]*Psim[1][-1]<0):
			
			E1=Em
		else:
			E0 = Em	


	Energias.append([E0,E1])													  #Guarda el intervalo de energías cuando la diferencia entre
																				  #los valores de E0 y E1 es menor a la tolerancia deseada.


E=[]                                                                              #Finalmente, se guarda el valor medio de cada intervalo en 
																				  #'Energías'. 
for i in range(len(Energias)):
	E.append((Energias[i][0]+Energias[i][1])/2)


#------------------------------------------------------------------------------------------------#
#--------------------------------GRÁFICA DE LAS FUNCIONES DE ONDA--------------------------------#
#------------------------------------------------------------------------------------------------#

	'''
	A partir de las energías encontradas (almacenadas en 'E'), se encuentra la solución de la ecuación
	y se grafica en el mismo gráfico del potencial.
	'''

for n in range(len(E)):

	En=E[n]


	x=list(np.arange(0,L,dx))
	Psi_in=[x[0:2],[0.0,0.1]]

	Psi=Sch_equation(x,Psi_in,dx, En)

	Psi=normalizar(Psi)

	plt.plot(Psi[0], [I**1 + En for I in Psi[1]],'b-')

#------------------------------------------------------------------------------------------------#
#--------------------------------GRÁFICA DE LOS NIVELES ENERGÉTICOS------------------------------#
#------------------------------------------------------------------------------------------------#

	'''
	Se grafican en el mismo gráfico los niveles energéticos mediante líneas horizontales, añadiendo
	un nuevo eje y en el lado derecho.
	'''
x=list(np.arange(-2,5,dx))

plt.twinx()
plt.axis(GraphAxis)
h=plt.ylabel("  $E'=\\frac{mE}{\hbar^2}$ ", color='red', size=20, y=1.09, x=2)
h.set_rotation(0)
h.set
plt.tick_params(axis='y', labelcolor='red')


for n in range(len(E)):

	y=[E[n] for i in range(len(x))]

	plt.plot(x,y,'r--')



plt.show()