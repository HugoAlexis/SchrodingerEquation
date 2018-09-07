import matplotlib.pyplot as plt
import numpy as np 


def v_p(xv):                                  #Define el potencial en el que se busca la solución a la ecuación de 
	l=0
	Z=1
	if xv>0:

		yv=-Z/xv+(l*(l+1))/(2*xv**2)                         #Schrodinger
	else:
		yv=200
	

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

	return Ps      

def normalizar(Psin):                         #Ya que la solución encontrada utiliza dos puntos aproximados a la solución,
                                              #pero sin llegar a serlo, quedará escalada por algun factor, por lo que debe
                                              #ser normalizada.

	x=Psin[0]
	y=[i**2 for i in Psin[1]]                 #Cuadrado de la función, que debe ser integrada.
	dxn=x[2]-x[1]                             #Define dx a partir de la diferencia de dos valores de x en Psin
	Int=0                                     #Sumará los pequeños diferenciales de la integal

	for i in range(len(y)-1):                 #Integración numérica
		Int += 0.5*(y[i]+y[i])*(dx)

	Psin[1] = [k/np.sqrt(-Int) for k in Psin[1]]  #Divide la solución entre la integral, para normalizarla.
	
	print(type(Psin))

	#Psin[0].reverse()
	#Psin[1].reverse()
	return Psin                               #Devuelve la misma función que se le pasa, pero normalizada.

dx=-0.1
Rmax=100
N=int(abs(Rmax/dx))
ValR=list(np.linspace(Rmax, 0, N))
GraphAxis=[0,30,-0.55,0.02] 
#GraphAxis=[0,20,-0.1,0.01] 

#---------------Gráfica Potencial--------------------------------------->
#plt.subplots(facecolor='white')
xp=np.arange(-10,50,-dx)
vp=[v_p(i) for i in xp]
xp2=[-3 for i in xp]
plt.plot(xp,vp,color='gray')
plt.fill_between(xp,xp2,vp,color='gray')
plt.axis(GraphAxis)
#plt.title(r'$-\frac{\hbar}{2m}\frac{\partial^2{\psi}}{\partial{x^2}}+U(x)\psi(x)=E\Psi(x)$', size=20, y=1.025)
plt.title('Niveles energéticos del\n átomo de hidrógeno', size=18, y=1.015)
plt.xlabel('$r(a_o)$', size=16)
h=plt.ylabel('$V_{ef}(r)$', size=16, y=1.05,)
h.set_rotation(0)


#----------------------------------------------------------------------->



E_psb=list(np.arange(-0.6,0,0.001))                #Define los valores de la energía con los que se encontrará la solución
                                               #a la ecuación de Schrodinger

x=list(np.linspace(Rmax, 0, N))                   #Valores de x donde se busca la solución
Eo=-0.6                                           #Comienza buscando la solución con esta energía.


E_vfun=[]                                      #Guarda los pares de valores de la energía que se aproximan más a la solución
                                               #Si E es una energía, entonces la función diverge hacia +infinito o -infinito
                                               #cuando la energía es más baja, pero cuando pasa a un valor más alto, la
                                               #solución diverge hacia el otro lado. Se utiliza este hecho para buscar los
                                               #valores de la energía donde la ecuación diverge hacia extremos distintos.


while Eo<0:                                   #Busca la solución de la ecuación para valores de energía entre 0 y 20 y
                                               #guarda el último punto (de la derecha) en el vector E_vfun, con los que 
                                               #se encuentra hacia dónde diverge la solución
	
	Psi_in=[x[0:2],[0.0,0.01]] 
	P=Sch_equation(x,Psi_in,dx,Eo,v_p)
#	P=normalizar(P)
	E_vfun.append(P[1][-1])                    #Pequeños pasos de aumento de la energía. 
	Eo += 0.001

E_par=[]                                       #Guarda los pares de puntos finales(a partir de E_vfun) donde este cambia 
                                               #de signo, ya que esto se ignifica que la solución con esas energías divergen
                                               #en direcciones contrarias.
for i in range(len(E_vfun)-1):

	if E_vfun[i]*E_vfun[i+1]<0:
		E_par.append([E_psb[i],E_psb[i+1]])

Ei_t=[]

for n in range(len(E_par)):

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

##-------------------------------------------------------------##





plt.twinx()
plt.axis(GraphAxis)
h=plt.ylabel(" $E(E_h=27.2 eV)$", color='red', size=16, x=1, y=1.08)
h.set_rotation(0)
h.set
plt.tick_params(axis='y', labelcolor='red')


x=list(np.arange(0,51,1)) 


for j in range(6):
	y1=[Ei_t[j] for i in range(len(x))]
	leyenda="$E = $ " + str(Ei_t[j]) 
	leyenda = leyenda[0:18]  + "$E_h$"
	plt.plot(x,y1,label=leyenda)

plt.legend(loc="lower right")
plt.show()
