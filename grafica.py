import numpy as np
import matplotlib.pylab as plt

def minimo(A):
    j=0
    mine=np.min(A[:,1])
    for i in range(len(A)):
        p=A[i,1]
        j=j+1
        if p==mine:
            return j,mine



# fig,ax=plt.subplots(1,1)
# plt.suptitle("Funciones no correspondientes al estado fundamental de la molécula de hidrógeno")
# A=np.loadtxt("datos2.txt")
# ax.plot(A[:,0],A[:,1],'.-r',label=r"$\phi_{2}$")
# ax.set_xlabel("Distancia \t" r"$(\AA$)")



# ax.grid()

# A=np.loadtxt("datos3.txt")
# ax.plot(A[:,0],A[:,1],'.-k',label=r"$\phi_{3}$")





# A=np.loadtxt("datos4.txt")
# ax.plot(A[:,0],A[:,1],'.-b',label=r"$\phi_{4}$")

# ax.set_ylabel("Energía (eV)")
# plt.legend()
# plt.show()

# fig,ax=plt.subplots(1)
# plt.suptitle("Función del estado fundamental de la molécula de Hidrógeno")
# A=np.loadtxt("datos_lcao_nj.txt")
# a,b=minimo(A)
# ax.plot(A[:,0],A[:,1],'-k')
# ax.plot([A[a-1,0],A[a-1,0]],[A[a-1,1],0],'k-')
# ax.plot([A[0,0],A[-1,0]],[0,0],'k-')
# ax.text(A[a,0],b-0.2,'(%f , %f )' %(A[a-1,0],b))
# ax.set_xlabel("Distancia \t" r"$(\AA$)")
# ax.set_ylabel("Energía (eV)")
# ax.set_title(r"$\phi_{1} $")
# ax.set_xlim(0,4)
# ax.set_ylim(0,35)
# ax.grid()
# plt.show()



# AHORA PLOTEAMOS TODO JUNTO
# fig,ax=plt.subplots(1)
# plt.suptitle("Distancia frente a Energía para diferentes funciones de onda")
# A=np.loadtxt("datos_lcao_nj.txt")
# ax.plot(A[:,0],A[:,1],'-k',label="MO-LCAO" r"($\phi_{1} $)")

# B=np.loadtxt("datos_lcao_jas.txt")
# ax.plot(B[:,0],B[:,1],'-r',label="MO-LCAO con Jastrow" r"($\phi_{1}\phi_{j} $)")

# C=np.loadtxt("dfock.txt")
# ax.plot(C[:,0],C[:,1],'-b',label="Hartree-Fock")

# D=np.loadtxt("datos_combinacion.txt")
# ax.plot(D[:,0],D[:,1],'-g',label="Combinación" r"($\phi_{1}+ \lambda \phi_{2} $)")




# ax.set_xlim(0.3,1.5)
# ax.set_xlabel("Distancia \t" r"$(\AA$)")
# ax.set_ylabel("Energía (eV)")
# plt.grid()
# plt.legend()
# plt.show()



# AQUI VAMOS A PONER LAS ENERGIAS DE CADA TERMINO DEL HAMILTONIANO
# CON LA DISTANCIA

# A=np.loadtxt("energias_jastrow.txt")
# plt.title("Energía de los diferentes términos del Hamiltoniano para la función del estado fundamental con término de Jastrow y una distancia de equilibrio R=0.74")
A=np.loadtxt("energias_njastrow.txt")
plt.title("Energía de los diferentes términos del Hamiltoniano para la función del estado fundamental sin término de Jastrow y una distancia de equilibrio R=0.74")
plt.plot(A[:,0],A[:,1],'.',label="Términos cinéticos")
# plt.plot(A[:,0],A[:,2],'.',label="electron1nucleoA")
# plt.plot(A[:,0],A[:,3],'.',label="electron1nucleoB")
# plt.plot(A[:,0],A[:,4],'.',label="electron2nucleoA")
# plt.plot(A[:,0],A[:,5],'.',label="electron2nucleoB")
plt.plot(A[:,0],A[:,6],'.',label="Término \t"  r"$\frac{1}{|r_{12}|}$")
# plt.plot(A[:,0],A[:,7],'.',label="repulsión nuclear")
plt.xlabel("Distancia \t" r"$(\AA$)")
plt.ylabel("Energía (eV)")
plt.legend()
print(A[0,:])
plt.show()
