import numpy as np
import matplotlib.pylab as plt
from numba import jit
# la funcion va a ser x*ln(x) de 0 a 1 que de resultado da -0.25

h= lambda z,u: (z**2)*np.log(u)

n=20000
b=1
a=0

x=np.zeros(n)
y=np.zeros(n)
resu=np.zeros(n)
resu[0]=0
acu=np.zeros(n)
acu[0]=0
el=0

for i in range(1,n):

    x[i]=np.random.uniform(a,b)
    y[i]=np.random.uniform(a,b)
    
    resu[i]=(h(x[i],y[i]))
    el=el+h(x[i],y[i])
    acu[i]=el/i
   
print((b-a)*np.mean(resu))


fig,ax=plt.subplots(1)
# plt.suptitle((b-a)*np.mean(resu))
# ax[0].plot([x[:],x[:],],[np.zeros(n),h(x[:])],'r-')
# ax[0].plot(np.linspace(a,b,100),h(np.linspace(a,b,100)))
# ax[0].grid()
# ax[1].plot((b-a)*acu)
# ax[1].plot([0,n],[(b-a)*np.mean(resu),(b-a)*np.mean(resu)],'r--')
# ax[1].grid()
# plt.show()


    