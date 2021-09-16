import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from scipy.ndimage.filters import gaussian_filter
from scipy import stats
from matplotlib import cm
from sklearn.linear_model import LinearRegression

def myplot(x, y,s, bins=1000):
    heatmap, xedges, yedges = np.histogram2d(x, y, bins=bins)
    heatmap = gaussian_filter(heatmap, sigma=s)

    extent = [xedges[0], xedges[-1], yedges[0], yedges[-1]]
    return heatmap.T, extent


# fig,ax=plt.subplots(1,2)
# A=np.loadtxt("el_jastrow.txt")
# ax[0].plot(A[:,1],A[:,0],'.k')
# ax[0].plot([0,np.max(A[:,1])],[np.sum(A[:,0])/len(A[:,0]),np.sum(A[:,0])/len(A[:,0])],'-g')
# ax[0].text(np.max(A[:,1]),np.sum(A[:,0])/len(A[:,0])-5,'%f'%(np.sum(A[:,0])/len(A[:,0])))
# ax[0].set_title("MO-LCAO con Jastrow")
# ax[0].set_xlabel("Iteraciones")
# ax[0].set_ylabel('Energía local (eV)')
# ax[0].grid()


# A=np.loadtxt("el_njastrow.txt")
# ax[1].plot(A[:,1],A[:,0],'.k')
# ax[1].plot([0,np.max(A[:,1])],[np.sum(A[:,0])/len(A[:,0]),np.sum(A[:,0])/len(A[:,0])],'-g')
# ax[1].text(np.max(A[:,1]),np.sum(A[:,0])/len(A[:,0])-10,'%f'%(np.sum(A[:,0])/len(A[:,0])))
# ax[1].set_title("MO-LCAO sin Jastrow")
# ax[1].set_xlabel("Iteraciones")
# ax[1].set_ylabel('Energía local (eV)')
# ax[1].grid()
# plt.show()


# fig,ax=plt.subplots(1,2)
# sigma = 0.1
# A=np.loadtxt("el_jastrow.txt")

# x = A[:,1]
# y = A[:,0]
# img, extent = myplot(x, y,sigma)
# img1=ax[0].imshow(img, extent=extent, vmin=0, vmax=5, origin='lower', cmap=cm.jet,aspect='auto')
# ax[0].set_xlabel("Iteraciones")
# ax[0].set_ylabel('Energía local (eV)')
# ax[0].set_title("MO-LCAO con Jastrow")
# B=np.loadtxt("el_njastrow.txt")

# x = B[:,1]
# y = B[:,0]
# img, extent = myplot(x, y,sigma)
# img2=ax[1].imshow(img, extent=extent, vmin=0, vmax=5, origin='lower', cmap=cm.jet,aspect='auto')
# ax[1].set_xlabel("Iteraciones")
# ax[1].set_ylabel('Energía local (eV)')
# ax[1].set_title("MO-LCAO sin Jastrow")

# plt.show()



'''DISTANCIAS'''

A=np.loadtxt("di_jastrow.txt")
plt.plot(A[:,1],A[:,0],'.k')
plt.plot([0,np.max(A[:,1])],[np.sum(A[:,0])/len(A[:,0]),np.sum(A[:,0])/len(A[:,0])],'-g')
plt.text(np.max(A[:,1]),np.sum(A[:,0])/len(A[:,0])-5,'%f'%(np.sum(A[:,0])/len(A[:,0])))
plt.title("MO-LCAO con Jastrow")
plt.xlabel(r"Distancia $r_{12}$ " r"$(\AA$)")
plt.ylabel('Energía local (eV)')
plt.grid()
plt.show()

A=np.loadtxt("di_njastrow.txt")
plt.plot(A[:,1],A[:,0],'.k')
plt.plot([0,np.max(A[:,1])],[np.sum(A[:,0])/len(A[:,0]),np.sum(A[:,0])/len(A[:,0])],'-g')
plt.text(np.max(A[:,1]),np.sum(A[:,0])/len(A[:,0])-10,'%f'%(np.sum(A[:,0])/len(A[:,0])))
plt.title("MO-LCAO sin Jastrow")
plt.xlabel(r"Distancia $r_{12}$ " r"$(\AA$)")
plt.ylabel('Energía local (eV)')
plt.grid()
plt.show()