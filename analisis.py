# -*- coding: utf-8 -*-
"""
Created on Wed Mar 24 12:28:01 2021

@author: yo
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from scipy.ndimage.filters import gaussian_filter
from scipy import stats
from matplotlib import cm



# A=np.loadtxt("electrones_0,65.txt")


# fig=plt.figure()
# ax=fig.add_subplot(1,1,1,projection='3d')
# ax.scatter(A[:,0],A[:,1],A[:,2],color='red',cmap='jet')
# ax.scatter(A[:,3],A[:,4],A[:,5],color='blue',cmap='jet')
# ax.plot([0,0,0],[0,0,0],'go')
# ax.plot([1,0,0],[0,0,0],'go')
# ax.set_xlabel('X')
# ax.set_ylabel('Y')
# ax.set_zlabel('Z')

# fig=plt.figure()
# ax=fig.add_subplot(1,1,1)
# ax.scatter(A[:,0],A[:,1],c=np.abs(A[:,2]),marker='.',cmap=cm.coolwarm)
# ax.scatter(A[:,3],A[:,4],c=np.abs(A[:,5]),marker='.',cmap=cm.coolwarm)
# ax.plot(0,0,'go')
# ax.plot(0.65,0,'go')
# ax.set_xlabel('X')
# ax.set_ylabel('Y')

# fig.savefig("probando.svg")







# def myplot(x, y, s, bins=1000):
    # heatmap, xedges, yedges = np.histogram2d(x, y, bins=bins)
    # heatmap = gaussian_filter(heatmap, sigma=s)

    # extent = [xedges[0], xedges[-1], yedges[0], yedges[-1]]
    # return heatmap.T, extent


# fig, axs = plt.subplots(2, 2)


# x = A[:,0]
# y = A[:,1]
# x2= A[:,3]
# y2= A[:,4]
# xf=np.concatenate((x,x2))
# yf=np.concatenate((y,y2))
# sigmas = [0, 16]

# for ax, s in zip(axs.flatten(), sigmas):
    # if s == 0:
        # ax.plot(xf, yf, 'k.', markersize=5)
        # ax.set_title("Scatter plot")
    # else:
        # img, extent = myplot(xf, yf, s)
        # ax.imshow(img, extent=extent, origin='lower', cmap=cm.jet)
        # ax.set_title("Smoothing with  $\sigma$ = %d" % s)

# plt.show()

# def myplot(x, y, s, bins=1000):
    # heatmap, xedges, yedges = np.histogram2d(x, y, bins=bins)
    # heatmap = gaussian_filter(heatmap, sigma=s)

    # extent = [xedges[0], xedges[-1], yedges[0], yedges[-1]]
    # return heatmap.T, extent

def myplot(x, y,s, bins=1000):
    heatmap, xedges, yedges = np.histogram2d(x, y, bins=bins)
    heatmap = gaussian_filter(heatmap, sigma=s)

    extent = [xedges[0], xedges[-1], yedges[0], yedges[-1]]
    return heatmap.T, extent



# fig,ax=plt.subplots(2,3)
# sigma = 16
# A=np.loadtxt("electrones_MOLCAO_0,5.txt")
# plt.suptitle("Forma de los orbitales para distintas distancias")
# x = A[:,0]
# y = A[:,1]
# x2= A[:,3]
# y2= A[:,4]
# xf=np.concatenate((x,x2))
# yf=np.concatenate((y,y2))
# img, extent = myplot(xf, yf,sigma)
# ax[0,0].imshow(img, extent=extent, origin='lower', cmap=cm.jet)
# ax[0,0].set_title("r=0,5")

# A=np.loadtxt("electrones_MOLCAO_0,65.txt")
# x = A[:,0]
# y = A[:,1]
# x2= A[:,3]
# y2= A[:,4]
# xf=np.concatenate((x,x2))
# yf=np.concatenate((y,y2))
# img, extent = myplot(xf, yf,sigma)
# ax[0,1].imshow(img, extent=extent, origin='lower', cmap=cm.jet)
# ax[0,1].set_title("r=0,65")

# A=np.loadtxt("electrones_MOLCAO_0,74.txt")
# x = A[:,0]
# y = A[:,1]
# x2= A[:,3]
# y2= A[:,4]
# xf=np.concatenate((x,x2))
# yf=np.concatenate((y,y2))
# img, extent = myplot(xf, yf,sigma)
# ax[0,2].imshow(img, extent=extent, origin='lower', cmap=cm.jet)
# ax[0,2].set_title("r=0,74")

# A=np.loadtxt("electrones_MOLCAO_0,9.txt")
# x = A[:,0]
# y = A[:,1]
# x2= A[:,3]
# y2= A[:,4]
# xf=np.concatenate((x,x2))
# yf=np.concatenate((y,y2))
# img, extent = myplot(xf, yf,sigma)
# ax[1,0].imshow(img, extent=extent, origin='lower', cmap=cm.jet)
# ax[1,0].set_title("r=0,9")

# A=np.loadtxt("electrones_MOLCAO_1,5.txt")
# x = A[:,0]
# y = A[:,1]
# x2= A[:,3]
# y2= A[:,4]
# xf=np.concatenate((x,x2))
# yf=np.concatenate((y,y2))
# img, extent = myplot(xf, yf,sigma)
# ax[1,1].imshow(img, extent=extent, origin='lower', cmap=cm.jet)
# ax[1,1].set_title("r=1,5")

# A=np.loadtxt("electrones_MOLCAO_2,3.txt")
# x = A[:,0]
# y = A[:,1]
# x2= A[:,3]
# y2= A[:,4]
# xf=np.concatenate((x,x2))
# yf=np.concatenate((y,y2))
# img, extent = myplot(xf, yf,sigma)
# ax[1,2].imshow(img, extent=extent, origin='lower', cmap=cm.jet)
# ax[1,2].set_title("r=2,3")

fig,ax=plt.subplots(1,2)
sigma=1
A=np.loadtxt("enlazante.txt")
plt.suptitle("Orbitales para las funciones de onda enlazante y antienlazante")
x = A[:,0]
y = A[:,2]
x2= A[:,3]
y2= A[:,5]
xf=np.concatenate((x,x2))
yf=np.concatenate((y,y2))
img, extent = myplot(xf, yf,sigma)
ax[0].imshow(img, extent=extent, origin='lower', cmap=cm.jet)
ax[0].set_title(r"$\phi_{1}$ r=0,74")
ax[0].set_xlim(-2,2)
ax[0].set_ylim(-2,2)
ax[0].axes.get_xaxis().set_visible(False)
ax[0].axes.get_yaxis().set_visible(False)


A=np.loadtxt("antienlazante.txt")
x = A[:,0]
y = A[:,2]
x2= A[:,3]
y2= A[:,5]
xf=np.concatenate((x,x2))
yf=np.concatenate((y,y2))
img, extent = myplot(xf, yf,sigma)
ax[1].imshow(img, extent=extent, origin='lower', cmap=cm.jet)
ax[1].set_title(r"$\phi_{2}$ r=0,74")
ax[1].set_xlim(-2,2)
ax[1].set_ylim(-2,2)
ax[1].axes.get_xaxis().set_visible(False)
ax[1].axes.get_yaxis().set_visible(False)



plt.show()

