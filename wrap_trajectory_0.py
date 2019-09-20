from scipy.integrate import odeint
#import sdf
import matplotlib
import matplotlib as mpl
#mpl.style.use('https://gist.github.com/danielballan/be066529de85e87a5fe7/raw')
#matplotlib.use('agg')
import matplotlib.pyplot as plt
import numpy as np
from numpy import ma
from matplotlib import colors, ticker, cm
from matplotlib.mlab import bivariate_normal
from optparse import OptionParser
import os
from mpl_toolkits.mplot3d import Axes3D
import random
from mpl_toolkits import mplot3d
from matplotlib import rc
import matplotlib.transforms as mtransforms
from matplotlib.colors import LinearSegmentedColormap
#rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
## for Palatino and other serif fonts use:
#rc('font',**{'family':'serif','serif':['Palatino']})
#rc('text', usetex=True)

font = {'family' : 'monospace',
        'color'  : 'black',
        'weight' : 'normal',
        'size'   : 25,
        }

font_size = 25
font_size2 = 20

# get colormap
ncolors = 256
color_array = plt.get_cmap('Spectral_r')(range(ncolors))
# change alpha values
color_array[:,-1] = abs(np.linspace(-1.0,1.0,ncolors))**1.5
# create a colormap object
map_object = LinearSegmentedColormap.from_list(name='bwr_alpha',colors=color_array)

def cal_chi(Ey,Bz,px,py):
    temp_g = (px**2+py**2+1.)**0.5
    a_s = 4.1e5
    return ((temp_g*Ey-px*Bz)**2+(py*Bz)**2-(Ey*py)**2)**0.5/a_s

#plt.rc('text', usetex=True)
#plt.rc('font', family='serif')

# function that returns dz/dt

# solve ODE

part_number = 1
from_path = './part_01_no/'
nsteps      = int(sum(1 for line in open(from_path+'t_0000.txt'))/part_number)

from_path = './part_01_no/'
t0  = np.loadtxt(from_path+'t_0000.txt')/2/np.pi
x0  = np.loadtxt(from_path+'x_0000.txt')/2/np.pi
y0  = np.loadtxt(from_path+'y_0000.txt')/2/np.pi
px0 = np.loadtxt(from_path+'px_0000.txt')
py0 = np.loadtxt(from_path+'py_0000.txt')
Ey0 = np.loadtxt(from_path+'e_part_0000.txt')
Bz0 = np.loadtxt(from_path+'b_part_0000.txt')
t0  = np.reshape(t0,(part_number,nsteps))
x0  = np.reshape(x0,(part_number,nsteps))
y0  = np.reshape(y0,(part_number,nsteps))
px0 = np.reshape(px0,(part_number,nsteps))
py0 = np.reshape(py0,(part_number,nsteps))
Ey0 = np.reshape(Ey0,(part_number,nsteps))
Bz0 = np.reshape(Bz0,(part_number,nsteps))
gg0 = (px0**2+py0**2+1)**0.5
xi_0 = t0-x0
R_0  = gg0-px0
C0 = (gg0-px0)[:,0]
chi_0 = cal_chi(Ey0,Bz0,px0,py0)

from_path = './part_01_rr/'
t1  = np.loadtxt(from_path+'t_0000.txt')/2/np.pi
x1  = np.loadtxt(from_path+'x_0000.txt')/2/np.pi
y1  = np.loadtxt(from_path+'y_0000.txt')/2/np.pi
px1 = np.loadtxt(from_path+'px_0000.txt')
py1 = np.loadtxt(from_path+'py_0000.txt')
Ey1 = np.loadtxt(from_path+'e_part_0000.txt')
Bz1 = np.loadtxt(from_path+'b_part_0000.txt')
radt1 = np.loadtxt(from_path+'radt_0000.txt')
rad_px1 = np.loadtxt(from_path+'rad_px_0000.txt')
t1  = np.reshape(t1,(part_number,nsteps))
x1  = np.reshape(x1,(part_number,nsteps))
y1  = np.reshape(y1,(part_number,nsteps))
px1 = np.reshape(px1,(part_number,nsteps))
py1 = np.reshape(py1,(part_number,nsteps))
Ey1 = np.reshape(Ey1,(part_number,nsteps))
Bz1 = np.reshape(Bz1,(part_number,nsteps))
radt1 = np.reshape(radt1,(part_number,nsteps))
rad_px1 = np.reshape(rad_px1,(part_number,nsteps))

gg1 = (px1**2+py1**2+1)**0.5
xi_1 = t1-x1
R_1  = gg1-px1
rad_R = radt1-rad_px1
print(rad_R)
C1_R = C0 - rad_R
print(C1_R)
chi_1 = cal_chi(Ey1,Bz1,px1,py1)


from_path = './part_01_qe/'
t2  = np.loadtxt(from_path+'t_0000.txt')/2/np.pi
x2  = np.loadtxt(from_path+'x_0000.txt')/2/np.pi
y2  = np.loadtxt(from_path+'y_0000.txt')/2/np.pi
px2 = np.loadtxt(from_path+'px_0000.txt')
py2 = np.loadtxt(from_path+'py_0000.txt')
Ey2 = np.loadtxt(from_path+'e_part_0.txt')
Bz2 = np.loadtxt(from_path+'b_part_0.txt')
radt2 = np.loadtxt(from_path+'radt_0000.txt')
radn2 = np.loadtxt(from_path+'radn_0000.txt')
rad_px2 = np.loadtxt(from_path+'rad_px_0000.txt')
part_number=10
t2  = np.reshape(t2,(part_number,nsteps))
x2  = np.reshape(x2,(part_number,nsteps))
y2  = np.reshape(y2,(part_number,nsteps))
px2 = np.reshape(px2,(part_number,nsteps))
py2 = np.reshape(py2,(part_number,nsteps))
Ey2 = np.reshape(Ey2,(part_number,nsteps))
Bz2 = np.reshape(Bz2,(part_number,nsteps))
radt2 = np.reshape(radt2,(part_number,nsteps))
radn2 = np.reshape(radn2,(part_number,nsteps))
rad_px2 = np.reshape(rad_px2,(part_number,nsteps))
gg2 = (px2**2+py2**2+1)**0.5
xi_2 = t2-x2
R_2  = gg2-px2
rad_R = radt2-rad_px2
C2_R = C0 - rad_R
chi_2 = cal_chi(Ey2,Bz2,px2,py2)


alpha =10**(2.2-3.5)
y_max = (C0/alpha)**0.5/(2*np.pi)
x_grid = np.linspace(0,30.0,601)
y_grid = np.linspace(-10,10,201)

[x_grid,y_grid] = np.meshgrid(x_grid,y_grid)
Ey = np.cos(2*np.pi*x_grid)

eee=np.max([-np.min(Ey),np.max(Ey)])
levels = np.linspace(-eee, eee, 64)

plt.subplot(1,1,1)
#### manifesting colorbar, changing label and axis properties ####
#cbar=plt.colorbar(ticks=[-eee, -eee/2, 0, eee/2, eee])
#cbar.set_label('Normalized electric field',fontdict=font)        
Bz = -2*alpha*y_grid*(alpha*C0)**0.5
eee=np.max([-np.min(Bz),np.max(Bz)])
levels = np.linspace(-eee, eee, 512)
plt.contourf(x_grid, y_grid, Ey, levels=np.linspace(-1., 1., 64), cmap='bwr',alpha=0.6,zorder=0)
plt.contourf(x_grid, y_grid, Bz, levels=levels, cmap=map_object,zorder=1)


x_grid_t = np.linspace(0,17.0,12)
y_grid_t = np.linspace(-8.5,8.5,6)
[x_grid_t,y_grid_t] = np.meshgrid(x_grid_t,y_grid_t)
b_color='k'
plt.plot(x_grid_t[y_grid_t>7], y_grid_t[y_grid_t>7], 'x', color=b_color, markersize=15, markeredgewidth=2,zorder=2)
plt.plot(x_grid_t[(y_grid_t>4)&(y_grid_t<7)], y_grid_t[(y_grid_t>4)&(y_grid_t<7)], 'x', color=b_color, markersize=11.5, markeredgewidth=1.8,zorder=2)
plt.plot(x_grid_t[(y_grid_t>0)&(y_grid_t<4)], y_grid_t[(y_grid_t>0)&(y_grid_t<4)], 'x', color=b_color, markersize=6.5, markeredgewidth=1.5,zorder=2)
plt.plot(x_grid_t[y_grid_t<-7], y_grid_t[y_grid_t<-7], '.', color=b_color, markersize=15, markeredgewidth=2,zorder=2)
plt.plot(x_grid_t[(y_grid_t<-4)&(y_grid_t>-7)], y_grid_t[(y_grid_t<-4)&(y_grid_t>-7)], '.', color=b_color, markersize=11.5, markeredgewidth=1.8,zorder=2)
plt.plot(x_grid_t[(y_grid_t<0)&(y_grid_t>-4)], y_grid_t[(y_grid_t<0)&(y_grid_t>-4)], '.', color=b_color, markersize=6.5, markeredgewidth=1.5,zorder=2)



#plt.scatter(xi_0[xi_0<30], y0[xi_0<30], c=R_0[xi_0<30], s=10, cmap='rainbow', edgecolors='None',zorder=3)
plt.plot(xi_0[xi_0<30], y0[xi_0<30],  linewidth=4.3, color='blue', zorder=3,label='No RR')
plt.plot(np.linspace(0,25,500), np.zeros_like(np.linspace(0,25,500))+y_max,':',color='blue',linewidth=3.5,zorder=3)
plt.plot(np.linspace(0,25,500), np.zeros_like(np.linspace(0,25,500))-y_max,':',color='blue',linewidth=3.5,zorder=3)

y_max = (C1_R/alpha)**0.5/(2*np.pi)
#plt.scatter(xi_1[xi_1<30], y1[xi_1<30], c=R_1[xi_1<30], s=10, cmap='magma', edgecolors='None',zorder=4)
plt.plot(xi_1[xi_1<30], y1[xi_1<30], linewidth=4.3, color='red', zorder=4,label='Classic RR')
plt.plot(xi_1[0,:], y_max[0,:],':',color='r',linewidth=3.5,zorder=4)
plt.plot(xi_1[0,:],-y_max[0,:],':',color='r',linewidth=3.5,zorder=4)
print(xi_1,y_max)

condition = ((radn2[0,1:]-radn2[0,:-1]) > 0) & ((radt2[0,1:]-radt2[0,:-1]) > 10.0 ) 
y_max = (C2_R/alpha)**0.5/(2*np.pi)
plt.plot(xi_2[0,:], y2[0,:], linewidth=4.3, color='lime',zorder=4,label='QED')
plt.scatter(xi_2[0,:-1][condition], y2[0,:-1][condition], color='lime', s=150, marker='*', edgecolors='k', alpha=1,zorder=5)
plt.plot(xi_2[0,:], y_max[0,:],':',color='lime',linewidth=3.5,zorder=4)
plt.plot(xi_2[0,:], -y_max[0,:],':',color='lime',linewidth=3.5,zorder=4)

#plt.plot((t[index,:])/2/np.pi,np.sqrt(px[index,:]**2+py[index,:]**2+1),'--k',linewidth=2.5,label='No RR')
#plt.legend(loc='upper right')
#cbar=plt.colorbar(ticks=np.linspace(np.min(gamma), np.max(gamma), 5))
#cbar=plt.colorbar(ticks=np.linspace(np.min(lgR), np.max(lgR), 5))
#cbar.set_label(r'$R$', fontdict=font)#plt.xlim(47,53)
plt.xlabel(r'$\xi\ [2\pi]$',fontdict=font)
plt.ylabel(r'$y\ [\mu m]$',fontdict=font)
plt.xticks([0,4,8,12],fontsize=font_size); 
plt.yticks([-8,-4,0,4,8],fontsize=font_size);
plt.ylim(-9.8,9.8)
plt.xlim(.0,13.25)
plt.legend(loc='lower right',fontsize=(font_size2-3),framealpha=0.5)
#plt.grid(color='k', linestyle='--', linewidth=1)


#plt.show()
plt.subplots_adjust(top=0.99, bottom=0.15, left=0.12, right=0.98, hspace=0.03, wspace=0.03)
fig = plt.gcf()
fig.set_size_inches(9, 6)
#fig.set_size_inches(5, 4.5)
fig.savefig('./wrap_trajectory_0.png',format='png',dpi=160)
plt.close("all")
