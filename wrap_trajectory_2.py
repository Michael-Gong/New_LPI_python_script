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
        'size'   : 23,
        }

font_2 = {'family' : 'monospace',
        'color'  : 'black',
        'weight' : 'normal',
        'size'   : 23,
        }

font_size = 23
font_size2 = 18

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

#plt.subplot2grid((6,2),(0,0),colspan=2,rowspan=2)

plt.subplot(2,1,1)
plt.scatter(xi_0[xi_0<30], gg0[xi_0<30]/1e3,  s=7, color='blue', zorder=3,label='No RR')
#plt.scatter(xi_1[xi_1<30], y1[xi_1<30], c=R_1[xi_1<30], s=10, cmap='magma', edgecolors='None',zorder=4)
plt.scatter(xi_1[xi_1<30], gg1[xi_1<30]/1e3,  s=7, color='red', zorder=4,label='Classic RR')

plt.scatter(xi_2[0,:], gg2[0,:]/1e3, s=7, color='lime',zorder=4,label='QED')
#plt.xlabel(r'$\xi\ [2\pi]$',fontdict=font)
plt.ylabel(r'$\gamma\ [10^3]$',fontdict=font_2)
plt.xticks([0,4,8,12],fontsize=0.001); 
plt.yticks([0,5,10],fontsize=23);
#plt.yticks([-8,-4,0,4,8],fontsize=font_size);
plt.xlim(.0,13.25)
plt.grid(color='k', linestyle='--', linewidth=0.5)


plt.subplot(2,1,2)
#plt.scatter(xi_0[xi_0<30], y0[xi_0<30], c=R_0[xi_0<30], s=10, cmap='rainbow', edgecolors='None',zorder=3)
plt.plot(xi_0[0,:], (1983.*chi_0**2)[0,:],  linewidth=3.3, color='blue', zorder=1,label='No RR',alpha=1)

plt.plot(xi_2[0,:], (1983.*chi_2**2)[0,:], linewidth=3.3, color='lime',zorder=2,label='QED',alpha=1)
#plt.scatter(xi_1[xi_1<30], y1[xi_1<30], c=R_1[xi_1<30], s=10, cmap='magma', edgecolors='None',zorder=4)
plt.plot(xi_1[0,:], (1983.*chi_1**2)[0,:],  linewidth=3.3, color='red', zorder=2,label='Classic RR',alpha=1)

plt.xlabel(r'$\xi\ [2\pi]$',fontdict=font)
plt.ylabel(r'$\mathcal{P}_{rad}\ [m_ec^2\omega]$',fontdict=font_2)
plt.xticks([0,4,8,12],fontsize=font_size); 
plt.yticks([0,5,10,15],fontsize=23);
#plt.yticks([-8,-4,0,4,8],fontsize=font_size);
plt.xlim(.0,13.25)
plt.grid(color='k', linestyle='--', linewidth=0.5)


plt.subplots_adjust(top=0.99, bottom=0.15, left=0.12, right=0.98, hspace=0.03, wspace=0.03)
fig = plt.gcf()
fig.set_size_inches(9, 6)
#fig.set_size_inches(5, 4.5)
fig.savefig('./wrap_trajectory_2.png',format='png',dpi=160)
plt.close("all")
