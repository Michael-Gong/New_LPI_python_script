import sdf
import matplotlib
matplotlib.use('agg')
#%matplotlib inline
import matplotlib.pyplot as plt
import numpy as np
#from numpy import ma
from matplotlib import colors, ticker, cm
from matplotlib.mlab import bivariate_normal
from optparse import OptionParser
import os
import matplotlib.colors as mcolors
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import matplotlib.gridspec as gridspec
from numpy import ma
from matplotlib import colors, ticker, cm
from matplotlib.mlab import bivariate_normal
import matplotlib.colors as mcolors
import scipy.ndimage as ndimage



######## Constant defined here ########
pi        =     3.1415926535897932384626
q0        =     1.602176565e-19 # C
m0        =     9.10938291e-31  # kg
v0        =     2.99792458e8    # m/s^2
kb        =     1.3806488e-23   # J/K
mu0       =     4.0e-7*pi       # N/A^2
epsilon0  =     8.8541878176203899e-12 # F/m
h_planck  =     6.62606957e-34  # J s
wavelength=     1.0e-6
frequency =     v0*2*pi/wavelength

exunit    =     m0*v0*frequency/q0
bxunit    =     m0*frequency/q0
denunit    =     frequency**2*epsilon0*m0/q0**2
print('electric field unit: '+str(exunit))
print('magnetic field unit: '+str(bxunit))
print('density unit nc: '+str(denunit))

font = {'family' : 'monospace',  
        'style'  : 'normal',
        'color'  : 'black',  
	    'weight' : 'normal',  
        'size'   : 25,  
       }  

font2 = {'family' : 'monospace',  
        'style'  : 'normal',
        'color'  : 'black',  
	    'weight' : 'normal',  
        'size'   : 25,  
       }

font_size =25

to_path   = './two/'
grid_data = np.load(to_path+'ey_vph.npy')
grid_time = np.load(to_path+'grid_time.npy')/1e-15
grid_x    = np.load(to_path+'grid_x.npy')

xx  =  np.load(to_path+'xx2d.npy')
yy  =  np.load(to_path+'yy2d.npy')
px  =  np.load(to_path+'px2d.npy')
py  =  np.load(to_path+'py2d.npy')
work_x = np.load(to_path+'workx2d.npy')
work_y = np.load(to_path+'worky2d.npy')
gg  =  (1+px**2+py**2)**0.5
  
ey  =  np.load(to_path+'ey2d.npy') 
ex  =  np.load(to_path+'ex2d.npy') 
bz  =  np.load(to_path+'bz2d.npy') 

for n in range(22):
    plt.subplot(2,3,1) 
    #axin1 = inset_axes(ax, width='15%', height='5%', loc='upper left')
    #axin2 = inset_axes(ax, width='15%', height='5%', loc='upper center')
    X, Y = np.meshgrid(grid_x, grid_time) 
    levels = np.linspace(-60.1, 60.1, 50)
    grid_data[grid_data < -60]=-60
    grid_data[grid_data >  60]= 60
    plt.contourf(30+X-v0*Y*1e-15/1e-6, Y, grid_data.T, levels=levels, norm=mcolors.Normalize(vmin=levels.min(), vmax=levels.max()), cmap=cm.bwr)
    #cbar=plt.colorbar(ticks=np.linspace(-60,60,3))#,orientation="horizontal")
    #cbar.set_label('$E_y$ [$m_ec\omega/e$]', fontdict=font2)
    #for n in range(7):
    #    plt.plot(xx[n,:]+30-v0*grid_time*1e-15/1e-6, grid_time,'-k')
    #plt.contourf(X, Y, grid_data.T, levels=levels, norm=mcolors.Normalize(vmin=levels.min(), vmax=levels.max()), cmap=cm.bwr)
    #n = 3
    plt.scatter(xx[n,:]+30-v0*grid_time*1e-15/1e-6, grid_time, c=gg[n,:], norm=colors.Normalize(vmin=0,vmax=1500), s=30, cmap='rainbow', edgecolors='None', alpha=1,zorder=3)
    #plt.xlabel('Energy [MeV]',fontdict=font)
    plt.ylabel('time [fs]',fontdict=font)
    plt.xlabel('$\Delta=x-ct+30$ [$\mu$m]',fontdict=font)
    plt.xticks(fontsize=20); plt.yticks(fontsize=20);
    #plt.yscale('log')
    plt.ylim(150,290)
    plt.xlim(-5,5)
    #plt.legend(loc='best',fontsize=20,framealpha=0.5)
    plt.grid()
    
    
    plt.subplot(2,3,2)
    plt.plot(py[n,:], grid_time, color='k', linewidth=3, linestyle='-')
    plt.scatter(py[n,:], grid_time, c=gg[n,:], norm=colors.Normalize(vmin=0,vmax=1500), s=30, cmap='rainbow', edgecolors='None', alpha=1,zorder=3)
    #  cbar=plt.colorbar( ticks=np.linspace(0, 500, 3) ,pad=0.005)
    #  cbar.ax.set_yticklabels(cbar.ax.get_yticklabels(), fontsize=font_size)
    #  cbar.set_label('$\gamma$',fontdict=font)
    #plt.ylabel('time [fs]',fontdict=font)
    plt.xlabel('$p_y$ [$m_ec$]',fontdict=font)
    plt.ylim(150,290)
    plt.xticks(fontsize=font_size); plt.yticks(fontsize=0.01);
    plt.grid()
    #plt.xscale('log')
    #plt.xlim(-200,14000)
    #plt.ylim(-800,800)
    
    plt.subplot(2,3,3)
    plt.plot(work_y[n,:], grid_time, color='k', linewidth=3, linestyle='-')
    plt.scatter(work_y[n,:], grid_time, c=gg[n,:], norm=colors.Normalize(vmin=0,vmax=1500), s=30, cmap='rainbow', edgecolors='None', alpha=1,zorder=3)
    #  cbar=plt.colorbar( ticks=np.linspace(0, 500, 3) ,pad=0.005)
    #  cbar.ax.set_yticklabels(cbar.ax.get_yticklabels(), fontsize=font_size)
    #  cbar.set_label('$\gamma$',fontdict=font)
    #plt.ylabel('time [fs]',fontdict=font)
    plt.xlabel('$W_y$ [$m_ec^2$]',fontdict=font)
    plt.ylim(150,290)
    plt.xticks(fontsize=font_size); plt.yticks(fontsize=0.01);
    plt.grid()
    #plt.xscale('log')
    #plt.xlim(-200,14000)
    #plt.ylim(-800,800)
    
    plt.subplot(2,3,4)
    plt.plot(ey[n,:], grid_time, color='k', linewidth=3, linestyle='-')
    plt.scatter(ey[n,:], grid_time, c=gg[n,:], norm=colors.Normalize(vmin=0,vmax=1500), s=30, cmap='rainbow', edgecolors='None', alpha=1,zorder=3)
    #  cbar=plt.colorbar( ticks=np.linspace(0, 500, 3) ,pad=0.005)
    #  cbar.ax.set_yticklabels(cbar.ax.get_yticklabels(), fontsize=font_size)
    #  cbar.set_label('$\gamma$',fontdict=font)
    plt.ylabel('time [fs]',fontdict=font)
    plt.xlabel('$E_y$ [$m_ec\omega/|e|$]',fontdict=font)
    plt.ylim(150,290)
    plt.xticks(fontsize=font_size); plt.yticks(fontsize=font_size);
    plt.grid()
    #plt.xscale('log')
    #plt.xlim(-200,14000)
    #plt.ylim(-800,800)
    
    plt.subplot(2,3,5)
    plt.plot(yy[n,:], grid_time, color='k', linewidth=3, linestyle='-')
    plt.scatter(yy[n,:], grid_time, c=gg[n,:], norm=colors.Normalize(vmin=0,vmax=1500), s=30, cmap='rainbow', edgecolors='None', alpha=1,zorder=3)
    #  cbar=plt.colorbar( ticks=np.linspace(0, 500, 3) ,pad=0.005)
    #  cbar.ax.set_yticklabels(cbar.ax.get_yticklabels(), fontsize=font_size)
    #  cbar.set_label('$\gamma$',fontdict=font)
    #plt.ylabel('time [fs]',fontdict=font)
    plt.xlabel('$y$ [$\mu m$]',fontdict=font)
    plt.ylim(150,290)
    plt.xticks(fontsize=font_size); plt.yticks(fontsize=0.01);
    plt.grid()
    #plt.xscale('log')
    #plt.xlim(-200,14000)
    #plt.ylim(-800,800)
    
    plt.subplot(2,3,6)
    plt.plot(xx[n,:], grid_time, color='k', linewidth=3, linestyle='-')
    plt.scatter(xx[n,:], grid_time, c=gg[n,:], norm=colors.Normalize(vmin=0,vmax=1500), s=30, cmap='rainbow', edgecolors='None', alpha=1,zorder=3)
    #  cbar=plt.colorbar( ticks=np.linspace(0, 500, 3) ,pad=0.005)
    #  cbar.ax.set_yticklabels(cbar.ax.get_yticklabels(), fontsize=font_size)
    #  cbar.set_label('$\gamma$',fontdict=font)
    #plt.ylabel('time [fs]',fontdict=font)
    plt.xlabel('$x$ [$\mu m$]',fontdict=font)
    plt.ylim(150,290)
    plt.xticks(fontsize=font_size); plt.yticks(fontsize=0.01);
    plt.grid()
    #plt.xscale('log')
    #plt.xlim(-200,14000)
    #plt.ylim(-800,800)
    
    
    plt.subplots_adjust(left=0.10, bottom=0.08, right=0.95, top=0.92,
                    wspace=0.05, hspace=0.15)
    
    fig = plt.gcf()
    fig.set_size_inches(30.2, 20.3)
    fig.savefig(to_path+'comb_fig3_zoomin_'+str(n).zfill(4)+'.png',format='png',dpi=160)
    plt.close("all")

