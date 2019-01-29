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

for n in range(1):
    n = 7
    dw_y = -py[n,:]/gg[n,:]*ey[n,:]*2*np.pi/3.333 
    ax = plt.subplot(4,1,1)
    plt.scatter(grid_time[dw_y>0],work_y[n,dw_y>0], color='lime',alpha=1,s=150)
    plt.plot(grid_time, work_y[n,:], color='k', linewidth=5, linestyle='-',label='W$_y$')
    plt.plot(grid_time, work_x[n,:], color='royalblue', linewidth=5, linestyle='-',label='W$_x$')
    plt.plot(grid_time, gg[n,:], color='r', linewidth=5, linestyle='--',label='$\gamma$')
#    plt.scatter(grid_time, work_y[n,:], c=gg[n,:], norm=colors.Normalize(vmin=0,vmax=1500), s=30, cmap='magma', edgecolors='None', alpha=1,zorder=3)
    #  cbar=plt.colorbar( ticks=np.linspace(0, 500, 3) ,pad=0.005)
    #  cbar.ax.set_yticklabels(cbar.ax.get_yticklabels(), fontsize=font_size)
    #  cbar.set_label('$\gamma$',fontdict=font)
    #plt.ylabel('time [fs]',fontdict=font)
    plt.ylabel('W$_{x,y}$ [m$_e$c$^2$]',fontdict=font)
    plt.xlim(125,290)
    plt.xticks([150,200,250],fontsize=0.01); plt.yticks([0,500,1000,1500],fontsize=font_size);
    plt.grid(linestyle='--')
    plt.legend(loc='upper left',fontsize=18)
    #plt.xscale('log')
    #plt.xlim(-200,14000)
    #plt.ylim(-800,800)
    
    ax = plt.subplot(4,1,2)
    y2 = np.zeros_like(grid_time)
    ax.fill_between(grid_time,dw_y,y2, where=dw_y>=y2, facecolor='lime',alpha=0.5,interpolate=True)
    ax.fill_between(grid_time,dw_y,y2, where=dw_y<=y2, facecolor='black',alpha=0.5,interpolate=True)
    plt.plot(grid_time,dw_y, color='k', linewidth=5, linestyle='-')
#    plt.scatter(grid_time[dw_y>0],dw_y[dw_y>0], color='limegreen',alpha=1,s=100)
#    plt.scatter(grid_time,-py[n,:]/gg[n,:]*ey[n,:]*2*np.pi/3.333, c=gg[n,:], norm=colors.Normalize(vmin=0,vmax=1200), s=30, cmap='jet', edgecolors='None', alpha=1,zorder=3)
    #  cbar=plt.colorbar( ticks=np.linspace(0, 500, 3) ,pad=0.005)
    #  cbar.ax.set_yticklabels(cbar.ax.get_yticklabels(), fontsize=font_size)
    #  cbar.set_label('$\gamma$',fontdict=font)
    #plt.ylabel('time [fs]',fontdict=font)
    plt.ylabel('dW$_y$/dt [m$_e$c$^2$/fs]',fontdict=font)
    plt.xlim(125,290)
    plt.ylim(-8,38)
    plt.xticks([150,200,250],fontsize=0.01); plt.yticks([0,10,20,30],fontsize=font_size);
    plt.grid(linestyle='--')
    #plt.xscale('log')
    #plt.xlim(-200,14000)
    #plt.ylim(-800,800)
    
    
    plt.subplot(4,1,3)
    grid_y = np.linspace(-2,2,401)
    X, Y = np.meshgrid(grid_time, grid_y)
    Z = np.zeros_like(X)
    Z = Z+ey[n,:] 
    eee  = 50. 
    levels = np.linspace(-eee, eee, 101)
    Z[Z < -49]=-49
    Z[Z >  49]= 49
    plt.contourf(X, Y, Z, levels=levels, norm=mcolors.Normalize(vmin=levels.min(), vmax=levels.max()), cmap=cm.bwr) 
    plt.plot(grid_time,yy[n,:], color='k', linewidth=5, linestyle='-')
    plt.scatter(grid_time[dw_y>0],yy[n,dw_y>0], color='lime',alpha=1,s=150)
#    plt.scatter(grid_time,yy[n,:], c=gg[n,:], norm=colors.Normalize(vmin=0,vmax=1200), s=30, cmap='jet', edgecolors='None', alpha=1,zorder=3)
    #  cbar=plt.colorbar( ticks=np.linspace(0, 500, 3) ,pad=0.005)
    #  cbar.ax.set_yticklabels(cbar.ax.get_yticklabels(), fontsize=font_size)
    #  cbar.set_label('$\gamma$',fontdict=font)
    #plt.ylabel('time [fs]',fontdict=font)
    plt.ylabel('$y$ [$\mu m$]',fontdict=font)
    plt.xlim(125,290)
    plt.ylim(-1.5,1.5)
    plt.xticks([150,200,250],fontsize=0.01); plt.yticks([-1,0,1],fontsize=font_size);
    plt.grid(linestyle='--')
    #plt.xscale('log')
    #plt.xlim(-200,14000)
    #plt.ylim(-800,800)
    
    plt.subplot(4,1,4) 
    #axin1 = inset_axes(ax, width='15%', height='5%', loc='upper left')
    #axin2 = inset_axes(ax, width='15%', height='5%', loc='upper center')
    X, Y = np.meshgrid(grid_x, grid_time)
    eee  = 50. 
    levels = np.linspace(-eee, eee, 101)
    grid_data[grid_data < -49]=-49
    grid_data[grid_data >  49]= 49
    plt.contourf(Y.T, (30+X-v0*Y*1e-15/1e-6).T, grid_data, levels=levels, norm=mcolors.Normalize(vmin=levels.min(), vmax=levels.max()), cmap=cm.bwr)
    #cbar=plt.colorbar(ticks=np.linspace(-60,60,3))#,orientation="horizontal")
    #cbar.set_label('$E_y$ [$m_ec\omega/e$]', fontdict=font2)
    #for n in range(7):
    #    plt.plot(xx[n,:]+30-v0*grid_time*1e-15/1e-6, grid_time,'-k')
    #plt.contourf(X, Y, grid_data.T, levels=levels, norm=mcolors.Normalize(vmin=levels.min(), vmax=levels.max()), cmap=cm.bwr)
    #n = 3
#    plt.scatter(grid_time,xx[n,:]+30-v0*grid_time*1e-15/1e-6, c=gg[n,:], norm=colors.Normalize(vmin=0,vmax=1200), s=30, cmap='jet', edgecolors='None', alpha=1,zorder=3)
    plt.plot(grid_time,xx[n,:]+30-v0*grid_time*1e-15/1e-6, color='k', linewidth=5, linestyle='-')
    plt.scatter(grid_time[dw_y>0],(xx[n,:]+30-v0*grid_time*1e-15/1e-6)[dw_y>0], color='lime',alpha=1,s=150)
    #plt.xlabel('Energy [MeV]',fontdict=font)
    plt.xlabel('t [fs]',fontdict=font)
#    plt.ylabel('$\Delta=(x-ct)/\lambda+30$ [$\mu$m]',fontdict=font)
    plt.ylabel('$\Delta$ [$\mu m$]',fontdict=font)
    plt.xticks([150,200,250],fontsize=font_size); plt.yticks([-4,-3,-2],fontsize=font_size);
    #plt.yscale('log')
    plt.xlim(125,290)
    plt.ylim(-4,-2)
    #plt.legend(loc='best',fontsize=20,framealpha=0.5)
    plt.grid(linestyle='--')
    
    
    plt.subplots_adjust(left=0.10, bottom=0.08, right=0.95, top=0.92,
                    wspace=0.05, hspace=0.08)
    
    fig = plt.gcf()
    fig.set_size_inches(16, 15)
    fig.savefig(to_path+'wrap_trace_'+str(n).zfill(4)+'.png',format='png',dpi=160)
    plt.close("all")
    print('ok')

