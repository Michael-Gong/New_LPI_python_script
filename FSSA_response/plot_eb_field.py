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
font_size = 25
######### Parameter you should set ###########

upper = matplotlib.cm.jet(np.arange(256))
lower = np.ones((int(256/4),4))
for i in range(3):
    lower[:,i] = np.linspace(1, upper[0,i], lower.shape[0])
cmap = np.vstack(( lower, upper ))
mycolor_jet = matplotlib.colors.ListedColormap(cmap, name='myColorMap', N=cmap.shape[0])


if __name__ == "__main__":
  to_path = './two/'
#  nsteps      = 79580 #sum(1 for line in open(from_path+'x_0000.txt'))/part_number

  px  =  np.load(to_path+'px2d.npy')
  py  =  np.load(to_path+'py2d.npy')
#  pz  =  np.load(to_path+'pz3d.npy')
  xx  =  np.load(to_path+'xx2d.npy')
  yy  =  np.load(to_path+'yy2d.npy')
#  zz  =  np.load(to_path+'zz3d.npy')
  work_x = np.load(to_path+'workx2d.npy')
  work_y = np.load(to_path+'worky2d.npy')
#  work_z = np.load(to_path+'workz3d.npy')
  #gg  =  (px**2+py**2+pz**2+1.0)**0.5
  gg  =  (px**2+py**2+1.0)**0.5

  part_number = len(px[:,0])
  nsteps      = len(px[0,:])

  ey  =  np.load(to_path+'ey2d.npy') 
  ex  =  np.load(to_path+'ex2d.npy') 
  bz  =  np.load(to_path+'bz2d.npy') 

  for n in range(22):
      plt.subplot(3,3,1)
      plt.plot(px[n,:], py[n,:], color='k', linewidth=3, linestyle='-')
      plt.scatter(px[n,:], py[n,:], c=gg[n,:], norm=colors.Normalize(vmin=0,vmax=1500), s=30, cmap='rainbow', edgecolors='None', alpha=1,zorder=3)
#  cbar=plt.colorbar( ticks=np.linspace(0, 500, 3) ,pad=0.005)
#  cbar.ax.set_yticklabels(cbar.ax.get_yticklabels(), fontsize=font_size)
#  cbar.set_label('$\gamma$',fontdict=font)
      plt.xlabel('$p_x$ [$m_ec$]',fontdict=font)
      plt.ylabel('$p_y$ [$m_ec$]',fontdict=font)
      plt.xticks(fontsize=font_size); plt.yticks(fontsize=font_size);
      #plt.xscale('log')
      #plt.xlim(-200,14000)
      #plt.ylim(-800,800)

      plt.subplot(3,3,2)
      plt.plot(xx[n,:], yy[n,:], color='k', linewidth=3, linestyle='-')
      plt.scatter(xx[n,:], yy[n,:], c=gg[n,:], norm=colors.Normalize(vmin=0,vmax=1500), s=30, cmap='rainbow', edgecolors='None', alpha=1,zorder=3)
#  cbar=plt.colorbar( ticks=np.linspace(0, 500, 3) ,pad=0.005)
#  cbar.ax.set_yticklabels(cbar.ax.get_yticklabels(), fontsize=font_size)
#  cbar.set_label('$\gamma$',fontdict=font)
      plt.xlabel('$x$ [$\mu m$]',fontdict=font)
      plt.ylabel('$y$ [$\mu m$]',fontdict=font)
      plt.xticks(fontsize=font_size); plt.yticks(fontsize=font_size);
      #plt.xscale('log')
      #plt.xlim(-200,14000)
      #plt.ylim(-800,800)

      plt.subplot(3,3,3)
      plt.plot(xx[n,:], work_x[n,:], color='salmon', linewidth=3, linestyle='--',label='$W_x$')
      plt.plot(xx[n,:], work_y[n,:], color='limegreen', linewidth=3, linestyle='--',label='$W_y$')
#  cbar=plt.colorbar( ticks=np.linspace(0, 500, 3) ,pad=0.005)
#  cbar.ax.set_yticklabels(cbar.ax.get_yticklabels(), fontsize=font_size)
#  cbar.set_label('$\gamma$',fontdict=font)
      plt.legend(loc='upper left',fontsize=18)
      plt.xlabel('$x$ [$\mu m$]',fontdict=font)
      plt.ylabel('$W_{x,y}$ [$m_ec^2$]',fontdict=font)
      plt.xticks(fontsize=font_size); plt.yticks(fontsize=font_size);
      #plt.xscale('log')
      #plt.xlim(-200,14000)
      #plt.ylim(-800,800)

      plt.subplot(3,3,4)
      plt.plot(xx[n,:], py[n,:], color='k', linewidth=3, linestyle='-')
      plt.scatter(xx[n,:], py[n,:], c=gg[n,:], norm=colors.Normalize(vmin=0,vmax=1500), s=10, cmap='rainbow', edgecolors='None', alpha=1,zorder=3)
#  cbar=plt.colorbar( ticks=np.linspace(0, 500, 3) ,pad=0.005)
#  cbar.ax.set_yticklabels(cbar.ax.get_yticklabels(), fontsize=font_size)
#  cbar.set_label('$\gamma$',fontdict=font)
      plt.xlabel('$x$ [$\mu m$]',fontdict=font)
      plt.ylabel('$p_y$ [$m_ec$]',fontdict=font)
      plt.xticks(fontsize=font_size); plt.yticks(fontsize=font_size);
      #plt.xscale('log')
      #plt.xlim(-200,14000)
      #plt.ylim(-800,800)

      plt.subplot(3,3,5)
      plt.plot(xx[n,:], (gg-px)[n,:], color='k', linewidth=3, linestyle='-')
      plt.scatter(xx[n,:], (gg-px)[n,:], c=gg[n,:], norm=colors.Normalize(vmin=0,vmax=1500), s=10, cmap='rainbow', edgecolors='None', alpha=1,zorder=3)
#  cbar=plt.colorbar( ticks=np.linspace(0, 500, 3) ,pad=0.005)
#  cbar.ax.set_yticklabels(cbar.ax.get_yticklabels(), fontsize=font_size)
#  cbar.set_label('$\gamma$',fontdict=font)
      plt.xlabel('$x$ [$\mu m$]',fontdict=font)
      plt.ylabel('$R$',fontdict=font)
      plt.xticks(fontsize=font_size); plt.yticks(fontsize=font_size);
      #plt.xscale('log')
      #plt.xlim(-200,14000)
      #plt.ylim(-800,800)

      plt.subplot(3,3,6)
      plt.plot(xx[n,:], gg[n,:], color='k', linewidth=3, linestyle='-')
      plt.scatter(xx[n,:], gg[n,:], c=gg[n,:], norm=colors.Normalize(vmin=0,vmax=1500), s=10, cmap='rainbow', edgecolors='None', alpha=1,zorder=3)
#  cbar=plt.colorbar( ticks=np.linspace(0, 500, 3) ,pad=0.005)
#  cbar.ax.set_yticklabels(cbar.ax.get_yticklabels(), fontsize=font_size)
#  cbar.set_label('$\gamma$',fontdict=font)
      plt.xlabel('$x$ [$\mu m$]',fontdict=font)
      plt.ylabel('$\gamma$',fontdict=font)
      plt.xticks(fontsize=font_size); plt.yticks(fontsize=font_size);
      #plt.xscale('log')
      #plt.xlim(-200,14000)
      #plt.ylim(-800,800)

      plt.subplot(3,3,7)
      
      print(np.size(ey[ey>1.0]),np.shape(ey))
      plt.plot(xx[n,:], ey[n,:], color='k', linewidth=3, linestyle='-')
      plt.scatter(xx[n,:], ey[n,:], c=gg[n,:], norm=colors.Normalize(vmin=0,vmax=1500), s=10, cmap='rainbow', edgecolors='None', alpha=1,zorder=3)
#  cbar=plt.colorbar( ticks=np.linspace(0, 500, 3) ,pad=0.005)
#  cbar.ax.set_yticklabels(cbar.ax.get_yticklabels(), fontsize=font_size)
#  cbar.set_label('$\gamma$',fontdict=font)
      plt.xlabel('$x$ [$\mu m$]',fontdict=font)
      plt.ylabel('$E_y$ [$m_ec\omega/|e|$]',fontdict=font)
      plt.xticks(fontsize=font_size); plt.yticks(fontsize=font_size);
      #plt.xscale('log')
      #plt.xlim(-200,14000)
      #plt.ylim(-800,800)

      plt.subplot(3,3,9)
      plt.plot(xx[n,:], (bz-ey)[n,:], color='k', linewidth=3, linestyle='-')
      plt.scatter(xx[n,:], (bz-ey)[n,:], c=gg[n,:], norm=colors.Normalize(vmin=0,vmax=1500), s=10, cmap='rainbow', edgecolors='None', alpha=1,zorder=3)
#  cbar=plt.colorbar( ticks=np.linspace(0, 500, 3) ,pad=0.005)
#  cbar.ax.set_yticklabels(cbar.ax.get_yticklabels(), fontsize=font_size)
#  cbar.set_label('$\gamma$',fontdict=font)
      plt.xlabel('$x$ [$\mu m$]',fontdict=font)
      plt.ylabel('$B_z-E_y$ [$m_ec\omega/|e|$]',fontdict=font)
      plt.xticks(fontsize=font_size); plt.yticks(fontsize=font_size);
      #plt.xscale('log')
      #plt.xlim(-200,14000)
      #plt.ylim(-800,800)

      plt.subplot(3,3,8)
      plt.plot(xx[n,:], bz[n,:], color='k', linewidth=3, linestyle='-')
      plt.scatter(xx[n,:], bz[n,:], c=gg[n,:], norm=colors.Normalize(vmin=0,vmax=1500), s=10, cmap='rainbow', edgecolors='None', alpha=1,zorder=3)
#  cbar=plt.colorbar( ticks=np.linspace(0, 500, 3) ,pad=0.005)
#  cbar.ax.set_yticklabels(cbar.ax.get_yticklabels(), fontsize=font_size)
#  cbar.set_label('$\gamma$',fontdict=font)
      plt.xlabel('$x$ [$\mu m$]',fontdict=font)
      plt.ylabel('$B_z$ [$m_e\omega/|e|$]',fontdict=font)
      plt.xticks(fontsize=font_size); plt.yticks(fontsize=font_size);
      #plt.xscale('log')
      #plt.xlim(-200,14000)
      #plt.ylim(-800,800)


      plt.subplots_adjust(left=0.10, bottom=0.08, right=0.95, top=0.92,
                wspace=0.25, hspace=0.35)

      #plt.show()
      #lt.figure(figsize=(100,100))
      fig = plt.gcf()
      fig.set_size_inches(30, 30)
      fig.savefig(to_path+'field_particle'+str(n).zfill(2)+'.png',format='png',dpi=160)
      plt.close("all")
      print('finished!'+str(n))
      print(np.max(gg[n,:]))
