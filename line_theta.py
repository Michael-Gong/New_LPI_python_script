#%matplotlib inline
#import sdf
import matplotlib
import matplotlib as mpl
#mpl.style.use('https://raw.githubusercontent.com/Michael-Gong/DLA_project/master/style')
matplotlib.use('agg')
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
import sys
#rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
## for Palatino and other serif fonts use:
#rc('font',**{'family':'serif','serif':['Palatino']})
#rc('text', usetex=True)
font = {'family' : 'monospace',
        'color'  : 'black',
        'weight' : 'normal',
        'size'   : 25,
       }

font_2 = {'family' : 'monospace',
        'color'  : 'black',
        'weight' : 'normal',
        'size'   : 20,
       }


font_size = 25
font_size2 = 20
#plt.rc('text', usetex=True)
#plt.rc('font', family='serif')

upper = matplotlib.cm.jet(np.arange(256))
lower = np.ones((int(256/4),4))
for i in range(3):
  lower[:,i] = np.linspace(1, upper[0,i], lower.shape[0])
cmap = np.vstack(( lower, upper ))
mycolor_jet = matplotlib.colors.ListedColormap(cmap, name='myColorMap', N=cmap.shape[0])

upper = matplotlib.cm.rainbow(np.arange(256))
lower = np.ones((int(256/4),4))
for i in range(3):
  lower[:,i] = np.linspace(1, upper[0,i], lower.shape[0])
cmap = np.vstack(( lower, upper ))
mycolor_rainbow = matplotlib.colors.ListedColormap(cmap, name='myColorMap', N=cmap.shape[0])

def pxpy_to_energy(gamma, weight):
    binsize = 200
    en_grid = np.linspace(50,19950,200)
    en_bin  = np.linspace(0,20000.0,201)
    en_value = np.zeros_like(en_grid) 
    for i in range(binsize):
#        if i == binsize-1:
#            en_value[i] = sum(weight[en_bin[i]<=gamma])
#        else:
            en_value[i] = sum(weight[ (en_bin[i]<=gamma) & (gamma<en_bin[i+1]) ])
    return (en_grid, en_value)

def theta_to_grid(theta, weight):
    binsize = 90
    theta_grid = np.linspace(0.5,89.5,90)
    theta_bin  = np.linspace(0,90,91)
    theta_value = np.zeros_like(theta_grid) 
    for i in range(binsize):
#        if i == binsize-1:
#            en_value[i] = sum(weight[en_bin[i]<=gamma])
#        else:
            theta_value[i] = sum(weight[ (theta_bin[i]<=theta) & (theta<theta_bin[i+1]) ])
    return (theta_grid, theta_value)


if __name__ == "__main__":
  part_number = 50000
  from_path = './p50000_rr_T150/'
  nsteps      = int(sum(1 for line in open(from_path+'t_tot_s.txt'))/part_number)

  from_path = './p50000_no_T150/'
  t0  = np.loadtxt(from_path+'t_tot_s.txt')/2/np.pi
  #x2  = np.loadtxt(from_path+'x_tot_s.txt')/2/np.pi
  #y2  = np.loadtxt(from_path+'y_tot_s.txt')/2/np.pi
  px0 = np.loadtxt(from_path+'px_tot_s.txt')
  py0 = np.loadtxt(from_path+'py_tot_s.txt')
  t0  = np.reshape(t0,(part_number,nsteps))
  #x2  = np.reshape(x1,(part_number,nsteps))
  #y2  = np.reshape(y1,(part_number,nsteps))
  px0 = np.reshape(px0,(part_number,nsteps))
  py0 = np.reshape(py0,(part_number,nsteps))
  gg0 = (px0**2+py0**2+1)**0.5
  ww0 = np.zeros_like(gg0)+gg0
  theta0 = np.arctan2(py0,px0)/np.pi*180


  from_path = './p50000_qe_T150/'
  t2  = np.loadtxt(from_path+'t_tot_s.txt')/2/np.pi
  #x2  = np.loadtxt(from_path+'x_tot_s.txt')/2/np.pi
  #y2  = np.loadtxt(from_path+'y_tot_s.txt')/2/np.pi
  px2 = np.loadtxt(from_path+'px_tot_s.txt')
  py2 = np.loadtxt(from_path+'py_tot_s.txt')
  t2  = np.reshape(t2,(part_number,nsteps))
  #x2  = np.reshape(x1,(part_number,nsteps))
  #y2  = np.reshape(y1,(part_number,nsteps))
  px2 = np.reshape(px2,(part_number,nsteps))
  py2 = np.reshape(py2,(part_number,nsteps))
  gg2 = (px2**2+py2**2+1)**0.5
  ww2 = np.zeros_like(gg2)+gg2
  theta2 = np.arctan2(py2,px2)/np.pi*180

  from_path = './p50000_rr_T150/'
  t1  = np.loadtxt(from_path+'t_tot_s.txt')/2/np.pi
  #x1  = np.loadtxt(from_path+'x_tot_s.txt')/2/np.pi
  #y1  = np.loadtxt(from_path+'y_tot_s.txt')/2/np.pi
  px1 = np.loadtxt(from_path+'px_tot_s.txt')
  py1 = np.loadtxt(from_path+'py_tot_s.txt')
  t1  = np.reshape(t1,(part_number,nsteps))
  #x1  = np.reshape(x1,(part_number,nsteps))
  #y1  = np.reshape(y1,(part_number,nsteps))
  px1 = np.reshape(px1,(part_number,nsteps))
  py1 = np.reshape(py1,(part_number,nsteps))
  gg1 = (px1**2+py1**2+1)**0.5
  ww1 = np.zeros_like(gg1)+gg1
  theta1 = np.arctan2(py1,px1)/np.pi*180
 
  for i in range(0,nsteps,2):
      print('start ploting: '+str(round(100.*i/nsteps,1))+'%')
      plt.subplot(1,1,1)        
      grid_x, grid_y = theta_to_grid(abs(theta0[:,i]),ww0[:,i])
      plt.plot(grid_x, grid_y*0.51e-3/(np.pi/180.), '-', color='royalblue', linewidth=3,label='No RR')
      grid_x, grid_y = theta_to_grid(abs(theta1[:,i]),ww1[:,i])
      plt.plot(grid_x, grid_y*0.51e-3/(np.pi/180.), '-', color='mediumseagreen', linewidth=3,label='Classic RR')
      grid_x, grid_y = theta_to_grid(abs(theta2[:,i]),ww2[:,i])
      plt.plot(grid_x, grid_y*0.51e-3/(np.pi/180.), '-', color='tomato', linewidth=3,label='QED RR')

#### manifesting colorbar, changing label and axis properties ####
      plt.xlabel(r'$|\theta|$'+' [$^\circ$]',fontdict=font)
      plt.ylabel('d$E_e$/d'+r'$\theta$'+' [GeV$\cdot$rad$^{-1}$]',fontdict=font)
      plt.xticks([10,30,50,70,90],fontsize=font_size); 
#  plt.title('dN/dE for QED RR', fontsize=font_size)
      plt.xlim(1,91)
      plt.yscale('log')
      plt.ylim(1e0,2e6)
      plt.yticks([1e0,1e2,1e4,1e6],fontsize=font_size);
      plt.legend(loc='best',fontsize=font_size2)
#  plt.ylim(0,500)
      plt.grid(b=True, which='major', color='grey', linestyle='--')
      plt.grid(b=True, which='minor', color='grey', linestyle=':')
      plt.title('time at t='+str(round(t1[0,i]*3.333,1)),fontdict=font)
#  plt.text(3.,2,'t = 500 fs',fontdict=font)

      plt.subplots_adjust(top=0.95, bottom=0.12, left=0.15, right=0.98, hspace=0.20, wspace=0.23)

      fig = plt.gcf()
      fig.set_size_inches(8, 7.5)
#fig.set_size_inches(5, 4.5)
      fig.savefig('./wrap_percent_'+str(i).zfill(4)+'.png',format='png',dpi=160)
      plt.close("all")

  
          
          
