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

font_size = 25

#plt.rc('text', usetex=True)
#plt.rc('font', family='serif')

upper = matplotlib.cm.jet(np.arange(256))
lower = np.ones((int(256/4),4))
for i in range(3):
  lower[:,i] = np.linspace(1, upper[0,i], lower.shape[0])
cmap = np.vstack(( lower, upper ))
mycolor_jet = matplotlib.colors.ListedColormap(cmap, name='myColorMap', N=cmap.shape[0])

def pxpy_to_energy(gamma, weight):
    binsize = 15
    en_grid = np.linspace(400,11600,15)
    en_bin  = np.linspace(0,12000.0,16)
    en_value = np.zeros_like(en_grid) 
    for i in range(binsize):
#        if i == binsize-1:
#            en_value[i] = sum(weight[en_bin[i]<=gamma])
#        else:
            en_value[i] = sum(weight[ (en_bin[i]<=gamma) & (gamma<en_bin[i+1]) ])
    return (en_grid, en_value)



if __name__ == "__main__":
  part_number = 50000
  from_path = './p50000_rr_T150/'
  nsteps      = int(sum(1 for line in open(from_path+'t_tot_s.txt'))/part_number)


  from_path = './p50000_rr_T150/'
  #to_path   = './Data_no/'
  t0  = np.loadtxt(from_path+'t_tot_s.txt')/2/np.pi
  #x0  = np.loadtxt(from_path+'x_tot_s.txt')/2/np.pi
  #y0  = np.loadtxt(from_path+'y_tot_s.txt')/2/np.pi
  px0 = np.loadtxt(from_path+'px_tot_s.txt')
  py0 = np.loadtxt(from_path+'py_tot_s.txt')
  t0  = np.reshape(t0,(part_number,nsteps))
  #x0  = np.reshape(x0,(part_number,nsteps))
  #y0  = np.reshape(y0,(part_number,nsteps))
  px0 = np.reshape(px0,(part_number,nsteps))
  py0 = np.reshape(py0,(part_number,nsteps))
  gg0 = (px0**2+py0**2+1)**0.5
  ww0 = np.zeros_like(gg0)+1


  axis_time = np.linspace(0,nsteps*3.333,nsteps)  
  axis_gg   = np.linspace(400,11600,15)*0.51*1e-3
  data_rr   = np.zeros([15,nsteps])

  for i in range(nsteps):
      axis_temp, data_rr[:,i] = pxpy_to_energy(gg0[:,i], ww0[:,i]) 
      print(np.max(data_rr[:,i]))

  fig = plt.figure()
  ax=fig.add_subplot(111, projection='3d')
  axis_gg, axis_time = np.meshgrid(axis_gg,axis_time[::11])
  GG   = axis_gg.flatten()
  Time = axis_time.flatten()
  Zi   = np.zeros((data_rr[:,::11].T).size)

  dx = .24 * np.ones((data_rr[:,::11].T).size)
  dy = 20. * np.ones((data_rr[:,::11].T).size)
  dz = (data_rr[:,::11].T).flatten()
  ax.bar3d(GG, Time, Zi, dx, dy, dz, color='deepskyblue')

  for t in ax.xaxis.get_major_ticks(): t.label.set_fontsize(font_size)
  for t in ax.yaxis.get_major_ticks(): t.label.set_fontsize(font_size)
  for t in ax.zaxis.get_major_ticks(): t.label.set_fontsize(font_size)

  #ax.set_xlim([x_start/20-5,x_stop/20-5])
  #ax.set_ylim([0,500])
  ax.set_xlabel('\n\n'+r'$\varepsilon_e$'+' [GeV]',fontdict=font)
  ax.set_ylabel('\n\n'+'t [fs]',fontdict=font)
  ax.set_zlabel('\n\n'+'dN/d'+r'$\varepsilon_e$'+' [A.U.]',fontdict=font)
  #ax.zaxis.set_scale('log')
  #ax.set_zlim([1,50000])
  #ax.zaxis.set_scale('log')
  print(np.shape(GG),np.shape(Time),np.shape(Zi),np.shape(dx),np.shape(dy),np.shape(dz))

  #ax.grid(False)
  #ax.xaxis.pane.set_edgecolor('black')
  #ax.yaxis.pane.set_edgecolor('black')
  #ax.zaxis.pane.set_edgecolor('black')
  #ax.xaxis.pane.fill = False
  #ax.yaxis.pane.fill = False
  #ax.zaxis.pane.fill = False
  #ax.grid(linestyle='None', linewidth='0.5', color='white')
  plt.show()
  plt.subplots_adjust(top=0.93, bottom=0.16, left=0.1, right=0.95, hspace=0.10, wspace=0.05)

  fig = plt.gcf()
  fig.set_size_inches(10, 8.5)
  #fig.set_size_inches(5, 4.5)
  ax.view_init(elev=30, azim=135) #Reproduce view
  fig.savefig('./3d_bar_rr_T150.png',format='png',dpi=160)
  plt.close("all")
    
          
          
