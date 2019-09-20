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



if __name__ == "__main__":
  part_number = 50000
  from_path = './p50000_no_T150/'
  nsteps      = int(sum(1 for line in open(from_path+'t_tot_s.txt'))/part_number)


  from_path = './p50000_no_T150/'
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
  ww2 = np.zeros_like(gg2)+1



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
  ww1 = np.zeros_like(gg1)+1



  axis_time = np.linspace(0,nsteps*3.333,nsteps)  
  axis_gg   = np.linspace(50,19950,200)  
  data_no   = np.zeros([200,nsteps])
  data_rr   = np.zeros([200,nsteps])
  data_qe   = np.zeros([200,nsteps])

  for i in range(nsteps):
      axis_temp, data_no[:,i] = pxpy_to_energy(gg0[:,i], ww0[:,i]) 
      axis_temp, data_qe[:,i] = pxpy_to_energy(gg2[:,i], ww2[:,i]) 
      axis_temp, data_rr[:,i] = pxpy_to_energy(gg1[:,i], ww1[:,i]) 
      print(np.max(data_rr[:,i]))

  plt.subplot(2,2,3)        
  x,y=np.meshgrid(axis_gg,axis_time)
  levels = np.logspace(0, np.log10(1e4), 41)
  plt.pcolormesh(x*0.51*1e-3, y, data_no.T, norm=colors.LogNorm(vmin=1, vmax=1e4), cmap='viridis')
#plt.pcolormesh(x, y, ex.T, norm=mpl.colors.Normalize(vmin=0,vmax=100,clip=True), cmap=cm.cubehelix_r)
#  plt.axis([x.min(), x.max(), y.min(), y.max()])
#### manifesting colorbar, changing label and axis properties ####
#  cbar=plt.colorbar(pad=0.01)#ticks=[np.min(ex), -eee/2, 0, eee/2, np.min()])
#  cbar.ax.set_yticklabels(cbar.ax.get_yticklabels(), fontsize=font_size)
#  cbar.set_label('dN/dE [A.U.]',fontdict=font)
#  a0=200.0
#  alpha=np.linspace(-3.5,0.5,501)
  plt.xlabel(r'$\epsilon_e$'+' [GeV]',fontdict=font)
  plt.ylabel('t [fs]',fontdict=font)
  plt.xticks(fontsize=font_size); 
  plt.yticks(fontsize=font_size);
#  plt.title('dN/dE for no RR', fontsize=font_size)
  plt.xlim(0,5.9)
  plt.ylim(0,500)
  plt.text(3.5,50,'No RR',fontdict=font)

  plt.subplot(2,2,1)        
  x,y=np.meshgrid(axis_gg,axis_time)
  levels = np.logspace(0, np.log10(1e4), 41)
  plt.pcolormesh(x*0.51*1e-3, y, data_rr.T, norm=colors.LogNorm(vmin=1, vmax=1e4), cmap='viridis')
#plt.pcolormesh(x, y, ex.T, norm=mpl.colors.Normalize(vmin=0,vmax=100,clip=True), cmap=cm.cubehelix_r)
#  plt.axis([x.min(), x.max(), y.min(), y.max()])
#### manifesting colorbar, changing label and axis properties ####
#  cbar=plt.colorbar(pad=0.01)#ticks=[np.min(ex), -eee/2, 0, eee/2, np.min()])
#  cbar.ax.set_yticklabels(cbar.ax.get_yticklabels(), fontsize=font_size)
#  cbar.set_label('dN/dE [A.U.]',fontdict=font)
#  a0=200.0
#  alpha=np.linspace(-3.5,0.5,501)
  plt.xlabel(r'$\epsilon_e$'+' [GeV]',fontdict=font)
  plt.ylabel('t [fs]',fontdict=font)
  plt.xticks(fontsize=font_size); 
  plt.yticks(fontsize=font_size);
#  plt.title('dN/dE for QED RR', fontsize=font_size)
  plt.xlim(0,5.9)
  plt.ylim(0,500)
  plt.text(3.5,50,'Classic RR',fontdict=font)

  plt.subplot(2,2,2)        
  x,y=np.meshgrid(axis_gg,axis_time)
  levels = np.logspace(0, np.log10(1e4), 41)
  plt.pcolormesh(x*0.51*1e-3, y, data_qe.T, norm=colors.LogNorm(vmin=1, vmax=1e4), cmap='viridis')
#plt.pcolormesh(x, y, ex.T, norm=mpl.colors.Normalize(vmin=0,vmax=100,clip=True), cmap=cm.cubehelix_r)
#  plt.axis([x.min(), x.max(), y.min(), y.max()])
#### manifesting colorbar, changing label and axis properties ####
#  cbar=plt.colorbar(pad=0.01)#ticks=[np.min(ex), -eee/2, 0, eee/2, np.min()])
#  cbar.ax.set_yticklabels(cbar.ax.get_yticklabels(), fontsize=font_size)
#  cbar.set_label('dN/dE [A.U.]',fontdict=font)
#  a0=200.0
#  alpha=np.linspace(-3.5,0.5,501)
  plt.xlabel(r'$\epsilon_e$'+' [GeV]',fontdict=font)
  plt.ylabel('t [fs]',fontdict=font)
  plt.xticks(fontsize=font_size); 
  plt.yticks(fontsize=font_size);
#  plt.title('dN/dE for QED RR', fontsize=font_size)
  plt.xlim(0,5.9)
  plt.ylim(0,500)
#plt.title('electron at y='+str(round(y[n,0]/2/np.pi,4)),fontdict=font)
  plt.text(3.5,50,'QED RR',fontdict=font)


  plt.subplot(2,2,4)        
  x,y=np.meshgrid(axis_gg,axis_time)
  levels = np.logspace(0, np.log10(1e4), 41)
  plt.plot((axis_temp*0.51*1e-3), data_no[:,-10], '-', color='royalblue', linewidth=3,label='No RR')
  plt.plot((axis_temp*0.51*1e-3), data_rr[:,-10], '-', color='mediumseagreen', linewidth=3,label='Classic RR')
  plt.plot((axis_temp*0.51*1e-3), data_qe[:,-10], '-', color='tomato', linewidth=3,label='QED RR')

#plt.pcolormesh(x, y, ex.T, norm=mpl.colors.Normalize(vmin=0,vmax=100,clip=True), cmap=cm.cubehelix_r)
#  plt.axis([x.min(), x.max(), y.min(), y.max()])
#### manifesting colorbar, changing label and axis properties ####
#  cbar=plt.colorbar(pad=0.01)#ticks=[np.min(ex), -eee/2, 0, eee/2, np.min()])
#  cbar.ax.set_yticklabels(cbar.ax.get_yticklabels(), fontsize=font_size)
#  cbar.set_label('dN/dE [A.U.]',fontdict=font)
#  a0=200.0
#  alpha=np.linspace(-3.5,0.5,501)
  plt.xlabel(r'$\epsilon_e$'+' [GeV]',fontdict=font)
  plt.ylabel('dN/d'+r'$\epsilon_e$'+' [A.U.]',fontdict=font)
  plt.xticks(fontsize=font_size); 
  plt.yticks(fontsize=font_size);
#  plt.title('dN/dE for QED RR', fontsize=font_size)
  plt.xlim(0,5.9)
  plt.yscale('log')
  plt.legend(loc='upper center',fontsize=font_size2)
#  plt.ylim(0,500)
  plt.grid()
#plt.title('electron at y='+str(round(y[n,0]/2/np.pi,4)),fontdict=font)
  plt.text(3.,2,'t = 500 fs',fontdict=font)


  plt.subplots_adjust(top=0.98, bottom=0.12, left=0.1, right=0.98, hspace=0.20, wspace=0.23)

  fig = plt.gcf()
  fig.set_size_inches(14, 12)
#fig.set_size_inches(5, 4.5)
  fig.savefig('./wrap_dN_dE_time_T100.png',format='png',dpi=160)
  plt.close("all")

  
          
          
