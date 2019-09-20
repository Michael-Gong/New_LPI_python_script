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
        'size'   : 28,
       }
font2 = {'family' : 'monospace',
        'color'  : 'black',
        'weight' : 'normal',
        'size'   : 15,
       }
font_size = 28
font_size_2 = 15

#plt.rc('text', usetex=True)
#plt.rc('font', family='serif')

upper = matplotlib.cm.jet(np.arange(256))
lower = np.ones((int(256/4),4))
for i in range(3):
  lower[:,i] = np.linspace(1, upper[0,i], lower.shape[0])
cmap = np.vstack(( lower, upper ))
mycolor_jet = matplotlib.colors.ListedColormap(cmap, name='myColorMap', N=cmap.shape[0])

upper = matplotlib.cm.viridis(np.arange(256))
lower = np.ones((int(256/4),4))
for i in range(3):
  lower[:,i] = np.linspace(1, upper[0,i], lower.shape[0])
cmap = np.vstack(( lower, upper ))
mycolor_viridis = matplotlib.colors.ListedColormap(cmap, name='myColorMap', N=cmap.shape[0])

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
    binsize = 240
    theta_grid = np.linspace(-119.5,119.5,240)
    theta_bin  = np.linspace(-120,120,241)
    theta_value = np.zeros_like(theta_grid) 
    for i in range(binsize):
#        if i == binsize-1:
#            en_value[i] = sum(weight[en_bin[i]<=gamma])
#        else:
            theta_value[i] = sum(weight[ (theta_bin[i]<=theta) & (theta<theta_bin[i+1]) ])
    return (theta_grid, theta_value)






if __name__ == "__main__":
  part_number = 50000
  from_path = './p50000_no_T150/'
  nsteps      = int(sum(1 for line in open(from_path+'t_tot_s.txt'))/part_number)

  ntheta = 270
  ngg    = 120

  from_path_list = ['./p50000_no_T150/','./p50000_rr_T150/','./p50000_qe_T150/']
  #from_path_list = ['./Data_qe_T500_p50000_try/']

  for i in range(np.size(from_path_list)):
      from_path = from_path_list[i] #'./Data_qe_T050_p50000/'
      to_path   = from_path
      t0  = np.loadtxt(from_path+'t_tot_s.txt')/2/np.pi
      px0 = np.loadtxt(from_path+'px_tot_s.txt')
      py0 = np.loadtxt(from_path+'py_tot_s.txt')
      t0  = np.reshape(t0,(part_number,nsteps))
      px0 = np.reshape(px0,(part_number,nsteps))
      py0 = np.reshape(py0,(part_number,nsteps))
      gg0 = (px0**2+py0**2+1)**0.5*0.51e-3
      ww0 = np.zeros_like(gg0)+1
      ww0 = np.zeros_like(gg0)+gg0
      theta0 = np.arctan2(py0,px0)
    
      theta_edges = np.linspace(-np.pi,np.pi,  ntheta +1)
      gg_edges    = np.linspace(0.1, 6,  ngg +1)
    
      theta_edges_1 = np.linspace(-np.pi,np.pi,ntheta)
      gg_edges_1    = np.linspace(0.1, 6,  ngg)
       
      for n in range(np.size(t0[0,:])):
          H, _, _ = np.histogram2d(gg0[:,n], theta0[:,n], [gg_edges, theta_edges], weights=gg0[:,n])
          print('Max H:',np.max(H))
          Theta, R = np.meshgrid(theta_edges_1,gg_edges_1)
          H_temp = np.sum(H[:,:]*R,0)
          print('averaged |theta|=',np.sum(H_temp*abs(theta_edges_1))/np.sum(H_temp)/np.pi*180)
          fig, ax = plt.subplots(1, 1, subplot_kw=dict(projection='polar'))
          ax.set_facecolor('whitesmoke')
          levels = np.logspace(1,5, 101)
          H[H<0.01] = np.nan
          img=ax.pcolormesh(Theta,  R,  H, norm=colors.LogNorm(vmin=0.01, vmax=1e3), cmap='viridis')
        #  cax = fig.add_axes([0.68,0.97,0.25,0.02])
        #  cbar=fig.colorbar(img,cax=cax, ticks=[1e3,1e5],orientation='horizontal')
        #  cbar.ax.set_xticklabels(cbar.ax.get_xticklabels(), fontsize=font_size_2)
        #  cbar.set_label(r'dI/d$\theta$dE [A.U.]',fontdict=font2)
        #  ax.tick_params(axis="y", pad=25)
          ax.tick_params(axis="x", pad=10)
        #  ax.set_xticks([])
          if (i%3 != 2):
              ax.set_xticklabels([])
              #ax.set_xlim(10,50)
              #ax.set_ylim(0.,1.)
          ax.set_xlabel(r'$\theta\ [^o]$',fontdict=font)
        #  ax.set_rlim(1e-1,1e3)
        #  ax.set_rmax(1e3)
          l_r = np.array([0,1,2,3])
          ax.set_rticks(l_r+1)
          ax.set_yticklabels([])
        #  ax.set_yticklabels(['$10^%d$' % x for x in (l_r+1)])
          ax.set_rlim(0, 6)
          ax.set_rlabel_position(90) 
        #  ax.set_rscale('log')
        #  ax.set_rscale('log')
        #  ax.set_thetamin(-90)
        #  ax.set_thetamax(90)
        #  ax.set_yticklabels([0.1,1,10,100,1000])
          ax.set_xticklabels([0,90,180,270])
        
          #ax.set_theta_zero_location('N')
            #  ax.set_ylabel(r'$\theta\ [^o]$',fontdict=font)
          ax.tick_params(axis='x',labelsize=font_size) 
          ax.tick_params(axis='y',labelsize=font_size_2)
          #ax.set_title('proton_angular_time='+str(time1), va='bottom', y=1., fontsize=20)
            #  plt.text(-100,650,' t = '++' fs',fontdict=font)
          ax.grid(True,linestyle='--',linewidth=1.5,color='grey')
        #plt.pcolormesh(x, y, ex.T, norm=mpl.colors.Normalize(vmin=0,vmax=100,clip=True), cmap=cm.cubehelix_r)
        #  plt.axis([x.min(), x.max(), y.min(), y.max()])
        #### manifesting colorbar, changing label and axis properties ####
        #  cbar=plt.colorbar(pad=0.01)#ticks=[np.min(ex), -eee/2, 0, eee/2, np.min()])
        #  cbar.ax.set_yticklabels(cbar.ax.get_yticklabels(), fontsize=font_size)
        #  cbar.set_label('dN/dE [A.U.]',fontdict=font)
        #  a0=200.0
        #  alpha=np.linspace(-3.5,0.5,501)
        #  plt.xlabel(r'$\theta$'+' [degree]',fontdict=font)
        #  plt.ylabel('time [fs]',fontdict=font)
         # plt.xticks([-135,-90,-45,0,45,90,135],fontsize=font_size); 
          #plt.yticks([0,500,1000,1500],fontsize=font_size);
        #  plt.title(r'$dN/d\theta$'+' for no RR', fontsize=font_size)
        #  plt.xlim(-120,120)
        #  plt.ylim(0,1650)
        #plt.title('electron at y='+str(round(y[n,0]/2/np.pi,4)),fontdict=font)
        
          plt.subplots_adjust(top=0.90, bottom=0.11, left=0.1, right=0.93, hspace=0.10, wspace=0.05)
        
          fig = plt.gcf()
          fig.set_size_inches(6., 6.)
        #fig.set_size_inches(5, 4.5)
          fig.savefig(to_path+'theta_en_dist_'+to_path[7:-1]+'_'+str(n).zfill(4)+'.png',format='png',dpi=160)
          plt.close("all")

      
          
          
