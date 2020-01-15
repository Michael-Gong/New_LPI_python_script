import sdf
#import matplotlib
#matplotlib.use('agg')
#import matplotlib.pyplot as plt
import numpy as np
#from numpy import ma
#from matplotlib import colors, ticker, cm
#from matplotlib.mlab import bivariate_normal
#from optparse import OptionParser
#import os
#from colour import Color

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
#print 'electric field unit: '+str(exunit)
#print 'magnetic field unit: '+str(bxunit)
#print 'density unit nc: '+str(denunit)

font = {'family' : 'monospace',  
        'color'  : 'black',  
    	'weight' : 'normal',  
        'size'   : 20,  
       }  

from_path = './PW_w020/'
to_path = './PW_w020/'
part_name='proton'
part_mass=1836.0*1

data = sdf.read(from_path+"0001.sdf",dict=True)
grid_x = data['Grid/Particles/'+part_name].data[0]/wavelength
grid_y = data['Grid/Particles/'+part_name].data[1]/wavelength
px = data['Particles/Px/'+part_name].data/(part_mass*m0*v0)
py = data['Particles/Py/'+part_name].data/(part_mass*m0*v0)
theta = np.arctan2(py,px)*180.0/np.pi
gg = (px**2+py**2+1)**0.5
Ek = (gg-1)*part_mass*0.51

part13_id = data['Particles/ID/'+part_name].data
#part13_id = part13_id[ (Ek>225) & (abs(grid_y)<8) & (Ek<245)]
part13_id = part13_id[abs(grid_y)<8]
choice = np.random.choice(range(part13_id.size), 100000, replace=False)
part13_id = part13_id[choice]
print('part13_id size is ',part13_id.size,' max ',np.max(part13_id),' min ',np.min(part13_id))

######### Parameter you should set ###########
start   =  1  # start time
stop    =  19  # end time
step    =  1  # the interval or step

#  if (os.path.isdir('jpg') == False):
#    os.mkdir('jpg')
######### Script code drawing figure ################
part_id = part13_id
for n in range(start,stop+step,step):
    #### header data ####
    data = sdf.read(from_path+str(n).zfill(4)+".sdf",dict=True)
    header=data['Header']
    time=header['time']
    part_id = np.intersect1d(data['Particles/ID/'+part_name].data, part_id)
    print('Particle_ID size is ',part_id.size,' max ',np.max(part_id),' min ',np.min(part_id))

print('After intersecting with all.sdf')
print('Particle_ID size is ',part_id.size,' max ',np.max(part_id),' min ',np.min(part_id))


######### Parameter you should set ###########
start   =  1  # start time
stop    =  19  # end time
step    =  1  # the interval or step

px_d = np.zeros([part_id.size,stop-start+1])
py_d = np.zeros([part_id.size,stop-start+1])
xx_d = np.zeros([part_id.size,stop-start+1])
yy_d = np.zeros([part_id.size,stop-start+1])
for n in range(start,stop+step,step):
    #### header data ####
    data = sdf.read(from_path+str(n).zfill(4)+".sdf",dict=True)
    px = data['Particles/Px/'+part_name].data/(part_mass*m0*v0)
    py = data['Particles/Py/'+part_name].data/(part_mass*m0*v0)
    grid_x = data['Grid/Particles/'+part_name].data[0]/wavelength
    grid_y = data['Grid/Particles/'+part_name].data[1]/wavelength
    temp_id =  data['Particles/ID/'+part_name].data

    px = px[np.in1d(temp_id,part_id)]
    py = py[np.in1d(temp_id,part_id)]
    grid_x = grid_x[np.in1d(temp_id,part_id)]
    grid_y = grid_y[np.in1d(temp_id,part_id)]
    temp_id = temp_id[np.in1d(temp_id,part_id)]

    for ie in range(part_id.size):
        px_d[ie,n-start] = px[temp_id==part_id[ie]]
        py_d[ie,n-start] = py[temp_id==part_id[ie]]
        xx_d[ie,n-start] = grid_x[temp_id==part_id[ie]]
        yy_d[ie,n-start] = grid_y[temp_id==part_id[ie]]
    print('finised '+part_name+' '+str(round(100.0*(n-start+step)/(stop-start+step),4))+'%')

np.savetxt(to_path+part_name+'_px.txt',px_d)
np.savetxt(to_path+part_name+'_py.txt',py_d)
np.savetxt(to_path+part_name+'_xx.txt',xx_d)
np.savetxt(to_path+part_name+'_yy.txt',yy_d)
