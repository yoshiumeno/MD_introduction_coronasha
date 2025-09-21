import numpy as np
from scipy import interpolate

fin  = open('EAM_Al_den.orig','r') # Input file
fout = open('EAM_Al_den.dat','w')  # Output file

# number of fine mesh grids to create
eam_mesh_fine = 10000
# read from fin
eam_x = []
eam_y = []
for line in fin:
    xy = line.split()
    if (xy[0] != "#"):
        eam_x += [float(xy[0])]
        eam_y += [float(xy[1])]
eam_mesh = len(eam_x)
eam_xmin = eam_x[0]
eam_xmax = eam_x[-1]
print('Number of mesh grids in file: ',eam_mesh)
print('Range of x: %f --- %f' % (eam_xmin, eam_xmax))
print('Number of fine mesh grids to create: '\
,eam_mesh_fine)
# get spline interpolation
f = interpolate.interp1d(eam_x, eam_y,kind="cubic")
x = np.linspace(eam_x[0],eam_x[-1]\
,num=eam_mesh_fine,endpoint=True)
y = f(x)
# calculate slope of y
dy = np.zeros(eam_mesh_fine)
dy[ 0] = (y[ 1]-y[ 0])/(x[ 1]-x[ 0])
dy[-1] = (y[-1]-y[-2])/(x[-1]-x[-2])
for i in range(1,eam_mesh_fine-1):
    dy[i] = (y[i+1]-y[i-1])/(x[i+1]-x[i-1])
# output to fout
print(x[0],x[-1],file=fout)
for i in range(eam_mesh_fine):
    print(y[i],dy[i],file=fout)