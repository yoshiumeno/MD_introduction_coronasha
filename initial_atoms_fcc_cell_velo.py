import math
import random

lo = 4.04 # Lattice constant of Al
# how many unit cells?
irepx = 3
irepy = 3
irepz = 3

kb = 1.380649e-23 # Boltzmann constant
temp = 300.0 # temperature
wm = 27.0 * 1.67e-27 # Al mass (kg)

def genv(t,m): # to generate random velocity (ang/s)
    fac = math.sqrt(kb * t / m)
    xx = -6.0
    for i in range(12):
        xx += random.random()
    return fac * xx * 1.0e10

f = open('conf_fcc_cell.d','w')

# Atoms in unit cell
nu = 4
rux = []
ruy = []
ruz = []
rux += [0.0]
ruy += [0.0]
ruz += [0.0]
rux += [lo / 2.0]
ruy += [lo / 2.0]
ruz += [0.0]
rux += [lo / 2.0]
ruy += [0.0]
ruz += [lo / 2.0]
rux += [0.0]
ruy += [lo / 2.0]
ruz += [lo / 2.0]

# Unit cell size
lx = lo * irepx
ly = lo * irepy
lz = lo * irepz

# prepare array
rx = []; ry = []; rz = []; vx = []; vy = []; vz = []

icnt = 0
for ix in range(irepx):
    for iy in range(irepy):
        for iz in range(irepz):
            for ii in range(nu):
                rxx = rux[(icnt % nu)] + lo * ix
                ryy = ruy[(icnt % nu)] + lo * iy
                rzz = ruz[(icnt % nu)] + lo * iz
                rx += [rxx]
                ry += [ryy]
                rz += [rzz]
                vx += [genv(temp,wm)]
                vy += [genv(temp,wm)]
                vz += [genv(temp,wm)]
                icnt += 1
natom = icnt # number of total atoms

print('# of atoms = ',natom)
print('%12.8f %12.8f %12.8f' % (lx,ly,lz), file=f)
for i in range(natom):
    print('%d %12.8f %12.8f %12.8f %12.8e %12.8e %12.8e' % (i,rx[i],ry[i],rz[i],vx[i],vy[i],vz[i]))
    print('%12.8f %12.8f %12.8f %12.8e %12.8e %12.8e' % (rx[i],ry[i],rz[i],vx[i],vy[i],vz[i]), file=f)
f.close()

