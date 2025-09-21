import math

ep = 0.2703; al = 1.1646; ro = 3.253
ev = 1.6021892e-19
ang = 1.0e-10

def v(r):
    rang = r/ang
    return ep*(math.exp(-2.0*al*(rang-ro))\
    -2.0*math.exp(-al*(rang-ro))) *ev # [J]

def vp(r):
    rang = r/ang
    return -2.0*al*ep*(math.exp(-2.0*al*(rang-ro))\
    -math.exp(-al*(rang-ro))) *ev/ang # [J/m]

def write_cfg(fname,n,rx,ry,rz,lx,ly,lz):
    fout = open(fname,'w')
    print('Number of particles = %d'%(n), file=fout)
    print('A = 1.0 Angstrom', file=fout)
    print('H0(1,1) = %f'%(lx/ang), file=fout)
    print('H0(1,2) = 0.0', file=fout)
    print('H0(1,3) = 0.0', file=fout)
    print('H0(2,1) = 0.0', file=fout)
    print('H0(2,2) = %f'%(ly/ang), file=fout)
    print('H0(2,3) = 0.0', file=fout)
    print('H0(3,1) = 0.0', file=fout)
    print('H0(3,2) = 0.0', file=fout)
    print('H0(3,3) = %f'%(lz/ang), file=fout)
    for i in range(n):
        print('1.0 Al %e %e %e %e %e %e' \
        % (rx[i]/lx,ry[i]/ly,rz[i]/lz,0.0,0.0,0.0), file=fout)


# declaration of arrays
rx = []; ry = []; rz = [] # m
vx = []; vy = []; vz = [] # m/s
fx = []; fy = []; fz = [] # N
epot = []
ekin = []
dt = 1.0e-15 # s
wm = 26.98 * 1.0e-3 / 6.022e23  #  Al mass (kg)
kb = 1.380649e-23 # Boltzmann constant
rc = 10.0*ang # cutoff radius (m)
nei  = [] # neighbor list
rcbk = rc + 1.0*ang # cutoff for neighbor list (bookkeep)
rcbk2 = rcbk**2
inbk = 10 # interval for bookkeep
temp_set = 300.0 # target temperature
vscale = 1.0 # velocity scaling factor
ifvscale = False # velocity scaling on/off

rc2 = rc**2
delt = 1.0e-3*ang # for avoiding self-interaction
del2 = delt**2

# read initial position and velocity
f = open('initial_3d_cell.d', 'r')
fo= open('ene.d','w')
n = 0

line = f.readline()
xy = line.split()
lx = float(xy[0])*ang
ly = float(xy[1])*ang
lz = float(xy[2])*ang
print("cell size: ",lx/ang,ly/ang,lz/ang)
for line in f:
    xy = line.split()
    rx += [float(xy[0])*ang]
    ry += [float(xy[1])*ang]
    rz += [float(xy[2])*ang]
    vx += [float(xy[3])*ang]
    vy += [float(xy[4])*ang]
    vz += [float(xy[5])*ang]
    fx += [0]; fy += [0]; fz += [0]
    epot += [0]
    ekin += [0]
    nei  += [[0]]
    n += 1 # number of total atoms
print("number of atoms = ",n)

step = 0
stepend = 1000

while step < stepend:
    # make neighbor list (bookkeep)
    if (step % inbk == 0):
        # Keep atoms in cell
        for i in range(n):
            if (rx[i]>lx): rx[i] -= lx
            if (ry[i]>ly): ry[i] -= ly
            if (rz[i]>lz): rz[i] -= lz
            if (rx[i]<0):  rx[i] += lx
            if (ry[i]<0):  ry[i] += ly
            if (rz[i]<0):  rz[i] += lz

        rep = [-1,0,1] # replica cell index
        for i in range(n):
            nei[i].clear()
            for j in range(n):
                for ix in rep:
                    for iy in rep:
                        for iz in rep:
                            rxj = rx[j] + lx * float(ix)
                            ryj = ry[j] + ly * float(iy)
                            rzj = rz[j] + lz * float(iz)
                            drx = rx[i] - rxj
                            dry = ry[i] - ryj
                            drz = rz[i] - rzj
                            rr2 = drx**2 + dry**2 + drz**2
                            if (rr2 < rcbk2 and rr2 > del2):
                                nei[i] += [[j, ix, iy, iz]]
    # Verlet(1)
    for i in range(n):
        rx[i] += dt * vx[i] + (dt*dt/2.0) * fx[i] / wm
        ry[i] += dt * vy[i] + (dt*dt/2.0) * fy[i] / wm
        rz[i] += dt * vz[i] + (dt*dt/2.0) * fz[i] / wm
        vx[i] += dt/2.0 * fx[i] / wm
        vy[i] += dt/2.0 * fy[i] / wm
        vz[i] += dt/2.0 * fz[i] / wm
    # Force and energy
    for i in range(n):
        fx[i] = 0; fy[i] = 0; fz[i] = 0
        epot[i] = 0
    for i in range(n):
        for jj in range(len(nei[i])):
            j  = nei[i][jj][0]
            ix = nei[i][jj][1]
            iy = nei[i][jj][2]
            iz = nei[i][jj][3]
            rxj = rx[j] + lx * float(ix)
            ryj = ry[j] + ly * float(iy)
            rzj = rz[j] + lz * float(iz)
            drx = rx[i] - rxj
            dry = ry[i] - ryj
            drz = rz[i] - rzj
            rr2 = drx**2 + dry**2 + drz**2
            if (rr2 < rc2 and rr2 > del2):
                rr = math.sqrt(rr2)
                fx[i] += -vp(rr)/rr*drx
                fy[i] += -vp(rr)/rr*dry
                fz[i] += -vp(rr)/rr*drz
                epot[i] += v(rr)/2.0
    # Verlet(2)
    for i in range(n):
        vx[i] += dt/2.0 * fx[i] / wm
        vy[i] += dt/2.0 * fy[i] / wm
        vz[i] += dt/2.0 * fz[i] / wm
        ekin[i] = 0.5 * wm * (vx[i]**2 + vy[i]**2 + vz[i]**2)

    # Calc pot/kin energy
    epot_t = 0.0
    ekin_t = 0.0
    for i in range(n):
        epot_t += epot[i]
        ekin_t += ekin[i]
    temp = ekin_t *2.0 / (3.0 * n * kb)
    
    step += 1
    # Temperature control by velocity scaling
    if (ifvscale):
        vscale = math.sqrt(temp_set/temp);
        for i in range(n):
            vx[i] *= vscale; vy[i] *= vscale; vz[i] *= vscale;

    # output to file
    print('%d %e %e %e' % (step,epot_t,ekin_t,epot_t+ekin_t), file=fo)
    cfg_interval = 10 # cfg file is written every XX steps
    if ((step % cfg_interval)==0):
        write_cfg('out'+'%04d'%(step//cfg_interval)\
        +'.cfg',n,rx,ry,rz,lx,ly,lz)
    # output to console
    print('step:%d E_p:%e E_k:%e E:%e  temp:%f' \
    % (step,epot_t,ekin_t,epot_t+ekin_t,temp))
