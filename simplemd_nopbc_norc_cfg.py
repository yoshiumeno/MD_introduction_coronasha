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


# read initial position and velocity
f = open('initial_3d.d', 'r')
fo= open('ene.d','w')
n = 0
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
    n += 1 # number of total atoms
print("number of atoms = ",n)
# cell size must be set for cfg format
lx = max(max(rx),10*ang)
ly = max(max(ry),10*ang)
lz = max(max(rz),10*ang)

step = 0
stepend = 1000

while step < stepend:
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
        for j in range(n):
            if (i != j):
                rr = math.sqrt((rx[i]-rx[j])**2 \
                + (ry[i]-ry[j])**2 + (rz[i]-rz[j])**2)
                drx = rx[i] - rx[j]
                dry = ry[i] - ry[j]
                drz = rz[i] - rz[j]
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
    # output to file
    print('%d %e %e %e' % (step,epot_t,ekin_t,epot_t+ekin_t), file=fo)
    cfg_interval = 10 # cfg file is written every XX steps
    if ((step % cfg_interval)==0):
        write_cfg('out'+'%04d'%(step//cfg_interval)\
        +'.cfg',n,rx,ry,rz,lx,ly,lz)
    # output to console
    print('step:%d E_p:%e E_k:%e E:%e  temp:%f' \
    % (step,epot_t,ekin_t,epot_t+ekin_t,temp))
