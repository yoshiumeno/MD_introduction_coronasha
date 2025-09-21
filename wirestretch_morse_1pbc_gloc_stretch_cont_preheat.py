import math
import numpy as np

ep = 0.2703; al = 1.1646; ro = 3.253
ev = 1.6021892e-19
ang = 1.0e-10
# for stress
mpa = 1.0e-6

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
def read_cfg(fname,rx,ry,rz,lxyz):
    fin = open(fname,'r')
    line = fin.readline(); xy = line.split()
    if (int(xy[4]) != n):
        print('Warning: different number of atoms found')
    line = fin.readline() # skip 2nd line
    for i in range(3):
        for j in range(3):
            line = fin.readline(); xy = line.split()
            if (i==0 and j==0): lxyz[0] = float(xy[2])*ang
            if (i==1 and j==1): lxyz[1] = float(xy[2])*ang
            if (i==2 and j==2): lxyz[2] = float(xy[2])*ang
    print('in',lxyz)
    for i in range(n):
        line = fin.readline(); xy = line.split()
        rx[i] = float(xy[2])*lxyz[0]
        ry[i] = float(xy[3])*lxyz[1]
        rz[i] = float(xy[4])*lxyz[2]

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
ifgloc = True # gloc on/off
# quasistatic deformation setting
fnl_eps = 0.30 # maximum strain
del_eps = 0.01 # strain increment
rlx_tol = 0.02 # force tolerance in relaxation (eV/A)
rlx_maxstep = 3000 # maximum iteration for relaxation
# bulk deformation setting
#ini_eps = -0.1 # initial strain
#fnl_eps =  0.1 # final strain
#del_eps = 0.01 # strain increment
# control perodic boundaries
# pbc -> [-1,0,1] / no pbc -> [0]
repx = [0] # no pbc in x
repy = [0] # no pbc in y
repz = [-1,0,1] # pbc in z
ifmd = True # move atoms on/off
ifcont_deform = True # T:continue, F:from initial
cont_cfg_fname = 'eps00230_unrlx.cfg'
# relax preheat (fluctuation by NVT) setting
ifrlx_preheat = True
rlx_preheat_temp = 100
rlx_preheat_step = 50
if (ifrlx_preheat):
    rlx_maxstep += rlx_preheat_step

rc2 = rc**2
delt = 1.0e-3*ang # for avoiding self-interaction
del2 = delt**2

# read initial position and velocity
f = open('initial_3d_cell.d', 'r')
fo= open('ene.d','w')
fs= open('stress.d','w')
fl= open('lat_ene_sig.d','w')
if (ifcont_deform):
    fd= open('ss.d','a') # output stress-strain (add)
else:
    fd= open('ss.d','w') # output stress-strain (new)

n = 0
line = f.readline()
xy = line.split()
lx = float(xy[0])*ang
ly = float(xy[1])*ang
lz = float(xy[2])*ang
print('initial cell size:',lx/ang,ly/ang,lz/ang)
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

# for atomic stress
satom = np.zeros((n,3,3))
sig = np.zeros((3,3))
vol = lx*ly*lz
vol_a = vol/n

step = 0
stepend = 1000

# setting before deformation
#eps = ini_eps
#lx0 = lx; ly0 = ly; lz0 = lz # L in file
# setting before stretch
eps = 0.0
lz0 = lz
lzmax = lz0 * (1.0 + fnl_eps)
# continue from cont_cfg_fname
if (ifcont_deform):
    lxyz = [0.0,0.0,0.0]
    read_cfg(cont_cfg_fname,rx,ry,rz,lxyz)
    lx = lxyz[0]; ly = lxyz[1]; lz = lxyz[2]
    vol = lx*ly*lz
    vol_a = vol/n
    eps = lz/lz0 - 1.0
    
#while step < stepend:
while lz <= lzmax:
#while eps <= fnl_eps:
#    eps = round(eps,10)
#    # cell deformation
#    lx_next = lx0 * (1.0 + eps)
#    ly_next = ly0 * (1.0 + eps)
#    lz_next = lz0 * (1.0 + eps)
#    for i in range(n):
#        rx[i] *= lx_next/lx
#        ry[i] *= ly_next/ly
#        rz[i] *= lz_next/lz
#    lx = lx_next
#    ly = ly_next
#    lz = lz_next
#    print('eps:%f cell size: %f %f %f'\
#    % (eps,lx/ang,ly/ang,lz/ang))
#    vol = lx*ly*lz
#    vol_a = vol/n
#    eps += del_eps
    # make neighbor list (bookkeep)
    if (step % inbk == 0):
        # Keep atoms in cell
        for i in range(n):
            if (len(repx)>1 and rx[i]>lx): rx[i] -= lx
            if (len(repy)>1 and ry[i]>ly): ry[i] -= ly
            if (len(repz)>1 and rz[i]>lz): rz[i] -= lz
            if (len(repx)>1 and rx[i]<0):  rx[i] += lx
            if (len(repy)>1 and ry[i]<0):  ry[i] += ly
            if (len(repz)>1 and rz[i]<0):  rz[i] += lz
        for i in range(n):
            nei[i].clear()
            for j in range(n):
                for ix in repx:
                    for iy in repy:
                        for iz in repz:
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
    if (ifmd):
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
        for ii in range(3):
            for jj in range(3):
                satom[i][ii][jj] = 0
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
                # stress calc
                ad1 = vp(rr)/rr/vol_a
                vec1 = [drx,dry,drz]
                for ii in range(3):
                    for jj in range(3):
                        satom[i][ii][jj] \
                        += ad1*vec1[ii]*vec1[jj]
                
    # Verlet(2)
    if (ifmd):
        for i in range(n):
            vx[i] += dt/2.0 * fx[i] / wm
            vy[i] += dt/2.0 * fy[i] / wm
            vz[i] += dt/2.0 * fz[i] / wm
            ekin[i] = 0.5 * wm * (vx[i]**2 + vy[i]**2 + vz[i]**2)
    else:
        for i in range(n):
            ekin[i] = 0.0
    # stress calc
    ad1 = wm/vol_a
    vec1 = [vx[i],vy[i],vz[i]]
    for ii in range(3):
        for jj in range(3):
            satom[i][ii][jj] \
            -= ad1*vec1[ii]*vec1[jj]
    # total stress calc
    for ii in range(3):
        for jj in range(3):
            sig[ii][jj] = 0
            for i in range(n):
                sig[ii][jj] += satom[i][ii][jj]/n

    # Calc pot/kin energy
    epot_t = 0.0
    ekin_t = 0.0
    for i in range(n):
        epot_t += epot[i]
        ekin_t += ekin[i]
    temp = ekin_t *2.0 / (3.0 * n * kb)
    
    if (not ifmd):
        print('%f %f %f %e %f %f %f'\
        % (lx/ang,ly/ang,lz/ang,epot_t\
        ,sig[0][0]*mpa,sig[1][1]*mpa,sig[2][2]*mpa)\
        ,file=fl)
        continue
    # preheat control
    if (ifrlx_preheat):
        if (step <= rlx_preheat_step):
            ifgloc = False
            ifvscale = True
        else:
            ifgloc = True
            ifvscale = False
    # Structural relaxation by gloc
    if (ifgloc):
        vf = 0.0
        for i in range(n):
            vf += vx[i]*fx[i] + vy[i]*fy[i] + vz[i]*fz[i]
        if (vf < 0.0):
            for i in range(n):
                vx[i] = 0.0
                vy[i] = 0.0
                vz[i] = 0.0
    # Temperature control by velocity scaling
    if (ifvscale):
        vscale = math.sqrt(temp_set/temp);
        for i in range(n):
            vx[i] *= vscale; vy[i] *= vscale; vz[i] *= vscale;
    
    # relaxation check and stretch
    fmax = 0.0
    for i in range(n):
        fnor = fx[i]**2 + fy[i]**2 + fz[i]**2
        if (fnor > fmax): fmax = fnor
    fmax = math.sqrt(fmax)/ev*ang
    if (fmax < rlx_tol): # relaxed
        print('\neps = %f relaxed by %d steps.' % (eps,step))
        print('sigma_z = %f MPa\n' % (sig[2][2]*mpa))
        print('%f %f' % (eps,sig[2][2]*mpa), file=fd)
        fd.flush()
        write_cfg('eps'+'%05d'%(round(eps/0.001))\
        +'.cfg',n,rx,ry,rz,lx,ly,lz)
        step = 0
        eps += del_eps
        lz_next = lz0 * (1.0 + eps)
        for i in range(n):
            rz[i] *= lz_next/lz
        lz = lz_next
        vol = lx*ly*lz
        vol_a = vol/n
    elif (step > rlx_maxstep):
        print('eps = %f needs over %d steps' % (eps,step))
        print('fmax = %f eV/A' % (fmax))
        write_cfg('eps'+'%05d'%(round(eps/0.001))\
        +'_unrlx.cfg',n,rx,ry,rz,lx,ly,lz)
        exit(0)

    step += 1
    # output to file
    #print('%d %e %e %e' \
    #% (step,epot_t,ekin_t,epot_t+ekin_t), file=fo)
    #print('%d %e %e %e' \
    #% (step,sig[0][0]*mpa,sig[1][1]*mpa,sig[2][2]*mpa), file=fs)
    #cfg_interval = 10 # cfg file is written every XX steps
    #if ((step % cfg_interval)==0):
    #    write_cfg('out'+'%04d'%(step//cfg_interval)\
    #    +'.cfg',n,rx,ry,rz,lx,ly,lz)
    # output to console (for quasi-static stretch)
    print('eps:%f step:%d temp:%f fmax:%f' \
    % (eps,step,temp,fmax))
