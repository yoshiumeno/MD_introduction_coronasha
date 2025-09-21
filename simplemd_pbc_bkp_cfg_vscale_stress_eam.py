import math
import numpy as np

ep = 0.2703; al = 1.1646; ro = 3.253
ev = 1.6021892e-19
ang = 1.0e-10
# for stress
mpa = 1.0e-6

# eam potential
eam_pair  = []; eam_paird = []
eam_den  = []; eam_dend = []
eam_embed  = []; eam_embedd = []
f_eam_pair  = open('EAM_Al_pair','r')
f_eam_den   = open('EAM_Al_den','r')
f_eam_embed = open('EAM_Al_embed','r')
# read pair
line = f_eam_pair.readline()
xy = line.split()
eam_pairmin = float(xy[0])
eam_pairmax = float(xy[1])
for line in f_eam_pair:
    xy = line.split()
    eam_pair  += [float(xy[0])]
    eam_paird += [float(xy[1])]
eam_pairmesh = len(eam_pair)
dpair =(eam_pairmax-eam_pairmin)/float(eam_pairmesh-1)
# read den
line = f_eam_den.readline()
xy = line.split()
eam_denmin = float(xy[0])
eam_denmax = float(xy[1])
for line in f_eam_den:
    xy = line.split()
    eam_den  += [float(xy[0])]
    eam_dend += [float(xy[1])]
eam_denmesh = len(eam_den)
dden =(eam_denmax-eam_denmin)/float(eam_denmesh-1)
# read embed
line = f_eam_embed.readline()
xy = line.split()
eam_embedmin = float(xy[0])
eam_embedmax = float(xy[1])
for line in f_eam_embed:
    xy = line.split()
    eam_embed  += [float(xy[0])]
    eam_embedd += [float(xy[1])]
eam_embedmesh = len(eam_embed)
dembed =(eam_embedmax-eam_embedmin)/float(eam_embedmesh-1)
# prepare coefficients for cubic interpolation
eam_apair = np.zeros(eam_pairmesh)
eam_bpair = np.zeros(eam_pairmesh)
eam_cpair = np.zeros(eam_pairmesh)
eam_dpair = np.zeros(eam_pairmesh)
eam_aden = np.zeros(eam_denmesh)
eam_bden = np.zeros(eam_denmesh)
eam_cden = np.zeros(eam_denmesh)
eam_dden = np.zeros(eam_denmesh)
eam_aembed = np.zeros(eam_embedmesh)
eam_bembed = np.zeros(eam_embedmesh)
eam_cembed = np.zeros(eam_embedmesh)
eam_dembed = np.zeros(eam_embedmesh)
for i in range(eam_pairmesh-1):
    eam_apair[i] = eam_pair[i] - eam_pair[i + 1] \
    + eam_paird[i + 1]*dpair
    eam_bpair[i] = -2 * (eam_pair[i] - eam_pair[i + 1]) \
    - (eam_paird[i] + eam_paird[i+1])*dpair
    eam_cpair[i] = eam_paird[i]*dpair
    eam_dpair[i] = eam_pair[i]
for i in range(eam_denmesh-1):
    eam_aden[i] = eam_den[i] - eam_den[i + 1] \
    + eam_dend[i + 1]*dden
    eam_bden[i] = -2 * (eam_den[i] - eam_den[i + 1]) \
    - (eam_dend[i] + eam_dend[i+1])*dden
    eam_cden[i] = eam_dend[i]*dden
    eam_dden[i] = eam_den[i]
for i in range(eam_embedmesh-1):
    eam_aembed[i] = eam_embed[i] - eam_embed[i + 1] \
    + eam_embedd[i + 1]*dembed
    eam_bembed[i] = -2 * (eam_embed[i] - eam_embed[i + 1]) \
    - (eam_embedd[i] + eam_embedd[i+1])*dembed
    eam_cembed[i] = eam_embedd[i]*dembed
    eam_dembed[i] = eam_embed[i]
# eam functions
def eam_func_pair(x):
    if (x > eam_pairmax*ang): x = eam_pairmax*ang -1e-14
    if (x < eam_pairmin*ang): x = eam_pairmin*ang
    xi = float(eam_pairmesh - 1)*(x/ang - eam_pairmin)\
    /(eam_pairmax - eam_pairmin) + 1.0
    ii = int(xi)
    xx = xi-float(ii)
    return ev*(eam_apair[ii] * xx * xx * xx + eam_bpair[ii]\
     * xx * xx + eam_cpair[ii] * xx + eam_dpair[ii])
def eam_func_pair_d(x):
    if (x > eam_pairmax*ang): x = eam_pairmax*ang -1e-14
    if (x < eam_pairmin*ang): x = eam_pairmin*ang
    xi = float(eam_pairmesh - 1)*(x/ang - eam_pairmin)\
    /(eam_pairmax - eam_pairmin) + 1.0
    ii = int(xi)
    xx = xi-float(ii)	
    return (eam_pairmesh - 1) * (3*eam_apair[ii] * xx * xx\
     + 2*eam_bpair[ii] * xx + eam_cpair[ii])\
     /(eam_pairmax - eam_pairmin);
def eam_func_den(x):
    if (x > eam_denmax*ang): x = eam_denmax*ang -1e-14
    if (x < eam_denmin*ang): x = eam_denmin*ang
    xi = float(eam_denmesh - 1)*(x/ang - eam_denmin)\
    /(eam_denmax - eam_denmin) + 1.0
    ii = int(xi)
    xx = xi-float(ii)
    return eam_aden[ii] * xx * xx * xx + eam_bden[ii]\
     * xx * xx + eam_cden[ii] * xx + eam_dden[ii]
def eam_func_den_d(x):
    if (x > eam_denmax*ang): x = eam_denmax*ang -1e-14
    if (x < eam_denmin*ang): x = eam_denmin*ang
    xi = float(eam_denmesh - 1)*(x/ang - eam_denmin)\
    /(eam_denmax - eam_denmin) + 1.0
    ii = int(xi)
    xx = xi-float(ii)
    return (eam_denmesh - 1) * (3*eam_aden[ii] * xx * xx\
     + 2 * eam_bden[ii] * xx + eam_cden[ii]) \
     / (eam_denmax - eam_denmin)
def eam_func_embed(x):
    if (x > eam_embedmax): x = eam_embedmax -1e-14
    if (x < eam_embedmin): x = eam_embedmin
    xi= float(eam_embedmesh - 1)*(x - eam_embedmin)\
    /(eam_embedmax - eam_embedmin) + 1.0
    ii = int(xi)
    xx = xi-float(ii)
    return ev*(eam_aembed[ii] * xx * xx * xx + eam_bembed[ii]\
     * xx * xx + eam_cembed[ii] * xx + eam_dembed[ii])
def eam_func_embed_d(x):
    if (x > eam_embedmax): x = eam_embedmax -1e-14
    if (x < eam_embedmin): x = eam_embedmin
    xi = float(eam_embedmesh - 1)*(x - eam_embedmin)\
    /(eam_embedmax - eam_embedmin)+1.0
    ii = int(xi)
    xx = xi-float(ii)
    return  (eam_embedmesh - 1) * (3*eam_aembed[ii] * xx * xx\
     + 2*eam_bembed[ii] * xx + eam_cembed[ii])\
     /(eam_embedmax - eam_embedmin)

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
rhob = []
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
fs= open('stress.d','w')
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
    rhob += [0]
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
    # Charge density
    for i in range(n):
        rhob[i] = 0
    for i in range(n):
        for jj in range(len(nei[i])):
            j = nei[i][jj][0];
            if (j>=i):
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
                    rho_ij = eam_func_den(rr)
                    rhob[i] += rho_ij
                    if (j!=i):
                        rhob[j] += rho_ij
    # Force and energy
    for i in range(n):
        fx[i] = 0; fy[i] = 0; fz[i] = 0
        epot[i] = eam_func_embed(rhob[i]) # embedded energy
        for ii in range(3):
            for jj in range(3):
                satom[i][ii][jj] = 0
    for i in range(n):
        for jj in range(len(nei[i])):            
            j  = nei[i][jj][0]
            if (j>=i):
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
                    vij  = eam_func_pair(rr)
                    dv   = eam_func_pair_d(rr)
                    df_i = eam_func_embed_d(rhob[i])
                    df_j = eam_func_embed_d(rhob[j])
                    dr   = eam_func_den_d(rr)
                    f = -(df_i*dr+df_j*dr+dv)*ev/ang
                    ki = f/rr
                    if (j==i): ki /= 2.0
                    kif = ki/2.0

                    fx[i] += ki*drx
                    fy[i] += ki*dry
                    fz[i] += ki*drz
                    fx[j] -= ki*drx
                    fy[j] -= ki*dry
                    fz[j] -= ki*drz
                    epot[i] += vij/2.0
                    if (j!=i): epot[j] += vij/2.0
                    # stress calc
                    vec1 = [drx,dry,drz]
                    for ii in range(3):
                        for jj in range(3):
                            satom[i][ii][jj] \
                            -= kif*vec1[ii]*vec1[jj]
                            satom[j][ii][jj] \
                            -= kif*vec1[ii]*vec1[jj]

    # Verlet(2)
    for i in range(n):
        vx[i] += dt/2.0 * fx[i] / wm
        vy[i] += dt/2.0 * fy[i] / wm
        vz[i] += dt/2.0 * fz[i] / wm
        ekin[i] = 0.5 * wm * (vx[i]**2 + vy[i]**2 + vz[i]**2)
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
    
    # Temperature control by velocity scaling
    if (ifvscale):
        vscale = math.sqrt(temp_set/temp);
        for i in range(n):
            vx[i] *= vscale; vy[i] *= vscale; vz[i] *= vscale;

    step += 1
    # output to file
    print('%d %e %e %e' % (step,epot_t,ekin_t,epot_t+ekin_t), file=fo)
    print('%d %e %e %e' % (step,sig[0][0]*mpa,sig[1][1]*mpa,sig[2][2]*mpa), file=fs)
    cfg_interval = 10 # cfg file is written every XX steps
    if ((step % cfg_interval)==0):
        write_cfg('out'+'%04d'%(step//cfg_interval)\
        +'.cfg',n,rx,ry,rz,lx,ly,lz)
    # output to console
    print('step:%d E_p:%e E_k:%e E:%e  temp:%f' \
    % (step,epot_t,ekin_t,epot_t+ekin_t,temp))
