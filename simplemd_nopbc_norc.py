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
    vx += [float(xy[3])]
    vy += [float(xy[4])]
    vz += [float(xy[5])]
    fx += [0]; fy += [0]; fz += [0]
    epot += [0]; ekin += [0]
    n += 1 # number of total atoms
print("number of atoms = ",n)

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
        ekin[i] = 0.5 * wm \
        * (vx[i]**2 + vy[i]**2 + vz[i]**2)

    # Calc pot/kin energy
    epot_t = 0.0
    ekin_t = 0.0
    for i in range(n):
        epot_t += epot[i]
        ekin_t += ekin[i]
    temp = ekin_t *2.0 / (3.0 * n * kb)

    step += 1
    # output energy to file
    print('%d %e %e %e' \
    % (step,epot_t,ekin_t,epot_t+ekin_t), file=fo)
    # output to console
    print('step:%d E_p:%e E_k:%e E:%e  temp:%f' \
    % (step,epot_t,ekin_t,epot_t+ekin_t,temp))
