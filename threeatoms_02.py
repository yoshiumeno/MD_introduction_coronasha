import math

ep = 0.2703
al = 1.1646
ro = 3.253
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
rx = [15.0*ang, 20.0*ang, 20.0*ang] # position [m]
ry = [10.0*ang, 10.0*ang, 15.0*ang]
vx = [0.0, 0.0, 0.0] # velocity [m/s]
vy = [0.0, 0.0, 0.0]
fx = [0.0, 0.0, 0.0] # force [N]
fy = [0.0, 0.0, 0.0]
epot = [0.0, 0.0, 0.0] # potential energy
epot_tot = 0.0 # total potential energy
n = 3 # number of atoms
dt = 1.0e-15 # s
wm = 26.98 * 1.0e-3 / 6.022e23  #  Al mass (kg)

# open file for output (trajectory)
fout = open('trajectory.d','w')

step = 0
stepend = 1000

while step < stepend:
    # Verlet(1)
    for i in range(n):
        rx[i] += dt * vx[i] + (dt*dt/2.0) * fx[i] / wm
        ry[i] += dt * vy[i] + (dt*dt/2.0) * fy[i] / wm
        vx[i] += dt/2.0 * fx[i] / wm
        vy[i] += dt/2.0 * fy[i] / wm
    # Force and energy
    for i in range(n):
        fx[i] = 0
        fy[i] = 0
        epot[i] = 0
    for i in range(n):
        for j in range(n):
            if (i != j):
                rr = math.sqrt((rx[i]-rx[j])**2 \
                + (ry[i]-ry[j])**2)
                drx = rx[i] - rx[j]
                dry = ry[i] - ry[j]
                fx[i] += -vp(rr)/rr*drx
                fy[i] += -vp(rr)/rr*dry
                epot[i] += v(rr)/2.0
    # Verlet(2)
    for i in range(n):
        vx[i] += dt/2.0 * fx[i] / wm
        vy[i] += dt/2.0 * fy[i] / wm
    step += 1
    # output
    print(step,rx[0]/ang,ry[0]/ang,rx[1]/ang\
    ,ry[1]/ang,rx[2]/ang,ry[2]/ang\
    ,file=fout)
    print('step=%d r1:(%9.5f, %9.5f) \
r2:(%9.5f, %9.5f) r3:(%9.5f, %9.5f)'\
     % (step,rx[0]/ang,ry[0]/ang,rx[1]/ang\
     ,ry[1]/ang,rx[2]/ang,ry[2]/ang))
