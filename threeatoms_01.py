import math

ep = 0.2703
al = 1.1646
ro = 3.253
ev = 1.6021892e-19
ang = 1.0e-10

def v(r): # Morse function
    rang = r*ang
    return ep*(math.exp(-2.0*al*(rang-ro))\
    -2.0*math.exp(-al*(rang-ro))) *ev # J

def vp(r): # Derivatve of Morse function
    rang = r*ang
    return -2.0*al*ep*(math.exp(-2.0*al*(rang-ro))\
    -math.exp(-al*(rang-ro))) *ev/ang # J/m

# declaration of arrays
rx = [15.0*ang, 20.0*ang, 20.0*ang] # position [m]
ry = [10.0*ang, 10.0*ang, 15.0*ang]
fx = [0.0, 0.0, 0.0] # force [N]
fy = [0.0, 0.0, 0.0]
epot = [0.0, 0.0, 0.0] # potential energy
epot_tot = 0.0 # total potential energy
n = 3 # number of atoms

# Force and energy
for i in range(n):
    for j in range(n):
        if (i != j):
            rr = math.sqrt((rx[i]-rx[j])**2 + (ry[i]-ry[j])**2)
            drx = rx[i] - rx[j]
            dry = ry[i] - ry[j]
            fx[i] += -vp(rr)/rr*drx
            fy[i] += -vp(rr)/rr*dry
            epot[i] += v(rr)/2.0
    epot_tot += epot[i]

# Output
print('Epot_tot [J] = ',epot_tot)
for i in range(n):
    print('i = ',i,'Position [ang]:',rx[i]/ang,ry[i]/ang)
    print('       Force [N]:',fx[i],fy[i])
    print('       Epot [J]:',epot[i])
