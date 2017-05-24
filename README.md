import numpy as np
import math
import matplotlib.pyplot as plt



def marsinit(filename):
    marstable = np.loadtxt(filename,skiprows=2)
    return marstable

def marsatm(h,marstable):
    for i in range(1,len(marstable)):
        if h==marstable[i][0]:
            return marstable[i][1:4]
        elif h>marstable[i-1][0] and h<marstable[i][0]:
            k =  (h - marstable[i-1][0])/(marstable[i][0]-marstable[i-1][0])
            a = (1-k) * marstable[i-1][1] + k *marstable[i][1]
            b = (1-k) * marstable[i-1][2] + k *marstable[i][2]
            c = (1-k) * marstable[i-1][3] + k *marstable[i][3]
            return a,b,c   
        
        
mt = marsinit("marsatm.txt")
speed = 262
h = 20.0



R = 191.84
vref = -2
throttle_factor = 0.05        
cs = 4.92
ve = 4400
dry_mass = 699
fluid_mass = 201 #assumption
m = dry_mass + fluid_mass
g0 = 3.711
cutoff = 0.3
alt_thrustle = 1500 #assumption
v_start = 262 #initital speed
H= np.array([0.,20000.]) # array of heights
angle = -20 #gamma angle       
V=np.array([np.cos(angle),np.sin(angle)])*v_start #array of velocity components
results = {}
results={'time':[],'mfuel':[],'me':[],'Hx':[],'Hy':[],'V':[],'Vx':[],'Vy':[],'angle':[]} 
dt = 0.01
time = 0    
vy = 0
vx = 0
m_limit= 5

while H[1]>0:
    density, temperature, speed = marsatm(h,mt)
    pressure  = density * temperature * R
    Fdrag = 1/2*cs*density*V
    x = Fdrag * math.cos(angle)
    y = Fdrag * math.sin(angle)
    Fgrav = m*g0
    F = Fdrag + Fgrav
    a=F/m 
    angle=math.degrees(-(math.atan2(V[0],V[1])-np.pi/2))
    if m>dry_mass and H[1]<=alt_thrustle and H[1]>cutoff:
        me = m*g0/ve +throttle_factor*(vref-V[1])
        me = math.min(me,m_limit)
    else:
        me = 0
    m -= me*dt                                      
    V += a*dt                                        
    H += V*dt
    h += V[1]*dt
    time = time + dt;                                        
    results['me'].append(me)
    results['time'].append(time)
    results['mfuel'].append(m-dry_mass)
    results['Hx'].append(H[0])
    results['Hy'].append(H[1])
    results['V'].append((V[0]**2+V[1]**2)**0.5)
    results['Vx'].append(V[0])
    results['Vy'].append(V[1])
    results['angle'].append(angle)

fig = plt.figure()
fig.subplots_adjust(bottom=0.2)
ax1 = fig.add_subplot(211)
ax2 = ax1.twinx()
line1 = ax1.plot(list(results['H']),'bo-',label='Alt')
line2 = ax2.plot(list(results['time']),'yo-',label='Time')
ax1.set_ylim(0,500)
ax2.set_ylim(0,500)
plt.show()
