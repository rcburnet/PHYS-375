import matplotlib.pyplot as plt
import numpy as np

#Part a)

fid = open('hipparcos.txt','r')

llines = fid.readlines()

fid.close()

lines_list = []

for i in range(len(llines)):
    lines_list.append([])
    lines_list[i] = llines[i].split()

parallax = []
V = []
B = []
I = []
D = []

for i in range(len(lines_list)):
    parallax.append(float(lines_list[i][0]))
    V.append(float(lines_list[i][1]))
    B.append(float(lines_list[i][2]))
    I.append(float(lines_list[i][3]))
    D.append(float(1/float(lines_list[i][0])*1000))

def absmag(x,y):
    return x-5*np.log10(y/10)

abs_Vmag = []

for i in range(len(D)):
    abs_Vmag.append([])
    abs_Vmag[i] = absmag(V[i],D[i])

B_V = []

for i in range(len(B)):
    B_V.append(B[i] - V[i])

def plota(x,y):
    plt.plot(x,y, 'ro')
    plt.gca().invert_yaxis()
    plt.grid(True)
    plt.title('2a) M$_V$ Magnitude vs. B-V Colour of %s stars'%(len(B_V)), fontsize=10)
    plt.xlabel('B-V Colour',fontsize=10)
    plt.ylabel('M$_V$ Magnitude', fontsize=10)
    plt.savefig('a1q2a.png', format='png')
    plt.close()

plota(B_V, abs_Vmag)


#Part b)

def T(x):
    return 9000/(0.93+x)

temp_list = []

for i in range(len(B_V)):
    temp_list.append(T(B_V[i]))

temp_list_log = map(np.log10, temp_list)

log_lumin = []

def loglum(v_mag):
    D_sun = 0.000004848
    #m_sun = -27
    m_sun = 4.83
    #return (v_mag -  m_sun + 5*np.log10(D_sun/10))/(-2.5)
    return (v_mag - m_sun)/(-2.5)

log_lumin = map(loglum, abs_Vmag)

def plotb(x,y):
    plt.plot(x,y, 'ro')
    plt.grid(True)
    plt.title('2b) log(L$_V$/L$_{\odot}$) vs log(T) of %s stars'%(len(B_V)))
    plt.xlabel('log(T)')
    plt.ylabel('log(L$_V$/L$_{\odot}$)')
    #plt.savefig('a1q2b.png', format='png')
    #plt.close()

#plotb(temp_list_log, log_lumin)

def stefan(T,R):
    sigma = 5.670367*10**(-8)
    factor = 6.995*10**8
    return np.log10(4*np.pi*(R*factor)**2*sigma)+4*T

def plotc():
    plotb(temp_list_log, log_lumin)
    radii = (1, 0.2, 5)
    solar_lum = 3.846*10**26
    x = np.arange(3162,11200,10)
    x = map(np.log10,x)
    y1 = []
    y2 = []
    y3 = []
    for i in range(len(x)):
        y1.append([])
        y2.append([])
        y3.append([])
        y1[i] = stefan(x[i], radii[0])-np.log10(solar_lum)
        y2[i] = stefan(x[i], radii[1])-np.log10(solar_lum)
        y3[i] = stefan(x[i], radii[2])-np.log10(solar_lum)
    #print x[0], radii[0], stefan(x[0],radii[0]), y1[0]
    #print y2[0]
    #y1 = map(np.log10, y1)
    #y2 = map(np.log10, y2)
    #y3 = map(np.log10, y3)
    #x = map(np.log10, x)
    plt.plot(x,y1, label='1 R$_\odot$')
    plt.plot(x,y2, label='0.2 R$_\odot$')
    plt.plot(x,y3, label='5 R$_\odot$')
    plt.gca().set_aspect('equal', adjustable='box')
    plt.yticks(fontsize=8)
    plt.xticks(fontsize=8)
    plt.axis('tight')
    plt.title('2b,c) log(L$_V$/L$_{\odot}$) vs log(T) of 200 stars and log(L/L$_\odot$) vs log(T) for stars of various Solar Radii (Stefan-Boltzmann Law)', fontsize=10)
    plt.xlabel('log(T)', fontsize=10)
    plt.ylabel('log(L/L$_\odot$), log(L$_V$/L$_{\odot}$)', fontsize=10)
    lgd = plt.legend(bbox_to_anchor=(1.02, 1.02), loc=2,title='log(L/L$_\odot$) vs log(T) of various Solar Radii')
    lgd.get_title().set_fontsize('8')
    plt.savefig('a1q2bc.png', format='png', bbox_extra_artists=(lgd,),bbox_inches='tight')
    plt.close()

plotc()


#Question 3

fid = open('W16_assignment1_orbit.dat','r')

llines = fid.readlines()

fid.close()

lines_list = []

for i in range(len(llines)):
    lines_list.append([])
    lines_list[i] = llines[i].split()

velocity1 = []
velocity2 = []
mag = []

for i in range(len(lines_list)):
    velocity1.append([])
    velocity2.append([])
    mag.append([])
    velocity1[i] = float(lines_list[i][1])
    velocity2[i] = float(lines_list[i][2])
    mag[i] = float(lines_list[i][3])
    
def plot3a(v1,v2):
    t = np.arange(0,50-0.04985,0.04985)
    plt.plot(t,v1, label='Star A')
    plt.plot(t,v2, label='Star B')
    plt.grid(True)
    #plt.gca().set_aspect('equal', adjustable='box')
    plt.title('3a) Radial velocities of two stars in a binary system over time', fontsize=10)
    plt.axis('tight')
    plt.xlabel('Time (days)', fontsize=10)
    plt.ylabel('Radial Velocities (km/s)', fontsize=10)
    plt.yticks(fontsize=8)
    plt.xticks(fontsize=8)
    lgd = plt.legend(bbox_to_anchor=(1.02, 1.02), loc=2,title='Radial Velocities of Star:')
    lgd.get_title().set_fontsize('8')
    plt.savefig('a1q3a.png', format='png', bbox_extra_artists=(lgd,),bbox_inches='tight')
    plt.close()

plot3a(velocity1,velocity2)

vbin1 = np.average(velocity1)

vbin2 = np.average(velocity2)

vbin = np.average([vbin1, vbin2])

def subtract(x):
    minmag = min(mag)
    return minmag - x

mag = map(subtract, mag)

def log_lum(x):
    return 0.4*x

lum = map(log_lum, mag)

def plot3c(x):
    t = np.arange(0,50-0.04985,0.04985)
    plt.plot(t,x)
    plt.grid(True)
    plt.title('3c) log(L/L$_o$) vs Time for eclipsing stars in binary system', fontsize=10)
    plt.axis('tight')
    plt.xlabel('Time (days)', fontsize=10)
    plt.ylabel('log(L/L$_o$)', fontsize=10)
    plt.yticks(fontsize=8)
    plt.xticks(fontsize=8)
    plt.savefig('a1q3c.png', format='png')
    plt.show()
    plt.close()

plot3c(lum)
