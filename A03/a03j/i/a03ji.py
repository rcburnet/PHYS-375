import numpy as np
import scipy as sp
import matplotlib.pyplot as plt

def f(n):
    x = np.linspace(0,7,1000)
    theta = 1
    dtheta = 0
    ddtheta = -1
    dx = x[1]-x[0]
    theta_list = [1]
    for i in range(len(x)-1):
        i += 1
        ddtheta = - (theta**n + ((2.0/x[i]) * dtheta))
        dtheta += ddtheta * dx
        theta += dtheta * dx
        theta_list.append(theta)
    return theta_list

f0=f(0)
f1=f(1)
f2=f(2)
f3=f(3)
f4=f(4)
f5=f(5)
x = np.linspace(0,7,1000)

l0, = plt.plot(x,f0,label='n=0')
l1, = plt.plot(x,f1,label='n=1')
l2, = plt.plot(x,f2,label='n=2')
l3, = plt.plot(x,f3,label='n=3')
l4, = plt.plot(x,f4,label='n=4')
l5, = plt.plot(x,f5,label='n=5')
l0.set_dashes([4,4])
l1.set_dashes([8,4,2,4])
l2.set_dashes([8,4,2,4,2,4])
l4.set_dashes([8,4,2,4,2,4,2,4])
l5.set_dashes([8,4,2,4,2,4,2,4,2,4])
plt.grid(True)
axes = plt.gca()
plt.title('3e) $\\theta$ vs x', fontsize=10)
plt.axis('tight')
plt.xlabel('x', fontsize=10)
plt.ylabel('$\\theta$', fontsize=10)
plt.yticks(fontsize=8)
plt.xticks(fontsize=8)
axes.set_ylim([0,1])
lgd = plt.legend(bbox_to_anchor=(1.02, 1.02), loc=2, numpoints = 1)
plt.savefig('a3q3e.png', format='png', bbox_extra_artists=(lgd,),bbox_inches='tight')
#plt.show()
plt.close()

rho_c = 1610

for i in range(len(f3)):
    f3[i] = f3[i]**3*rho_c
    #if f3[i] < 0:
        #f3[i]=0
        
rho = np.array(f3)
for i in range(len(rho)):
    if f3[i] < 0:
        f3[i] = 0
x = x/7

plt.plot(x,f3)
plt.grid(True)
axes = plt.gca()
plt.title('3h) $\\rho$ vs $\\frac{r}{R_{star}}$ for n=3 star with mass and radius of the sun', fontsize=10)
plt.axis('tight')
plt.xlabel('$\\frac{r}{R_{star}}$')
plt.ylabel('$\\rho (kg/m)$')
plt.yticks(fontsize=8)
plt.xticks(fontsize=8)
#axes.set_ylim([0,1])
plt.savefig('a3q3h_rho.png', format='png')
#plt.show()
plt.close()

K = 2.84764614e10

for i in range(len(rho)):
    f3[i] = K*(f3[i]**(4.0/3.0))

P = np.array(f3)/1e16

plt.plot(x,P)
plt.grid(True)
axes = plt.gca()
plt.title('3h) $P$ vs $\\frac{r}{R_{star}}$ for n=3 star with mass and radius of the sun', fontsize=10)
plt.axis('tight')
plt.xlabel('$\\frac{r}{R_{star}}$')
plt.ylabel('P (N/m$^2$) $\\times$ 1e16')
plt.yticks(fontsize=8)
plt.xticks(fontsize=8)
#axes.set_ylim([0,1])
plt.savefig('a3q3h_P.png', format='png')
#plt.show()
plt.close()

mu = 0.70175
mp = 1.672622e-27
kb = 1.38064852e-23

for i in range(len(f3)):
    try:
        f3[i] = f3[i]*mu*mp/(rho[i]*kb)
    except:
        f3[i]=0

T = np.array(f3)/1e7

plt.plot(x,T)
plt.grid(True)
axes = plt.gca()
plt.title('3h) $T$ vs $\\frac{r}{R_{star}}$ for n=3 star with mass and radius of the sun', fontsize=10)
plt.axis('tight')
plt.xlabel('$\\frac{r}{R_{star}}$')
plt.ylabel('T (K) $\\times$ 1e7')
plt.yticks(fontsize=8)
plt.xticks(fontsize=8)
#axes.set_ylim([0,1])
plt.savefig('a3q3h_T.png', format='png')
#plt.show()
plt.close()

epp = np.zeros([len(f3)])
ecno = np.zeros([len(f3)])
dL = np.zeros([len(f3)])


for i in range(len(f3)):
    epp[i] = rho[i]*0.55*0.55*(T[i]*10)**4*1e-7
    ecno[i] = rho[i]*0.55*0.55*0.03*(T[i]*10)**19.9*1e-7
    dL[i] = 4*np.pi*(x[i]*10*696000000)**2*rho[i]*((1.07e-5)*epp[i]+(8.24e-24)*ecno[i])#/(2.68e123)

#x = x*5.0/3.0
dL = np.array(dL)

plt.plot(x,epp)
plt.grid(True)
axes = plt.gca()
plt.title('3i) $\\frac{\\epsilon_{pp}}{C_{pp}}$ vs $\\frac{r}{R_{star}}$ for n=3 star with mass and radius of the sun', fontsize=10)
plt.axis('tight')
plt.xlabel('$\\frac{r}{R_{star}}$')
plt.ylabel('$\\frac{\\epsilon_{pp}}{C_{pp}}$')
plt.yticks(fontsize=8)
plt.xticks(fontsize=8)
axes.set_xlim([0,1])
plt.savefig('a3q3iepp.png', format='png')
#plt.show()
plt.close()

plt.plot(x,ecno)
plt.grid(True)
axes = plt.gca()
plt.title('3i) $\\frac{\\epsilon_{CNO}}{C_{CNO}}$ vs $\\frac{r}{R_{star}}$ for n=3 star with mass and radius of the sun', fontsize=10)
plt.axis('tight')
plt.xlabel('$\\frac{r}{R_{star}}$')
plt.ylabel('$\\frac{\\epsilon_{CNO}}{C_{CNO}}$')
plt.yticks(fontsize=8)
plt.xticks(fontsize=8)
axes.set_xlim([0,1])
plt.savefig('a3q3ieCNO.png', format='png')
#plt.show()
plt.close()

plt.plot(x,dL/1e18)
plt.grid(True)
axes = plt.gca()
plt.title('3i) $\\frac{dL}{dr}$ vs $\\frac{r}{R_{star}}$ for n=3 star with mass and radius of the sun', fontsize=10)
plt.axis('tight')
plt.xlabel('$\\frac{r}{R_{star}}$')
plt.ylabel('$\\frac{dL}{dr}$ $\\times$ 1e18')
plt.yticks(fontsize=8)
plt.xticks(fontsize=8)
axes.set_xlim([0,1])
plt.savefig('a3q3i.png', format='png')
#plt.show()
plt.close()

L = 0

for i in range(len(dL)):
    L += dL[i] * (x[1]-x[0]) * 696000000
