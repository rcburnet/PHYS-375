import numpy as np
import scipy as sp
from scipy import constants as con
import matplotlib.pyplot as plt

def rho(n):
    x = np.linspace(0,5,10000)
    theta = 1
    dtheta = 0
    ddtheta = -1
    dx = x[1]-x[0]
    theta_list = [1]
    M_list = [0]
    M = 0
    Ps_list = [0]
    x2 = np.linspace(0,5,10000)
    #MJ = 0
    MJ_list = [0]
    x3 = np.linspace(0,4000000000000000,10000)
    for i in range(len(x)-1):
        i += 1
        ddtheta = - (4*np.pi*theta**2 + ((2.0/x[i]) - (1.0/theta) * dtheta)*dtheta)
        dtheta += ddtheta * dx
        theta += dtheta * dx
        theta_list.append(theta)
        dM = 4*np.pi*x[i]**2*theta
        M += dM*dx
        M_list.append(M)
        Ps = theta*M**2*(con.k*10/(2.4*con.m_p))**4*(1.0/con.G)**3*(1.0/(1.989*10**30))**2
        x2[i] = x2[i]*(con.G*1.989*10**30*2.4*con.m_p)/(con.k*10*M*1.496*10**11)
        Ps_list.append(Ps)
        MJ = 0.2*1.989*10**30*((theta*(1.989*10**30/M)*(con.k*10/(con.G*1.989*10**15*2.4*con.m_p))**3)/(3.0*10**15))**(-1.0/2.0)
        MJ_list.append(MJ)
        x3[i] = x3[i]/(1.496*10**11)
    return (theta_list,M_list,Ps_list,x2,MJ_list,x3)

rho0=rho(0)

f0 = rho0[0]
x = np.linspace(0,5,10000)
l0, = plt.plot(x,f0,label='n=0')
plt.grid(True)
axes = plt.gca()
plt.title('2b) $\\tilde{\\rho}$ vs $\\tilde{r}$', fontsize=10)
plt.axis('tight')
plt.xlabel('$\\tilde{r}$', fontsize=10)
plt.ylabel('$\\tilde{\\rho}$', fontsize=10)
plt.yticks(fontsize=8)
plt.xticks(fontsize=8)
axes.set_ylim([0,1])
plt.savefig('a4q2b1.png', format='png')
#plt.show()
plt.close()

M0 = rho0[1]
l0, = plt.plot(x,M0,label='n=0')
plt.grid(True)
axes = plt.gca()
plt.title('2b) $\\tilde{M}$ vs $\\tilde{r}$', fontsize=10)
plt.axis('tight')
plt.xlabel('$\\tilde{r}$', fontsize=10)
plt.ylabel('$\\tilde{M}$', fontsize=10)
plt.yticks(fontsize=8)
plt.xticks(fontsize=8)
plt.savefig('a4q2b2.png', format='png')
#plt.show()
plt.close()

P0 = rho0[2]
#P0 = P0[::-1]
bob=rho0[3]
#bob[0]=8.2*10**7
#bob = bob[::-1]
l0, = plt.plot(bob,P0,label='n=0')
plt.grid(True)
axes = plt.gca()
axes.set_xlim([0,100000])
plt.title('2c) P$_s$ vs r', fontsize=10)
#plt.axis('tight')
plt.xlabel('r (AU)', fontsize=10)
plt.ylabel('P$_s$ (Pa)', fontsize=10)
plt.yticks(fontsize=8)
plt.xticks(fontsize=8)
plt.savefig('a4q2c.png', format='png',)
#plt.show()
plt.close()

P0 = rho0[4]
#P0 = P0[::-1]
bob=rho0[3]
#bob.pop()
#bob[0]=8.2*10**7
#bob = bob[::-1]
l0, = plt.plot(bob,P0,label='n=0')
plt.grid(True)
axes = plt.gca()
axes.set_xlim([0,100000])
plt.title('2d) M$_J$ vs r', fontsize=10)
#plt.axis('tight')
plt.xlabel('r (AU)', fontsize=10)
plt.ylabel('M$_J$ (kg)', fontsize=10)
plt.yticks(fontsize=8)
plt.xticks(fontsize=8)
plt.savefig('a4q2d.png', format='png',)
#plt.show()
plt.close()
