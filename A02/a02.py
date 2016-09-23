import matplotlib.pyplot as plt
import numpy as np
import scipy as sp
from scipy import constants as con

def b_lambda(l):
    t=5500
    return (2*con.h*con.c**2/l**5)*(np.e**(con.h*con.c/(l*con.k*t))-1)**(-1)

def b_f(nu):
    t=5500
    return (2*con.h*nu**3/con.c**2)*(np.e**(con.h*nu/(con.k*t))-1)**(-1)

def lam_blam1(l):
    t=3000
    return (2*con.h*con.c**2/l**4)*(np.e**(con.h*con.c/(l*con.k*t))-1)**(-1)

def lam_blam2(l):
    t=5500
    return (2*con.h*con.c**2/l**4)*(np.e**(con.h*con.c/(l*con.k*t))-1)**(-1)

def lam_blam3(l):
    t=30000
    return (2*con.h*con.c**2/l**4)*(np.e**(con.h*con.c/(l*con.k*t))-1)**(-1)

def T(s):
    Te = 10000
    p = 10**(-6)
    K = 3.0
    return ((3.0/4)*Te**4*(p*K*s+2.0/3))**(1.0/4)

def F2(t):
    mp = 1.67e-27
    me = 9.11e-31
    k = 1.38e-23
    hbar = 6.626e-34/(2.*np.pi)
    C = (mp/(1*10**(-6)))*(me * k * t / (2.*np.pi*hbar**2))**(3./2) * np.exp(-13.6*1.602e-19/(k*t))
    NI2_NI1 = 4.0*np.e**(-13.6*(3.0/4)/(8.6173324*10**(-5)*t))
    NI3_NI1 = 9.0*np.e**(-13.6*(8.0/9)/(8.6173324*10**(-5)*t))
    NI4_NI1 = 16.0*np.e**(-13.6*(15.0/16)/(8.6173324*10**(-5)*t))
    NI5_NI1 = 25.0*np.e**(-13.6*(24.0/25)/(8.6173324*10**(-5)*t))
    NI6_NI1 = 36.0*np.e**(-13.6*(35.0/36)/(8.6173324*10**(-5)*t))
    NII_Nt = (-1.0*C/2.0)+(C/2.0)*(np.sqrt(1.0+4.0/C))
    NI_Nt = 1.0-NII_Nt
    N2I_NI = NI2_NI1/(1+NI2_NI1+NI3_NI1++NI4_NI1+NI5_NI1+NI6_NI1)
    return N2I_NI*NI_Nt

def tau(s):
    p = 1*10**(-6)
    k = 3.0
    return p*k*s

def tau_balm(f,s):
    p = 1*10**(-6)
    k = 3.5*10**5
    return f*p*k*s

s = np.arange(0,1*10**6,1000)
T1 = map(T, s)
f2 = map(F2, T1)

l = np.arange(1*10**(-7),2*10**(-6),1*10**(-8))
b_lam = map(b_lambda, l)
b_lam1 = map(lam_blam1, l)
b_lam2 = map(lam_blam2, l)
b_lam3 = map(lam_blam3, l)

nu = np.arange(0.15*10**14,30*10**14,0.05*10**14)
b_nu = map(b_f, nu)

l1 = np.arange(4*10**(-7),8*10**(-7),0.05*10**(-7))
b_lam4 = map(lam_blam1, l1)
b_lam5 = map(lam_blam2, l1)
b_lam6 = map(lam_blam3, l1)

l1= map(np.log10, l1)
b_lam4 = map(np.log10, b_lam4)
b_lam5 = map(np.log10, b_lam5)
b_lam6 = map(np.log10, b_lam6)

def plot2a(t,x):
    plt.plot(t,x)
    plt.grid(True)
    plt.plot(t[x.index(max(x))],max(x),'xk')
    plt.text(t[x.index(max(x))]+0.000000015,max(x)+0.01*10**(13),'(%5.2e,%5.2e)'%(t[x.index(max(x))],max(x)),fontsize=10)
    plt.title('2a) Plank\'s Function vs Wavelength for wavelengths of 100nm to 2$\mu$m', fontsize=10)
    #plt.axis('tight')
    plt.xlabel('$\lambda$ (m)', fontsize=10)
    plt.ylabel('B$_\lambda$ (W/sr/m$^3)$', fontsize=10)
    plt.yticks(fontsize=8)
    plt.xticks(fontsize=8)
    plt.savefig('a2q2a.png', format='png')
    #plt.show()
    plt.close()

def plot2b(t,x):
    plt.plot(t,x)
    plt.grid(True)
    plt.plot(t[x.index(max(x))],max(x),'xk')
    plt.text(t[x.index(max(x))]+0.03*10**15,max(x),'(%5.2e,%5.2e)'%(t[x.index(max(x))],max(x)),fontsize=10)
    plt.title('2b) Plank\'s Function vs Frequency for frequencies of (0.15-30)$\\times$10$^{14}$Hz', fontsize=10)
    #plt.axis('tight')
    plt.xlabel('$\\nu$ (Hz)', fontsize=10)
    plt.ylabel('B$_\\nu$ (W/sr/m$^2$/Hz)', fontsize=10)
    plt.yticks(fontsize=8)
    plt.xticks(fontsize=8)
    plt.savefig('a2q2b.png', format='png')
    #plt.show()
    plt.close()
    
def plot2d(t,x,y,z):
    line1, = plt.plot(t,x, 'r', label='T=3000K')
    line1.set_dashes([8, 4])
    plt.plot(t,y, 'g', label='T=5500K')
    line2, = plt.plot(t,z, 'b', label='T=30000K')
    line2.set_dashes([8, 4, 2, 4, 2, 4])
    plt.grid(True)
    plt.title('2d.1) log($\lambda$B$_\lambda$) vs log($\lambda$) for wavelengths of 100nm to 2$\mu$m', fontsize=10)
    plt.axis('tight')
    plt.xlabel('log($\lambda$)', fontsize=10)
    plt.ylabel('log($\lambda$B$_\lambda$)', fontsize=10)
    plt.yticks(fontsize=8)
    plt.xticks(fontsize=8)
    lgd = plt.legend(bbox_to_anchor=(1.02, 1.02), loc=2,title='Temperature')
    plt.savefig('a2q2d1.png', format='png', bbox_extra_artists=(lgd,),bbox_inches='tight')
    #plt.show()
    plt.close()

def plot2d2(t,x,y,z):
    line1, = plt.plot(t,x, 'r', label='T=3000K')
    line1.set_dashes([8, 4])
    plt.plot(t,y, 'g', label='T=5500K')
    line2, = plt.plot(t,z, 'b', label='T=30000K')
    line2.set_dashes([8, 4, 2, 4, 2, 4])
    plt.grid(True)
    plt.plot(t[y.index(max(y))],max(y),'xk')
    plt.text(t[y.index(max(y))]+0.001,max(y)+0.05,'(%5.2f,%5.2f)'%(t[y.index(max(y))],max(y)),fontsize=10)
    plt.gca().set_aspect('equal', adjustable='box')
    plt.title('2d.2) log($\lambda$B$_\lambda$) vs log($\lambda$) for wavelengths of 400nm to 800nm', fontsize=10)
    plt.axis('tight')
    plt.xlabel('log($\lambda$)', fontsize=10)
    plt.ylabel('log($\lambda$B$_\lambda$)', fontsize=10)
    plt.yticks(fontsize=8)
    plt.xticks(fontsize=8)
    lgd = plt.legend(bbox_to_anchor=(1.02, 1.02), loc=2,title='Temperature')
    plt.savefig('a2q2d2.png', format='png', bbox_extra_artists=(lgd,),bbox_inches='tight')
    #plt.show()
    plt.close()

def plot4a(t,x):
    plt.plot(t,x)
    plt.grid(True)
    plt.title('4a) Temperature of typical A0V star vs its Depth', fontsize=10)
    plt.axis('tight')
    plt.xlabel('Depth (m)', fontsize=10)
    plt.ylabel('Temperature (K)', fontsize=10)
    plt.yticks(fontsize=8)
    plt.xticks(fontsize=8)
    plt.savefig('a2q4a.png', format='png')
    #plt.show()
    plt.close()

def plot4b(t,x):
    plt.plot(t,x)
    plt.grid(True)
    plt.plot(t[x.index(max(x))],max(x),'xk')
    plt.text(t[x.index(max(x))]+10000,max(x)+0.0000001,'(%5.2e,%5.2e)'%(t[x.index(max(x))],max(x)),fontsize=10)
    plt.title('4b) Fraction of hydrogen atoms in first excited state of typical A0V star vs its Depth', fontsize=10)
    #plt.axis('tight')
    plt.xlabel('Depth (m)', fontsize=10)
    plt.ylabel('Fraction of hydrogen atoms in first excited state', fontsize=10)
    plt.yticks(fontsize=8)
    plt.xticks(fontsize=8)
    plt.savefig('a2q4b.png', format='png')
    #plt.show()
    plt.close()

def plot4c(t,x,y):
    plt.plot(t,x,label='$\\tau$')
    line, = plt.plot(t,y,label='$\\tau_{balmer}$')
    line.set_dashes([8, 4])
    plt.grid(True)
    axes = plt.gca()
    axes.set_ylim([0,1])
    s1 = np.linspace(0,350000,10000)
    s2 = []
    for i in range(len(s1)):
        s2.append(2.0/3.0)
    line2, = plt.plot(s1,s2)
    line2.set_dashes([2,2])
    plt.plot(113113.0,2.0/3.0,'xk')
    plt.text(113113.0+5000,2.0/3.0+0.01,'(%5.2e,$\\frac{2}{3}$)'%(113113.0),fontsize=10)
    plt.plot(222222,2.0/3.0,'xk')
    plt.text(222222+5000,2.0/3.0-0.03,'(%5.2e,$\\frac{2}{3}$)'%(222222),fontsize=10)
    plt.title('4c) Optical Depth of star considering the Balmer series vs its Depth', fontsize=10)
    #plt.axis('tight')
    lgd = plt.legend(bbox_to_anchor=(1.02, 1.02), loc=2,title='Optical Depths')
    plt.xlabel('Depth (m)', fontsize=10)
    plt.ylabel('$\\tau$', fontsize=10)
    plt.yticks(fontsize=8)
    plt.xticks(fontsize=8)
    plt.savefig('a2q4c.png', format='png', bbox_extra_artists=(lgd,),bbox_inches='tight')
    #plt.show()
    plt.close()

plot2a(l,b_lam)
plot2b(nu,b_nu)

l = map(np.log10, l)
b_lam1 = map(np.log10, b_lam1)
b_lam2 = map(np.log10, b_lam2)
b_lam3 = map(np.log10, b_lam3)

plot2d(l,b_lam1,b_lam2,b_lam3)
plot2d2(l1,b_lam4,b_lam5,b_lam6)
plot4a(s, T1)
plot4b(s, f2)

s = np.linspace(0,333333,1000)
T1 = map(T, s)
f2 = map(F2, T1)

tau_list = map(tau,s)
tau_balm_list=[]
for i in range(len(s)):
    tau_balm_list.append(tau_balm(f2[i],s[i])+tau_list[i])

plot4c(s, tau_list, tau_balm_list)
s = np.linspace(0,1000000000,100000)
T1 = map(T, s)
f2 = map(F2, T1)
area = sum(f2)
#area = np.trapz(f2,dx=10000)
