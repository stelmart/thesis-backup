from scipy import *
from scipy.integrate import quad
from numpy import *
from scipy.special import jv
import matplotlib.pyplot as plt

def A1(x, H, h, z):
    sxH = sinh(x*H)
    sxh = sinh(x*h)
    sxHh = sinh(x*(H-h))
    cothxH = cosh(x*H)/sxH

    ddxHh = ((H-h)*cosh(x*(H-h)) - H*sxHh*cothxH)/sxH
    ddxh =  (h*cosh(x*h) - H*sxh*cothxH)/sxH

    Aa = x*h*H*z*sinh(x*(H-z))*sxHh
    Ab = z*(h*sxH*cosh(x*(H-h-z)) - H*sxh*cosh(x*z))
    Ac = x*H*z*sxH*cosh(x*(H-z))*(ddxHh)
    Ad = -H*sinh(x*z)*(sxH*ddxh + x*H*ddxHh)

    A1fin = (sxH**2 - (x*H)**2)**(-1)*(Aa + Ab + Ac + Ad)
    return A1fin

def integral1(rho, z, H, h):
    func1 = lambda x: x*jv(1,rho*x)*A1(x,H,h,z)
    return quad(func1,0,inf)[0]

def integral2(rho, z, H, h):
    func2 = lambda x: x*x*(jv(0,rho*x)-jv(2,rho*x))*A1(x,H,h,z)/2
    return quad(func2,0,inf)[0]


def W(alpha, beta, x, y, z, H, h):
    xsq = x*x
    ysq = y*y
    rhosq = xsq + ysq
    rho = sqrt(rhosq)

    if alpha==1 and beta==1:
        return -(1-(xsq/rhosq))/rho*integral1(rho, z, H, h) - xsq/rhosq*integral2(rho, z, H, h) 

    elif (alpha==1 and beta==2) or (alpha==2 and beta==1):
        return x*y/rhosq*(integral1(rho, z, H, h) - integral2(rho, z, H, h))  

    elif alpha==2 and beta==2:
        return -(1-(ysq/rhosq))/rho*integral1(rho, z, H, h) - ysq/rhosq*integral2(rho, z, H, h)  

    else:
        print("alpha or beta out of range")
        return 0


H = 1.0
h = 0.5
z = 0.0

#leaving off the 1/4 pi mu

xlist = arange(0.2,16*H,0.2)
divs = len(xlist)
Wlist = zeros(divs)

for i in range(divs):
    Wlist[i] = W(2,2,xlist[i],xlist[i],z,H,h)

plt.title('w11 for H=1, h=0.5, x3 = 0.0')
plt.xlabel('x1')
plt.ylabel('w')
plt.plot(xlist, Wlist)
plt.show()
###############################

#plot correlation functions
#plt.title('x(t)')
#plt.xlabel('t')
#plt.ylabel('x')
#plt.plot(t, psoln[:,0])
#plt.plot(t, psoln2[:,0])
#plt.show()
