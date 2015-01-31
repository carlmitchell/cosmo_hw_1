import numpy as np
import matplotlib.pyplot as plt
import scipy
from scipy import integrate
from scipy import optimize

H0 = 70/(3.08567758*10**19) #Hubble constant today in s^-1 (70 km/s/Mpc)
gyr = 3.15569*10**16 #Number of seconds in a Gyr

def adot(a,om_rad,om_mat,om_de,w=-1):
    """Calculates the time derivative of the scale factor when the scale
    factor is 'a'

    """
    return H0*np.sqrt( om_rad/a**2 + om_mat/a + om_de/a**(1+3*w) + (1 - om_rad - om_de - om_mat) )

def time(a,om_rad,om_mat,om_de,w=-1):
    """Calculates the time a universe would have when the scale factor is
    'a'

    """
    return scipy.integrate.quad(lambda x: 1/adot(x,om_rad,om_mat,om_de,w),1,a)[0]

def timediff(a,t,om_rad,om_mat,om_de,w=-1):
    """
    """
    return abs(t-time(a,om_rad,om_mat,om_de,w))

def scale(t,om_rad,om_mat,om_de,w=-1):
    """Calculates the scale factor a universe would have when the time is
    't' (in seconds).

    """
    print t
    return scipy.optimize.minimize(timediff,1,args=(t,om_rad,om_mat,om_de,w),
                                   method='powell').x

def hubble(t,om_rad,om_mat,om_de,w=-1):
    """Calculates the hubble 'constant' for a universe as a function of
    time (in seconds).

    """
    a = scale(t,om_rad,om_mat,om_de,w)
    return adot(a,om_rad,om_mat,om_de,w)/a

def density(t,om_rad,om_mat,om_de,w=-1):
    """Calculates the various components of the density (in units of
    critical density at a=1 as a function of time (in seconds).

    """
    a = scale(t,om_rad,om_mat,om_de,w)
    rho_rad = om_rad/a**4
    rho_mat = om_mat/a**3
    rho_de = om_de/a**(3+3*w)
    return rho_rad, rho_mat, rho_de, rho_rad+rho_mat+rho_de

def max_scale(om_rad,om_mat,om_de,w=-1):
    """Determines when/if a maximum scale factor will occur, and what that
    scale factor will be.

    """
    a_max = scipy.optimize.minimize(adot,1,args=(om_rad,om_mat,om_de,w),method='powell').x
    
    return time(a_max,om_rad,om_mat,om_de,w)

if __name__ == '__main__':
    
    om_rad = 0
    om_mat = 1
    om_de = 0
    w = -1
    
    #print scale(time(0,0,1,0),0,1,0)

    #max_time = max_scale(om_rad,om_mat,om_de,w)

    times = np.linspace(time(0,om_rad,om_mat,om_de,w),-time(0,om_rad,om_mat,om_de,w),10)
    #times = np.linspace(time(0,om_rad,om_mat,om_de,w),max_time,1000)
    
    As = np.array([scale(t,om_rad,om_mat,om_de,w) for t in times])
    plt.plot(times[1:]/gyr,As[1:])
    plt.scatter([0],[1])
    plt.show()
    
    #hubbles = np.array([hubble(t,om_rad,om_mat,om_de,w) for t in times])*70/H0
    #plt.semilogy(times[1:]/gyr,hubbles[1:])
    #plt.scatter([0],[70])
    #plt.show()
    
    #densitys = np.array([density(t,om_rad,om_mat,om_de,w) for t in times])
    #plt.semilogy(times[1:]/gyr,densitys[1:,3])
    #plt.scatter([0],[1])
    #plt.show()
    
