import numpy as np
import matplotlib.pyplot as plt
from scipy import integrate
from scipy import optimize

H0 = 70/(3.08567758*10**19) #Hubble constant today in s^-1 (70 km/s/Mpc)
gyr = 3.15569*10**16 #Number of seconds in a Gyr
suffix = '.png'
#rad_den = 0.0001
rad_den = 0

plt.rc('font', family='serif') #Changes all plotting fonts.

def adot(a,om_rad,om_mat,om_de,w):
    """Calculates the time derivative of the scale factor when the scale
    factor is 'a'

    """
    return H0*np.sqrt( om_rad/a**2 + om_mat/a + om_de/a**(1+3*w) + (1-(om_rad+om_de+om_mat)) )

def time(a,om_rad,om_mat,om_de,w):
    """Calculates the time a universe would have when the scale factor is
    'a'

    """
    return integrate.quad(lambda x: 1/adot(x,om_rad,om_mat,om_de,w),1,a)[0]

def scale(t,om_rad,om_mat,om_de,w): #Currently obsolete
    """Calculates the scale factor a universe would have when the time is
    't' (in seconds).

    """
    return optimize.minimize(lambda x: (time(x,om_rad,om_mat,om_de,w)-t)**2, 1, method='powell').x

def hubble(a,om_rad,om_mat,om_de,w):
    """Calculates the hubble 'constant' for a universe as a function of
    'a'.

    """
    return adot(a,om_rad,om_mat,om_de,w)/a

def density(a,om_rad,om_mat,om_de,w):
    """Calculates the various components of the density (in units of
    critical density at a=1 as a function of 'a'.

    """
    rho_rad = om_rad/a**4
    rho_mat = om_mat/a**3
    rho_de = om_de/a**(3+3*w)
    return np.array([rho_rad, rho_mat, rho_de, rho_rad+rho_mat+rho_de])

def t_start(om_rad,om_mat,om_de,w): #Currently obsolete
    """Returns the time (or *A* time) that the scale factor is equal to zero.
    
    """
    return time(0,om_rad,om_mat,om_de,w)

if __name__ == '__main__':
    
    om_rad = np.array([0, 0, 0, 0, 0, 0])
    om_mat = np.array([1, 2, 0.3, 0.3, 0.3, 0.3])
    om_de = np.array([0, 0, 0, 0.7, 0.7, 0.7])
    w = np.array([-1, -1, -1, -1, -2/3., -4/3.])
    uni = ['EdS','Closed','Open','LCDM','Quint','Phant']
    npoints = 10000

    if rad_den != 0:
        om_rad = om_rad+rad_den
        om_mat = om_mat-rad_den
        for i in range(len(uni)):
            uni[i] = uni[i]+'+R'
    
    for i in range(len(uni)):
        
        #Determine if the universe is closed or not.
        if om_rad[i]+om_mat[i]+om_de[i] > 1:
            
            #Universe is closed
            #Find the value of the maximum scale factor
            max_scale = optimize.minimize(lambda x: adot(x,om_rad[i],om_mat[i],om_de[i],w[i]),1,method='Nelder-Mead',options={'xtol': 10**-10}).x[0]
            
            #Calculate the first half of the universe
            As = np.linspace(0,max_scale,npoints/2)
            times = np.array([time(a,om_rad[i],om_mat[i],om_de[i],w[i]) for a in As])
            hubbles = hubble(As,om_rad[i],om_mat[i],om_de[i],w[i])
            densitys = density(As,om_rad[i],om_mat[i],om_de[i],w[i])
            
            #Put my thing down, flip it, and reverse it
            As = np.hstack((As,As[::-1]))
            times = np.hstack((times,-times[::-1]+2*np.max(times)))
            hubbles = np.hstack((hubbles,-hubbles[::-1]))
            densitys = np.hstack((densitys,densitys[:,::-1]))
            
        else:
            
            #Universe is flat or open
            As = np.linspace(0,20,npoints)
            times = np.array([time(a,om_rad[i],om_mat[i],om_de[i],w[i]) for a in As])
            hubbles = hubble(As,om_rad[i],om_mat[i],om_de[i],w[i])
            densitys = density(As,om_rad[i],om_mat[i],om_de[i],w[i])

        fig = plt.figure()
        ax1 = plt.subplot2grid((2,2), (0,0))
        ax2 = plt.subplot2grid((2,2), (0,1))
        ax3 = plt.subplot2grid((2,2), (1,0))
        ax4 = plt.subplot2grid((2,2), (1,1))

        ax1.plot(times/gyr,As,label=uni[i],color='black')
        ax1.scatter(0,1,label='Today',color='black')
        ax1.set_title('Scale Factor vs. Time')
        ax1.set_xlabel(r'$t$ [Gyr]')
        ax1.set_ylabel(r'$a$')
        ax1.set_xlim([np.min(times)/gyr,np.max(times)/gyr])
        ax1.set_ylim([0,np.max(As)])
        ax1.legend(loc='best', scatterpoints=1,prop={'size':8}, framealpha=0)
        
        ax2.plot(times/gyr,np.arcsinh(hubbles*70/H0),label=uni[i],color='black')
        ax2.scatter(0,np.arcsinh(70),label='Today',color='black')
        ax2.set_title('Hubble "Constant" vs. Time')
        ax2.set_xlabel(r'$t$ [Gyr]')
        ax2.set_ylabel(r'$H$ [km/s/Mpc]')
        newticks = np.sinh(ax2.get_yticks())
        newticks = ["%.1e" % tick for tick in newticks]
        for j in range(len(newticks)):
            [mant,exp] = newticks[j].split('e')
            exp = str(int(exp))
            newticks[j] = r'${!s}\times 10^{!s}$'.format(mant,exp)
        ax2.set_yticklabels(newticks)
        ax2.set_xlim([np.min(times)/gyr,np.max(times)/gyr])
        ax2.legend(loc='best', scatterpoints=1,prop={'size':8}, framealpha=0)
        
        if om_rad[i]>0:
            ax3.semilogy(times/gyr,densitys[0,:],label=uni[i]+': Radiation',color='red')
            ax4.semilogy(times/gyr,densitys[0,:]/(hubbles/H0)**2,label=uni[i]+': Radiation',color='red')
        if om_mat[i]>0:
            ax3.semilogy(times/gyr,densitys[1,:],label=uni[i]+': Matter',color='blue')
            if i==0: #Special case for EdS universe because log plot is bad
                ax4.plot(times/gyr,densitys[1,:]/(hubbles/H0)**2,label=uni[i]+': Matter',color='blue')
                ax4.set_ylim([0,2])
            else:
                ax4.semilogy(times/gyr,densitys[1,:]/(hubbles/H0)**2,label=uni[i]+': Matter',color='blue')
        if om_de[i]>0:
            ax3.semilogy(times/gyr,densitys[2,:],label=uni[i]+': Dark Energy',color='green')
            ax4.semilogy(times/gyr,densitys[2,:]/(hubbles/H0)**2,label=uni[i]+': Dark Energy',color='green')
        ax3.semilogy(times/gyr,densitys[3,:],label=uni[i]+': Total',color='black')
        if i==0: #Special case for EdS universe because log plot is bad
            ax4.plot(times/gyr,densitys[3,:]/(hubbles/H0)**2,label=uni[i]+': Total',color='black')
        else:
            ax4.semilogy(times/gyr,densitys[3,:]/(hubbles/H0)**2,label=uni[i]+': Total',color='black')
        ax3.scatter([0,0,0,0],[om_rad[i],om_mat[i],om_de[i],om_rad[i]+om_mat[i]+om_de[i]],label='Today',color='black')
        if i==0: #Special case for EdS universe because log plot is bad
            ax4.scatter([0],[om_mat[i]],label='Today',color='black')
        else:
            ax4.scatter([0,0,0,0],[om_rad[i],om_mat[i],om_de[i],om_rad[i]+om_mat[i]+om_de[i]],label='Today',color='black')
            
        ax3.set_title('Density vs. Time')
        ax3.set_xlabel(r'$t$ [Gyr]')
        ax3.set_ylabel(r'$\rho/\rho_{c,0}$')
        ax3.set_xlim([np.min(times)/gyr,np.max(times)/gyr])
        ax3.legend(loc='best', scatterpoints=1,prop={'size':8}, framealpha=0)
        
        ax4.set_title('Density vs. Time')
        ax4.set_xlabel(r'$t$ [Gyr]')
        ax4.set_ylabel(r'$\rho/\rho_c(t)$')
        ax4.set_xlim([np.min(times)/gyr,np.max(times)/gyr])
        ax4.legend(loc='best', scatterpoints=1,prop={'size':8}, framealpha=0)
        
        #Save the figure to file and close
        plt.tight_layout()
        plt.savefig(uni[i]+suffix)
        plt.close()
        