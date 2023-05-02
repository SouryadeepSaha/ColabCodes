'''
Created on 08.05.2018

@author: rall
'''
import numpy as np
import matplotlib.pyplot as plt

class GAUSSLET(object):
    '''
    classdocs
    '''


    def __init__(self, slit = 1, gridp = 501, wavelength = 1, nrGausslet = 91, propagationDistance = 1 , plot_domain_factor = 5, steps_to_plot = 5):
        '''
        Constructor
        '''
    #%gaussian beam initlal after prpagataing a distance zR

        self.waist = slit/2;
        
        
        self.plot_domain = plot_domain_factor*self.waist;
        self.steps_to_plot = steps_to_plot
        self.gridp = gridp;
        self.wavelength = wavelength;
        self.k = 2*np.pi/self.wavelength;
        self.nrGausslet = nrGausslet;
        self.waist_gausslet = self.waist / self.nrGausslet*2.5;
        self.gausslet = np.zeros([self.nrGausslet, self.gridp],dtype=complex);
        self.gausslet_sum = np.zeros(self.gridp,dtype=complex);
        self.z = propagationDistance;
        self.rayleigh = np.pi*self.waist_gausslet*self.waist_gausslet/self.wavelength;
        self.waistofsmallgausslet = self.waistAt(self.waist_gausslet,175e-3, self.wavelength);

        self.zR =175e-3;#pi*self.waist*self.waist/self.wavelength;
        self.subplotindex = 1;
        self.radius = np.linspace(-self.plot_domain,self.plot_domain,self.gridp);
        self.shift = np.linspace(-self.waist, self.waist, self.nrGausslet);
        self.scale = np.exp(-2*self.shift**2/self.waist**2)*0+0.1;
        self.distance_gausslets = np.abs(self.shift[1]-self.shift[2]);
        self.ripple_ratio = self.waist_gausslet / self.distance_gausslets;
        self.IrrGausslet =  np.zeros(self.gridp,dtype=complex);
        print(self.ripple_ratio)
        
    def calcGausslet(self):
        print(['fresnel number: ', self.waist*self.waist/self.wavelength/self.z])
        
        for z in np.linspace(0,self.z,self.steps_to_plot):    
            #self.gauss = (self.waist**2/ self.waistAt(self.waist,z, self.wavelength)**2)*np.exp(-2*self.radius**2/self.waistAt(self.waist,z, self.wavelength)**2);
            
            for i in range(0,self.nrGausslet):
                #self.gausslet[i,:] = self.scale[i]*(self.waist_gausslet**2/waistAt(self.waist_gausslet,z, self.wavelength)**2)*np.exp \
                #    (-2*(radius-shift[i])**2/waistAt(self.waist_gausslet,z, self.wavelength)**2);
            
                self.gausslet[i,:] = self.scale[i]*(self.waist_gausslet/self.waistAt(self.waist_gausslet,z, self.wavelength))*np.exp \
                    (-((self.radius-self.shift[i])/self.waistAt(self.waist_gausslet,z, self.wavelength))**2 \
                     -1j*(self.k*z+self.k*(self.radius-self.shift[i])**2/(2*self.curvature(self.waist_gausslet, z, self.wavelength))));
                     #%gausslet_sum = gausslet_sum + gausslet(i,:);
            #self.plotResults()         
    
    def plotResults(self):

        for i in range(self.nrGausslet):
            plt.plot(self.radius,np.abs(self.gausslet[i,:]))
            #print(self.gausslet[i,:])

        plt.show()

        plt.subplot(1,3,self.subplotindex)
        self.IrrGausslet = ((np.sum(self.gausslet[:,:],0))*np.conj(np.sum(self.gausslet[:,:],0)))
        plt.plot(self.radius,self.IrrGausslet)
        plt.legend("intensity")        
        self.subplotindex = self.subplotindex +1;
     
        plt.subplot(1,3,self.subplotindex)          
        self.subplotindex = self.subplotindex +1;
        plt.plot(self.radius,np.real(np.sum(self.gausslet[:,:],0)))
        plt.legend("efield")
        
        plt.subplot(1,3,self.subplotindex)        
        self.subplotindex = self.subplotindex +1;
        plt.plot(self.radius,np.angle(np.sum(self.gausslet[:,:],0)))
        plt.legend("phase")
        
        #if(self.subplotindex == 3*self.steps_to_plot+1):
        plt.show()
        
    def waistAt(self, waist, z , wavelength):
        w = waist *np.sqrt(1+z**2/(np.pi*waist**2/wavelength)**2);
        return w

    def curvature(self, waist, z, wavelength):
        if z==0:
            R = float("inf")
        else:
            R = z*(1+((np.pi*waist**2/wavelength)/z)**2)
        return R
