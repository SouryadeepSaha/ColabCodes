'''
Created on 08.05.2018

@author: rall
'''
import numpy as np
import matplotlib.pyplot as plt
class DIFFRACTION_FRESNEL(object):
    '''
    classdocs
    '''

    def __init__(self, slit = 5, dx = 0.01, wavelength = 1, propDistance = 5 , plot_range = 5, steps_to_plot = 5):
        '''
        Constructor
        '''
        #fresnelratio = 100
        self.amplitude    = 1     # Volt / sqrt(micron)
        self.slitWidth    = slit     # microns        
        self.propDistance = propDistance # microns (= 10 mm)      
        self.steps_to_plot = steps_to_plot
        self.dx = dx
        self.x = np.arange(-plot_range, plot_range, self.dx) #microns
        #print(self.x)
        self.wavelength = wavelength;
        self.k = 2*np.pi/self.wavelength;
        self.field = np.zeros(self.x.size, dtype = complex)
        self.diffractedField = np.zeros(self.x.size, dtype = complex)

        #self.subplotindex = 1;
        #self.radius = np.linspace(-self.plot_domain,self.plot_domain,self.gridp);
        self.fS = 1 / self.dx      # Spatial sampling frequency, units are inverse microns
        #self.f  = (self.fS /(plot_range*2)) * np.arange(0, plot_range*2, step = self.dx)
        self.f  = (self.fS / self.x.size) * np.arange(0, self.x.size, step = 1) # inverse microns
         # inverse microns
        #self.xPrime   = np.hstack((self.f[-(self.f.size/2):] - self.fS, self.f[0:self.f.size/2])) * self.wavelength * self.propDistance
        
        self.xPrime = np.array(self.f-np.max(self.f)/2) * self.wavelength * self.propDistance
        
        #print(self.xPrime)
        self.IrradTheory = 0
        self.IrradFFT = 0
        
    def setupField(self): 
        
        
        #setting slit of width "self.slitWidth"
        self.field[np.logical_and(self.x > -self.slitWidth / 2, self.x <= self.slitWidth / 2)] = self.amplitude + 0j              
        # w0 = 2.5;
        # zR = np.pi * w0 * w0 / self.wavelength;
        # Rz = self.propDistance * (1 + np.square(zR / self.propDistance));
        # #self.propDistance = 0.000001;
        # wz = w0 * np.sqrt(1 + np.square(self.propDistance / zR));
        print('FRESNEL APPROXIMATION')
        print(['fresnel number: ', 0.25*self.slitWidth*self.slitWidth/self.wavelength/self.propDistance])
        print('self.slitWidth ',self.slitWidth)
        print('self.wavelength ',self.wavelength)
        print('self.propDistance ',self.propDistance)
        #self.field = np.exp(-self.x * self.x / wz / wz) * np.exp(-1j * self.k * self.x * self.x / 2.0 / Rz);
        #set two slits
        #self.field[np.logical_and(self.x > -self.slitWidth / 4, self.x <= self.slitWidth / 4)] = 0 + 0j


    def analytical(self):       
        # self.IrradTheory = self.amplitude / (self.wavelength * self.propDistance) * \
        # (self.slitWidth * np.sinc(self.xPrime * self.slitWidth / self.wavelength / self.propDistance))**2
        # plt.plot(self.xPrime, self.IrradTheory)
        # plt.show()
        
        self.IrradTheory = self.amplitude / (self.wavelength * self.propDistance) * \
        (self.slitWidth * np.sinc(self.xPrime * self.slitWidth / self.wavelength / self.propDistance))**2
    def fresnel(self):
        
        self.diffractedField = self.dx * np.fft.fft((self.field*np.exp(1j*self.k*np.square(self.x)/2/self.propDistance)))
        self.diffractedField = np.fft.fftshift(self.diffractedField)
        self.diffractedField = self.diffractedField / (1j* self.wavelength * self.propDistance) * np.exp(1j*self.k*self.propDistance)*np.exp(1j*self.k / 2 /self.propDistance * np.square(self.x ) )
        
        self.IrradFFT    = (self.diffractedField * np.conj(self.diffractedField)) / self.wavelength / self.propDistance
   