'''
Created on 08.05.2018

@author: rall
'''
import numpy as np
import matplotlib.pyplot as plt
class DIFFRACTION_FRAUNHOFER(object):
    '''
    classdocs
    '''

    def __init__(self, slit = 5, dx = 0.1, wavelength = 1, propDistance = 1000 , plot_range = 10000, steps_to_plot = 5):
        '''
        Constructor
        '''
        self.amplitude    = 1     # Volt / sqrt(micron)
        self.slitWidth    = slit     # microns        
        self.propDistance = propDistance # microns (= 10 mm)      
        self.steps_to_plot = steps_to_plot
        self.dx = dx
        self.x = np.arange(-plot_range, plot_range, self.dx) #microns
        self.wavelength = wavelength;
        self.k = 2*np.pi/self.wavelength;
        #self.dx = self.x[1]-self.x[0]
        self.field = np.zeros(self.x.size, dtype = complex)
        self.diffractedField = np.zeros(self.x.size, dtype = complex)

        #self.zR =175e-3;
        #self.subplotindex = 1;
        #self.radius = np.linspace(-self.plot_domain,self.plot_domain,self.gridp);
        #self.dx = self.x[1] - self.x[0] # Spatial sampling period, microns
        self.fS = 1 / self.dx      # Spatial sampling frequency, units are inverse microns
        self.f  = (self.fS / self.x.size) * np.arange(0, self.x.size, step = 1) # inverse microns
        
        #self.xPrime   = np.hstack((self.f[-(self.f.size/2):] - self.fS, self.f[0:self.f.size/2])) * self.wavelength * self.propDistance
        
        self.xPrime = np.array(self.f-np.max(self.f)/2) * self.wavelength * self.propDistance
        self.IrradTheory = 0
        self.IrradFFT = 0
        
    def setupField(self): 
        
        self.field[np.logical_and(self.x >= -self.slitWidth / 2, self.x <= self.slitWidth / 2)] = self.amplitude + 0j
        #print(self.field)
        #plt.plot(self.x, self.field)
        #plt.show()
        print('FRAUNHOFER APPROXIMATION')
        print(['fresnel number: ', 0.25*self.slitWidth*self.slitWidth/self.wavelength/self.propDistance])
        print('self.slitWidth ',self.slitWidth)
        print('self.wavelength ',self.wavelength)
        print('self.propDistance ',self.propDistance)

    def analytical(self):       
        self.IrradTheory = self.amplitude / (self.wavelength * self.propDistance) * \
        (self.slitWidth * np.sinc(self.xPrime * self.slitWidth / self.wavelength / self.propDistance))**2
        dxFourierPlane = self.xPrime[1] - self.xPrime[0];
        #print('power SINC', np.trapz(self.IrradTheory, x=None, dx = dxFourierPlane))
        #plt.plot(self.xPrime, self.IrradTheory)
        #plt.show()
    def fraunhofer(self):
        #self.field_fraunhofer = 0



        self.diffractedField = self.dx * np.fft.fftshift(np.fft.fft(self.field)) / (1j* self.wavelength * self.propDistance) * np.exp(1j*self.k*self.propDistance)*np.exp(1j*self.k / 2 /self.propDistance * np.square(self.x ) )
        #plt.plot(self.x, np.real(self.diffractedField))
        #plt.show()
        self.IrradFFT    = (self.diffractedField * np.conj(self.diffractedField)) / (self.wavelength * self.propDistance)

        dxFourierPlane = self.xPrime[1]-self.xPrime[0];
        #print('power FOUR', np.trapz(np.abs(self.IrradFFT),x = None, dx = np.abs(dxFourierPlane)))
        #print('power SLIT', np.trapz(self.field,x = None, dx =  self.dx))
        #print('numbr GRID', np.size(self.x))
        
#         plt.plot(self.xPrime, np.abs(self.IrradFFT), '.', label = 'FFT')
#         plt.plot(self.xPrime, self.IrradTheory, label = 'Theory')
#         plt.xlim((-70, 70))
#         plt.xlabel(r'x-position, $\mu m$')
#         plt.ylabel(r'Power density, $V^2 / \mu m$')
#         plt.grid(True)
#         plt.legend()
#         plt.show()


    def fresnel(self):
        
        self.diffractedField = self.dx * np.fft.fft((self.field*np.exp(1j*self.k*np.square(self.x)/2/self.propDistance)))
        self.diffractedField = np.fft.fftshift(self.diffractedField)
        self.diffractedField = self.diffractedField / (1j* self.wavelength * self.propDistance) * np.exp(1j*self.k*self.propDistance)*np.exp(1j*self.k / 2 /self.propDistance * np.square(self.x ) )
        
        self.IrradFFT    = (self.diffractedField * np.conj(self.diffractedField)) / self.wavelength / self.propDistance