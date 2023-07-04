import numpy as np
import matplotlib.pyplot as plt
class Diffraction(object):
    
    #classdocs
    

    def __init__(self, slit = 5, n = 1024, wavelength = 1, propDistance = 1000 , domainSize = 10000):
        
        #Constructor
        
        # self. : variables, which are part of class Diffraction
        self.amplitude    = 1     # Volt / sqrt(micron)
        self.slitWidth    = slit           
        self.propDistance = propDistance 
        self.coordinate = np.linspace(-domainSize, domainSize, n)
        print(self.coordinate)
        self.dx = self.coordinate[1]-self.coordinate[0]
        #self.coordinate = np.arange(-domainSize, domainSize+1e-10, self.dx) #microns
        self.wavelength = wavelength;
        self.wavenumber = 2*np.pi/self.wavelength;
        self.field = np.zeros(self.coordinate.size, dtype = complex)
        
        self.fS = 1 / self.dx
        self.f = (self.fS / self.coordinate.size) * np.arange(0, self.coordinate.size, step = 1)
        self.f -= np.max(self.f)/2
        
        self.phasefactorArg = np.zeros(self.coordinate.size, dtype=float)
        self.phasefactor = np.zeros(self.coordinate.size, dtype=complex)
        

        
        self.xprime = self.wavelength * self.propDistance * self.f
        #self.xprime -= np.max(self.xprime)/2
        #self.xprimeflip = np.flip(self.xprime)
        
        self.diffractedFieldFraunhofer = np.zeros(self.coordinate.size, dtype=complex)
        self.irradFraunhofer = np.zeros(self.coordinate.size, dtype=float)
        self.diffractedFieldFresnel = np.zeros(self.coordinate.size, dtype=complex)
        self.irradFresnel = np.zeros(self.coordinate.size, dtype=float)
        self.irradAnalyticalFarField = np.zeros(self.coordinate.size, dtype=float)
        self.angularSpectrumField = np.zeros(self.coordinate.size, dtype=complex)
        
    
    def calcPhasefactor(self):
        self.phasefactorArg = +(self.wavenumber*self.wavenumber) - 4*np.pi*np.pi*(self.f*self.f)
        #self.phasefactorArg = +(1  - 4*np.pi*np.pi*(self.f*self.f))
        
        #plt.plot(self.coordinate, self.phasefactorArg)
        self.phasefactorArg[self.phasefactorArg<0] = 0 #setzt alle neg. auf 0
        #self.phasefactor = 1j*np.sqrt(self.phasefactorArg)*self.propDistance * self.wavenumber
        self.phasefactor = 1j*np.sqrt(self.phasefactorArg)*self.propDistance
        print('pasehfac ', np.max(np.abs(self.phasefactor)))
        print('phasefactorArg ', np.max(np.abs(self.phasefactorArg)))
#         print('Wurzel: ', self.phasefactorArg)
#         print('phasefac :', self.phasefactor)
#         print('phasefac real :', np.real(self.phasefactor))
#         print('phasefac imag :', np.imag(self.phasefactor))
#        plt.plot(self.coordinate, self.phasefactorArg)
#         plt.plot(self.coordinate, np.imag(self.phasefactorArg))
#         plt.plot(self.coordinate, np.real(self.phasefactorArg))
         
    def setupField(self): 
        for i in np.arange(np.size(self.coordinate)):
            if np.abs(self.coordinate[i])<(self.slitWidth/2):
                 self.field[i] = self.amplitude          
            else:
                 self.field[i] = 0
       # self.ampiltude in slit, 0 otherwise
    
    def angularSpectrum(self):
        self.angularSpectrumField = np.fft.fft(self.field)
        
        self.angularSpectrumField = np.fft.fftshift(self.angularSpectrumField)
        self.angularSpectrumField = self.angularSpectrumField * np.exp(self.phasefactor)
        self.angularSpectrumField = np.fft.ifft(self.angularSpectrumField)
        #self.angularSpectrumField /= np.max(np.abs(self.angularSpectrumField))
        plt.plot(self.coordinate, np.real(self.angularSpectrumField))
        plt.plot(self.coordinate, np.imag(self.angularSpectrumField))
        plt.plot(self.coordinate, np.abs(self.angularSpectrumField))
        #self.angularSpectrumField = np.fft.fftshift(self.angularSpectrumField)
        plt.xlim([-2*self.slitWidth, 2*self.slitWidth])
        
        
    def plotAngularSpectrum(self):
        print(self.angularSpectrumField)
        #plt.plot(self.coordinate, np.abs(self.angularSpectrumField)**2)
        plt.xlim([-2*self.slitWidth, 2*self.slitWidth])

        
    
    def fresnelNumber(self):
        return (self.slitWidth/2)**2/(self.wavelength*self.propDistance) 
                
    def analyticalSolution(self):
        #print(self.coordinate)
        #print(self.xprime)
        print(self.xprime.size)
        self.irradAnalyticalFarField = (self.slitWidth * np.sinc((self.slitWidth*self.xprime)/(self.wavelength*self.propDistance)))**2
        #print(self.irradAnalyticalFarField)
        
    def plotAnalyticalSolution(self):
        plt.plot(self.coordinate, self.irradAnalyticalFarField)
        plt.xlim([-2*self.slitWidth, 2*self.slitWidth])
                 
    def fraunhofer(self):   
        self.diffractedFieldFraunhofer = np.exp((1j*self.wavenumber*self.propDistance)/(self.wavelength*self.propDistance)) * np.exp(1j*(self.wavenumber/(2*self.propDistance))*(self.xprime**2)) * np.fft.fft(self.field) * self.dx 
        self.diffractedFieldFraunhofer = np.fft.fftshift(self.diffractedFieldFraunhofer)
        self.diffractedFieldFraunhofer /= np.linalg.norm(self.diffractedFieldFraunhofer)
        self.irradFraunhofer = np.abs(self.diffractedFieldFraunhofer)**2
       
    def plotFraunhofer(self):
        plt.plot(self.xprime, self.irradFraunhofer)
        
    def plotField(self):
        plt.plot(self.coordinate, np.abs(self.field))  
        
    def fresnel(self):
        print(self.coordinate)
        print(self.xprime)
        self.diffractedFieldFresnel = np.exp((1j*self.wavenumber*self.propDistance))/(1j*self.wavelength*self.propDistance) * np.exp(1j*(self.wavenumber/(2*self.propDistance))*(self.xprime**2)) * np.fft.fft(self.field * np.exp(1j*self.wavenumber/(2*self.propDistance)*(self.coordinate**2))) 
        self.diffractedFieldFresnel = np.fft.fftshift(self.diffractedFieldFresnel)
        self.diffractedFieldFresnel /= np.max(np.abs(self.diffractedFieldFresnel))
        self.irradFresnel = np.abs(self.diffractedFieldFresnel)**2
        
    def plotFresnel(self):
        plt.plot(self.xprime, self.irradFresnel)
        plt.xlim([-2*self.slitWidth, 2*self.slitWidth])
        
       