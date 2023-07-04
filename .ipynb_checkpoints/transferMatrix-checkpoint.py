'''
Created on 29.05.2018

@author: rall
'''

import numpy as np
import matplotlib.pyplot as plt
from numpy import number


class TRANSFERMATRIX(object):
    '''
    classdocs
    '''

    def __init__(self, numberLayers = 101, refractiveIndex1 = 1.7445, refractiveIndex2 = 1.8100, fMin = 1e+12, fMax = 1000e+12, plotPoints = 101, layerThickness1 = 4e-7, layerThickness2 = 4e-7 ):
        '''
        Constructor
        '''
        self.numberLayers = numberLayers
        self.refractiveIndex1 = refractiveIndex1
        self.refractiveIndex2 = refractiveIndex2
        self.fMin = fMin
        self.fMax = fMax
        self.plotPoints = plotPoints
        self.layerThickness1 = layerThickness1
        self.layerThickness2 = layerThickness2
        #self.transferMatrix = np.zeros(2,2,dType = complex)
        self.scatteringMatrix = np.zeros([2,2],dtype = complex)
        self.frequency = np.linspace(fMin, fMax, self.plotPoints)
        
        self.reflectivity = np.zeros(self.plotPoints)
        self.transmission = np.zeros(self.plotPoints)
        self.layerThickness = np.zeros(numberLayers+2)
        self.refractiveIndex = np.zeros(numberLayers+2)
        self.c0 = 3e+8
        self.wavelengthPlot = self.c0/self.frequency
        
              
    def setUpLayers(self):    
        for i in range(0,self.numberLayers+2):
            if(i%2 == 0):
                
                self.refractiveIndex[i] = self.refractiveIndex2 
            else:
               
                self.refractiveIndex[i] = self.refractiveIndex1
        for i in range(1,self.numberLayers+2):
            if(i%2 == 0):
                self.layerThickness[i] = self.layerThickness[i-1]+ self.layerThickness2
            else:
                self.layerThickness[i] = self.layerThickness[i-1]+ self.layerThickness1              
        self.layerThickness[0] = 0.0      
        self.layerThickness[-1] = self.layerThickness[-2]   
        self.refractiveIndex[0] = 1.0
        self.refractiveIndex[-1] = 1.0
        #print(self.layerThickness)
        #print(self.refractiveIndex)

    def runSimulation(self):
        for i in range(0, self.plotPoints):
            #print(["frequency ", i, " of ", self.plotPoints, " calculated"])
            self.setupMatrix(i)
            #self.setupMatrixNew(self.frequency[i])
        # transfer.setupMatrixNew(transfer.frequency[1])
        #fig, ax = plt.plot(self.frequency, self.reflectivity)
        fig, ax = plt.subplots()
        ax.plot(self.frequency, self.reflectivity, self.frequency, self.transmission)
        ax.legend(["R","T"])

        ax.set_xlabel('Frequency [1/s]')
        ax.set_ylabel('Reflectivity')
        ax.set_title('HR-Coating with $\lambda /4$ layers')

        # Tweak spacing to prevent clipping of ylabel
        fig.tight_layout()
        plt.show()





    def setupMatrix(self, frequencyIndex):
        frequency = self.frequency[frequencyIndex]
        k = 2 * np.pi * frequency  / self.c0  
        transferMatrix = np.zeros([2,2],dtype = complex)
        tempMatrix1 = np.zeros([2,2],dtype = complex)
        tempMatrix2 = np.zeros([2,2],dtype = complex)
        transferMatrix[0,0] = k * (self.refractiveIndex[1]+self.refractiveIndex[0])/(2 * k * self.refractiveIndex[1])
        transferMatrix[1,0] = k * (self.refractiveIndex[1]-self.refractiveIndex[0])/(2 * k * self.refractiveIndex[1])
        transferMatrix[0,1] = k * (self.refractiveIndex[1]-self.refractiveIndex[0])/(2 * k * self.refractiveIndex[1])
        transferMatrix[1,1] = k * (self.refractiveIndex[1]+self.refractiveIndex[0])/(2 * k * self.refractiveIndex[1])
        #tempMatrix[0,0]
        #tempMatrix[0,0]
                
        for i in range(0,self.numberLayers+1):
            tempMatrix1[0,0] = k * (self.refractiveIndex[i+1]+self.refractiveIndex[i])/(2 * k * self.refractiveIndex[i+1])
            tempMatrix1[1,0] = k * (self.refractiveIndex[i+1]-self.refractiveIndex[i])/(2 * k * self.refractiveIndex[i+1])
            tempMatrix1[0,1] = k * (self.refractiveIndex[i+1]-self.refractiveIndex[i])/(2 * k * self.refractiveIndex[i+1])
            tempMatrix1[1,1] = k * (self.refractiveIndex[i+1]+self.refractiveIndex[i])/(2 * k * self.refractiveIndex[i+1])
            #print(["tempmatrix ",np.linalg.det(tempMatrix1)])
            tempMatrix2[0,0] = np.exp(-1j* k * self.refractiveIndex[i] * (self.layerThickness[i]-self.layerThickness[i-1]))
            #print(["tempmatrix ",abs(np.exp(-1j* k * self.refractiveIndex[i] * (self.layerThickness[i]-self.layerThickness[i-1]))))
            tempMatrix2[1,1] = np.exp(1j* k * self.refractiveIndex[i] * (self.layerThickness[i]-self.layerThickness[i-1]))
            #print(["tempmatrix ",abs(np.exp(1j* k * self.refractiveIndex[i] * (self.layerThickness[i]-self.layerThickness[i-1]))))
            transferMatrix = np.dot( np.dot(tempMatrix1 , tempMatrix2) , transferMatrix)
            
            #print("determinant for layer num ", i , transferMatrix[0,0]*transferMatrix[1,1]-transferMatrix[1,0]*transferMatrix[0,1])      
            #print("transmission matrix  for layer num ", i , transferMatrix)
            
            #print(["transfermatrix det ",np.linalg.det(transferMatrix)])
        tempMatrix1[0,0] = k * (self.refractiveIndex[-1]+self.refractiveIndex[-2])/(2 * k * self.refractiveIndex[-1])
        tempMatrix1[1,0] = k * (self.refractiveIndex[-1]-self.refractiveIndex[-2])/(2 * k * self.refractiveIndex[-1])
        tempMatrix1[0,1] = k * (self.refractiveIndex[-1]-self.refractiveIndex[-2])/(2 * k * self.refractiveIndex[-1])
        tempMatrix1[1,1] = k * (self.refractiveIndex[-1]+self.refractiveIndex[-2])/(2 * k * self.refractiveIndex[-1])
        transferMatrix = np.dot  (tempMatrix1, transferMatrix)
        
        
        
        self.scatteringMatrix[0,0] = -transferMatrix[0,1]/transferMatrix[1,1]
        self.scatteringMatrix[1,0] = 1.0 / transferMatrix[1,1]
        self.scatteringMatrix[0,1] = 1.0 / transferMatrix[1,1]
        self.scatteringMatrix[1,1] = transferMatrix[0,1]/transferMatrix[1,1]
        
        
        
        #print("determinant ", transferMatrix[0,0]*transferMatrix[1,1]-transferMatrix[1,0]*transferMatrix[0,1])      
        #print("transmission matrix ", transferMatrix)
        
        a = (self.scatteringMatrix[0,0])
        b = (self.scatteringMatrix[1,0])
        ref = self.scatteringMatrix[0,0]
        trn = self.scatteringMatrix[0,1]
        print(["S11^2+S21^2 ",a*np.conjugate(a)+b*np.conjugate(b)])    
        print(["R+T ",ref*np.conjugate(ref)+trn * np.conj(trn)])    
        self.reflectivity[frequencyIndex] = ref*np.conjugate(ref);
        self.transmission[frequencyIndex] = trn * np.conj(trn);
            
    def setupMatrixNew(self, frequency):
        k = 2 * np.pi * frequency  / self.c0  
        transferMatrix = np.zeros([2,2],dtype = complex)
        tempMatrix1 = np.zeros([2,2],dtype = complex)
        tempMatrix2 = np.zeros([2,2],dtype = complex)
        transferMatrix[0,0] = k * (self.refractiveIndex[1]+self.refractiveIndex[0])/(k * (self.refractiveIndex[1]+self.refractiveIndex[2]))
        transferMatrix[1,0] = k * (self.refractiveIndex[1]-self.refractiveIndex[0])/(k * (self.refractiveIndex[1]+self.refractiveIndex[2]))
        transferMatrix[0,1] = k * (self.refractiveIndex[1]-self.refractiveIndex[0])/(k * (self.refractiveIndex[1]+self.refractiveIndex[2]))
        transferMatrix[1,1] = k * (self.refractiveIndex[1]+self.refractiveIndex[0])/(k * (self.refractiveIndex[1]+self.refractiveIndex[2]))
        #tempMatrix[0,0]
        #tempMatrix[0,0]
                
        for i in range(1,self.numberLayers+1):
            tempMatrix1[0,0] = k * (self.refractiveIndex[i+1]+self.refractiveIndex[i])/(k * (self.refractiveIndex[i+1]+self.refractiveIndex[i]))
            tempMatrix1[1,0] = k * (self.refractiveIndex[i+1]-self.refractiveIndex[i])/(k * (self.refractiveIndex[i+1]+self.refractiveIndex[i]))
            tempMatrix1[0,1] = k * (self.refractiveIndex[i+1]-self.refractiveIndex[i])/(k * (self.refractiveIndex[i+1]+self.refractiveIndex[i]))
            tempMatrix1[1,1] = k * (self.refractiveIndex[i+1]+self.refractiveIndex[i])/(k * (self.refractiveIndex[i+1]+self.refractiveIndex[i]))
            tempMatrix2[0,0] = np.exp(-1j* k * self.refractiveIndex[i] * (self.layerThickness[i]-self.layerThickness[i-1]))
            tempMatrix2[1,1] = np.exp(1j* k * self.refractiveIndex[i] * (self.layerThickness[i]-self.layerThickness[i-1]))
            transferMatrix = np.dot( np.dot(tempMatrix1 , tempMatrix2) , transferMatrix)
        tempMatrix1[0,0] = k * (self.refractiveIndex[-1]+self.refractiveIndex[-2])/(k * (self.refractiveIndex[-1]+self.refractiveIndex[-2]))
        tempMatrix1[1,0] = k * (self.refractiveIndex[-1]-self.refractiveIndex[-2])/(k * (self.refractiveIndex[-1]+self.refractiveIndex[-2]))
        tempMatrix1[0,1] = k * (self.refractiveIndex[-1]-self.refractiveIndex[-2])/(k * (self.refractiveIndex[-1]+self.refractiveIndex[-2]))
        tempMatrix1[1,1] = k * (self.refractiveIndex[-1]+self.refractiveIndex[-2])/(k * (self.refractiveIndex[-1]+self.refractiveIndex[-2]))
        transferMatrix = np.dot  (tempMatrix1, transferMatrix)
        self.scatteringMatrix[0,0] = -transferMatrix[0,1]/transferMatrix[1,1]
        self.scatteringMatrix[1,0] = 1.0 / transferMatrix[1,1]
        self.scatteringMatrix[0,1] = 1.0 / transferMatrix[1,1]
        self.scatteringMatrix[1,1] = transferMatrix[0,1]/transferMatrix[1,1]      
        
        #print(["determinante ", np.linalg.det(transferMatrix)])
        
        a = (self.scatteringMatrix[0,0])
        b = (self.scatteringMatrix[1,0])

        #print(["S11^2+S21^2 ",a*np.conjugate(a)+b*np.conjugate(b)]) 
          
    def test(self):
        k = 2 * np.pi * self.fMin
        tempMatrix1 = np.zeros([2,2],dtype = complex)
        tempMatrix2 = np.zeros([2,2],dtype = complex)
        tempMatrix1[0,0] = k * (self.refractiveIndex[2]+self.refractiveIndex[1])/(k * (self.refractiveIndex[1]+self.refractiveIndex[2]))
        tempMatrix1[1,0] = k * (self.refractiveIndex[2]-self.refractiveIndex[1])/(k * (self.refractiveIndex[1]+self.refractiveIndex[2]))
        tempMatrix1[0,1] = k * (self.refractiveIndex[2]-self.refractiveIndex[1])/(k * (self.refractiveIndex[1]+self.refractiveIndex[2]))
        tempMatrix1[1,1] = k * (self.refractiveIndex[2]+self.refractiveIndex[1])/(k * (self.refractiveIndex[1]+self.refractiveIndex[2]))
        tempMatrix2[0,0] = np.exp(-1j* k * self.refractiveIndex[1] * (self.layerThickness[1]-self.layerThickness[0]))
        tempMatrix2[1,1] = np.exp(1j* k * self.refractiveIndex[1] * (self.layerThickness[1]-self.layerThickness[0]))  
        transferMatrix = np.dot(tempMatrix1, tempMatrix2)
        #print(["transfer matrix ", transferMatrix])      
        #print(["determinante ", np.linalg.det(transferMatrix)])
        self.scatteringMatrix[0,0] = -transferMatrix[0,1]/transferMatrix[1,1]
        self.scatteringMatrix[1,0] = 1.0 / transferMatrix[1,1]
        self.scatteringMatrix[0,1] = 1.0 / transferMatrix[1,1]
        self.scatteringMatrix[1,1] = transferMatrix[0,1]/transferMatrix[1,1]
        a = (self.scatteringMatrix[0,0])
        b = (self.scatteringMatrix[1,0])

        #print(["S11^2+S21^2 ",a*np.conjugate(a)+b*np.conjugate(b)])