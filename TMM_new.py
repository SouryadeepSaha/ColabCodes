import numpy as np
import matplotlib.pyplot as plt


class TransferMatrix(object):
    '''
    Transfer matrices : T and S calculations of a multilayered system.
    '''

    def __init__(self, numberOfLayers=1, freqMin=1e12, freqMax=1e15, gridpoints=101, **RIandLT):

        self.c0 = 299792458.0
        self.numberOfLayers = numberOfLayers
        self.gridPoints = gridpoints
        # if self.gridPoints % 2 == 0:
        #     self.gridPoints += 1
        self.RI = np.ones([self.numberOfLayers + 2, self.gridPoints])  # refractive indices
        self.RIandLT = RIandLT
        self.layerThickness = np.zeros([self.numberOfLayers + 2, self.gridPoints])
        self.freq = np.linspace(freqMin, freqMax, self.gridPoints)
        self.wavelength = self.c0 / self.freq
        self.k = 2 * np.pi / self.wavelength
        self.reflectance = np.zeros(self.gridPoints)  # S_11*conj(S_11)
        self.transmittance = np.zeros(self.gridPoints)  # S_21*conj(S_21)
        self.TMatrix = np.zeros([2, 2, self.gridPoints], dtype='complex')
        self.SMatrix = np.zeros([2, 2, self.gridPoints], dtype='complex')
        self.M = np.zeros([2, 2, self.numberOfLayers + 1, self.gridPoints], dtype='complex')
        self.KLayers = np.ones([self.numberOfLayers + 2, self.gridPoints])

    def setupLayers(self):
        RI1_keys = ['RI1', 'RefractiveIndex1', 'refractiveIndex1', 'ri1', 'refractiveindex1']
        RI2_keys = ['RI2', 'RefractiveIndex2', 'refractiveIndex2', 'ri2', 'refractiveindex2']
        LT1_keys = ['LT1', 'LayerThickness1', 'layerThickness1', 'lt1', 'layerthickness1', 'layer1_thickness',
                    'layerthickness']
        LT2_keys = ['LT2', 'LayerThickness2', 'layerThickness2', 'lt2', 'layerthickness2', 'layer2_thickness',
                    'layerthickness']

        for keys in self.RIandLT:
            if keys in RI1_keys:
                self.RI[1, :] = self.RIandLT[keys]
            if keys in RI2_keys:
                self.RI[2, :] = self.RIandLT[keys]
            if keys in LT1_keys:
                self.layerThickness[1, :] = self.RIandLT[keys]
            if keys in LT2_keys:
                self.layerThickness[2, :] = self.RIandLT[keys]
        for i in range(0, self.numberOfLayers + 2):
            if i % 2 == 1:
                self.RI[i] = self.RI[1]
                self.layerThickness[i, :] = self.layerThickness[1, :]
            if i % 2 == 0:
                self.RI[i] = self.RI[2]
                self.layerThickness[i, :] = self.layerThickness[2, :]
        self.layerThickness[0, :] = 0
        self.RI[0, :] = 1.0
        self.RI[-1, :] = 1.0
        self.layerThickness[-1, :] = 0
        self.KLayers = self.RI * self.k

    def m_Matrices(self):
        self.M[0, 0, :, :] = (self.KLayers[1:, :] + self.KLayers[:-1, :]) / (2 * self.KLayers[1:, :]) \
                             * np.exp(-1j * self.KLayers[:-1, :] * self.layerThickness[:-1, :])
        self.M[0, 1, :, :] = (self.KLayers[1:, :] - self.KLayers[:-1, :]) / (2 * self.KLayers[1:, :]) \
                             * np.exp(1j * self.KLayers[:-1, :] * self.layerThickness[:-1, :])
        self.M[1, 0, :, :] = (self.KLayers[1:, :] - self.KLayers[:-1, :]) / (2 * self.KLayers[1:, :]) \
                             * np.exp(-1j * self.KLayers[:-1, :] * self.layerThickness[:-1, :])
        self.M[1, 1, :, :] = (self.KLayers[1:, :] + self.KLayers[:-1, :]) / (2 * self.KLayers[1:, :]) \
                             * np.exp(1j * self.KLayers[:-1, :] * self.layerThickness[:-1, :])

    def transferMatrices(self):
        self.m_Matrices()
        for i in range(self.gridPoints):
            temp_T = self.M[:, :, 0, i]
            for layer in range(self.numberOfLayers):
                temp_T = np.matmul(temp_T, self.M[:, :, layer + 1, i])
            self.TMatrix[:, :, i] = temp_T
            self.SMatrix[0, 0, i] = -temp_T[0, 1]/temp_T[1, 1]
            self.SMatrix[0, 1, i] = 1/temp_T[1, 1]
            self.SMatrix[1, 0, i] = 1 / temp_T[1, 1]
            self.SMatrix[1, 1, i] = temp_T[0, 1] / temp_T[1, 1]

        # self.SMatrix[0, 0, :] = -self.TMatrix[0, 1, :] / self.TMatrix[1, 1, :]#-self.TMatrix[1, 0, :] / self.TMatrix[1, 1, :]
        # self.SMatrix[0, 1, :] = 1 / self.TMatrix[1, 1, :]
        # self.SMatrix[1, 0, :] = 1 / self.TMatrix[1, 1, :]#self.TMatrix[0, 0, :] - self.TMatrix[0, 1, :] * self.TMatrix[1, 0, :] \
        #                         #/ self.TMatrix[1, 1, :]
        # self.SMatrix[1, 1, :] = self.TMatrix[0, 1, :] / self.TMatrix[1, 1, :]

    def run(self):
        self.setupLayers()
        self.transferMatrices()
        for i in range(self.gridPoints):
            self.reflectance[i] = np.real(self.SMatrix[0, 0, i] * np.conj(self.SMatrix[0, 0, i]))
            self.transmittance[i] = np.real(self.SMatrix[0, 1, i] * np.conj(self.SMatrix[0, 1, i]))

        fig, ax = plt.subplots()
        ax.plot(self.freq, self.reflectance, label='R')
        ax.plot(self.freq, self.transmittance, label='T')
        ax.plot(self.freq, self.reflectance + self.transmittance, label='R+T')
        ax.legend()
        ax.set_xlabel('Frequency [Hz]')
        ax.set_ylabel('Transmission Coeffs')
        # ax.set_title('HR-Coating with $\lambda /4$ layers')

        # Tweak spacing to prevent clipping of ylabel
        fig.tight_layout()
        plt.show()


if __name__ == '__main__':
    wavelength = 980e-9
    frequency = 3e+8 / wavelength
    nSiN = 2.0
    nSiO2 = 1.5
    dSiN = wavelength / (4 * nSiN)
    dSiO2 = wavelength / (4 * nSiO2)
    TMM = TransferMatrix(numberOfLayers=101, freqMin=1e12, freqMax=1e15, gridpoints=1001,
                         RI1=1.7445, RI2=1.81, LT1=4e-7, LT2=4e-7)
    TMM.run()
