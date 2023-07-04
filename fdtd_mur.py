"""
Created on 24.04.2018

@author: rall
"""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from datetime import datetime


class FDTD_MUR(object):
    """
    classdocs
    """

    def __init__(self, Lx=10, Ly=5, timesteps=50, ppw=15, cfl=0.5, n=1.5, k=0):
        """
        Constructor
        """
        self.mu0 = (np.pi * 4.0e-7)  # [N / (A^2)];
        self.eps0 = 8.854187817e-12  # [(As)/(Vm)]
        self.c0 = 299792458.0  # [m/s]
        self.sig = None
        self.eps = None
        self.H_y = None
        self.H_x = None
        self.E_z = None
        self.E_z_old = None

        # input variables
        self.wavelength = 500e-9
        self.ppw = ppw
        self.cfl = cfl
        self.timesteps = timesteps
        self.actualStep = 0

        self.timestart = 0
        self.timeend = 0

        # vacuum
        self.n = 1.0
        self.k = 0.0
        self.rho = 1.0

        # "artificial" material for our obstacle
        self.n_slit = n
        self.k_slit = k

        c = self.c0 / self.n

        self.obstacle = False

        self.Lx = Lx * self.wavelength
        self.Ly = Ly * self.wavelength

        self.dx = self.wavelength / self.ppw
        self.dy = self.wavelength / self.ppw

        self.Nx = int(self.Lx / self.dx)
        self.Ny = int(self.Ly / self.dy)
        # self.source_y = int(self.Ny/2)
        # self.source_x = int(self.Nx/2)
        self.source_y = np.arange(int(self.Ny * 0.45), int(self.Ny * 0.65))
        # self.source_y = int(self.Ny/2)
        self.source_y_2 = int(self.Ny * 0.25)
        self.source_x = int(self.Nx * 0.25)
        self.source_x_2 = np.arange(int(self.Nx * 0.25), int(self.Nx * 0.55))

        self.omega = 2.0 * np.pi / self.wavelength * self.c0
        # print(self.omega)
        self.mu = self.mu0 * 1.0
        self.epsilon = self.eps0 * (self.n * self.n - self.k * self.k)
        self.sigma = 4.0 * np.pi * np.sqrt(self.eps0 / self.mu0) * (self.n * self.k) / self.wavelength

        '''What is CFL?'''
        self.dt = self.cfl * np.sqrt(self.mu * self.epsilon) / np.sqrt(
            1.0 / (self.dx * self.dx) + 1.0 / (self.dy * self.dy))
        #print(self.dt)
        # self.cfl = 0.9
        #self.dt = self.cfl * np.sqrt(self.mu * self.epsilon) * self.dx / np.sqrt(2) * self.n_slit
        #print(self.dt)
        #self.dt = self.cfl / c * self.dx / np.sqrt(2) / 1.0
        #print(self.dt)

        print('cfl-condition: ', c * self.dt * np.sqrt(1.0 / (self.dx * self.dx) + 1.0 / (self.dy * self.dy)))

        self.tau = self.wavelength / 3 / self.c0 / (np.sqrt(n * n + k * k))


        self.mur_const = (self.c0 * self.dt - self.dx) / (self.c0 * self.dt + self.dx)
        self.mur = True
        self.periodic = False
        self.reflective = False

        print(['time', self.dt])

    def createStaggeredGrid(self):
        self.E_z = np.zeros([self.Ny + 1, self.Nx + 1])
        self.E_z_old = np.zeros([self.Ny + 1, self.Nx + 1])
        self.H_x = np.zeros([self.Ny, self.Nx + 1])
        self.H_y = np.zeros([self.Ny + 1, self.Nx])
        self.eps = np.ones([self.Ny + 1, self.Nx + 1]) * self.epsilon
        self.sig = np.ones([self.Ny + 1, self.Nx + 1]) * self.sigma

    def setBoundaryConditionMur(self):
        self.mur = True
        self.periodic = False
        self.reflective = False

    def setBoundaryConditionPeriodic(self):
        self.mur = False
        self.periodic = True
        self.reflective = False

    def setBoundaryConditionReflective(self):
        self.mur = False
        self.periodic = False
        self.reflective = True

    def updateEzEfficient(self):
        self.E_z_old = np.copy(self.E_z)
        Ez_const = 1.0 / (1.0 + ((self.dt / 2. * self.sig[1:-1, 1:-1]) / self.eps[1:-1, 1:-1]))

        #         Hy_E = self.H_y[1:, 1:-1]
        #         Hy_W = self.H_y[:-1, 1:-1]
        #         Hx_S = self.H_x[1:-1, :-1]
        #         Hx_N = self.H_x[1:-1, 1:]

        #         self.E_z[1:-1, 1:-1] = self.dt /1.0 / self.eps[1:-1, 1:-1] * ((((Hx_S) - (Hx_N)) / self.dy \
        #                              - ((Hy_W) - (Hy_E)) / self.dx) - self.sig[1:-1, 1:-1] *self.E_z[1:-1, 1:-1])+ self.E_z[1:-1, 1:-1];

        Hy_E = self.H_y[1:-1, 1:]
        Hy_W = self.H_y[1:-1, :-1]
        Hx_S = self.H_x[:-1, 1:-1]
        Hx_N = self.H_x[1:, 1:-1]
        #         print(np.shape(Hy_E))
        #         print(np.shape(Hy_W))
        #         print(np.shape(Hx_N))
        #         print(np.shape(Hx_S))
        #         print(np.shape( self.E_z[1:-1, 1:-1]))

        self.E_z[1:-1, 1:-1] = Ez_const*(self.E_z_old[1:-1, 1:-1] +
                                         self.dt / 1.0 / self.eps[1:-1, 1:-1] * (
                                        ((Hx_S - Hx_N) / self.dy + (Hy_E - Hy_W) / self.dx)))

        # self.E_z[1:-1, 1:-1] = self.dt / 1.0 / self.eps[1:-1, 1:-1] * (
        #         ((Hx_S - Hx_N) / self.dy + (Hy_E - Hy_W) / self.dx)
        #         - self.sig[1:-1, 1:-1] * self.E_z[1:-1, 1:-1]) + self.E_z[1:-1, 1:-1]

        # print('Ez new')
        # print(self.E_z[y_iter_Ez, x_iter_Ez])

    def updateHyEfficient(self):

        #         self.H_y = self.H_y + (self.dt / self.mu) \
        #                                          * ((self.E_z[1:,:]) - (self.E_z[:-1,:])) / self.dx;

        self.H_y = self.H_y + (self.dt / self.mu) * (self.E_z[:, 1:] - self.E_z[:, :-1]) / self.dx

    def updateHxEfficient(self):

        #         self.H_x = self.H_x - (self.dt / self.mu) \
        #                                          * ((self.E_z[:, 1:]) - (self.E_z[:, :-1])) / self.dy;

        self.H_x = self.H_x - (self.dt / self.mu) * (self.E_z[1:, :] - self.E_z[:-1, :]) / self.dy



    def setObstacle(self, xMin=0.55, xMax=0.85, yMin=0.15, yMax=0.75):
        self.obstacle = True
        self.eps[int(self.Ny * yMin):int(self.Ny * yMax), int(self.Nx * xMin):int(self.Nx * xMax)] \
            = self.eps0 * (self.n_slit * self.n_slit - self.k_slit * self.k_slit)
        self.sig[int(self.Ny * yMin):int(self.Ny * yMax), int(self.Nx * xMin):int(self.Nx * xMax)] \
            = 4.0 * np.pi * np.sqrt(self.eps0 / self.mu0) * (self.n_slit * self.k_slit) / self.wavelength

    def applyMurBoundaryCondition(self):
        # print('using mur bc')
        self.E_z[0, :] = self.E_z_old[1, :] + self.mur_const * (self.E_z[1, :] - self.E_z_old[0, :])
        self.E_z[-1, :] = self.E_z_old[-2, :] + self.mur_const * (self.E_z[-2, :] - self.E_z_old[-1, :])
        self.E_z[:, 0] = self.E_z_old[:, 1] + self.mur_const * (self.E_z[:, 1] - self.E_z_old[:, 0])
        self.E_z[:, -1] = self.E_z_old[:, -2] + self.mur_const * (self.E_z[:, -2] - self.E_z_old[:, -1])

    def updateEzPeriodic(self, x_iter_Ez, y_iter_Ez):
        Ez_const = 1.0 / (1.0 + ((self.dt / 2. * self.sig[y_iter_Ez, x_iter_Ez]) / self.eps[y_iter_Ez, x_iter_Ez])) #1/(1+dt*sig/eps)
        Hx_N = self.H_y[y_iter_Ez % self.Ny, x_iter_Ez % (self.Nx - 1)]
        Hx_S = self.H_y[y_iter_Ez - 1 % self.Ny, x_iter_Ez % (self.Nx - 1)]
        Hy_W = self.H_x[y_iter_Ez % (self.Ny - 1), x_iter_Ez - 1 % self.Nx]
        Hy_E = self.H_x[y_iter_Ez % (self.Ny - 1), x_iter_Ez % self.Nx]

        self.E_z[y_iter_Ez, x_iter_Ez] = Ez_const * (self.E_z[y_iter_Ez, x_iter_Ez]
                                                     + (self.dt / 2. / self.eps[y_iter_Ez, x_iter_Ez])
                                                     * ((Hx_S - Hx_N) / self.dy + (Hy_E - Hy_W) / self.dx))

    def applyPeriodicBoundaryCondition(self):
        for x_iter_Ez in range(0, self.Nx + 1):
            for y_iter_Ez in range(0, self.Ny + 1):
                if y_iter_Ez == 0:
                    self.updateEzPeriodic(x_iter_Ez, y_iter_Ez)
                    self.E_z[self.Ny, x_iter_Ez] = self.E_z[y_iter_Ez, x_iter_Ez]
                if x_iter_Ez == 0:
                    self.updateEzPeriodic(x_iter_Ez, y_iter_Ez)
                    self.E_z[y_iter_Ez, self.Nx] = self.E_z[y_iter_Ez, x_iter_Ez]

    def updateFieldEfficient(self):
        self.actualStep += 1
        t = self.actualStep * self.dt
        '''
        source
        in python: array(y,x) !!! first arg y, second arg x: row major.
        '''
        # gauss = np.exp(-((t-6*self.tau)/self.tau)**2)
        self.E_z[self.source_y, self.source_x] = (1.0 * np.sin(self.omega * t))
        # self.E_z[ self.source_y , self.source_x_2 ] = ( 1.0 * np.sin( self.omega * t ) )

        '''
        Ez-Field update
        '''

        self.updateEzEfficient()

        '''
        Ez-Field boundary: mur, periodic or reflective
        '''
        if self.mur:
            self.applyMurBoundaryCondition()

        elif self.periodic:
            self.applyPeriodicBoundaryCondition()

        '''
        For reflective, nothing has to be changed. boundary is and stays zero!
        '''

        '''
        Hx-Field update
        '''

        self.updateHxEfficient()

        '''
        Hy-Field update
        '''
        self.updateHyEfficient()
        # self.E_z_old = self.E_z

    def updateEz(self, x_iter_Ez, y_iter_Ez):
        self.E_z_old[y_iter_Ez, x_iter_Ez] = self.E_z[y_iter_Ez, x_iter_Ez]

        Ez_const = 1.0 / (1.0 + ((self.dt / 2. * self.sig[y_iter_Ez, x_iter_Ez]) / self.eps[y_iter_Ez, x_iter_Ez]))

        Hx_N = self.H_y[y_iter_Ez, x_iter_Ez]
        Hx_S = self.H_y[y_iter_Ez - 1, x_iter_Ez]
        Hy_W = self.H_x[y_iter_Ez, x_iter_Ez - 1]
        Hy_E = self.H_x[y_iter_Ez, x_iter_Ez]

        self.E_z[y_iter_Ez, x_iter_Ez] = Ez_const * (self.E_z[y_iter_Ez, x_iter_Ez]
                                                     + (self.dt / 2. / self.eps[y_iter_Ez, x_iter_Ez])
                                                     * ((Hx_S - Hx_N) / self.dy + (Hy_E - Hy_W) / self.dx))
        # print('Ez old')
        # print(self.E_z[y_iter_Ez, x_iter_Ez])
        self.E_z[y_iter_Ez, x_iter_Ez] = self.dt / 1.0 / self.eps[y_iter_Ez, x_iter_Ez] * (
                ((Hx_S - Hx_N) / self.dy + (Hy_E - Hy_W) / self.dx)
                - self.sig[y_iter_Ez, x_iter_Ez] * self.E_z[y_iter_Ez, x_iter_Ez]) + self.E_z[y_iter_Ez, x_iter_Ez]
        # print('Ez new')
        # print(self.E_z[y_iter_Ez, x_iter_Ez])

    def updateHy(self, x_iter_Hy, y_iter_Hy):
        self.H_y[y_iter_Hy, x_iter_Hy] = self.H_y[y_iter_Hy, x_iter_Hy] \
                                         + (self.dt / self.mu) * ((self.E_z[y_iter_Hy, x_iter_Hy])
                                                                  - (self.E_z[y_iter_Hy + 1, x_iter_Hy])) / self.dx

    def updateHx(self, x_iter_Hx, y_iter_Hx):
        self.H_x[y_iter_Hx, x_iter_Hx] = self.H_x[y_iter_Hx, x_iter_Hx] - (self.dt / self.mu) \
                                         * ((self.E_z[y_iter_Hx, x_iter_Hx]) - (
            self.E_z[y_iter_Hx, x_iter_Hx + 1])) / self.dy

    def updateField(self):
        self.actualStep += 1
        t = self.actualStep * self.dt
        '''
        source
        in python: array(y,x) !!! first arg y, second arg x
        '''
        # gauss = np.exp(-((t-6*self.tau)/self.tau)**2)
        self.E_z[self.source_y, self.source_x] = (1.0 * np.sin(self.omega * t))

        '''
        Ez-Field update
        '''
        for x_iter_Ez in range(1, self.Nx):
            for y_iter_Ez in range(1, self.Ny):
                self.updateEz(x_iter_Ez, y_iter_Ez)

        '''
        Ez-Field boundary: mur, periodic or reflective
        '''
        if self.mur:
            self.applyMurBoundaryCondition()
            # for x_iter_Ez in range(0, self.Nx + 1):
            #     for y_iter_Ez in range(0, self.Ny + 1):
            #         if y_iter_Ez == 0:
            #             self.E_z[y_iter_Ez, x_iter_Ez] = self.E_z_old[y_iter_Ez + 1, x_iter_Ez] + self.mur_const * (
            #                     self.E_z[y_iter_Ez + 1, x_iter_Ez] - self.E_z[y_iter_Ez, x_iter_Ez])
            #
            #         if y_iter_Ez == self.Ny:
            #             self.E_z[y_iter_Ez, x_iter_Ez] = self.E_z_old[y_iter_Ez - 1, x_iter_Ez] + self.mur_const * (
            #                     self.E_z[y_iter_Ez - 1, x_iter_Ez] - self.E_z[y_iter_Ez, x_iter_Ez])
            #
            #         if x_iter_Ez == 0:
            #             self.E_z[y_iter_Ez, x_iter_Ez] = self.E_z_old[y_iter_Ez, x_iter_Ez + 1] + self.mur_const * (
            #                     self.E_z[y_iter_Ez, x_iter_Ez + 1] - self.E_z[y_iter_Ez, x_iter_Ez])
            #
            #         if x_iter_Ez == self.Nx:
            #             self.E_z[y_iter_Ez, x_iter_Ez] = self.E_z_old[y_iter_Ez, x_iter_Ez - 1] + self.mur_const * (
            #                     self.E_z[y_iter_Ez, x_iter_Ez - 1] - self.E_z[y_iter_Ez, x_iter_Ez])

        elif self.periodic:
            self.applyPeriodicBoundaryCondition()
            # for x_iter_Ez in range(0, self.Nx + 1):
            #     for y_iter_Ez in range(0, self.Ny + 1):
            #         if y_iter_Ez == 0:
            #             self.updateEzPeriodic(x_iter_Ez, y_iter_Ez)
            #             self.E_z[self.Ny, x_iter_Ez] = self.E_z[y_iter_Ez, x_iter_Ez]
            #         if x_iter_Ez == 0:
            #             self.updateEzPeriodic(x_iter_Ez, y_iter_Ez)
            #             self.E_z[y_iter_Ez, self.Nx] = self.E_z[y_iter_Ez, x_iter_Ez]

        # for reflective, nothing has to be changed. boundary is and stays zero!
        # elif( self.reflective):

        '''
        Hx-Field update
        '''
        for x_iter_Hx in range(0, self.Nx):
            for y_iter_Hx in range(0, self.Ny + 1):
                self.updateHx(x_iter_Hx, y_iter_Hx)

        '''
        Hy-Field update
        '''
        for x_iter_Hy in range(0, self.Nx + 1):
            for y_iter_Hy in range(0, self.Ny):
                self.updateHy(x_iter_Hy, y_iter_Hy)

    def runSimulation(self, plotSteps=10):
        self.timestart = datetime.now()

        for t in range(0, self.timesteps):
            self.updateFieldEfficient()
            if t % plotSteps == 0:
                print(['timestep ', t])
                plt.subplot(1, 1, 1)
                plt.title('E_z')
                plt.pcolormesh(self.E_z, vmin=-0.5, vmax=0.5, cmap='jet')
                # print([ 'Ez max min',np.max(((fdtd.E_z))),np.min(((fdtd.E_z)))])

                #                 plt.subplot(3, 1, 2)
                #                 plt.title('H_x')
                #                 plt.pcolormesh((self.H_x))  # , vmin = -1, vmax = 1)
                #                 # print(['Hx max min',np.max(((fdtd.H_x))),np.min(((fdtd.H_x)))])

                #                 plt.subplot(3, 1, 3)
                #                 plt.title('H_y')
                #                 plt.pcolormesh((self.H_y))  # , vmin = -1, vmax = 1)
                #                 # print(['Hy max min',np.max(((fdtd.H_y))),np.min(((fdtd.H_y)))])

                plt.pause(0.05)

                # cbar.remove()
                # cbar = plt.colorbar(plot)
            # if (t == 0):
            # plt.pcolor(X, Y, f(data), cmap=cm, vmin=-4, vmax=4)
            # plt.colorbar()
            # plt.clim(-4,4)
            # self.updateField()


        self.timeend = datetime.now()

        plt.show()

    def printTime(self):
        print('time for simulation', (self.timeend - self.timestart).total_seconds())
