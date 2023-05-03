'''
Created on 11.04.2018

@author: er96apow, aq28ebyp
'''
import numpy as np
import matplotlib.pyplot as plt

class GaussSeidel(object):
    '''
    classdocs
    '''


    def __init__(self, h_grid = 0.01, x_len = 1.0, y_len = 1.0, i_iter = 100):
        '''
        Constructor
        '''
        self.x_len = x_len
        self.y_len = y_len
        self.i_iter = i_iter
        self.h_x = h_grid 
        self.h_y = h_grid
        self.n_x = int(self.x_len / self.h_x)
        self.n_y = int(self.y_len / self.h_y)
        # assert self.h_x==self.h_y, f'Grid sizes do not match. $h_x$={self.h_x}, $h_y$={self.h_y}'
        
   
    def run(self):     
        n_x = self.n_x
        n_y = self.n_y        
        h_x = self.h_x
        h_y = self.h_y
        # assert h_x==h_y, f'Grid sizes do not match. $h_x$={h_x}, $h_y$={h_y}'
        
        self.u_exact = np.zeros((n_y+1)*(n_x+1))
        self.u_h     = np.zeros((n_y+1)*(n_x+1))
        self.f       = np.zeros((n_y+1)*(n_x+1))
        
        # 1. Schritt exakte Loesung setzen
        for y in range(0,n_y+1):
            for x in range(0,n_x+1):  
                self.u_exact[x + (n_x+1)*y] = np.sin(np.pi*float(x)*h_x) * np.sinh(np.pi*float(y)*h_y) / np.sinh(np.pi)
        
        
        # 2. Randwerte setzen
        for x in range(0,n_x+1):
            self.u_h[x + (n_x+1)*n_y] = self.u_exact[x + (n_x+1)*n_y]
            self.f[x + (n_x+1)*n_y]   = self.u_exact[x + (n_x+1)*n_y]
        
        errorVector = np.zeros([self.i_iter,1])

        error = 0
        # 3. Loesen
        #python indexing starts at 0!
        for i in range(0,self.i_iter):
            #print('l_2 norm : ', error)
            # print('iteration step : ', i)
            for y in range(1,n_y):
                # for x in range(1,n_x):            
                #     self.u_h[x + (n_x+1)*y] = self.f[x + (n_x+1)*n_y] * 0.5 * h_x * h_x +0.5 * (self.u_h[x   + (n_x+1)*(y-1)] \
                #                                           +self.u_h[x + (n_x+1)*(y+1)])
            
                for x in range(1,n_x):
                    self.u_h[x + (n_x+1)*y] =  0.25 * (self.f[x + (n_x+1)*y]*h_x**2\
                                                       +self.u_h[x-1 + (n_x+1)* y] \
                                                       +self.u_h[x+1 + (n_x+1)* y] \
                                                       +self.u_h[x   + (n_x+1)*(y-1)] \
                                                       +self.u_h[x   + (n_x+1)*(y+1)])
            # 4. l_2 norm
            error = np.sqrt(np.sum(np.square(self.u_h-self.u_exact))/(float(n_y+1)*float(n_x+1)))
            errorVector[i] = error


        plt.title('Log error vs iterations')
        plt.plot( np.arange(self.i_iter), np.log10(errorVector))
        plt.show()


    def runSlice(self):  
        n_x = self.n_x
        n_y = self.n_y        
        h_x = self.h_x
        h_y = self.h_y
        # assert h_x==h_y, f'Grid sizes do not match. h_x={h_x}, h_y={h_y}'
        self.u_exact = np.zeros((n_y+1)*(n_x+1)).reshape((n_y+1),(n_x+1))      
        self.u_h     = np.zeros((n_y+1)*(n_x+1)).reshape((n_y+1),(n_x+1)) 
        self.f       = np.zeros((n_y+1)*(n_x+1)).reshape((n_y+1),(n_x+1)) 
        errorVector = np.zeros([self.i_iter,1])
        x = np.linspace(0, 1, self.n_x+1)
        y = np.linspace(0, 1, self.n_x+1)
        xv, yv = np.meshgrid(x, y)
        
        # 1. Schritt exakte Loesung setzen  
        self.u_exact = np.sin(np.pi*xv) * np.sinh(np.pi*yv) / np.sinh(np.pi)
                #self.u_exact[y][x] = y
        
        
        # 2. Randwerte setzen      
        self.u_h[n_y][:] = self.u_exact[n_y][:]    
        #self.u_h[n_y][:] = np.sin(np.pi*xv[n_y][:]) /np.pi/np.pi      
        self.f[n_y][:]   = self.u_exact[n_y][:]      
        #self.f[n_y][:]   = np.sin(np.pi*xv[n_y][:])

        # 3. Loesen
        #python indexing starts at 0!
        for i in range(0,self.i_iter):
            #print('l_2 norm : ', error)
            #print('iteration step : ', i) 
            #Stencil implementation without for-loops
            self.u_h[1:-1, 1:-1] = (self.u_h[2:, 1:-1] + self.u_h[:-2, 1:-1] + self.u_h[1:-1, 2:] + self.u_h[1:-1, :-2])/4.0   
            
            # 4. l_2 norm
            errorVector[i] = np.sqrt(np.divide(np.sum(np.square(self.u_h-self.u_exact)),(float(n_y+1)*float(n_x+1))))
            
            # plt.pcolormesh(self.u_h,edgecolors='k', linewidth=1)
            # plt.pause(0.01)
        plt.title('Log error vs iterations')
        plt.plot( np.arange(self.i_iter), np.log10(errorVector))
        plt.show()
        
    def plotLaplace(self):          
        #self.u_laplace = np.zeros((self.n_y-1)*(self.n_x-1)).reshape((self.n_y-1),(self.n_x-1)) 
        self.u_laplace[:, :] = (self.u_h[2:, 1:-1] + self.u_h[:-2, 1:-1] + self.u_h[1:-1, 2:] + self.u_h[1:-1, :-2] - 4.0 *self.u_h[1:-1, 1:-1])/self.h_x / self.h_y 
        self.u_laplace = np.zeros((self.n_y-1)*(self.n_x+1)).reshape((self.n_y-1),(self.n_x+1)) 
        # self.u_laplace[:, :] = (self.u_h[:-2,:] + self.u_h[2:,:] - 2.0 *self.u_h[1:-1,:]) / self.h_y 
        # self.u_laplace[:, :] = (self.u_h[-1:,:] + self.u_h[1:,:] + self.u_h[:,-1:] + self.u_h[:,1:] - 4.0*self.u_h[:,:])/self.h_x**2
        
        fig, (ax1, ax2) = plt.subplots(1,2)

        
        ax1.set_title('$\\nabla^2u$')
        im1 = ax1.imshow(self.u_laplace) 
        fig.colorbar(im1, ax=ax1)
        ax2.set_title('$f(x,y)$')
        im2 = ax2.imshow(self.f)
        fig.colorbar(im2, ax=ax2)
               
        plt.show()            
                      
        
    def run2d(self):     
        n_x = self.n_x      
        h_x = self.x_len / self.n_x
        self.u_exact = np.zeros((n_x+1))
        self.u_h     = np.zeros((n_x+1))
        
        # 1. Schritt exakte Loesung setzen
    
        for x in range(0,n_x+1):
            self.u_exact[x] = np.sin(np.pi*x*h_x) *np.sinh(np.pi)
                    #print x,y
        
        
        # 2. Randwerte setzen
        for x in range(0,n_x+1):
            self.u_h[x ] = self.u_exact[x ]
        
    
        
        error = 0
        # 3. Loesen
        for i in range(0,self.i_iter):
            error = np.sqrt(np.sum(np.square(self.u_h-self.u_exact))/(float(n_x+1)))
            print('l_2 norm : ', error)
            print('iteration step : ', i)
            for x in range(1,n_x):
                self.u_h[x] =  0.5 * (self.u_h[x-1]+self.u_h[x+1 ]) 
                print(x-1)
                print(x+1)
                     
        
        # 4. l_2 norm
        
        error = 0
        for yx in range(0,n_x+1):
            error += np.square(self.u_exact[x]-self.u_h[x]) 

        error =np.sqrt(error/((n_x+1)))
                         
                        
    def plot(self):         
        fig, (ax1, ax2, ax3) = plt.subplots(1,3, figsize=(12.6929,8.26772))

        self.u_h = np.reshape(self.u_h, (self.n_y+1, self.n_x+1))
        self.u_exact = np.reshape(self.u_exact,(self.n_y+1, self.n_x+1))
        
        
        im1 = ax1.imshow(self.u_h)  
        ax1.set_title('$u_h$')
        fig.colorbar(im1, ax=ax1,fraction=0.046, pad=0.04)
        im2 = ax2.imshow(self.u_exact)
        fig.colorbar(im2, ax=ax2,fraction=0.046, pad=0.04)
        ax2.set_title('$u_{exact}$')
        im3 = ax3.imshow(self.u_h - self.u_exact)
        ax3.set_title('$u_h-u_{exact}$')
        fig.colorbar(im3, ax=ax3,fraction=0.046, pad=0.04)
        fig.tight_layout()   
        plt.show()            
                                                