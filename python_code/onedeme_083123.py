#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 20 15:32:52 2023

@author: martingaric
"""

##### Code to study coevolution of virus on an antigenic space
##### Works with n identical populations

import numpy as np
import matplotlib.pyplot as plt
plt.style.use('ggplot')
import utils as ut
import matplotlib.animation as animation 
from tqdm import tqdm


class Population():
    def __init__(self, alpha, beta, gamma, r0, D, M, Nh, N0, Nc, Nstoch, nt, nx, burn_time, dx, dt):
        self.alpha = alpha;                     # Virulence
        self.beta = beta;                       # Viral transmission rate
        self.gamma = gamma;                     # Recovery rate
        self.r0 = r0;                           # Receptor cross reactivity
        self.D = D;                             # Mutation diffusion coef.
        self.M = M;                             # Protections per hosts
        self.Nh = Nh;                           # Number of hosts
        self.N0 = N0;                           # Initial number of virus
        self.Nc = Nc;
        self.Nstoch = Nstoch;                   # Stochastic behaviour threshold
        self.R0 = beta/(alpha + gamma);         # Reproductive ratio
        self.Fmax = beta - alpha - gamma;       # Maximal growth rate
        self.rho = r0 / (self.R0**(1/M)-1) ;    # Serves for stochastic calculation
        
        self.nt = nt;                           # Number of time steps
        self.nx = nx;                           # Size of antigenic space
        self.burn_time = burn_time;             # Burn time before reaching steady state
        self.dx = dx;                           # Lattice spacing
        self.dt = dt;                           # Time interval between time steps
        self.x_end = nx * dx;                   # Last point in x
        self.t_end = nt * dt;                   # Last time point

    def Check_Parameters(self):
        prec = dx/r0
        if prec > 0.1:
            print('WARNING: dx is not small compared to r0. dx/r0=', prec)
        CFL = D * dt / dx**2
        if CFL > 0.1:
            print('WARNING: CFL is not small, choose a smaller dt. CFL=', CFL)
        if prec <=0.1 and CFL<=0.1:
            print("Current parameters will not create numerical errors.\n")
        return
        
    def Init_Arrays(self):
        x_end = self.x_end;
        t_end = self.t_end;
        nx = self.nx;
        nt = self.nt;
        
        self.X = np.linspace(0, x_end, nx, endpoint=False);                # Space lattice array
        self.T = np.linspace(0, t_end, nt, endpoint=False);                # Time array
        self.N = np.zeros((nt, nx));                                      # Viral population array
        self.Ntot = np.zeros(nt);                                          # Total viral population array
        self.C = np.zeros((nt, nx));                                      # Receptor coverage array
        self.h = np.zeros((nt, nx));                                      # Density of immune receptors array
        
        ## Not sure yet if I put this in the initialization step...
        self.N[0, :] = ut.skewed_gaussian(x=self.X, sk=4, sigma=1,mu=0)  # Creating correctly normalized N
        Kn = self.Nh/100 * 1 / (self.N[0, :].sum()*self.dx)
        self.N[0, :] = self.N[0, :] * Kn
        
        self.h[0, :] = np.exp(-np.abs(self.X)/self.rho)                  # Creating correctly normalized h 
        Kh = 1 / (self.h[0, :].sum()*self.dx)
        self.h[0, :] = self.h[0, :] * Kh
        return

    def Evolution(self, i): #i is the time step
        # N:
        Dterm = self.D * (np.roll(self.N[i,:],1) - 2*self.N[i,:] + np.roll(self.N[i,:],-1)) * self.dt/(self.dx**2)
        self.C[i,:] = ut.cxtfun(self.h[i,:], self.r0, self.dx, self.nx, self.Nh)
        Ski = ut.S(self.C[i],self.M)
        Fki = ut.F(self.beta, Ski, self.alpha, self.gamma)
        Fterm = self.N[i,:] * Fki * self.dt
        self.N[i+1,:] = self.N[i, :] + Dterm + Fterm*np.heaviside(self.N[i,:]-self.Nc, 1)
        
        # Stochastic step
        p = 1 - np.exp(-2*self.D*self.dt/(self.dx**2))
        pos_small = np.where(self.N[i,:] < self.Nstoch)
        Ndt2 = self.N[i,:] + (1/self.dx) * (np.random.binomial((np.roll(self.N[i,:],1)*self.dx).astype(int), p/2) + \
                                    np.random.binomial((np.roll(self.N[i,:],-1)*self.dx).astype(int), p/2) - \
                                    2*np.random.binomial((self.N[i,:]*self.dx).astype(int), p/2))
        Nbar = (1 + Fki*self.dt) * Ndt2 * self.dx
        Nbar[np.where(Nbar<0)[0]]=0
        Poiss_evol = np.random.poisson(lam=Nbar)/self.dx
        self.N[i+1, pos_small[0]] = Poiss_evol[pos_small[0]]

        # h:
        for j in range(self.nx):
            self.h[i+1, j] = self.h[i,j] + 1/(self.Nh*self.M) * (self.N[i,j] - self.Ntot[i]*self.h[i,j]) * self.dt
    
    def Ntot_evo(self, i):
        self.Ntot[i] = self.N[i, :].sum() * self.dx 

        
        
#%%    

# Constants:
alpha = 0.3;
beta = 2;
gamma = 0.7;
r0 = 2;
D = 3e-4;
M = 5;
Nh = 10**6.;
N0 = 10**10.;
Npops = 2;
Nc = 1;
Nstoch = 1e6;

nt = 10000;
burn_time = 1000;
dx = .2;
dt = .4;
nx = 750;

jump_probas = [1e-8]

pop1 = Population(alpha, beta, gamma, r0, D, M, Nh, N0, Nc, Nstoch, nt, nx, burn_time, dx, dt)
pop1.Check_Parameters()
pop1.Init_Arrays()


for i in tqdm(range(pop1.nt-1)):
    pop1.Ntot_evo(i)
    
    # Growth step
    pop1.Evolution(i)
    
    # Jump step
    # ut.Jump_step(populations, jump_probas, i)


















# END #