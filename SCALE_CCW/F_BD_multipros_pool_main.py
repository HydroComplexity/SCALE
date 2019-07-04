#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 28 15:47:58 2017
Python MPI Practise
@author: qina
"""
#http://docs.sympy.org/latest/modules/solvers/solvers.html
#http://docs.sympy.org/latest/modules/solvers/solvers.html

from parameters import *

import numpy as np
#from sympy import Symbol,solve,nsolve,re
#import scipy.optimize
#import multiprocessing

from CN_layers_depth import ini_z_dz
import timeit
import matplotlib.pyplot as plt
import F_BD_multipros_pool_fn 
#

nrows = 2
ncols = 2
z_matrix, dz_matrix = ini_z_dz(nrows, ncols)
z_matrix[:,1]=z_matrix[:,1]*2.0
start               = timeit.default_timer()   

Cm = F_BD_multipros_pool_fn.output(z_matrix,dz_matrix, nrows,ncols)

stop = timeit.default_timer()-start


fig, ax = plt.subplots(figsize=(4,6)) 
ax.plot(Cm[:,0]*100,z_matrix[:,0],'o-', label='Cm')
#ax.plot(Cv[:,0]*1000.0,z_matrix[:,0], 'o-',label='Cv')
ax.set_xlabel(r'C [$KgC m^{-3}$]')
ax.set_ylabel(r'Soil Depth [m]')
ax.invert_yaxis()
ax.legend()
plt.show()

    
    #fig, ax = plt.subplots(figsize=(4,6)) 
    #ax.plot(BD[:,0]/1000.0,z_matrix[:,0],'o-', label='Cm')
    #ax.set_xlabel(r'C [$KgC m^{-3}$]')
    #ax.set_ylabel(r'Soil Depth [m]')
    #ax.invert_yaxis()
    #ax.legend()