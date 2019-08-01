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
from sympy import Symbol,solve,nsolve,re
import scipy.optimize


from CN_layers_depth import ini_z_dz, z_dz_update,C_interpolate,C_interpolate_simple
import timeit
import matplotlib.pyplot as plt

nrows = 1
ncols = 1
z_matrix, dz_matrix = ini_z_dz(nrows, ncols)

Ch_all_ini  = 36000.0*np.exp(-3.0*z_matrix) # 18000
Cb_all_ini  = 300.0*np.exp(-1.5*z_matrix) #300
Cl_all_ini  = 5600.0*np.exp(-1.0*z_matrix) #2800


Cl         =Ch_all_ini+Cb_all_ini+Cl_all_ini #np.array([[ 36919.1187679, 36605.7340711 , 34042.24123639, 31006.29236706, 27306.71422354,22542.72175682, 15350.17572777]])#39000.0*np.exp(-1.0*z_matrix) # 18000
#Cl = np.transpose(Cl)
Cm = np.zeros((soil_layer_num,ncols*nrows))
ClayPerct = 70.0 #
SiltPerct = 17.8#
Cv = Cl*1e-6

start               = timeit.default_timer()   

for cell in xrange(0,(ncols)*(nrows)):
    Zr = z_matrix[:,cell]
    aConst = 1.69-0.013*theta_R+0.0079*ClayPerct-0.0007*SiltPerct+ 0.00014*Zr*100.0

    for lanum in xrange(0,soil_layer_num):
        x = Symbol("x") 
        ff = (aConst[lanum]*x-0.2*(1000.0*x)**0.5*x-Cv[lanum,cell])
        s= nsolve(ff, x,0.0)
        Cm[lanum,cell] = re(s)

stop = timeit.default_timer()-start
print stop

#BD = 1.69-0.013*theta_R+0.0079*ClayPerct-0.0007*SiltPerct+ 0.00014*z_matrix*100.0-0.2*(1000*Cm)**0.5
#Cl2 = BD*Cm
#print (Cl2-Cv) # should be zero

fig, ax = plt.subplots(figsize=(4,6)) 
#ax.plot(Cm[:,0]*100.0,z_matrix[:,0],'o-', label='Cm')
ax.plot(Cv[:,0]*1000.0,z_matrix[:,0], 'o-',label='Cv')
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