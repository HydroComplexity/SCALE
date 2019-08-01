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
import multiprocessing
from multiprocessing import Process, Lock, Value, Array


from CN_layers_depth import ini_z_dz, z_dz_update,C_interpolate,C_interpolate_simple
import timeit
import matplotlib.pyplot as plt


def f(cell, Zr, Cmmm):
    Cl         = 40000.0*np.exp(-1.0*Zr) # 18000
    ClayPerct = 70.0 #
    SiltPerct = 17.8#
    Cv = Cl*1e-6

    aConst = 1.69-0.013*theta_R+0.0079*ClayPerct-0.0007*SiltPerct+ 0.00014*Zr*100.0
    Cmm = np.zeros(soil_layer_num)
    for lanum in xrange(0,soil_layer_num):
        x = Symbol("x") 
        ff = (aConst[lanum]*x-0.2*(1000.0*x)**0.5*x-Cv[lanum])
        s= nsolve(ff, x,0.0)
        Cmm[lanum] = re(s)

    for i in range(soil_layer_num):
        Cmmm[cell*soil_layer_num+i] = Cmm[i]


def output(z_matrix, grids):

    #    save_2d = np.zeros((100,100))
    arr = Array('d', range(soil_layer_num*grids))

    for num in range(grids):
#        arr = Array('d', range(100))
        p = Process(target=f, args=(num, z_matrix[:,num],arr))
        p.start()

    p.join()
#    p.terminate()
#        pl = p
#    for pl:
#        p.join()
    T = np.asarray(arr[:])   
    T = T.reshape((grids,soil_layer_num))
    T = np.transpose(T)
    return T
 

#nrows = 2
#ncols = 2
#z_matrix, dz_matrix = ini_z_dz(nrows, ncols)
##z_matrix[:,1]=z_matrix[:,1]*2.0
#Cl         = 40000.0*np.exp(-1.0*z_matrix) # 18000
#ClayPerct = 70.0 #
#SiltPerct = 17.8#
#Cv = Cl*1e-6
#
#
#start               = timeit.default_timer()   
#
#Cm = output(nrows,ncols)
#
#stop = timeit.default_timer()-start
#
#print stop



#BD = 1.69-0.013*theta_R+0.0079*ClayPerct-0.0007*SiltPerct+ 0.00014*z_matrix*100.0-0.2*(1000*Cm)**0.5
#Cl2 = BD*Cm
#print (Cl2-Cv) # should be zero

#fig, ax = plt.subplots(figsize=(4,6)) 
#ax.plot(Cm[:,0]*100,z_matrix[:,0],'o-', label='Cm')
##ax.plot(Cv[:,0]*1000.0,z_matrix[:,0], 'o-',label='Cv')
#ax.set_xlabel(r'C [$KgC m^{-3}$]')
#ax.set_ylabel(r'Soil Depth [m]')
#ax.invert_yaxis()
#ax.legend()
#plt.show()
#    



