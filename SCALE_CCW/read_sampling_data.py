#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 20 22:31:30 2017

@author: qina
"""

import numpy as np
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt

from parameters import *
from CN_biogeochemical_pool_parallel import aveg_k

ClayPerct = 70.0 #
SiltPerct = 17.8#


# core #12
Cw12 = 10.0*np.array([2.4,2.399,2.254,2.184,2.244,1.751,1.369,1.23,1.142,1.088,0.976,0.874,0.778,0.708,0.65,0.564,0.553,0.619,0.595,0.624,0.591,0.636])
z12  = 0.5-0.01*np.array([2.5,7.5,12.5,17.5,22.5,27.5,32.5,37.5,42.5,47.5,52,56,60.5,65.5,70.5,79.5,84.5,94.5,99.5,104.5,109.5,114])
BD12 = 1.69-0.013*theta_R+0.0079*ClayPerct-0.0007*SiltPerct+ 0.00014*z12 -0.198*(Cw12)**0.5
Cv12 = Cw12*BD12
# core #11
Cw11 = 10.0*np.array([2.334813619,2.396640917,2.295153614,1.733681785,1.355335877,1.511039998,3.658472134,2.970083409,2.485586733,2.090665552,1.879475761,1.640683289,1.434020341,1.320679241,1.203702052,0.960498734,0.871943156,0.715620704,0.701808517,0.660839588,0.652810187,0.603719973,0.623887767])
z11  = 0.5-0.01*np.array([2.5,7,11.5,21.5,26,30,34.5,39.5,44,48.5,53.5,58.5,63.5,67.5,71.5,76.5,83.5,88.5,91.5,96.5,101.5,106.5,113.5])
BD11 = 1.69-0.013*theta_R+0.0079*ClayPerct-0.0007*SiltPerct+ 0.00014*z11 -0.198*(Cw11)**0.5
Cv11 = Cw11*BD11
    
# core #8
Cw8  = 10.0*np.array([2.084208623,2.038287566,1.887615121,1.565343997,1.551546306,1.685875551,1.933095599,2.313340242,2.423447071,3.378953122,3.048604062,2.737810469,2.400531524,2.267791104,1.960423104,1.532604451,1.275584329,1.109133087,1.040599785,0.874593428,0.762099021,0.692337579])
z8   = 0.5-0.01*np.array([2.5,7.5,12.5,24,30.5,35.5,40.5,45.5,50.5,55,60,65.5,70.3,75.5,80.5,85.5,90.5,95.5,101,106.5,111,115.5])
BD8  = 1.69-0.013*theta_R+0.0079*ClayPerct-0.0007*SiltPerct+ 0.00014*z8 -0.198*(Cw8)**0.5
Cv8  = Cw8*BD8
    
## core #7
Cw7  = 10.0*np.array([1.964268892,1.917928939,1.874821409,1.288167598,1.516715735,1.381454461,1.573562723,1.641422157,1.784593494,2.122220705,2.77166079,3.284496077,3.611431957,3.228559866,2.756112033,2.310921873,1.914294135,1.616272527,1.378807909,1.200271594,1.108539361,1.007035187])
z7   = 0.5-0.01*np.array([2.5,7.5,12.5,17.7,22.5,27.5,32.5,38,44,49.5,54,58,62.5,67.5,72.5,77.5,82.5,87.5,92.5,97.5,102.5,107.5])
BD7  = 1.69-0.013*theta_R+0.0079*ClayPerct-0.0007*SiltPerct+ 0.00014*z7 -0.198*(Cw7)**0.5
Cv7  = Cw7*BD7
# core #6
Cw6   = 10.0*np.array([1.95,1.78,1.81,1.52,1.52,1.6,1.79,2.36,2.67,2.81,2.62,2.36,1.96,1.63,1.4,1.14,1.13,1.11])  # convert from 100*g/g --> g/kg
z6    = 0.5-0.01*np.array([3,9,15,20.5,25,36,41.5,47,53,58.5,64,69.5,74.5,79.5,84.5,94.5,99.5,104.5])
BD6   = 1.69-0.013*theta_R+0.0079*ClayPerct-0.0007*SiltPerct+ 0.00014*z6 -0.198*(Cw6)**0.5
Cv6   = Cw6*BD6

#cemeterey   #Cw_c[0] = 9.0103, z_c[0] = 1.0,
Cw_c = np.array([4.4949,3.1273,2.4737,1.8589,1.6998,1.614,1.5702,1.4874,1.4647,1.2802,1.0801,0.9341,0.7594,0.698,0.6107,0.4834,0.36,0.3358,0.315,0.2591,0.2627,0.2526,0.2349,0.2418])*10.0
z_c  = 0.5-0.01*np.array([4.5,9.5,14.5,19.5,24.5,29.5,34.5,39.5,44.5,49.5,54.5,59.5,64.5,69.5,74.5,79.5,84.5,89.5,94.5,99.5,104.5,109.5,114.5,119.5])
BDc  = 1.69-0.013*theta_R+0.0079*ClayPerct-0.0007*SiltPerct+ 0.00014*z_c -0.198*(Cw_c)**0.5
Cvc  = Cw_c*BDc

#========================================================================
# plot initial values with sampling data
#========================================================================

z_matrix_ini  =  np.array([0.05509642, 0.11937557, 0.19651056, 0.29292929, 0.4214876 ,
       0.61432507, 1.  ])
    
Ch_all_ini  = 43000.0*np.exp(-2.5*z_matrix_ini) # 45000
Cb_all_ini  = 3260.0*np.exp(-2.5*z_matrix_ini) #320
Cl_all_ini  = 7000.0*np.exp(-2.5*z_matrix_ini) #2800
C_all_ini   = (Ch_all_ini+Cb_all_ini+Cl_all_ini)/1000.0

fig, ax = plt.subplots(figsize=(2.5,5), dpi=150) 
ax.plot(Cv6,z6,'o--',label='core 6')
ax.plot(Cv7,z7+0.06,'o--',label='core 7')
ax.plot(Cv8,z8+0.02,'o--',label='core 8')
ax.plot(Cv11,z11-0.2,'o--',label='core 11')
ax.plot(Cv12,z12-0.52,'o--',label='core 12')
ax.plot(Cvc,z_c-0.5,'grey', linestyle='-.',)
ax.plot(C_all_ini,-z_matrix_ini,'k')

#ax.plot(C/1000.0,z7_final,'ok-',label='final results')
ax.legend(loc='lower right')
ax.set_xlabel(r'[$KgC \  m^{-3}$]')
ax.set_ylabel(r'Soil Depth [m]')
fig.tight_layout()
#fig.savefig('sampling_sites_withInitial.png') 
#ax.invert_yaxis()
           
#========================================================================
# converting back and forth to get the C value on Z_matrix_ini locations. 
#========================================================================
#
#z7_final    = 0.75 -np.array([0.        ,  0.11466046,  0.25225301,  0.42424369,  0.65356461,
#        0.99754598,  1.68550873])     #([ 0.05509642,0.11937557,0.19651056,0.29292929,0.4214876,0.61432507,1.])
##z10_need   = -np.array([ 0.04032524,0.08489524,0.13470876,0.1911641,0.25630487,0.33328941,0.42738163,0.54835734,0.71772334,1.])
#z = z6
#Cv = Cv6
#
#Cv7_final  = interp1d(z, Cv,fill_value='extrapolate')(z7_final)  #,fill_value='extrapolate'
##Cv10_need = interp1d(z, Cv)(z10_need)  #,fill_value='extrapolate'
#Cvc7_final = interp1d(z_c, Cvc,fill_value='extrapolate')(z7_final)  #,fill_value='extrapolate'
#C         = Cv7_final*1000.0
#
##============plot
## core 6
#fig, ax = plt.subplots(figsize=(4,6)) 
#ax.plot(Cv6,z6,'o--',label='core 6')
#ax.plot(Cv7_final, z7_final,'o--',label='7 final')
##ax.plot(Cv10_need,z10_need,'o--',label='10 layers')
##ax.plot(Cv_,z_ero[soil_layer_num-1] - z_ero, 'o-',label='ero')
##ax.plot(C_depo/1000.0,z_depo[soil_layer_num-1] - z_depo,'o-', label='depo')
#ax.legend(loc='lower right')
#ax.set_xlabel(r'[$KgC m^{-3}$]')
#ax.set_ylabel(r'Soil Depth [m]')
##ax.invert_yaxis()
#
## core cemetery
#fig, ax = plt.subplots(figsize=(4,6)) 
#ax.plot(Cvc,z_c,'o--',label='cemetery')
#ax.plot(Cvc7_final, z7_final,'o--',label='7 layers')
##ax.plot(Cv10_need,z10_need,'o--',label='10 layers')
##ax.plot(Cv_,z_ero[soil_layer_num-1] - z_ero, 'o-',label='ero')
##ax.plot(C_depo/1000.0,z_depo[soil_layer_num-1] - z_depo,'o-', label='depo')
#ax.legend(loc='lower right')
#ax.set_xlabel(r'[$KgC m^{-3}$]')
#ax.set_ylabel(r'Soil Depth [m]')
##ax.invert_yaxis()
#
##========================================================================
## find kl, kd, kh
##========================================================================
#### case2: ZERO water ponding is preserved all the time
#
#ADD = np.array([ 0.04725739,0.00428391,0.00672107,0.00676813,0.00491948,0.00283496,0.00125288])
#rh  = np.array([ 0.25,  0.25,  0.25,  0.25,  0.25,  0.25,  0.25])
#
#theta_mean_1col = np.array([ 0.25585338,  0.2628716 ,  0.26605831,  0.26727228,  0.26710455,
#    0.26546586,  0.26187034]) # 7 layers
#
## the exponential decreasing value
#Cb = np.array([ 276.20348168,  250.8158774 ,  223.41178417,  193.32805935,
#        159.42140939,  119.37799307,   66.93904804])
#
#Cl = np.array([ 5299.80580334,  4969.85678511,  4600.91886862,  4178.01939297,
#        3673.99268009,  3029.63311754,  2060.12487056])
#
##C[0:4] = C[3]
#C[0] = 20.0*1000.0
#C[1] = 21.0*1000.0
#C[2] = 23.0*1000.0
#C[3] = 27.0*1000.0
#C[4] = 33.0*1000.0
#C[5] = 17.0*1000.0
#C[6] = 10.0*1000.0
##C[4] = 30.0*1000.0
#Ch = C-Cb-Cl
#
#
## the final values of the 100 run with SWnLEnSCT
#
##Cl = 1000.0*np.array([5.16636439,  5.62313561,  5.80703598,  5.75180061,  5.28373701,        4.35678292,  2.91836437])
##Ch = 1000.0*np.array([20.25728709,  23.1711233 ,  25.3071049 ,  24.65176094,        19.68011562,  11.65604157,   3.74754419])
##Cb = 1000.0*np.array([0.48620445,  0.52109774,  0.54624091,  0.53861508,  0.47605081,        0.35886057,  0.19775208])
#fig, ax = plt.subplots(figsize=(4,6)) 
#ax.plot(Cw6,z6,'o--',label='core 6')
#ax.plot(Cw7,z7+0.06,'o--',label='core 7')
#ax.plot(Cw8,z8+0.02,'o--',label='core 8')
#ax.plot(Cw11,z11-0.2,'o--',label='core 11')
#ax.plot(Cw12,z12-0.52,'o--',label='core 12')
##ax.plot(Cvc,z_c-0.5,'grey',label='Cemetery')
##ax.plot(C/1000.0,z7_final,'ok-',label='final results')
#ax.legend(loc='lower right')
#ax.set_xlabel(r'[$KgC /Kg$]')
#ax.set_ylabel(r'Soil Depth [m]')
#
#fig, ax = plt.subplots(figsize=(2.5,5), dpi=150) 
#ax.plot(Cv6,z6,'o--',label='core 6')
#ax.plot(Cv7,z7+0.06,'o--',label='core 7')
#ax.plot(Cv8,z8+0.02,'o--',label='core 8')
#ax.plot(Cv11,z11-0.2,'o--',label='core 11')
#ax.plot(Cv12,z12-0.52,'o--',label='core 12')
#ax.plot(Cvc,z_c-0.5,'grey', linestyle='-.',)
##ax.plot(C/1000.0,z7_final,'ok-',label='final results')
#ax.legend(loc='lower right')
#ax.set_xlabel(r'[$KgC \  m^{-3}$]')
#ax.set_ylabel(r'Soil Depth [m]')
#fig.tight_layout()
##fig.savefig('sampling_sites_withInitial.png') 
##ax.invert_yaxis()
#
#Kl          = np.zeros((soil_layer_num))
#Kh          = np.zeros((soil_layer_num))
#Kd          = np.zeros((soil_layer_num))
#
#
#for ll in xrange(0,soil_layer_num):
#    Kl[ll],Kh[ll],Kd[ll] = aveg_k(theta_mean_1col[ll]/poros,np.mean(Cl[ll]), np.mean(Ch[ll]), np.mean(Cb[ll]),ADD[ll],rh[ll] )  
