#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 27 15:21:43 2017
This is to general the initial inputs of z_matrix, dz_matrix, Cl, Ch, Cb, Nl,
@author: qina
"""
from parameters import *
import numpy as np
from LEM_initial_landscape import ini_landscape, loadDEM
from F_Soil_Moisture_main_pool_parallel import vanGenuchten
from overland_Q import overland_setup
from F_save_netCDF import save2D_results
from CN_layers_depth import ini_z_dz
from F_root_fraction import root_fr
from CN_biogeochemical_pool_parallel import aveg_k
import os
from netCDF4 import Dataset
#------------------------------------------------------------------------
# LOAD INITIAL ELEVATION
#------------------------------------------------------------------------

def creating_ini_values(CCW,sin,litter_input,inputDatabase):
    #------------------------------------------------------------------------
    # INITIAL elevation, soil layer thickness
    #------------------------------------------------------------------------
    if CCW:
        X,Y,eta,mask_grid   = loadDEM()                               # loadDEM() or ini_landscape()
        mask_vector         = mask_grid.flatten()
        mask_ID = np.where(mask_vector==0)
    if sin:
        X,Y,eta             = ini_landscape()
        
    eta_vector_ini          = np.array(eta).flatten()                       # Convert into 1D array
    eta_vector              = eta_vector_ini*1.0    
    nrows, ncols            = eta.shape
    z_matrix_ini, dz_matrix_ini = ini_z_dz(nrows, ncols)
    Ch_all_ini  = 43000.0*np.exp(-2.5*z_matrix_ini) # 45000
    Cb_all_ini  = 180.0*np.exp(-2.5*z_matrix_ini) #320
    Cl_all_ini  = 6900.0*np.exp(-2.5*z_matrix_ini) #2800
    
    if inputDatabase:
        nc_f = "wd373.nc"  #test_wd500_noT
        nc_fid = Dataset(nc_f, 'r')
        
        soil_depth_change = np.array(nc_fid.variables['d_LE'])
        eta_vector          = eta_vector+soil_depth_change_d.flatten()+soil_depth_change_s.flatten()
        Cl_all = np.array(nc_fid.variables['Cl'])
        Cb_all = np.array(nc_fid.variables['Cb'])
        Ch_all = np.array(nc_fid.variables['Ch'])
        z_matrix            = np.array(nc_fid.variables['soil_depth'])
        dz_matrix           = np.zeros((z_matrix.shape))
        dz_matrix[0,:]      = z_matrix[0,:]*1.0
        for i in xrange(6): 
            dz_matrix[i+1,:]= z_matrix[i+1,:]*1.0-z_matrix[i,:]*1.0
        theta = np.array(nc_fid.variables['theta'])
    else:
        z_matrix        = z_matrix_ini*1.0
        dz_matrix       = dz_matrix_ini*1.0
        # np.save('z_matrix_ini',z_matrix_ini)
        # np.save('dz_matrix_ini',dz_matrix_ini)   
        #------------------------------------------------------------------------
        # INITIAL CARBON CONCENTRATION, CL, CB, CH, NL, HB, NH, N+, N-
        #------------------------------------------------------------------------

        Ch_all      = Ch_all_ini*1.0
        Cb_all      = Cb_all_ini*1.0
        Cl_all      = Cl_all_ini*1.0
        
    
        # np.save('Ch_all_ini',Ch_all_ini)
        # np.save('Cb_all_ini',Cb_all_ini)
        # np.save('Cl_all_ini',Cl_all_ini)
        
        #C_matrix_ini   = Ch_all_ini+Cb_all_ini+Cl_all_ini
        #C_matrix       = C_matrix_ini*1.0
        
    C1_vector   = np.array([Cl_all[0,:],Ch_all[0,:], Cb_all[0,:]])*1.0 #+Cb_all[0,:]+Ch_all[0,:]
    
    #CN_l_all    = CN_l_ini*np.ones((soil_layer_num,nrows*ncols)) #same as np.ones(z_matrix.shape)
    #N_plus_all  = 0.002*np.ones((soil_layer_num,nrows*ncols)) #0.002*np.exp(-1.0*z_matrix)
    #N_minus_all = 0.6*np.ones((soil_layer_num,nrows*ncols))  #0.6*np.exp(-1.0*z_matrix)
    
    #------------------------------------------------------------------------
    # INITIAL water deoth (wd) and water surface elevation (we)
    #------------------------------------------------------------------------
    wd_vector, we_vector    = overland_setup(eta_vector)                          # vector means reorganize 2D to 1D
    mann_bs                 = manning_bs * np.ones(nrows*ncols)          # load the Manning's coeff
    mann_vg                 = manning_vg * np.ones(nrows*ncols)
    #------------------------------------------------------------------------
    # INITIAL soil mositure and soil pressure head
    #------------------------------------------------------------------------
    Psi_ini_all  = -2.0*z_matrix-2.0 #-1.0*np.ones((soil_layer_num,nrows*ncols)) #np.linspace(2.0,2.0,nz)
    Psi_n_all    = Psi_ini_all*1.0
    
    C_ini_all,K_ini_all,theta_ini_all = vanGenuchten(Psi_ini_all,alpha, theta_S, theta_R, n, m, Ksat)
     
    theta_n_all    = theta_ini_all*1.0
    
#    theta_mean_1col = np.array([ 0.30887538, 0.31472939, 0.31682785, 0.31680671, 0.31508827, 0.31146861, 0.30576922])
    theta_mean_1col = np.array([0.35480344, 0.36059945, 0.36269936, 0.36273304, 0.3611157 ,0.3574786 , 0.35072803])
    #theta_mean_1col = np.array([0.255932  ,  0.26280422,  0.26647157,  0.26852458,  0.26961957, \
    #                            0.2700011 ,  0.26972128,  0.26868466,  0.26671475,  0.2637314]) # 10 layers
       
    #------------------------------------------------------------------------
    # INITIAL water deoth (wd) and water surface elevation (we)
    #------------------------------------------------------------------------
    wd_vector, we_vector = overland_setup(eta_vector)                          # vector means reorganize 2D to 1D
    

    #------------------------------------------------------------------------
    # Kl, Kh, Kd, and ADD
    #------------------------------------------------------------------------
    
    Kl          = np.zeros((soil_layer_num))
    Kh          = np.zeros((soil_layer_num))
    Kd          = np.zeros((soil_layer_num))
    
    Zr=z_matrix[:,1]  
    dZr=dz_matrix[:,1]
    Frooti = root_fr(Zr,dZr)
    
    litter_mean = np.mean(litter_input,1)
    Addb = Frooti*(Add_b+litter_mean[1]) #+litter_mean[1]
    ADD = Addb*1.0 
    ADD[0] = Addb[0]+ADD_s+litter_mean[0]  #+litter_mean[0]

    rh_2col = np.column_stack((rh_ini*np.ones(soil_layer_num), np.ones(soil_layer_num)*CN_h_ini/CN_l_ini))
    rh=np.min(rh_2col,axis =1)
        
    for ll in xrange(0,soil_layer_num):
        Kl[ll],Kh[ll],Kd[ll] = aveg_k(theta_mean_1col[ll]/poros,np.mean(Cl_all[ll,:]), np.mean(Ch_all[ll,:]), np.mean(Cb_all[ll,:]),ADD[ll],rh[ll] )  
    
    Kl = np.ones((soil_layer_num,nrows*ncols))*Kl[:,None]
    Kh = np.ones((soil_layer_num,nrows*ncols))*Kh[:,None]
    Kd = np.ones((soil_layer_num,nrows*ncols))*Kd[:,None]

    #------------------------------------------------------------------------
    #     SAVE VARIABLES 
    #------------------------------------------------------------------------
    water_depth_saved    = np.zeros(nrows*ncols)
#    dh_STs               = np.zeros(nrows*ncols)
#    dh_STd               = np.zeros(nrows*ncols)
#    dh_STr               = np.zeros(nrows*ncols)
    dh_STsave            = np.zeros(nrows*ncols)
#    qC_save_d            = np.zeros((3, nrows*ncols))
#    qC_save_s            = np.zeros((3, nrows*ncols))
#    qC_save_r            = np.zeros((3, nrows*ncols))
    dC_BG                = np.zeros((soil_layer_num,nrows*ncols))
    qC_LE              = np.zeros((nrows*ncols))
    #------------------------------------------------------------------------
    #     SAVE DATA 
    #------------------------------------------------------------------------
    folder =  'Output_data' #os.getcwd() #
    
    if not os.path.exists(folder):
        os.makedirs(folder)

    save2D_results(dh_STsave.reshape((nrows,ncols)), 
                   qC_LE.reshape((nrows,ncols)), 
                   Cl_all+Ch_all+Cb_all,
                   dC_BG,
                   theta_n_all,
                   z_matrix,
                    'd_LE', 'qC_LE', 'Clhb','dC_BG','theta','z_matrix', 
                   folder, "wd", 0) #"d_soildepth_s" 'qC_s' 
    
    return eta_vector,eta_vector_ini, nrows, ncols, z_matrix,z_matrix_ini,  dz_matrix,dz_matrix_ini, theta_n_all, Psi_n_all, wd_vector, we_vector, Cl_all, Ch_all, Cb_all,Cl_all_ini, Ch_all_ini, Cb_all_ini, C1_vector, Kl, Kh, Kd, Frooti, mann_bs, mann_vg, mask_ID

def forcing():
        
    #------------------------------------------------------------------------
    # ADD ATMOSPHERIC FORCING
    #------------------------------------------------------------------------
    #...................Rainfall
    rain = np.load('rainfall_100yr3_daily.npy') # mm/(day)
    rain = rain/1000.0
    
    #...................litter (2 * (365*2))
    # the first row is above ground, the second row is below ground
    litter_input = (np.load('litter_input.npy'))/dt
    #...................Rainfall
    #RF  = .02*1e-3           #Rainfall rate in unit [m/hr]
    #PPT = RF * dt                 # Rainfall depth per time step. [m]

    
    
    return (rain, litter_input)
    