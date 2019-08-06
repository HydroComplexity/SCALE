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
import gdal
from netCDF4 import Dataset
from functools import reduce
#------------------------------------------------------------------------
# LOAD forcing data
#------------------------------------------------------------------------
def forcing(unt, tre):
        
    #------------------------------------------------------------------------
    # ADD ATMOSPHERIC FORCING
    #------------------------------------------------------------------------
    #...................Rainfall
    rain = np.load('rainfall_50yr_daily.npy') # mm/(day)
    rain = rain/1000.0
    
    #...................litter (2 * (365*2))
    # the first row is above ground, the second row is below ground
    litter_input_unt  = (np.load('litter_input_untreated.npy'))/dt 
    litter_input_tre  = (np.load('litter_input_treated.npy'))/dt
    litter_input_treC = (np.load('litter_input_treatedCrop.npy'))/dt 
    if unt:
        litter_input  = litter_input_unt*1.0
        litter_inputC = litter_input_unt*1.0
    if tre:
        litter_input  = litter_input_tre*1.0
        litter_inputC = litter_input_treC*1.0
    #...................Rainfall
    #RF  = .02*1e-3           #Rainfall rate in unit [m/hr]
    #PPT = RF * dt                 # Rainfall depth per time step. [m]

    return (rain,litter_input,litter_inputC)

def creating_ini_values(CLP,sin,litter_input,litter_inputC,inputDatabase,tre,unt):
    #------------------------------------------------------------------------
    # INITIAL elevation, soil layer thickness, manning's coeff
    #------------------------------------------------------------------------
    if CLP:
        X,Y,eta,mann_vg         = loadDEM(tre,unt)                               # loadDEM() or ini_landscape()
        mann_vg                 = mann_vg.flatten()
    if sin:
        X,Y,eta                 = ini_landscape()
        mann_vg                 = manning_vg * np.ones(nrows*ncols)          # load the Manning's coeff
    
    ind_veg                     = np.where((mann_vg.flatten())<15)[0]
    eta_vector_ini              = np.array(eta).flatten()                       # Convert into 1D array
    eta_vector                  = eta_vector_ini*1.0    
    ind_dem                     = np.where(eta_vector>0.0)[0] # the ind_in is the index inside of the watershed boundary but within the big rectangula
    eta_vector[eta_vector<0.0]  = np.max(eta_vector)
#    eta_vector[eta_vector<0.0] = np.nan
    nrows, ncols                = eta.shape
    z_matrix_ini, dz_matrix_ini = ini_z_dz(nrows, ncols)
    mann_bs                     = manning_bs * np.ones(nrows*ncols)
    
    if sin:
        Ch_all_ini              = 43000.0*np.exp(-2.5*z_matrix_ini) # 45000
        Cl_all_ini              = 7000.0*np.exp(-2.5*z_matrix_ini) #2800
        Cb_all_ini              = 3260.0*np.exp(-2.5*z_matrix_ini) #320
        ind_out                 = np.where(eta_vector<0.0)[0]
    if CLP:
        current_dir     = os.getcwd()
        os.chdir(current_dir+"/"+'/Topo_SOCini_Landcover/')
        if tre:
            Manndata      = gdal.Open("landuse_tre_2m_clip.tif")
            Mann_vg       = Manndata.ReadAsArray()
            ind_crop      = np.where(Mann_vg.flatten()==7)[0]
            Clhb_ini      = np.load('a_surface_GLC.npy')
            Clhb_ini      = np.array(Clhb_ini).flatten()*1.0
            Clhb_profiles = (Clhb_ini*np.exp(-b_f_soc*z_matrix_ini) + c_f_soc)*1000.0
            Clhb_profiles[:,ind_crop] = (Clhb_ini[ind_crop]*np.exp(-b_c_soc*z_matrix_ini[:,ind_crop]) + c_c_soc)*1000.0
        if unt:  
            ind_crop      = np.array([0,1])
            Clhb_ini      = np.load('a_surface_Reference.npy')
            Clhb_ini      = np.array(Clhb_ini).flatten()*1.0 
            Clhb_profiles = (Clhb_ini*np.exp(-b_f_soc*z_matrix_ini) + c_f_soc)*1000.0
        
        ind_aSOC          = np.where(Clhb_ini>0.0)[0]
        os.chdir(current_dir)

        Cl_all_ini    = 0.1235*Clhb_profiles
        Ch_all_ini    = 0.8752*Clhb_profiles
        Cb_all_ini    = 0.0104*Clhb_profiles
        
        ind_in = reduce(np.intersect1d, (ind_dem, ind_veg, ind_aSOC))
    
    if inputDatabase:
        nc_f = "test_wd500_noT.nc"
        nc_fid = Dataset(nc_f, 'r')
        
        soil_depth_change_d = np.array(nc_fid.variables['d_soildepth_d'])
        soil_depth_change_s = np.array(nc_fid.variables['d_soildepth_s'])
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
        np.save('z_matrix_ini',z_matrix_ini)
        np.save('dz_matrix_ini',dz_matrix_ini)   
        #------------------------------------------------------------------------
        # INITIAL CARBON CONCENTRATION, CL, CB, CH, NL, HB, NH, N+, N-
        #------------------------------------------------------------------------

        Ch_all      = Ch_all_ini*1.0
        Cb_all      = Cb_all_ini*1.0
        Cl_all      = Cl_all_ini*1.0
        
#        np.save('Ch_all_ini',Ch_all_ini)
#        np.save('Cb_all_ini',Cb_all_ini)
#        np.save('Cl_all_ini',Cl_all_ini)
        
        
    C1_vector   = np.array([Cl_all[0,:],Ch_all[0,:], Cb_all[0,:]])*1.0 #+Cb_all[0,:]+Ch_all[0,:]
    
    #CN_l_all    = CN_l_ini*np.ones((soil_layer_num,nrows*ncols)) #same as np.ones(z_matrix.shape)
    #N_plus_all  = 0.002*np.ones((soil_layer_num,nrows*ncols)) #0.002*np.exp(-1.0*z_matrix)
    #N_minus_all = 0.6*np.ones((soil_layer_num,nrows*ncols))  #0.6*np.exp(-1.0*z_matrix)
    
    #------------------------------------------------------------------------
    # INITIAL water deoth (wd) and water surface elevation (we)
    #------------------------------------------------------------------------
    wd_vector, we_vector    = overland_setup(eta_vector)                          # vector means reorganize 2D to 1D

    #------------------------------------------------------------------------
    # INITIAL soil mositure and soil pressure head
    #------------------------------------------------------------------------
    Psi_ini_all  = -2.0*z_matrix-2.0 #-1.0*np.ones((soil_layer_num,nrows*ncols)) #np.linspace(2.0,2.0,nz)
    Psi_n_all    = Psi_ini_all*1.0
    
    C_ini_all,K_ini_all,theta_ini_all = vanGenuchten(Psi_ini_all,alpha, theta_S, theta_R, n, m, Ksat)
     
    theta_n_all    = theta_ini_all*1.0
    
    theta_mean_1col = np.array([0.2819225,0.28597336,0.28728486,0.28643568,0.28358753,0.27879374,0.27252435])

#    theta_mean_1col = np.array([0.27988195, 0.28893737, 0.29351663, 0.29522106, 0.29452995, 0.29169496, 0.28769183])
#    theta_mean_1col = np.array([ 0.30887538, 0.31472939, 0.31682785, 0.31680671, 0.31508827,
#       0.31146861, 0.30576922])

    #------------------------------------------------------------------------
    # INITIAL water deoth (wd) and water surface elevation (we)
    #------------------------------------------------------------------------
    wd_vector, we_vector = overland_setup(eta_vector)                          # vector means reorganize 2D to 1D
    
    #------------------------------------------------------------------------
    # Kl, Kh, Kd, and ADD
    #------------------------------------------------------------------------
    
    Kl          = np.zeros((soil_layer_num,nrows*ncols))
    Kh          = np.zeros((soil_layer_num,nrows*ncols))
    Kd          = np.zeros((soil_layer_num,nrows*ncols))
    
    Zr=z_matrix[:,1]  
    dZr=dz_matrix[:,1]
    Frooti = root_fr(Zr,dZr)
    
    rh_2col = np.column_stack((rh_ini*np.ones(soil_layer_num), np.ones(soil_layer_num)*CN_h_ini/CN_l_ini))
    rh=np.min(rh_2col,axis =1)
    

    litter             = np.mean(litter_input,1)[:,None]*np.ones((2,nrows*ncols))
    litterC            = np.mean(litter_inputC,1)[:,None]*np.ones((2,len(ind_crop)))
    litter[:,ind_crop] = litterC
        
    Addb     = Frooti[:,None]*(Add_b+litter[1,:])
    ADD      = Addb*1.0 
    ADD[0,:] = Addb[0]+ADD_s+litter[0,:]
    
    Cl_stable           = Cl_all_ini*0.865
    Ch_stable           = Ch_all_ini*0.865
    Cb_stable           = Cb_all_ini*0.865
    
#    Cl_stable[:,ind_crop]  = Cl_all_ini[:,ind_crop]*0.86 #0.1235*(4800.0*np.exp(-2.0*z_matrix_ini[:,ind_crop] )+ 1.5)
#    Ch_stable[:,ind_crop]  = Ch_all_ini[:,ind_crop]*0.86 #0.8752*(4800.0*np.exp(-2.0*z_matrix_ini[:,ind_crop] )+ 1.5)
#    Cb_stable[:,ind_crop]  = Cb_all_ini[:,ind_crop]*0.86 #0.0104*(4800.0*np.exp(-2.0*z_matrix_ini[:,ind_crop] )+ 1.5)
    
    
    for i in  ind_in:
        for ll in xrange(0,soil_layer_num):
            Kl[ll,i],Kh[ll,i],Kd[ll,i] = aveg_k(theta_mean_1col[ll]/poros,(Cl_stable[ll,i]), (Ch_stable[ll,i]), (Cb_stable[ll,i]),ADD[ll,i],rh[ll])  
    

#        ind_crop
#    Kl = np.ones((soil_layer_num,nrows*ncols))*Kl[:,None]
#    Kh = np.ones((soil_layer_num,nrows*ncols))*Kh[:,None]
#    Kd = np.ones((soil_layer_num,nrows*ncols))*Kd[:,None]

    #------------------------------------------------------------------------
    #     SAVE VARIABLES 
    #------------------------------------------------------------------------
    water_depth_saved    = np.zeros(nrows*ncols)
#    dh_STs               = np.zeros(nrows*ncols)
#    dh_STd               = np.zeros(nrows*ncols)
#    dh_STr               = np.zeros(nrows*ncols)
#    qC_save_d            = np.zeros((3, nrows*ncols))
#    qC_save_s            = np.zeros((3, nrows*ncols))
#    qC_save_r            = np.zeros((3, nrows*ncols))
    dh_STsave            = np.zeros(nrows*ncols)
    qC_LE                = np.zeros((nrows*ncols))
    dC_BG                = np.zeros((soil_layer_num,nrows*ncols))
    #------------------------------------------------------------------------
    #     SAVE DATA 
    #------------------------------------------------------------------------
    folder =  'Output_data'
    
    if not os.path.exists(folder):
        os.makedirs(folder)
    
    
    save2D_results( dh_STsave.reshape((nrows,ncols)), 
                   qC_LE.reshape((nrows,ncols)), 
                   Cl_all+Ch_all+Cb_all,
                   dC_BG,
                   theta_n_all,
                   z_matrix,
                    'd_LE', 'qC_LE', 'Clhb','dC_BG','theta','z_matrix', 
                   folder, "wd", 0) #"d_soildepth_s" 'qC_s' ,

    eta_vector_ini_in          = -1.0*np.ones(nrows*ncols)
    eta_vector_ini_in[ind_in]  = eta_vector_ini[ind_in]*1.0
    eta_2d_ini                 = eta_vector_ini_in.reshape(nrows,ncols)
    eta_2d_ini[eta_2d_ini<0.0] = np.nan
    np.save('eta_ini',eta_2d_ini)
    
    return eta_vector,eta_vector_ini, nrows, ncols, z_matrix,z_matrix_ini,  dz_matrix,dz_matrix_ini, theta_n_all, Psi_n_all, wd_vector, we_vector, Cl_all, Ch_all, Cb_all,Cl_all_ini, Ch_all_ini, Cb_all_ini, C1_vector, Kl, Kh, Kd, Frooti, mann_bs, mann_vg, ind_in,ind_crop


    