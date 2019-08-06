# -*- coding: utf-8 -*-
"""
Copyright (C) 2017, Qina Corporation
All rights reserved.

Distributed Hydrologicc and Regional Analysis (DHARA) Model
DHARA model is made available as a restricted, non-exclusive, 
non-transferable license for education and research purpose only, 
and not for commercial use. See the LICENSE.txt for more details.

Author: qinayan2@illinois.edu (Qina Yan)
"""

from parameters import *
import numpy as np
from matplotlib import cm
#from scipy.sparse import diags
import timeit
import csv
import os
from netCDF4 import Dataset
import matplotlib.pyplot as plt
#from LEM_initial_landscape import ini_landscape, loadDEM
from overland_Q import overland_setup,matrix_for_rolling, EstimateKOverland, ImplicitSolver_2nd_BC,ImplicitSolver_1st_BC,FlowVelocity
from LEM_sediment_transport import linear_diffusion,Overland_SedimentTransport, rill_erosion
from F_save_netCDF import save2D_results,save2D_results_physical_trans
from F_root_fraction import root_fr
from F_Plot_vertical import SOC_ver_profile, SOC_ver_profile2
from LEM_slope_n_direction import slope_direction
from initial_variables import creating_ini_values, forcing
# pool parallel
#from CN_layers_depth import ini_z_dz, z_dz_update,Klhd_interpolate,C_reslice
from CN_layers_depth_pool_parallel import ini_z_dz, z_dz_update,Klhd_interpolate,C_reslice

from F_Soil_Moisture_main_pool_parallel import vanGenuchten, Soil_Moisture
from CN_biogeochemical_pool_parallel import multi_layer_CN,aveg_k
from CN_mechanical_mix import C_mechamical_mix
#from CN_mechanical_mix_pool_parallel import C_mechamical_mix

#------------------------------------------------------------------------
# Decide which module to operate
#------------------------------------------------------------------------
# topography input: run 'sin' wave topography or the sub-catchment in 'CCW'
CLP = True # True or False

tre = False
unt = True
sin = False
litter_ratio = 1.0

#active different modules
# Soil moisutre
SM    = True
# surface water runoff 
SW    = True
# landscape evolution, diffusion
LED   = True
# landscape evolution, transport-limited
LEAtl = True # must have SW
## do not use this in CLP study landscape evolution, detachment-limited
LEAdl = False
# surface carbon transport
SCT   = True # must have LED or LEA
# Biogeochemical transformation
BG    = True # must have SM
# tillage (mechanical mixing only)
tillage = False
#obtain surface flow velocity 
obtain_velocity = False
# a faster erosion or slow erosion
slower_erosion = False
if slower_erosion:
    level = 2
#    nc_f = "Sine3_50yr_noTillage.nc" #Sine3_50yr_noTillage, Sine6_50yr_fastTillage
else:
    level = 1
#    nc_f = "wd0.nc"

rain,litter_input,litter_inputC = forcing(unt, tre)
inputDatabase                   = False

eta_vector,eta_vector_ini, nrows, ncols, z_matrix,z_matrix_ini,  dz_matrix,\
dz_matrix_ini, theta_n_all, Psi_n_all, wd_vector, we_vector, Cl_all, Ch_all,\
 Cb_all,Cl_all_ini, Ch_all_ini, Cb_all_ini, C1_vector, Kl, Kh, Kd, Frooti, \
 mann_bs, mann_vg,ind_in,ind_crop \
 = creating_ini_values(CLP,sin,litter_input,litter_inputC,inputDatabase,tre,unt)

#------------------------------------------------------------------------
# BUILD ROLLING MATRIX TO OBTAIN H_i+1,j, Hi,j+1 etc. 
#------------------------------------------------------------------------
A_p1,A_m1,A_pN,A_mN, A_p1N,A_m1N,A_pNm1,A_mNp1,BC_we  = matrix_for_rolling(nrows,ncols)

#------------------------------------------------------------------------
# UPDATE NEW ELEVATION ETA & TOP CARBON FROM DIFFUSION PROCESS
#------------------------------------------------------------------------
#matrix_Ad, matrix_BC,matrix_Adx,matrix_Ady = linear_diffusion(nrows,ncols,level)

matrix_Ad, matrix_BC,\
matrix_Adxp1, matrix_Adxm1, matrix_Adyp1, matrix_Adym1,\
BCx_p1, BCx_m1, BCy_p1,BCy_m1 = linear_diffusion(nrows,ncols,level)


"""
Time loop is created by Qina
This loop includes several submodels:
    
    - LEM: Surface soil+carbon transport
        + 
    - BIOGEOCHEM: C-N Transformation
        +
    - HYDROLOGY: Soil moisture
"""

Total_period =dt*365*50 #*365*50
time_steps = int(Total_period/dt)
#time_steps = 100
#mean_depth = np.zeros(time_steps)
print('dt = %4.2f ' % dt)
print('time_steps = %4.0f ' % time_steps)

#------------------------------------------------------------------------
#     SAVE VARIABLES 
#------------------------------------------------------------------------
i_save               = 0 #save the idex in NetCDF file
water_depth_saved    = np.zeros(nrows*ncols)
rain_depth           = np.zeros(nrows*ncols)
dh_ST                = np.zeros(nrows*ncols)
dh_STsave            = np.zeros(nrows*ncols)
dC_BG                = np.zeros((soil_layer_num, nrows*ncols)) 
qC_LE                = np.zeros((nrows*ncols))
Cl_all_update        = np.zeros((soil_layer_num, nrows*ncols)) 
Ch_all_update        = np.zeros((soil_layer_num, nrows*ncols)) 
Cb_all_update        = np.zeros((soil_layer_num, nrows*ncols)) 
dz1                  = dz_matrix[0,:]
dC_BG_monthly        = np.zeros((soil_layer_num, nrows*ncols)) 

#------------------------------------------------------------------------
#     temporal save variables
#------------------------------------------------------------------------
dCl_all              = np.zeros((soil_layer_num, nrows*ncols)) 
dCh_all              = np.zeros((soil_layer_num, nrows*ncols)) 
dCb_all              = np.zeros((soil_layer_num, nrows*ncols)) 
Cl1d_all             = np.zeros((soil_layer_num, nrows*ncols)) 
Ch1d_all             = np.zeros((soil_layer_num, nrows*ncols)) 
Cb1d_all             = np.zeros((soil_layer_num, nrows*ncols)) 

#eta_change_save1      = np.zeros(time_steps)
#eta_change_save2      = np.zeros(time_steps)

folder = 'Output_data'
if not os.path.exists(folder):
    os.makedirs(folder)
#/projects/users/qinayan2

start = timeit.default_timer()                        # trace the total time

for i in xrange(time_steps): 
    # 18252 is the maximum day can be put in
#    print i
    #::::::::::::::::::::::::::::::::::::::::::::::::::::::
    # SOIL MOISTURE 
    #::::::::::::::::::::::::::::::::::::::::::::::::::::::
    if SM:                                              #rain[i] or PPT
#        theta_n_all,Psi_n_all,rain_depth,Ts_n_all,Leak_n_all = \
#            Soil_Moisture(rain[i],theta_n_all, Psi_n_all, z_matrix,dz_matrix,wd_vector,nrows,ncols)
        
        theta_n_all_ins,Psi_n_all_ins,rain_depth_ins,Ts_n_all,Leak_n_all = \
            Soil_Moisture(rain[i],theta_n_all[:,ind_in], Psi_n_all[:,ind_in], z_matrix[:,ind_in],dz_matrix[:,ind_in],wd_vector[ind_in],len(ind_in)) 
         
        theta_n_all[:,ind_in] =  theta_n_all_ins
        Psi_n_all[:,ind_in]   = Psi_n_all_ins
        rain_depth[ind_in]    = rain_depth_ins  
            
    else:
        theta_n_all = 0
        Psi_n_all   = 0
        rain_depth  = rain[i] #rain[i] or PPT
        
        
    #::::::::::::::::::::::::::::::::::::::::::::::::::::::
    # SURFACE water flow and water depth, wd_vector
    # This modules models
    #::::::::::::::::::::::::::::::::::::::::::::::::::::::
    if SW:
        wd_vector, we_vector = overland_setup(eta_vector) 
        
        if  i>0: #( i%(365*2)> 163 and i%(365*2)< 263   ):
            mann = mann_vg
        else: 
            mann = mann_bs    
        
        for flowtime in xrange(2):
            #ke, kw, kn, ks, khe, khw, khn, khs
            ke,kw,kn,ks,Sse,Ssw,Ssn,Sss = EstimateKOverland(A_p1,A_m1,A_pN,A_mN, A_p1N,A_m1N,A_pNm1,A_mNp1, wd_vector, we_vector,mann)
    
            ##rain[i],PPT
            we_vector_new, wd_vector_new = ImplicitSolver_2nd_BC(we_vector,eta_vector,rain_depth/8.0,ke,kw,kn,ks,nrows, ncols,BC_we) 
            #. . . Update variables
    #    
            wd_vector   = wd_vector_new
            we_vector   = wd_vector+ eta_vector
            
        if obtain_velocity:
            we_p1 =  A_p1.dot(we_vector) 
            we_m1 =  A_m1.dot(we_vector) 
            we_mN =  A_mN.dot(we_vector) 
            we_pN =  A_pN.dot(we_vector)
            u_vector,v_vector   = FlowVelocity(we_p1,we_m1,we_mN,we_pN,wd_vector,ke,kw,kn,ks)


    #::::::::::::::::::::::::::::::::::::::::::::::::::::::
    # C-N TRANSFORMATION
    # This module simulates the C-N
    #::::::::::::::::::::::::::::::::::::::::::::::::::::::

    if BG:
        if LED or LEAtl:
            Klt = Kl[0,:]*1.0
            Kht = Kh[0,:]*1.0
            Kdt = Kd[0,:]*1.0
    
            dh_pID = np.where(dh_ST>0.0) # positive index of dh_ST
            Kl[0,dh_pID] = (Klt[dh_pID]*1.10*dh_ST[dh_pID]+ Klt[dh_pID]*dz1[dh_pID])/(dz1[dh_pID]+dh_ST[dh_pID])
            Kh[0,dh_pID] = (Kht[dh_pID]*1.11*dh_ST[dh_pID]+ Kht[dh_pID]*dz1[dh_pID])/(dz1[dh_pID]+dh_ST[dh_pID])
            Kd[0,dh_pID] = (Kdt[dh_pID]*1.06*dh_ST[dh_pID]+ Kdt[dh_pID]*dz1[dh_pID])/(dz1[dh_pID]+dh_ST[dh_pID])
        
        litter_above           = ADD_s + litter_input[0,i%(365)]*np.ones(nrows*ncols) #litter_input + ADD_s
        litter_above[ind_crop] = litter_inputC[0,i%(365)]
        
        litter_below           = (litter_input[1,i%(365)]*np.ones(nrows*ncols) + Add_b)*Frooti[:,None]
        litter_belowC          = Frooti[:,None]*litter_inputC[1,i%(365)]
        litter_below[:,ind_crop] = litter_belowC
#        litter_below = np.ones((soil_layer_num,nrows*ncols))*litter_below[:,None]

#        Cl_all,Ch_all,Cb_all \
#             = multi_layer_CN(Cl_all,Ch_all,Cb_all,theta_n_all/poros,z_matrix,dz_matrix,Kl,Kh,Kd,litter_above,litter_below,nrows*ncols)
        
        dCl_all[:,ind_in],dCh_all[:,ind_in],dCb_all[:,ind_in], Cl1d_all[:,ind_in], Ch1d_all[:,ind_in], Cb1d_all[:,ind_in] \
             = multi_layer_CN(Cl_all[:,ind_in],Ch_all[:,ind_in],Cb_all[:,ind_in],\
                              theta_n_all[:,ind_in]/poros,z_matrix[:,ind_in],\
                              dz_matrix[:,ind_in],Kl[:,ind_in],Kh[:,ind_in],\
                              Kd[:,ind_in],litter_above[ind_in]*litter_ratio,litter_below[:,ind_in]*litter_ratio,len(ind_in))

#        Cl_all[:,ind_in] = Cl_all_ins
#        Ch_all[:,ind_in] = Ch_all_ins
#        Cb_all[:,ind_in] = Cb_all_ins
    else:
        dCl_all  = 0.0
        dCh_all  = 0.0
        dCb_all  = 0.0
        Cl1d_all = 0.0
        Ch1d_all = 0.0
        Cb1d_all = 0.0
        
    #::::::::::::::::::::::::::::::::::::::::::::::::::::::
    # SURFACE SOIL TRANSPORT - LANDSCAPE EVOLUTION MODEL
    # This modules models
    #::::::::::::::::::::::::::::::::::::::::::::::::::::::
    if LED:
        eta_d_vector              = matrix_Ad.dot(eta_vector) + matrix_BC
        eta_d_vector[ind_crop]    = 1.1*eta_d_vector[ind_crop] 
#        eta_d_vector[ind_out]     = 0.0
        eta_d_vector[np.abs(eta_d_vector)>1.2e-6] = 0.0
    else:
        eta_d_vector              = 0.0

    if LEAtl or LEAdl:
        ### . . . Slope, flow direction, and drainage area
        slope,direction           = slope_direction(eta_vector,nrows,ncols)
    if LEAtl:
        if  i>0: #( i%(365*2)> 163 and i%(365*2)< 263   ):
            mann = mann_vg
        else: 
            mann = mann_bs   
#        eta_s_vector, qsCa_vector = Overland_SedimentTransport(wd_vector,Sse,Ssw,Ssn,Sss,A_p1,A_m1,A_pN,A_mN,C1_vector,nrows,ncols,level) 
        eta_s_vector, qsCa_vector = Overland_SedimentTransport(wd_vector,slope,direction,C1_vector,nrows,ncols,level,mann) 
        eta_s_vector[ind_crop]    = 1.1*eta_s_vector[ind_crop] 
        eta_s_vector[np.abs(eta_s_vector)>3.0e-4] = 0.0
        qsCa_vector[np.abs(qsCa_vector)>2.0] = 0.0
        eta_s_vector[np.abs(eta_s_vector)>2.2e-4] = 2.2e-4
        
    else:                                                      
        eta_s_vector = 0.0
        qsCa_vector  = 0.0
        
    if LEAdl:
        eta_r_vector,qrCa_vector = rill_erosion(wd_vector,slope,C1_vector,nrows,ncols)
    else:                                                      
        eta_r_vector = 0.0
        qrCa_vector  = 0.0
        
    dh_ST                     = ( eta_d_vector + eta_s_vector)/(1.0-lamdaP)+eta_r_vector  #eta_s_vector

    eta_vector_update         = dh_ST + uplift*dt+ eta_vector
    
#    dh_STd += eta_d_vector/(1.0-lamdaP)
#    dh_STs += eta_s_vector/(1.0-lamdaP)
#    dh_STr += eta_r_vector
    
    
    #::::::::::::::::::::::::::::::::::::::::::::::::::::::
    # SURFACE CARBON TRANSPORT 
    # This modules models
    #::::::::::::::::::::::::::::::::::::::::::::::::::::::
    dz1           = dz_matrix[0,:]

    if LED:
        qdCa_vector =  0.5*C1_vector*eta_d_vector \
        + np.transpose(A_p1.dot(C1_vector.T)) * (matrix_Adxp1.dot(eta_vector)+BCx_p1) \
        - np.transpose(A_m1.dot(C1_vector.T)) * (matrix_Adxm1.dot(eta_vector)+BCx_m1) \
        + np.transpose(A_pN.dot(C1_vector.T)) * (matrix_Adyp1.dot(eta_vector)+BCy_p1) \
        - np.transpose(A_mN.dot(C1_vector.T)) * (matrix_Adym1.dot(eta_vector)+BCy_m1)
        qdCa_vector[np.abs(qdCa_vector)>8e-3] = 0.0
    else:
        qdCa_vector = qdCa_vector  = 0.0

    if SCT:
#        qdCa_vector   = C1_vector*eta_d_vector #+ np.transpose(matrix_Adx.dot(C1_vector.T)) * matrix_Adx.dot(eta_vector) \
                    #+ np.transpose(matrix_Ady.dot(C1_vector.T))* matrix_Ady.dot(eta_vector)       

        qCa_vector    = (qdCa_vector+qsCa_vector)*Ksoc/(1.0-lamdaP)+qrCa_vector*Ksoc # qdCa_vector+qsCa_vector 

#        C1_vector_new = (qCa_vector+ C1_vector * dz1)/(dz1+dh_ST)
#        
#        C1_vector_new = np.maximum(C1_vector_new,np.array([Cl_all_ini[soil_layer_num-1,:],Ch_all_ini[soil_layer_num-1,:],Cb_all_ini[soil_layer_num-1,:]]))
#
##        print C1_vector_new[2,891]
#        qC_save_d     += qdCa_vector*Ksoc/(1.0-lamdaP)
#        qC_save_s     += qsCa_vector*Ksoc/(1.0-lamdaP)
#        qC_save_r     += qrCa_vector*Ksoc 
#        Cl_all[0,:]   = C1_vector_new[0,:]*1.0    #*prec_Cl*1.0 # *dhC+C_matrix[0,:]*dz1)/(dhC+dz1)
#        Ch_all[0,:]   = C1_vector_new[1,:]*1.0    #*prec_Ch*1.0
#        Cb_all[0,:]   = C1_vector_new[2,:]*1.0    #*prec_Cb*1.0
        
    else:
         qCa_vector    = np.zeros((3,nrows*ncols))
#        C1_vector_new = C1_vector*1.0
        
    
    #::::::::::::::::::::::::::::::::::::::::::::::::::::::
    # GLOBAL PROCESSING
    #::::::::::::::::::::::::::::::::::::::::::::::::::::::
    #. . . Update variables

    eta_vector         = eta_vector_update
    water_depth_saved += wd_vector
    dh_STsave         += dh_ST
    
    if BG:
        #dC_BG             += dCl_all + dCh_all + dCb_all
        dC_BG             = dCl_all + dCh_all + dCb_all - litter_below*dt
        dC_BG[0,:]        = dC_BG[0,:] - litter_above*dt
        dC_BG_monthly     = dC_BG_monthly + dC_BG
    # save accumulated SOC change by LE
    qC_LE             += np.sum(qCa_vector,0)
    
    # . . . update Cl_all, Ch_all, Cb_all caused by BG and bioturbation
    Cl_all_update = Cl_all  + dCl_all + Cl1d_all
    Ch_all_update = Ch_all  + dCh_all + Ch1d_all
    Cb_all_update = Cb_all  + dCb_all + Cb1d_all
    
    # SOC surface extra cost by surface transport:
    Cl_all_update[0,:] = (Cl_all_update[0,:]*dz1+ qCa_vector[0,:])/(dz1+dh_ST)
    Ch_all_update[0,:] = (Ch_all_update[0,:]*dz1+ qCa_vector[1,:])/(dz1+dh_ST)
    Cb_all_update[0,:] = (Cb_all_update[0,:]*dz1+ qCa_vector[2,:])/(dz1+dh_ST)
    
    Cl_all = Cl_all_update*1.0
    Ch_all = Ch_all_update*1.0
    Cb_all = Cb_all_update*1.0
    
#    eta_change_save1[i] = dh_STs[19*ncols+12]
#    eta_change_save2[i] = dh_STs[12*ncols+24]
    # . . . Soil depth and soil layer thickness    
    if SCT:
        z_matrix_update, dz_matrix = z_dz_update(z_matrix,dh_ST,nrows, ncols) #z_matrix_update, dz_matrix_update
    
    # . . . C interpolate
        Cl_all,Ch_all,Cb_all = C_reslice(Cl_all,Ch_all,Cb_all, z_matrix,z_matrix_update,dh_ST,nrows, ncols)
                               #C_interpolate(Cl_all,Ch_all,Cb_all,z_matrix,z_matrix_update,dh_ST,nrows,ncols)

    # . . . Kl, Kh, and Kd interpolate
#    if BG:
#        Kl, Kh, Kd = Klhd_interpolate(Kl, Kh, Kd, z_matrix,z_matrix_update,dh_ST,nrows, ncols)
        
    # . . . Mechanical mixing
    if tillage: # DOY = 105 every year
        if (i+1)%365 == 105: #(i+1)%(365)== 0:  (i+1)%(time_steps) == 0:  (i+1)%365 == 105: 
            Cl_all,Ch_all,Cb_all = C_mechamical_mix(Cl_all,Ch_all,Cb_all, z_matrix,dz_matrix,ind_crop, level)
    
    if SCT:
        z_matrix   = z_matrix_update*1.0
    
        C1_vector = np.array([Cl_all[0,:],Ch_all[0,:], Cb_all[0,:]])*1.0


    #::::::::::::::::::::::::::::::::::::::::::::::::::::::
    # plot within the loop
    #::::::::::::::::::::::::::::::::::::::::::::::::::::::
    if (i+1-365*(i//365))%366 in DOY_at_month: 
#    if (i+1)%(time_steps) == 0: 
    #------------------------------------------------------------------------
    #     SAVE DATA 
    #------------------------------------------------------------------------
        i_save += 1
        save2D_results( dh_STsave.reshape((nrows,ncols)), 
                       qC_LE.reshape((nrows,ncols)), 
                       Cl_all+Ch_all+Cb_all,
                       dC_BG_monthly,
                       theta_n_all,
                       z_matrix,
                        'd_LE', 'qC_LE', 'Clhb','dC_BG','theta','z_matrix', 
                       folder, "wd", i_save) #"d_soildepth_s" 'qC_s' ,


        dC_BG_monthly  = np.zeros((soil_layer_num, nrows*ncols))

        ### plot the flow direction
#        if obtain_velocity:
#            Y = np.arange(0, nrows*dy, dy)
#            X = np.arange(0, (ncols)*dx, dx) 
#            X, Y = np.meshgrid(X, Y)
#            u = u_vector.reshape((nrows,ncols))
#            v = v_vector.reshape((nrows,ncols))
#            plt.figure()  
#            Q = plt.quiver(X[::3, ::3],Y[::3, ::3],u[::3, ::3],v[::3, ::3],color='r', units='x',linewidths=(0.5,), edgecolors=('k'), headaxislength=5)
#            qk = plt.quiverkey(Q, 0.5, 0.03, 1, r'$1 \frac{m}{s}$', fontproperties={'weight': 'bold'})
#
#        ###.... plot the water depth  
#        if SW:
#            wd       = water_depth_saved.reshape((nrows,ncols))#qdC_vector.reshape((nrows,ncols))
#            fig1, ax1 = plt.subplots(figsize=(5.3,2.5 ), dpi=180)
#            plt.rc("font", size=10)
#            flip_plot = False
#            if flip_plot:
#                figplot  = ax1.matshow(np.fliplr(wd), extent=[0,ncols*dx,nrows*dy,0],origin="lower",vmin=0,vmax=1.0)  
#            else:
#                figplot  = ax1.matshow(wd, extent=[0,ncols*dx,nrows*dy,0]) #,vmin=0,vmax=0.032 #pcolor , matshow, imshow
#            cbar = fig1.colorbar(figplot,shrink=0.9, pad = 0.1)
#            cbar.set_label('Water Depth [m]')
#    #        ax1.set_title('Water Depth [m]', fontsize=20)
#            ax1.set_xlabel('x [m]', fontsize=10)
#            ax1.set_ylabel('y [m]', fontsize=10)
#        
#        ## save the water depth   
##        current_dir = os.getcwd()
##        os.chdir(current_dir+"/"+'/Figures/')
##        fig1.savefig('water_depth_'+str(round(i*dt,3))+'hr.eps') 
#        
#        ###... plot soil moisture
#        if SM:
#            S_ave = np.sum(theta_n_all*dz_matrix, 0)/z_matrix[soil_layer_num-1,:]#C1_vector_new*(dz1)
#            sd           = S_ave.reshape((nrows,ncols))#qdC_vector.reshape((nrows,ncols))
#            fig2, ax2 = plt.subplots(figsize=(5.3,2.5 ), dpi=90)     
#            plt.rc("font", size=14)
#            figplot  = ax2.matshow(sd, extent=[0,ncols*dx,nrows*dy,0],cmap=cm.terrain) #,vmin=0,vmax=0.032 #pcolor , matshow, imshow
#            cbar = fig2.colorbar(figplot,shrink=0.7, pad = 0.1)  #if shrink, shrink=0.9,
#            cbar.set_label(r'soil moisture')
#        #        ax2.set_title(r'$\Delta$ Soil Depth [m]', fontsize=20)
#            ax2.set_xlabel('x [m]', fontsize=18)
#            ax2.set_ylabel('y [m]', fontsize=18)
#    
#    
#        ###... plot the soil depth change 
#        if LEAdl or LEAtl or LED:
#
#            sd_vector    = eta_vector-eta_vector_ini
#            sd           = dh_STsave.reshape((nrows,ncols))#qdC_vector.reshape((nrows,ncols))
#            fig2, ax2 = plt.subplots(figsize=(5.3,2.5 ), dpi=180)     
#            plt.rc("font", size=14)
#            figplot  = ax2.matshow(sd, extent=[0,ncols*dx,nrows*dy,0],cmap=cm.terrain) #,vmin=0,vmax=0.032 #pcolor , matshow, imshow
#            cbar = fig2.colorbar(figplot,shrink=0.7, pad = 0.1)  #if shrink, shrink=0.9,
#            cbar.set_label(r'$\Delta$ Soil Depth [m]')
#    #        ax2.set_title(r'$\Delta$ Soil Depth [m]', fontsize=20)
#            ax2.set_xlabel('x [m]', fontsize=18)
#            ax2.set_ylabel('y [m]', fontsize=18)
#        
##        if LED:
##            sd           = dh_ST.reshape((nrows,ncols))#qdC_vector.reshape((nrows,ncols))
##            fig2, ax2 = plt.subplots(figsize=(10.3,4.5 ), dpi=120)     
##            plt.rc("font", size=12)
##            figplot  = ax2.matshow(sd, extent=[0,ncols*dx,nrows*dy,0],cmap=cm.terrain) #,vmin=0,vmax=0.032 #pcolor , matshow, imshow
##            cbar = fig2.colorbar(figplot,shrink=0.7, pad = 0.1)  #if shrink, shrink=0.9,
##            cbar.set_label(r'$\Delta$ Soil Depth [m]')
##    #        ax2.set_title(r'$\Delta$ Soil Depth [m]', fontsize=20)
##            ax2.set_xlabel('x [m]', fontsize=12)
##            ax2.set_ylabel('y [m]', fontsize=12)
#            
#
#        ###... plot the top soil carbon concentration 
#        if BG or SCT:
#            C_matrix             = Cl_all+Ch_all+Cb_all
#            C_matrix_ini         = Cl_all_ini+Ch_all_ini+Cb_all_ini
#            C_sum = np.sum(C_matrix*z_matrix,0)-np.sum(C_matrix_ini*z_matrix_ini, 0)#C1_vector_new*(dz1)
#            sd           = C_sum.reshape((nrows,ncols))#qdC_vector.reshape((nrows,ncols))
#            fig2, ax2 = plt.subplots(figsize=(10.3,4.5 ), dpi=120)     
#            plt.rc("font", size=12)
#            figplot  = ax2.matshow(sd, extent=[0,ncols*dx,nrows*dy,0],cmap=cm.terrain) #,vmin=0,vmax=0.032 #pcolor , matshow, imshow
#            cbar = fig2.colorbar(figplot,shrink=0.7, pad = 0.1)  #if shrink, shrink=0.9,
#            cbar.set_label(r'$SOC [g \ C/m^3]$')
#        #        ax2.set_title(r'$\Delta$ Soil Depth [m]', fontsize=20)
#            ax2.set_xlabel('x [m]', fontsize=12)
#            ax2.set_ylabel('y [m]', fontsize=12)
#
#
#		##... vertical plot of C concentration 
#            cell_ero_id = np.argmin(dh_ST)
#            cell_depo_id =np.argmax(dh_ST) #955, 915, 917, 918
#            C_matrix = Cl_all+Ch_all+Cb_all
#            SOC_ver_profile2(cell_ero_id, cell_depo_id, C_matrix, C_matrix_ini, z_matrix, z_matrix_ini)

############################################################ 
# END OF TIME LOOP 
############################################################
stop = timeit.default_timer()-start
print stop


















