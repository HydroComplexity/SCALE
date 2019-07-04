# -*- coding: utf-8 -*-
"""
Created on Thu Feb 23 15:47:54 2017

@author: Qina and Phong

This code is to simulate soil  moisture in 1D soil column

The current version has:
1. uniform grid size
2. uniform soil properties, e.g. builk density \rho_s, Ksat, 
3. BC: Constant head at two ends

Key state variables:
    Psi   : pressure head [L]
    theta : soil moisture [-]

Other variables:
    C    : Specific moisture capacity [1/L]
    K    : Hydraulic conductivity [L/T]

Notation:
    _ini    : Initial condition (time 0)
    _n      : At time step n
    _np1    : At time n+1 (known)
    _np1m   : At time n+1, iteration m (known)
    _np1mp1 : At time n+1, iteration m+1 (unknown)
    
"""

from parameters import *
import numpy as np
from scipy.sparse import diags
from scipy.sparse import dia_matrix
from scipy.sparse.linalg import spsolve
import matplotlib.pyplot as plt
import timeit
from F_root_fraction import root_fr

from CN_layers_depth import ini_z_dz
import matplotlib.pyplot as plt
from matplotlib import cm

def vanGenuchten(Psi,alpha, theta_S, theta_R, n, m, Ksat):

    #Convert unit of h from [m] to [cm] to match with alpha
    Psi = Psi*100.0
    
    #Compute the volumetric moisture content [eqn 21]
    theta = theta_R + (theta_S - theta_R)/(1.0 + (alpha*np.abs(Psi)*1.0)**n)**m
    
    theta[Psi>=0.0] = theta_S
    
    #Compute the effective saturation [eqn 2]
    Se = ((theta - theta_R)/(theta_S - theta_R))
     
    #Compute the hydraulic conductivity [eqn 8]
    K = Ksat*Se**(0.5)*(1.0 - (1.0 - Se**(1.0/m))**m)**2.0
     
    #Compute the specific moisture storage: 
    #C = d(theta)/dh
    C = alpha*n*(1.0/n - 1.0)*(alpha*np.abs(Psi))**(n - 1.0)*(theta_R - theta_S)*((alpha*np.abs(Psi))**n + 1.0)**(1.0/n - 2.0)*100.0
#    C = -100.0/Psi*m*(theta_S - theta_R)*Se*(1-Se**(1.0/m))/(1.0-m)
    C[Psi>=0.0] = 0.0

    return C,K,theta
    
def Khalf(nz,K):
    
    # K_i+1/2, Kp1
    diag_p = np.ones(nz,dtype = 'int')
    diag_p[nz-1] = 2
    A_p = diags( [diag_p,1], [0,1],shape=(nz, nz))
    K_ph = A_p.dot(K)/2.0 
    
    # K_i-1/2, Km1
    diag_m = np.ones(nz,dtype = 'int')
    diag_m[0] = 2
    A_m = diags([diag_m,1], [0,-1],shape=(nz,nz))
    K_mh = A_m.dot(K)/2.0

    return K_ph,K_mh

def ImplicitSolver_Richards(Ts, Psi_n, Psi_np1m, theta_n, theta_np1m, C_np1m, K_np1m_ph, K_np1m_mh, \
                   dz, topbound, Psi_top, q_top, botbound, Psi_bot, q_bot):
    
    # Build the diagal and up,low diagnal for matrix A, where Ax = b
    Mdia = Ss*theta_np1m/poros/dt + C_np1m/dt + (K_np1m_ph + K_np1m_mh)/dz**2
    Udia = -K_np1m_ph/dz**2
    Ldia = -K_np1m_mh/dz**2
    
    # known values on the RHS
    RHS = Ss*theta_np1m/poros*Psi_n/dt + (theta_n-theta_np1m)/dt + C_np1m*Psi_np1m/dt - (K_np1m_ph-K_np1m_mh)/dz - Ts
    
    Ldia_r  = np.roll(Ldia,-1)
    Udia_r  = np.roll(Udia,1)
    
    ### build the matrix A, Ax = b
    A1d = dia_matrix( ([Ldia_r,Mdia,Udia_r],[-1,0,1]),shape=(nz, nz))
    A1d = A1d.tocsr()
    
    # TOP BOUNDARY CONDITION
    if topbound == 0:         # Dirichlet BCs
        A1d[0,0] = 1.0
        A1d[0,1] = 0.0
        RHS[0] = Psi_top 
    else:
#        A1d[0,0] = A1d[0,0] + Ldia[0] #- K_np1m_mh[0]/dz[0]**2
#        RHS[0] = RHS[0] + (dz[0] + q_top*dz[0]/K_np1m_mh[0])*(-Ldia[0])  #/dz[0]**2 * K_np1m_ph[0]
        A1d[0,1] = A1d[0,1] + Ldia[0]
        RHS[0] = RHS[0] + (2.0*dz[0]+2.0*dz[0]*q_top/((K_np1m_ph[0]+K_np1m_mh[0])/2.0))*(-Ldia[0])
    # BOTTOM BOUNDARY CONDITION            
    if botbound == 0:         # Dirichlet BCs
        A1d[-1,-1] = 1.0
        A1d[-1,-2] = 0.0
        RHS[-1] = Psi_bot
#        RHS[-1] = RHS[-1]-Udia_r[-1]*Psi_np1m[-2]
    else:
#        A1d[-1,-1] = A1d[-1,-1] + Udia[-1] #- K_np1m_ph[-1]/dz[-1]**2
#        RHS[-1] = RHS[-1] + (-dz[-1] - q_bot*dz[-1]/K_np1m_ph[-1])*(-Udia[-1])  #/dz[-1]**2 * K_np1m_ph[-1]
        A1d[-1,-2] = A1d[-1,-2] + Udia[-1]
        RHS[-1] = RHS[-1] + (-2.0*dz[-1]-2.0*dz[-1]*q_bot/((K_np1m_ph[-1]+K_np1m_mh[-1])/2.0))*(-Udia[-1])

    ### Solve x, Ax = b
    Psi_np1mp1 = spsolve(A1d,RHS)
    
    return Psi_np1mp1

#------------------------------------------------------------------------
# Initial Condition and grid size
#------------------------------------------------------------------------
def Soil_Moisture(rainfall, theta_n_all, Psi_n_all, z_matrix,dz_matrix,wd_vector,nrows,ncols):
#......... Model Stopping tolerance for convegence 
    
    stop_tol = 0.2
    max_iter = 4
    soil_layer_nums = soil_layer_num+1
    theta_n_all_new = np.zeros((soil_layer_num,nrows*ncols))
    Psi_n_all_new = np.zeros((soil_layer_num,nrows*ncols))
    Leak_n_all_new = np.zeros((soil_layer_num,nrows*ncols))
    Ts_n_all_new = np.zeros((soil_layer_num,nrows*ncols))
    rain_depth = np.zeros(nrows*ncols)

    #------------------------------------------------------------------------
    # Boundary Conditions (constant head at both end)
    #------------------------------------------------------------------------
    """
    Boundary Conditions 
        0: Dirichlet Type (const head), 1: Neumann (flux) Type
    """
    
    # Bottom Boundary
#    botbound = 0
#    psi_lower = -2.0#*np.ones(time_steps)   # Lower boundary 
#    q_lower = -0.005#*np.ones(time_steps)    # Lower flux boundary 
     
    botbound = 1
    psi_lower = -1.0 #-np.random.rand(1)-1.7 #*np.ones(time_steps)   # Lower boundary 
#    q_lower = 0.2#Ksat #K_np1m #-0.005#*np.ones(time_steps)    # Lower flux boundary 
    
    Lp1 = np.arange(soil_layer_num)
    Lp1 = np.roll(Lp1,-1)
    
    for k in xrange(nrows*ncols):
#        print k
        z = z_matrix[:,k]
        dz = dz_matrix[:,k]
        theta_n = theta_n_all[:,k]
        Psi_n = Psi_n_all[:,k]
        
        # Tranpiration
        Frooti = root_fr(z,dz)
        Tsm = Frooti*T_max
        theta_Ts = (theta_n/theta_S-sw)/(s_star-sw)
        theta_Ts[theta_Ts<0.0]=0.0
        theta_Ts[theta_Ts>1.0]=1.0
        Ts = theta_Ts*Tsm

        #Evaporation
        theta_E =  (theta_n[0]/theta_S-sh)/(sw-sh)
        theta_E = np.max([0.0,theta_E]) 
        theta_E = np.min([1.0,theta_E]) 
        Evaps = theta_E*Evap
        
        # Top Boundary Condition
        min_cap = np.min([Infil_max, (theta_S-theta_n_all[0,k])*dz_matrix[0,k]])

        if ((wd_vector[k]>0.0) and (rainfall>0.0) ):
            topbound = 1
            rain_depth[k] = rainfall
            psi_upper = rainfall+wd_vector[k]
            q_upper= min_cap #0.0
        elif ((wd_vector[k]>0.0) and (rainfall==0.0) ):
            topbound = 1
            psi_upper = rainfall+wd_vector[k]-Evaps
            q_upper= min_cap-Evaps #0.0
        elif (rainfall>min_cap and wd_vector[k]==0.0 ) :
            topbound = 1
            rain_depth[k] = rainfall-min_cap
            psi_upper =  rainfall-min_cap #  # Upper boundary
            q_upper = min_cap #0.0
        elif (rainfall<=min_cap and wd_vector[k]==0.0 and rainfall>0.0 ) :  
            topbound = 1
            q_upper = rainfall #*np.ones(time_steps)    # Upper flux boundary
            psi_upper = 0.0
        elif (wd_vector[k]==0.0 and rainfall==0.0 ) :
            topbound = 1
            q_upper = -Evaps
            psi_upper = -1.0
#        print topbound
#        print q_upper
        #------------------------------------------------------------------------
        # The Main loop
        #------------------------------------------------------------------------
        
        diff_mp1=np.zeros((nz,max_iter))
#        
#        psi_store =np.zeros((nz,time_steps))
#        theta_store =np.zeros((nz,time_steps))
#        
        Psi_np1m = Psi_n*1.0
        loop_num = 0
        deltam   = 10000.0
        while deltam > stop_tol and loop_num<max_iter:
            loop_num = loop_num+1
    #        print loop_num
            C_np1m,K_np1m,theta_np1m = vanGenuchten(Psi_np1m,alpha, theta_S, theta_R, n, m, Ksat)
            K_np1m_ph,K_np1m_mh =  Khalf(nz,K_np1m)

#            start               = timeit.default_timer()
        
            Psi_np1mp1 = ImplicitSolver_Richards(Ts, Psi_n, Psi_np1m, theta_n, theta_np1m, \
                                             C_np1m, K_np1m_ph, K_np1m_mh, dz,\
                                             topbound, psi_upper, q_upper, \
                                             botbound, psi_lower , K_np1m[soil_layer_num-1] ) #botbound psi_lower  q_upper q_lower  K_np1m[soil_layer_num-1]
#            stop = timeit.default_timer()-start
#            print stop
            diff_mp1[:,loop_num-1] = (Psi_np1mp1 - Psi_np1m)
#            print Psi_np1m[0],Psi_np1mp1[0]
            deltam = np.max(np.abs(Psi_np1mp1 - Psi_np1m))
    
            Psi_np1m = Psi_np1mp1*1.0
    #        C_np1mp1, K_np1mp1, theta_np1mp1 = vanGenuchten(Psi_np1mp1, alpha, theta_S, theta_R, n, m, Ksat)
    #        plt.plot(Psi_np1m, -dz*np.arange(nz))
    
#        print loop_num
        Psi_n = Psi_np1mp1*1.0
        
        # Initialize the Picard iteration solver
        C_n,K_n,theta_n = vanGenuchten(Psi_n,alpha, theta_S, theta_R, n, m, Ksat)
        
        Leak     = -(K_n*((Psi_n[Lp1]-Psi_n)/dz+1))
        Leak[-1] = K_np1m[soil_layer_num-1]
        
        theta_n_all_new[:,k] = theta_n
        Psi_n_all_new[:,k]   = Psi_n
        Ts_n_all_new[:,k]  = Ts
        Leak_n_all_new[:,k]  = Leak
#        plt.plot(Psi_n_all_new[:,0], -z_matrix[:,0])

        
    return theta_n_all_new,Psi_n_all_new,rain_depth, Ts_n_all_new,Leak_n_all_new #ponding_depth


if __name__ == '__main__':
    #...................Rainfall
    rain = np.load('rainfall_100yr3_daily.npy') # mm/(day)
    rain = rain/1000.0
#    rain50 = rain62[0:50*365]/1000.0
#    rain = np.hstack([rain50,rain50[::-1]]) 
    
    
    time_steps = 365*100
    rp = 0
    #...................initial values
    theta_save = np.zeros((soil_layer_num,time_steps))
    Psi_save = np.zeros((soil_layer_num,time_steps))
    rain_depth_save = np.zeros(time_steps)
    z_matrix, dz_matrix = ini_z_dz(1, 1)
    
    ratio = 1.0 #1, 2, 0.6
    
    Psi_ini_all  = -2.0*z_matrix-2.0 #-1.0*np.ones((soil_layer_num,nrows*ncols)) #np.linspace(2.0,2.0,nz)
    Psi_n_all    = Psi_ini_all*1.0
    C_ini_all,K_ini_all,theta_ini_all = vanGenuchten(Psi_ini_all,alpha, theta_S, theta_R, n, m, Ksat)
    K_n_all       = K_ini_all*1.0
    theta_n_all    = theta_ini_all*1.0
    rain_depth = np.zeros(1*1)
     
    for i in xrange(time_steps): 
#        print i
        theta_n_all,Psi_n_all, rain_depth, Ts_n_all,Leak_n_all = \
            Soil_Moisture(rain[i+rp],theta_n_all, Psi_n_all, z_matrix/ratio,dz_matrix/ratio,rain_depth,1,1) #rain_depth

        theta_save[:,i] = theta_n_all[:,0]
        Psi_save[:,i]   = Psi_n_all[:,0]
        rain_depth_save[i] = rain_depth[0]
#        plt.plot(Psi_n_all[:,0],z_matrix[:,0])
#        
#    plt.gca().invert_yaxis()
#    
    theta = np.sum(theta_save*dz_matrix, axis=0)/np.sum(dz_matrix)
    print ratio 
    print np.mean(theta_save,1)
    
###### bar plot for rainfall
#    fig, ax = plt.subplots(figsize=(8,1.5 ), dpi=90)     
#    plt.rc("font", size=10)
#    plt.bar(np.arange(time_steps),1000*rain[rp:time_steps+rp])
#    ax.set_xlim([0,time_steps])
##    ax.set_ylim([0,50]) 
#    ax.set_xlabel('DOY', fontsize=12)
#    ax.set_ylabel('Rainfall [mm/day]', fontsize=12)
    
#    methods = [None, 'none', 'nearest', 'bilinear', 'bicubic', 'spline16', 'spline36', 'hanning', 'hamming', \
#           'hermite', 'kaiser', 'quadric', 'catrom', 'gaussian', 'bessel', 'mitchell', 'sinc', 'lanczos']
#   , interpolation='bicubic'
####### switch red and blue color for soil moisture
    fig, ax = plt.subplots(figsize=(11,2 ), dpi=90)     
    plt.rc("font", size=10)
    figplot  = ax.imshow(theta_save, extent=[0,6,1.0,0], interpolation='bicubic', cmap=cm.jet_r) #,vmin=0,vmax=0.032 #pcolor , matshow, imshow
    cbar = fig.colorbar(figplot,shrink=0.7, pad = 0.1)  #if shrink, shrink=0.9,
    cbar.set_label(r'Soil Moisture')
    ax.set_xticklabels([0, 50,100,150,200,2500,300])
    ax.set_xlabel('DOY', fontsize=12)
    ax.set_ylabel('Soil Depth [m]', fontsize=12)
    #    #ax.set_title(r'$\Delta$ Soil Depth [m]', fontsize=20)
    
    fig, ax = plt.subplots(figsize=(8,1.5 ), dpi=90)     
    plt.rc("font", size=10)
    plt.plot(np.arange(time_steps),theta)
    ax.set_ylim([0.08,poros])
    ax.set_xlim([0,time_steps])
    ax.set_xlabel('DOY', fontsize=12)
    ax.set_ylabel('Soil Moisture', fontsize=12)
    
    fig, ax = plt.subplots(figsize=(8,1.5 ), dpi=90)     
    plt.rc("font", size=10)
    plt.plot(np.arange(time_steps),theta_save[0,:])
    ax.set_ylim([0.08,poros])
    ax.set_xlim([0,time_steps])
    ax.set_xlabel('DOY', fontsize=12)
    ax.set_ylabel('Soil moisture', fontsize=12)
    
#    
#    plt.figure()
#    plt.plot(np.arange(time_steps),(theta_save[0,:]))
#    
#    plt.figure()
#    plt.plot(np.arange(time_steps),(theta_save[1,:]))
#    
#    plt.figure()
#    plt.plot(np.arange(time_steps),(theta_save[2,:]))
# 
#    plt.figure()
#    plt.plot(np.arange(time_steps),(theta_save[3,:]))
#    
#    plt.figure()
#    plt.plot(np.arange(time_steps),(theta_save[4,:]))
#    
#    plt.figure()
#    plt.plot(np.arange(time_steps),(theta_save[5,:]))
#
#    plt.figure()
#    plt.plot(np.arange(time_steps),(theta_save[6,:]))
