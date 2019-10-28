# -*- coding: utf-8 -*-
"""
The code follows the paper named 'Hydrological controls on soil carbon and nitrogen cycles. I & II in 2003'
"""


##=============================================##
##           C-N cycle Bucekt model            ##
##=============================================##

from __future__ import division
from parameters import *

import timeit
import matplotlib
import matplotlib.cm as cm
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import animation 
import pylab
from F_root_fraction import root_fr
from LEM_sediment_transport import linear_diffusion1D
import scipy.ndimage
import multiprocessing

from F_Soil_Moisture_main_pool_parallel import vanGenuchten, Soil_Moisture
from CN_layers_depth import ini_z_dz

def fd(s):     

    fds      = s/sfc
    ind      = np.where(s>sfc)

    fds[ind] = sfc/s[ind]
    return fds

def fn(s):
    fns      = s/sfc
    ind      = np.where(s>sfc)
    fns[ind] = (1.0-s[ind])/(1.0-sfc)
    return fns

#def T(s):      
#    Ts = (s-sw)/(s_star-sw)*T_max
#    ind1 = np.where(s<=sw)
#    ind2 = np.where(s>s_star)
#    Ts[ind1] = 0.0    
#    Ts[ind2] = T_max   
#    return Ts

#def L(s):
#    Ls = Ksat*(np.exp(beta*(s-sfc))-1.)/(np.exp(beta*(1.-sfc))-1.)
#    ind = np.where(s<sfc)
#    Ls[ind]=Ksat*s[ind]**(2*b+3)
#    return Ls
    
def UPa(DEM,UPp,ku,N):

    UPaN = DEM-UPp
    ind1 = np.where(DEM<=UPp)
    ind2 = np.where((ku*N)<=(DEM-UPp))
    UPaN[ind1] = 0.0
    UPaN[ind2] = ku[ind2]*N[ind2]
    return UPaN         
    
def phi_decide(PHI_big,IMM_max):
    phi       = np.ones(soil_layer_num)
    IMM       = np.zeros(soil_layer_num)
    MIN       = np.zeros(soil_layer_num)
    
    ind_p     = np.where(PHI_big >= 0.0)
    ind_n1    = np.where((PHI_big < 0.0) & (PHI_big >= -IMM_max))
    ind_n2    = np.where(PHI_big < -IMM_max)
        
    phi[ind_n2] = -IMM_max[ind_n2]/PHI_big[ind_n2]
    IMM[ind_n1] = -PHI_big[ind_n1]
    IMM[ind_n2] = IMM_max[ind_n2]
    MIN[ind_p]  = PHI_big[ind_p]
    
    return (phi,IMM,MIN)

def aveg_k(s,Cl,Ch,Cb,ADD,rh):
        
    s = np.array([s])
    A=np.zeros((3,3))
    B=np.zeros(3)
    A[0,0]=fd(s)*Cb*Cl
    A[0,1]=0.0
    A[0,2]=-Cb
    A[1][0]=rh*Cl*fd(s)*Cb
    A[1][1]=-Ch*fd(s)*Cb
    A[1][2]=0.0
    A[2][0]=(1.0-rh-rr)*fd(s)*Cl*Cb
    A[2][1]=(1.0-rr)*fd(s)*Ch*Cb
    A[2][2]=-Cb
    B[0]=ADD
    B[1]=0.0
    B[2]=0.0
    K=np.linalg.solve(A,B)
    
    return (K[0],K[1],K[2]) #(kl,kh,kd)


#==================the loop starts=================   

def f_SOM_parallel(params):
    
#    CN_h = CN_h_ini    
#    CN_b = CN_b_ini
  
    cell,Cl,Ch,Cb,s,Zr,dZr,Kl_cell,Kh_cell,Kd_cell,litter_above_cell,litter_below_cell = params

#    Frooti  = root_fr(Zr,dZr)
#    Addb    = Frooti*Add_b
    ##============================================##
    ##            bioturbation process          ##
    ##============================================## 

#    A1d_C = linear_diffusion1D(Zr,dZr)
    Cl1d   = np.zeros((soil_layer_num)) #A1d_C.dot(Cl)
    Ch1d   = np.zeros((soil_layer_num)) #A1d_C.dot(Ch)
    Cb1d   = np.zeros((soil_layer_num)) #A1d_C.dot(Cb)
    
#    Nl = Cl/CN_l
#    Nh = Ch/CN_h
#    Nb = Cb/CN_b
#
#    Nld = 0.0 #A1d_C.dot(Nl)
#    Nhd = 0.0 #A1d_C.dot(Nh)
#    Nbd = 0.0 #A1d_C.dot(Nb)
#    
    ##============================================##
    ##            biogeochemical process          ##
    ##============================================##   
    dCh=np.zeros((soil_layer_num))
    dCb=np.zeros((soil_layer_num))
    dCl=np.zeros((soil_layer_num))  
#    Nh_update=np.zeros((soil_layer_num))
#    Nb_update=np.zeros((soil_layer_num))
#    Nl_update=np.zeros((soil_layer_num)) 
#    N_plus_update = np.zeros((soil_layer_num)) 
#    N_minus_update = np.zeros((soil_layer_num))


#    rh_2col = np.column_stack((rh_ini*np.ones(soil_layer_num), CN_h/CN_l))
#    rh=np.min(rh_2col,axis =1)
    rh = 0.25*np.ones(soil_layer_num)
#        print rh_2col
    ##============================================##
    ##          To figure out phi & IMM           ##
    ##============================================##
    phi = 1.0
    fds = fd(s)
    DECl = phi*fds*Kl_cell*Cb*Cl
    DECh = phi*fds*Kh_cell*Cb*Ch
    
#    IMM_max=(ki_plus*N_plus+ki_minus*N_minus)*fds*Cb # paper typo, so delete Cb 
#    PHI_big=DECh*(1.0/CN_h-(1.0-rr)/CN_b)+DECl*(1.0/CN_l-rh/CN_h-(1.0-rr-rh)/CN_b)  
#    
#    phi,IMM,MIN=phi_decide(PHI_big,IMM_max)


#    IMM_plus  = ki_plus*N_plus*IMM/(ki_plus*N_plus+ki_minus*N_minus)
#    IMM_minus = ki_minus*N_minus*IMM/(ki_plus*N_plus+ki_minus*N_minus)

    ##============================================##
    ##              Litter Pool                   ##
    ##============================================##

    dCl = (litter_below_cell + Kd_cell*Cb - DECl)*dt #+Cl
#    Nl_update = (Addb/CN_addb + Kd_cell*Cb/CN_b-DECl/CN_l)*dt+Nl+ Nld
    dCl[0] = litter_above_cell*dt + dCl[0] #+Cl_update[0]
#    Nl_update[0]= ADD_s/CN_adds*dt+Nl_update[0]

    ##============================================##
    ##               Humus Pool                   ##
    ##============================================##

    dCh = (rh*DECl-DECh)*dt #Ch +
#    Nh_update = (rh*DECl/CN_h-DECh/CN_h)*dt+Nh + Nhd
    
    ##============================================##
    ##         Microbial Biomass Pool             ##
    ##============================================##
    dCb = ((1.-rh-rr)*DECl+(1.-rr)*DECh-Kd_cell*Cb)*dt  #Cb+
    Cb_update = dCb+Cb
    Cb_update[Cb_update<10.0] = 10.0
    dCb = Cb_update-Cb
#    Nb_update = (Cb_update)/CN_b + Nbd
    #Nb_update[x]=((1-rh*CN_l/CN_h)*DECl/CN_l+DECh/CN_h-kd*Cbs_update[x]/CN_b-PHI_big)*dt+Nb[x]

    
    ##============================================##
    ##     Mineralization & Immobilization        ##
    ##============================================##
#    ku_plus=a_plus*F*s**dd/s/poros/dZr
#    ku_minus=a_minus*F*s**dd/s/poros/dZr
#    LE_plus=a_plus*Ls/s/poros/dZr*N_plus
#    LE_minus=a_minus*Ls/s/poros/dZr*N_minus
#    UPp_plus=a_plus*Ts/s/poros/dZr*N_plus 
#    UPp_minus=a_minus*Ts/s/poros/dZr*N_minus
#    UPa_plus=UPa(DEM_plus,UPp_plus,ku_plus,N_plus)
#    UPa_minus=UPa(DEM_minus,UPp_minus,ku_minus,N_minus)
#    UP_plus=UPp_plus+UPa_plus
#    UP_minus=UPp_minus+UPa_minus
#    NIT=fn(s)*kn*N_plus*Cb

#    N_plus_update=(MIN-IMM_plus-NIT-LE_plus-UP_plus)*dt+N_plus
#
#    NegID_Np = np.where(N_plus_update<=0.001)
#    coeff_nu               = N_plus[NegID_Np]-0.001+(MIN[NegID_Np]-IMM_plus[NegID_Np])*dt
#    coeff_de               = (NIT[NegID_Np]+LE_plus[NegID_Np]+UP_plus[NegID_Np])*dt
#    NIT[NegID_Np]          = NIT[NegID_Np]*coeff_nu/coeff_de
#    LE_plus[NegID_Np]      = LE_plus[NegID_Np]*coeff_nu/coeff_de
#    UP_plus[NegID_Np]      = UP_plus[NegID_Np]*coeff_nu/coeff_de
#    N_plus_update[NegID_Np]= 0.001 # or same as (MIN-IMM_plus-NIT-LE_plus-UP_plus)*dt+N_plus[x] # should be zero
    
#    N_minus_update=(NIT-IMM_minus-LE_minus-UP_minus)*dt+N_minus  

#    NegID_Nm                 = np.where(N_minus_update<0.1)
#    coeff_nu                 = N_minus[NegID_Nm]-0.1- IMM_minus[NegID_Nm]*dt
#    coeff_de                 = (LE_minus[NegID_Nm]+UP_minus[NegID_Nm]-NIT[NegID_Nm])*dt
#    LE_minus[NegID_Nm]       = LE_minus[NegID_Nm]*coeff_nu/coeff_de
#    UP_minus[NegID_Nm]       = UP_minus[NegID_Nm]*coeff_nu/coeff_de
#    NIT[NegID_Nm]            = NIT[NegID_Nm]*coeff_nu/coeff_de
#    N_minus_update[NegID_Nm] = 0.1
#   
    
#    NegID_Nml                 = np.where(N_minus_update>1.2)
#    coeff_nu                 = N_minus[NegID_Nml]-1.2- IMM_minus[NegID_Nml]*dt
#    coeff_de                 = (LE_minus[NegID_Nml]+UP_minus[NegID_Nml]-NIT[NegID_Nml])*dt
#    LE_minus[NegID_Nml]       = LE_minus[NegID_Nml]*coeff_nu/coeff_de
#    UP_minus[NegID_Nml]       = UP_minus[NegID_Nml]*coeff_nu/coeff_de
#    NIT[NegID_Nml]            = NIT[NegID_Nml]*coeff_nu/coeff_de
#    N_minus_update[NegID_Nml] = 1.2
#    #N_minus_update = np.minimum(N_minus_update,np.ones(soil_layer_num)*1.4)
    
    ##============================================##
    ##             carbon mass balance            ##
    ##============================================##   

#        C_should_be_0=(Cl_update+Ch_update+Cb_update-Cl-Ch-Cb)/dt-(Addb-rr*DECl-rr*DECh) #-Cls_update-Chs_update-Cbs_update
#        C_should_be_0[0] = C_should_be_0[0] - ADD_s
#        print C_should_be_0
   
    ##============================================##
    ##           Nitrogen mass balance            ##
    ##============================================##   
 
#        N_should_be_0=(Nl_update+Nh_update+Nb_update+N_plus_update+N_minus_update-\
#                                (Nl+Nh+Nb+N_plus+N_minus))/dt-\
#                                (Addb/CN_addb-LE_plus-UP_plus-LE_minus-UP_minus)
#        print N_should_be_0
#==============================================================================
#         Save to plot
#==============================================================================
    
#    Cl_all[:,cell] = Cl_update
#    Ch_all[:,cell] = Ch_update
#    Cb_all[:,cell] = Cb_update
#    Nl_all[:,cell] = Nl_update
#    Nh_all[:,cell] = Nh_update
#    Nb_all[:,cell] = Nb_update
#    N_plus_all[:,cell] = N_plus_update
#    N_minus_all[:,cell] = N_minus_update
#    CN_l_all = Cl_all/Nl_all
    return dCl,dCh,dCb, Cl1d, Ch1d, Cb1d

#def multi_layer_CN(Cl_all,Ch_all,Cb_all,N_plus_all,N_minus_all,\
#                   CN_l_all,sm,z_matrix,dz_matrix,Kl,Kh,Kd,Ts_n_all,Leak_n_all,nrows,ncols): 

def multi_layer_CN(Cl_all,Ch_all,Cb_all,sm,z_matrix,dz_matrix,Kl,Kh,Kd,litter_above,litter_below,nrowsncols): 

    Cl_all_T    =  np.transpose(Cl_all)
    Ch_all_T   =  np.transpose(Ch_all)
    Cb_all_T =  np.transpose(Cb_all)

#    N_plus_all_T   =  np.transpose(N_plus_all)
#    N_minus_all_T  =  np.transpose(N_minus_all)
    
#    CN_l_all_T  =  np.transpose(CN_l_all)
    sm_T  =  np.transpose(sm)
    z_matrix_T    =  np.transpose(z_matrix)
    dz_matrix_T   =  np.transpose(dz_matrix)
    
#    Kl_2d = np.ones((soil_layer_num,nrows*ncols))*Kl[:,None]
#    Kh_2d = np.ones((soil_layer_num,nrows*ncols))*Kh[:,None]
#    Kd_2d = np.ones((soil_layer_num,nrows*ncols))*Kd[:,None]
    Kl_2d_T    =  np.transpose(Kl)
    Kh_2d_T   =  np.transpose(Kh)
    Kd_2d_T   =  np.transpose(Kd)

#    Ts_n_all_T    =  np.transpose(Ts_n_all)
#    Leak_n_all_T   =  np.transpose(Leak_n_all)

    litter_above_2dT = litter_above*1.0 
    litter_below_2dT = np.transpose(litter_below)
    params = zip(np.arange(nrowsncols),Cl_all_T,Ch_all_T,Cb_all_T,\
               sm_T,z_matrix_T,dz_matrix_T,Kl_2d_T,Kh_2d_T,Kd_2d_T,litter_above_2dT,litter_below_2dT)

    cores = multiprocessing.cpu_count()
    pool = multiprocessing.Pool(processes=cores)
    
    l = pool.map(f_SOM_parallel, params) 
    
    pool.close()
    pool.terminate()
    pool.join()
    l =  map(list, zip(*l))

    dCl_all          = np.asarray(l[0])
    dCh_all          = np.asarray(l[1])
    dCb_all          = np.asarray(l[2])
    Cl1d             = np.asarray(l[3])
    Ch1d             = np.asarray(l[4])
    Cb1d             = np.asarray(l[5])
#    N_plus_update_new    = np.asarray(l[3])
#    N_minus_update_new   = np.asarray(l[4])
#    CN_l_all_new         = np.asarray(l[5])
    
    dCl_all  = np.transpose(dCl_all)
    dCh_all  = np.transpose(dCh_all)
    dCb_all  = np.transpose(dCb_all)
    
    Cl1d_all = np.transpose(Cl1d)
    Ch1d_all = np.transpose(Ch1d)
    Cb1d_all = np.transpose(Cb1d)
    
#    N_plus_update_new  = np.transpose(N_plus_update_new)
#    N_minus_update_new  = np.transpose(N_minus_update_new)
#    CN_l_all_new  = np.transpose(CN_l_all_new)
    
    return dCl_all,dCh_all,dCb_all, Cl1d_all, Ch1d_all, Cb1d_all
   
if __name__ == '__main__':
    
    #...................Rainfall
    rain = np.load('rainfall_100yr3_daily.npy') # mm/(day)
    rain = rain/1000.0
    
    litter_input = (np.load('litter_input.npy'))/dt
    
    nrows = 1
    ncols = 1
    
    time_steps = 365*6
    rp = 0
    #........initial values, Soil moisture
    theta_save = np.zeros((soil_layer_num,time_steps))
    z_matrix, dz_matrix = ini_z_dz(nrows, ncols)
    dz_save = np.zeros((soil_layer_num,time_steps))
    
    Psi_ini_all  = -2.0*z_matrix-2.0 #-1.0*np.ones((soil_layer_num,nrows*ncols)) #np.linspace(2.0,2.0,nz)
    Psi_n_all    = Psi_ini_all*1.0
    C_ini_all,K_ini_all,theta_ini_all = vanGenuchten(Psi_ini_all,alpha, theta_S, theta_R, n, m, Ksat)
    
    ### case1: water ponding is preserved all the time
#    theta_mean_1col = np.array([ 0.29850789,  0.30059439,  0.30084879,  0.29989796,  0.29787837,
#        0.29459997,  0.28959179])
    
    ### case2: ZERO water ponding is preserved all the time
#    theta_mean_1col = np.array([ 0.25585338,  0.2628716 ,  0.26605831,  0.26727228,  0.26710455,
#        0.26546586,  0.26187034])
    theta_mean_1col = np.array([ 0.35480344, 0.36059945, 0.36269936, 0.36273304, 0.3611157 ,0.3574786 , 0.35072803])

    ### average all and ZERO water ponding is preserved all the time


    K_n_all       = K_ini_all*1.0
    theta_n_all    = theta_ini_all*1.0
    rain_depth = np.zeros(nrows*ncols)
     
    #...................initial values, Soil organic carbon
    Cl_save = np.zeros((soil_layer_num,time_steps))
    Ch_save = np.zeros((soil_layer_num,time_steps))
    Cb_save = np.zeros((soil_layer_num,time_steps))
    
    Ch_all_ini  = 43000.0*np.exp(-3.0*z_matrix) # 18000
    Cb_all_ini  = 180.0*np.exp(-3.0*z_matrix) #300
    Cl_all_ini  = 6900.0*np.exp(-3.0*z_matrix) #2800
    
    Ch_all      = Ch_all_ini*1.0
    Cb_all      = Cb_all_ini*1.0
    Cl_all      = Cl_all_ini*1.0
#    CN_l_all    = CN_l_ini*np.ones((soil_layer_num,1*1))
#    Nl_all      = Cl_all/CN_l_ini
#    Nh_all      = Ch_all/CN_h_ini
#    Nb_all      = Cb_all/CN_b_ini
#    N_plus_all  = 0.002*np.exp(-1.0*z_matrix)
#    N_minus_all = 0.6*np.exp(-1.0*z_matrix)
    Kl          = np.zeros((soil_layer_num))
    Kh          = np.zeros((soil_layer_num))
    Kd          = np.zeros((soil_layer_num))

    Zr=z_matrix[:,0]  
    dZr=dz_matrix[:,0]
    
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

    start = timeit.default_timer()
    for i in xrange(time_steps):
#        print i
#        print rain[i] ,rain_depth
        #rain[i]
        theta_n_all,Psi_n_all,rain_depth,Ts_n_all,Leak_n_all = \
        Soil_Moisture(rain[i+rp],theta_n_all, Psi_n_all, z_matrix,dz_matrix,rain_depth*0.0,nrows, ncols)


        theta_save[:,i]= theta_n_all[:,0]
        dz_save[:,i] = dz_matrix[:,0]
               
        litter_above = ADD_s + litter_input[0,i%(365*2)] #litter_input + ADD_s
        litter_below =(litter_input[1,i%(365*2)] + Add_b)*Frooti
        litter_below = np.ones((soil_layer_num,nrows*ncols))*litter_below[:,None]


        dCl_all,dCh_all,dCb_all, Cl1d_all, Ch1d_all, Cb1d_all \
             = multi_layer_CN(Cl_all,Ch_all,Cb_all,theta_n_all/poros,z_matrix,dz_matrix,Kl,Kh,Kd,litter_above,litter_below,nrows,ncols)
        
        Cl_all = Cl_all*1.0 + dCl_all
        Ch_all = Ch_all*1.0 + dCh_all
        Cb_all = Cb_all*1.0 + dCb_all
        
        Cl_save[:,i]= Cl_all[:,0]
        Ch_save[:,i]= Ch_all[:,0]
        Cb_save[:,i]= Cb_all[:,0]
    
    stop = timeit.default_timer()-start
    print stop
    
#        plt.plot(Psi_n_all[:,0],z_matrix[:,0])
#        
#    plt.gca().invert_yaxis()
#    
    theta = np.sum(theta_save*dz_save, axis=0)/np.sum(dz_matrix[:,0])
    Ccl   = np.sum(Cl_save*dz_save, axis=0)/np.sum(dz_matrix[:,0])
    Cch   = np.sum(Ch_save*dz_save, axis=0)/np.sum(dz_matrix[:,0])
    Ccb   = np.sum(Cb_save*dz_save, axis=0)/np.sum(dz_matrix[:,0])
    
    soil_moisture = False
    if soil_moisture:
        fig, ax = plt.subplots(figsize=(8,1.5 ), dpi=90)     
        plt.rc("font", size=6)
        plt.bar(np.arange(time_steps),1000*rain[rp:time_steps+rp])
        ax.set_xlim([0,time_steps])
#        ax.set_ylim([0,0.08]) 
        ax.set_xlabel('DOY', fontsize=12)
        ax.set_ylabel('Rainfall [mm/day]', fontsize=12)
        
    #    methods = [None, 'none', 'nearest', 'bilinear', 'bicubic', 'spline16', 'spline36', 'hanning', 'hamming', \
    #           'hermite', 'kaiser', 'quadric', 'catrom', 'gaussian', 'bessel', 'mitchell', 'sinc', 'lanczos']
    
        fig, ax = plt.subplots(figsize=(11,2 ), dpi=90)     
        plt.rc("font", size=10)
        figplot  = ax.imshow(theta_save, extent=[0,7,1.0,0], interpolation='bicubic', cmap=cm.jet_r) #,vmin=0,vmax=0.032 #pcolor , matshow, imshow
        cbar = fig.colorbar(figplot,shrink=0.7, pad = 0.1)  #if shrink, shrink=0.9,
        cbar.set_label(r'Soil moisture')
        ax.set_xticklabels([0, 50,100,150,200,250,300,350])
        ax.set_xlabel('DOY', fontsize=12)
        ax.set_ylabel('Soil Depth', fontsize=12)
        #    #ax.set_title(r'$\Delta$ Soil Depth [m]', fontsize=20)
        
        fig, ax = plt.subplots(figsize=(8,1.5 ), dpi=90)     
        plt.rc("font", size=6)
        plt.plot(np.arange(time_steps),theta)
        ax.set_ylim([0.08,0.4])
        ax.set_xlim([0,time_steps])
        ax.set_xlabel('DOY', fontsize=12)
        ax.set_ylabel('Soil moisture', fontsize=12)
    
    SOC = True
    if SOC:
        fig, ax = plt.subplots(figsize=(8,1.5 ), dpi=90)     
        plt.rc("font", size=12)
        ax.plot(np.arange(time_steps)/365,Ccl/1000.0)
        ax.set_ylabel(r'Cl [$Kg \ C m^{-3}$]', fontsize=12)
        ax.set_xlim([0,100])
#        fig.savefig("Cl_demo.png" )
        
        fig, ax = plt.subplots(figsize=(8,1.5 ), dpi=90)     
        plt.rc("font", size=12)
        ax.plot(np.arange(time_steps)/365,Cch/1000.0)
        ax.set_ylabel(r'Ch [$Kg \ C m^{-3}$]', fontsize=12)
        ax.set_xlim([0,100])
#        fig.savefig("Ch_demo.png" )
        
        fig, ax = plt.subplots(figsize=(8,1.5 ), dpi=90)     
        plt.rc("font", size=12)
        ax.plot(np.arange(time_steps)/365,Ccb/1000.0)
        ax.set_ylabel(r'Cb [$Kg \ C m^{-3}$]', fontsize=12)
        ax.set_xlim([0,100])
#        fig.savefig("Cb_demo.png" )
        
        fig, ax = plt.subplots(figsize=(8,1.5 ), dpi=90)     
        plt.rc("font", size=12)
        ax.plot(np.arange(time_steps)/365,(Ccb+Cch+Ccl)/1000.0)
        ax.set_ylabel(r'SOC [$Kg \ C m^{-3}$]', fontsize=12)
        ax.set_xlim([0,100])
#        fig.savefig("SOC_demo.png" )
        
#        fig, ax = plt.subplots(figsize=(8,1.5 ), dpi=90)     
#        plt.rc("font", size=6)
#        plt.plot(np.arange(time_steps),Cl_save[0,:])
#        ax.set_ylabel('Cl', fontsize=12)
        
    ### one column over time
    if False:
        fig, ax = plt.subplots(figsize=(11,2 ), dpi=90)     
        plt.rc("font", size=10)
        zi =  Cl_save/1000.0 #scipy.ndimage.zoom(Cl_save/1000.0, 1)
        xi = np.arange(0,time_steps)
        yi = z_matrix[:,0]
        xi, yi = np.meshgrid(xi, yi)
        CS = ax.contourf(xi,yi,zi,30, interpolation='bicubic', cmap=plt.cm.jet)
        fig.colorbar(CS, shrink=0.8 ) # draw colorbar #shrink=0.7, pad = 0.5
        ax.invert_yaxis()
        
#        figplot  = ax.imshow(Cl_save/1000.0, extent=[0,5,1.0,0], interpolation='bicubic', cmap=cm.jet) #,vmin=0,vmax=0.032 #pcolor , matshow, imshow
#        cbar = fig.colorbar(figplot,shrink=0.7, pad = 0.1)  #if shrink, shrink=0.9,
#        cbar.set_label(r'Cl [$\frac{kg \ C}{m^3}$]')
#        ax.set_xticklabels([0,20,40,60, 80,100])
#        ax.set_xlabel('Yr', fontsize=12)
#        ax.set_ylabel('Soil Depth', fontsize=12)
#        
#        fig, ax = plt.subplots(figsize=(8,1.5 ), dpi=90)     
#        plt.rc("font", size=6)
#        plt.plot(np.arange(time_steps),Cl_save[0,:])
#        ax.set_ylabel('Cl', fontsize=12)
        
    ### vertical plot
    if False:
        fig, ax = plt.subplots(figsize=(4,6)) 
        ax.plot(Cl_save[:,0]/1000.0,z_matrix[:,0],'o-', label='initial')
        ax.plot(Cl_save[:,36499]/1000.0,z_matrix[:,0], 'o-',label='final')
        ax.set_xlabel(r'Cl [$KgC m^{-3}$]')
        ax.set_ylabel(r'Soil Depth [m]')
        ax.invert_yaxis()
        ax.legend()

    
    

    
    
    
    
    
    
    