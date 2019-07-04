# -*- coding: utf-8 -*-
"""
The code follows the paper named 'Hydrological controls on soil carbon and nitrogen cycles. I & II in 2003'
"""


##=============================================##
##           C-N cycle Bucekt model            ##
##=============================================##

from __future__ import division
from parameters import *

import time
import matplotlib
import matplotlib.cm as cm
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import animation 
import pylab
from F_root_fraction import root_fr
from LEM_sediment_transport import linear_diffusion1D
import scipy.ndimage

from F_Soil_Moisture_main import vanGenuchten, Soil_Moisture
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

    
#def UPa(DEM,UPp,ku,N):
#
#    UPaN = DEM-UPp
#    ind1 = np.where(DEM<=UPp)
#    ind2 = np.where((ku*N)<=(DEM-UPp))
#    UPaN[ind1] = 0.0
#    UPaN[ind2] = ku[ind2]*N[ind2]
#    return UPaN         
#    
#def phi_decide(PHI_big,IMM_max):
#    phi       = np.ones(soil_layer_num)
#    IMM       = np.zeros(soil_layer_num)
#    MIN       = np.zeros(soil_layer_num)
#    
#    ind_p     = np.where(PHI_big >= 0.0)
#    ind_n1    = np.where((PHI_big < 0.0) & (PHI_big >= -IMM_max))
#    ind_n2    = np.where(PHI_big < -IMM_max)
#        
#    phi[ind_n2] = -IMM_max[ind_n2]/PHI_big[ind_n2]
#    IMM[ind_n1] = -PHI_big[ind_n1]
#    IMM[ind_n2] = IMM_max[ind_n2]
#    MIN[ind_p]  = PHI_big[ind_p]
#    
#    return (phi,IMM,MIN)

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

def multi_layer_CN(Cl_all,Ch_all,Cb_all,sm,z_matrix,dz_matrix,Kl,Kh,Kd,litter_above,litter_below,nrows,ncols): 

    dCl_all = np.zeros((soil_layer_num,(ncols)*(nrows)))
    dCh_all = np.zeros((soil_layer_num,(ncols)*(nrows)))
    dCb_all = np.zeros((soil_layer_num,(ncols)*(nrows)))
#==================the loop starts=================    
    for cell in xrange(0,(ncols)*(nrows)): #(ncols)*(nrows)
        Zr      = z_matrix[:,cell]  
        dZr     = dz_matrix[:,cell]
#        Frooti  = root_fr(Zr,dZr)
#        Addb    = Frooti*Add_b
        Cl      = Cl_all[:,cell]
        Ch      = Ch_all[:,cell]
        Cb      = Cb_all[:,cell]

        s       = sm[:,cell]

        ##============================================##
        ##           bioturbation diffusion           ##
        ##============================================##
#        A1d_C = linear_diffusion1D(Zr,dZr)

        Cl1d   = np.zeros((soil_layer_num)) #A1d_C.dot(Cl)
        Ch1d   = np.zeros((soil_layer_num)) #A1d_C.dot(Ch)
        Cb1d   = np.zeros((soil_layer_num)) #A1d_C.dot(Cb)

        ##============================================##
        ##            biogeochemical process          ##
        ##============================================##   
        dCh=np.zeros((soil_layer_num))
        dCb=np.zeros((soil_layer_num))
        dCl=np.zeros((soil_layer_num))  
        
        rh = 0.25*np.ones(soil_layer_num)
        ##============================================##
        ##          To figure out phi & IMM           ##
        ##============================================##
        phi = 1.0
        fds = fd(s)
        DECl = phi*fds*Kl*Cb*Cl
        DECh = phi*fds*Kh*Cb*Ch

        ##============================================##
        ##              Litter Pool                   ##
        ##============================================##

#        Cl_update = (litter_below + Kd*Cb - DECl)*dt
        dCl    = (litter_below + Kd*Cb - DECl)*dt 
#        Nl_update = (Addb/CN_addb + Kd*Cb/CN_b-DECl/CN_l)*dt+Nl+ Nld
#        Cl_update[0] = litter_above*dt+Cl_update[0]
        dCl[0] = litter_above*dt + dCl[0]
#        Nl_update[0]= ADD_s/CN_adds*dt+Nl_update[0]

        ##============================================##
        ##               Humus Pool                   ##
        ##============================================##

#        Ch_update = (rh*DECl-DECh)*dt+Ch + Chd
        dCh       = (rh*DECl-DECh)*dt
#        Nh_update = (rh*DECl/CN_h-DECh/CN_h)*dt+Nh + Nhd
        ##============================================##
        ##         Microbial Biomass Pool             ##
        ##============================================##
#        Cb_update = ((1.-rh-rr)*DECl+(1.-rr)*DECh-Kd*Cb)*dt+Cb+Cbd
#        Cb_update[Cb_update<10.0] = 10.0
        
        dCb       = ((1.-rh-rr)*DECl+(1.-rr)*DECh-Kd*Cb)*dt
        Cb_update = dCb+Cb
        Cb_update[Cb_update<10.0] = 10.0
        dCb = Cb_update-Cb
#        Cb_update[Cb_update>360.0] = 360.0
#        Nb_update = (Cb_update)/CN_b + Nbd
            
        ##============================================##
        ##             carbon mass balance            ##
        ##============================================##   

#        C_should_be_0=(Cl_update+Ch_update+Cb_update-Cl-Ch-Cb)/dt-(Addb-rr*DECl-rr*DECh) #-Cls_update-Chs_update-Cbs_update
#        C_should_be_0[0] = C_should_be_0[0] - ADD_s
#        print C_should_be_0

#==============================================================================
#         Save the who domain
#==============================================================================
        
        dCl_all[:,cell] = dCl
        dCh_all[:,cell] = dCh
        dCb_all[:,cell] = dCb

    return dCl_all,dCh_all,dCb_all, Cl1d, Ch1d, Cb1d
   
if __name__ == '__main__':
    
    #...................Rainfall
    rain = np.load('rainfall_100yr3_daily.npy') # mm/(day)
    rain = rain/1000.0
         
    litter_input = (np.load('litter_input.npy'))/dt
    
    time_steps = 1
    rp = 0
    #........initial values, Soil moisture
    theta_save = np.zeros((soil_layer_num,time_steps))
    z_matrix, dz_matrix = ini_z_dz(1, 1)
    
    Psi_ini_all  = -2.0*z_matrix-2.0 #-1.0*np.ones((soil_layer_num,nrows*ncols)) #np.linspace(2.0,2.0,nz)
    Psi_n_all    = Psi_ini_all*1.0
    C_ini_all,K_ini_all,theta_ini_all = vanGenuchten(Psi_ini_all,alpha, theta_S, theta_R, n, m, Ksat)
        
#    theta_mean_1col = np.array([0.30887538, 0.31472939, 0.31682785, 0.31680671, 0.31508827, 0.31146861, 0.30576922])
    #dry
    theta_mean_1col = np.array([0.35480344, 0.36059945, 0.36269936, 0.36273304, 0.3611157 ,0.3574786 , 0.35072803])
#    wet
#    theta_mean_1col = np.array([0.39126216, 0.39281667, 0.3923617,  0.39047139, 0.38715473, 0.38192914, 0.37363064])

    K_n_all       = K_ini_all*1.0
    theta_n_all    = theta_ini_all*1.0
    rain_depth = np.zeros(1*1)
    
    #...................initial values, Soil organic carbon
    Cl_save = np.zeros((soil_layer_num,time_steps))
    Ch_save = np.zeros((soil_layer_num,time_steps))
    Cb_save = np.zeros((soil_layer_num,time_steps))
    
    Ch_all_ini  = 43000.0*np.exp(-2.5*z_matrix) #45000, 43000 18000
    
    Cb_all_ini  = 180.0*np.exp(-2.5*z_matrix) #220  360
    
    Cl_all_ini  = 6900.0*np.exp(-2.5*z_matrix) # 6000 7000, 5600  2800

    Ch_all      = Ch_all_ini*1.0
    Cb_all      = Cb_all_ini*1.0
    Cl_all      = Cl_all_ini*1.0
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
    
#    Kl = np.array([5.50656191e-07,   2.00757751e-08,   3.86591174e-08,
#         4.97673891e-08,   4.98542976e-08,   4.62412079e-08,
#         5.28699819e-08])
#    
#    Kh = np.array([3.23695181e-08,   1.03736173e-09,   2.33081836e-09,
#         3.63872807e-09,   4.71380568e-09,   6.42975964e-09,
#         1.58989235e-08])
#    
#    Kd = np.array([0.00119767,   4.98163427e-05,   8.77443161e-05,
#         1.02108219e-04,   9.00033998e-05,   6.92643630e-05,
#         5.45903284e-05])

    for i in xrange(time_steps): 
#        print rain[i] ,rain_depth
        #rain[i]
        
        theta_n_all,Psi_n_all,rain_depth,Ts_n_all,Leak_n_all = \
        Soil_Moisture(rain[i+rp],theta_n_all, Psi_n_all, z_matrix,dz_matrix,rain_depth,1,1)

        theta_save[:,i]= theta_n_all[:,0]
        
        litter_above = ADD_s + litter_input[0,i%(365*2)] #litter_input + ADD_s
        litter_below =(litter_input[1,i%(365*2)] + Add_b)*Frooti
        
#        Cl_all,Ch_all,Cb_all \
        dCl_all,dCh_all,dCb_all, Cl1d_all, Ch1d_all, Cb1d_all \
             = multi_layer_CN(Cl_all,Ch_all,Cb_all,theta_n_all/poros,z_matrix,dz_matrix,Kl,Kh,Kd,litter_above,litter_below,1,1)
             
        Cl_all = Cl_all*1.0 + dCl_all
        Ch_all = Ch_all*1.0 + dCh_all
        Cb_all = Cb_all*1.0 + dCb_all
        
        Cl_save[:,i]= Cl_all[:,0]
        Ch_save[:,i]= Ch_all[:,0]
        Cb_save[:,i]= Cb_all[:,0]
        

        
#        plt.plot(Psi_n_all[:,0],z_matrix[:,0])
#        
#    plt.gca().invert_yaxis()
#    
    theta = np.sum(theta_save*dz_matrix, axis=0)/np.sum(dz_matrix)
    Ccl   = np.sum(Cl_save*dz_matrix, axis=0)/np.sum(dz_matrix)
    Cch   = np.sum(Ch_save*dz_matrix, axis=0)/np.sum(dz_matrix)
    Ccb   = np.sum(Cb_save*dz_matrix, axis=0)/np.sum(dz_matrix)
    
    soil_moisture = False
    if soil_moisture:
        fig, ax = plt.subplots(figsize=(8,1.5 ), dpi=90)     
        plt.rc("font", size=12)
        plt.bar(np.arange(time_steps),1000*rain[rp:time_steps+rp])
        ax.set_xlim([0,time_steps])
#        ax.set_ylim([0,0.08]) 
        ax.set_xlabel('DOY', fontsize=12)
        ax.set_ylabel('Rainfall [mm/day]', fontsize=12)
        
    #    methods = [None, 'none', 'nearest', 'bilinear', 'bicubic', 'spline16', 'spline36', 'hanning', 'hamming', \
    #           'hermite', 'kaiser', 'quadric', 'catrom', 'gaussian', 'bessel', 'mitchell', 'sinc', 'lanczos']
    
        fig, ax = plt.subplots(figsize=(11,2 ), dpi=90)     
        plt.rc("font", size=12)
        figplot  = ax.imshow(theta_save, extent=[0,7,1.0,0], interpolation='bicubic', cmap=cm.jet_r) #,vmin=0,vmax=0.032 #pcolor , matshow, imshow
        cbar = fig.colorbar(figplot,shrink=0.7, pad = 0.1)  #if shrink, shrink=0.9,
        cbar.set_label(r'Soil moisture')
        ax.set_xticklabels([0, 50,100,150,200,250,300,350])
        ax.set_xlabel('DOY', fontsize=12)
        ax.set_ylabel('Soil Depth', fontsize=12)
        #    #ax.set_title(r'$\Delta$ Soil Depth [m]', fontsize=20)
        
        fig, ax = plt.subplots(figsize=(8,1.5 ), dpi=90)     
        plt.rc("font", size=12)
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
        plt.rc("font", size=12)
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
    if True:
        fig, ax = plt.subplots(figsize=(4,6)) 
        plt.rc("font", size=12)
        ax.plot((Cl_save[:,0]+Ch_save[:,0]+Cb_save[:,0])/1000.0,z_matrix[:,0],'o-', label='initial')
        ax.plot((Cl_save[:,time_steps-1]+Ch_save[:,time_steps-1]+Cb_save[:,time_steps-1])/1000.0,z_matrix[:,0], 'o-',label='final')
        ax.set_xlabel(r'SOC [$KgC m^{-3}$]', fontsize=12)
        ax.set_ylabel(r'Soil Depth [m]', fontsize=12)
        ax.invert_yaxis()
        ax.legend()

    
    

    
    
    
#    fig, ax = plt.subplots(figsize=(4,6)) 
#    ax.plot((Cl_all_ini+Ch_all_ini+Cb_all_ini)/1000.0,z_matrix[:,0],'o-', label='initial')
#    ax.invert_yaxis()
    