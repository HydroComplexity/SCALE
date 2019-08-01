# -*- coding: utf-8 -*-
"""
Created on Mon Sep 26 21:53:32 2016

@author: a
"""

from parameters import *
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
import multiprocessing

# dz is the soil layer thickness
# z is soil depth for each layer
# h is the total soil depth, and dh is the soil depth change
# the logic is that: sum(dz)+La = z; max(z)=h

#==============================================================================
# initial dz andz
#==============================================================================
#soil_layer_num= 5
## initial dz and z
#m = np.linspace(0.3,0.9,soil_layer_num)
#dzz= -0.1/(m-1)
def ini_z_dz(nrows, ncols):
    dz = (initial_soil_depth-La)*dzz/np.sum(dzz) # delta z, soil layer thickness
    z = np.zeros(soil_layer_num) #soil_layer_depth
    z[0]=dz[0]+La
    for i in xrange(1,soil_layer_num):
        z[i]=dz[i]+z[i-1]

#plt.plot(z,dz,'o')

    z_matrix = np.zeros((soil_layer_num,nrows*ncols))
    dz_matrix = np.zeros((soil_layer_num,nrows*ncols))
    for i in xrange(nrows*ncols):
        z_matrix[:,i]=z
        dz_matrix[:,i]=dz
    return z_matrix,dz_matrix

def ini_z_dz_saved(nrows, ncols):
    
    zOS = np.arange(soil_layer_num)+1
    zOS = zOS.reshape(soil_layer_num,1)
    dz_matrix = dz*np.ones((soil_layer_num,nrows*ncols))
    z_matrix =  np.kron( np.ones(nrows*ncols) , dz*zOS)
    return dz_matrix,z_matrix
    
#==============================================================================
# update dz and z
#==============================================================================
def z_dz_update(z_matrix,dh_soil,nrows, ncols): #(i_t,eta_vector_update,eta_bedrock_vector,z_matrix, C_matrix,Ca_vector):

    soil_depth_vector = z_matrix[soil_layer_num-1,:] +dh_soil
#    print soil_depth_vector-z_matrix[soil_layer_num-1,:]
    z_matrix_update = np.zeros((soil_layer_num,nrows*ncols))
#    dz_matrix_update = np.zeros((soil_layer_num,nrows*ncols))
#    
    dzz_matrix = np.kron( np.ones(nrows*ncols) , dzz.reshape(soil_layer_num,1))

    dz_matrix_update = (soil_depth_vector-La)*dzz_matrix/np.sum(dzz)
#    
    z_matrix_update[0,:] = dz_matrix_update[0,:]+La
    for i in xrange(1,soil_layer_num):
        z_matrix_update[i,:] = dz_matrix_update[i,:] + z_matrix_update[i-1,:]
#        
##    print z_matrix[soil_layer_num-1,:]-z_matrix_update[soil_layer_num-1,:]
#    # interpolate C_matrix based on z_matrix and z_update
#    # use cubic spline interpolation method
#    z_matrix[soil_layer_num-1,:] = z_matrix_update[soil_layer_num-1,:]
#    z_matrix_interp = np.vstack((La*np.ones(nrows*ncols), z_matrix))
#    C_matrix_interp = np.vstack((Ca_vector, C_matrix))
#    
#    C_matrix_update = np.zeros((soil_layer_num,nrows*ncols))
#    for i in xrange(nrows*ncols):
#
#        C_matrix_update[:,i] = interp1d(z_matrix_interp[:,i], C_matrix_interp[:,i], kind='cubic',bounds_error=False)(z_matrix_update[:,i])
##        if np.min(C_matrix_update[:,i])<0:
##            print i

    return z_matrix_update, dz_matrix_update#,C_matrix_update
    

def z_dz_update_saved_not_finished(z_matrix,dh_soil,nrows, ncols): #(i_t,eta_vector_update,eta_bedrock_vector,z_matrix, C_matrix,Ca_vector):

    
    soil_depth_vector = z_matrix[soil_layer_num-1,:] +dh_soil
    
    dz_2D_vector = soil_depth_vector/soil_layer_num
   
    dz_matrix = dz_2D_vector
    z_matrix_update = dz_matrix

    return z_matrix_update, dz_matrix#,C_matrix_update
    

        
def C_inter_pool_parallel(params):
    
    
    cell,Z,Z_update,C_cell,dh_ST_cell = params
    C_update_cell = interp1d( Z+dh_ST_cell, C_cell,fill_value='extrapolate')(Z_update) # 

    return C_update_cell

def C_interpolate(Cl,Ch,Cb, z_matrix,z_matrix_update,dh_ST,nrows, ncols):
    
#    dh_ST = z_matrix_update[soil_layer_num-1,:] - z_matrix[soil_layer_num-1,:]

#    C_update = np.zeros((soil_layer_num,nrows*ncols))
    C = (Cl+Ch+Cb)*1.0
    ratio_l = Cl/C
    ratio_h = Ch/C
    ratio_b = Cb/C

    C_T = np.transpose(C)
    z_matrix_T = np.transpose(z_matrix)
    z_matrix_update_T = np.transpose(z_matrix_update)
    params = zip(np.arange(nrows*ncols),z_matrix_T,z_matrix_update_T,C_T,dh_ST)
    
    cores = multiprocessing.cpu_count()
    pool = multiprocessing.Pool(processes=cores)
    
    l = pool.map(C_inter_pool_parallel, params) 
    pool.close()
    pool.terminate()
    pool.join()
    
#    l =  map(list, zip(*l))

    C_update = np.asarray(l)
#    print C_update.shape
    C_update = np.transpose(C_update)

    Cl_update = ratio_l*C_update
    Ch_update = ratio_h*C_update
    Cb_update = ratio_b*C_update
    
    return Cl_update,Ch_update,Cb_update
        
        
def C_interpolate_simple(Cl,z_matrix,z_matrix_update,dh_ST,nrows, ncols):
    
#    dh_ST = z_matrix_update[soil_layer_num-1,:] - z_matrix[soil_layer_num-1,:]
    Cl_update = np.zeros((soil_layer_num,nrows*ncols))

    for i in xrange(nrows*ncols):


        Cl_update[:,i] = interp1d((z_matrix[:,i])+dh_ST[i], Cl[:,i],fill_value='extrapolate')(z_matrix_update[:,i]) #
        
    return Cl_update
        
        
        
        
def K_intep_pool_parallel(params):
    
    
    cell,Z,Z_update,Kl_cell,Kh_cell,Kd_cell,dh_ST_cell = params
    Kl_update_cell = interp1d( Z-Z[0], Kl_cell,fill_value='extrapolate')(Z_update-Z_update[0]) # 
    Kh_update_cell = interp1d( Z-Z[0], Kh_cell,fill_value='extrapolate')(Z_update-Z_update[0]) # 
    Kd_update_cell = interp1d( Z-Z[0], Kd_cell,fill_value='extrapolate')(Z_update-Z_update[0]) # 

    return Kl_update_cell,Kh_update_cell,Kd_update_cell

        
def Klhd_interpolate(Kl, Kh, Kd, z_matrix,z_matrix_update,dh_ST,nrows, ncols):

    z_matrix_T = np.transpose(z_matrix)
    z_matrix_update_T = np.transpose(z_matrix_update)
    Kl_T = np.transpose(Kl)
    Kh_T = np.transpose(Kh)
    Kd_T = np.transpose(Kd)
    
    params = zip(np.arange(nrows*ncols),z_matrix_T,z_matrix_update_T,Kl_T,Kh_T, Kd_T,dh_ST)
    
    cores = multiprocessing.cpu_count()
    pool = multiprocessing.Pool(processes=cores)
    
    l = pool.map(K_intep_pool_parallel, params) 
    pool.close()
    pool.terminate()
    pool.join()
#    Kl_update = np.zeros((soil_layer_num,nrows*ncols))
#    Kh_update = np.zeros((soil_layer_num,nrows*ncols))
#    Kd_update = np.zeros((soil_layer_num,nrows*ncols))

#    for i in xrange(nrows*ncols):
#        
#        Kl_update[:,i] = interp1d((z_matrix[:,i]-z_matrix[0,i]), Kl[:,i],fill_value='extrapolate')(z_matrix_update[:,i]-z_matrix_update[0,i]) #  (z_matrix_update-z_matrix_update[0])
#        Kh_update[:,i] = interp1d((z_matrix[:,i]-z_matrix[0,i]), Kh[:,i],fill_value='extrapolate')(z_matrix_update[:,i]-z_matrix_update[0,i]) # 
#        Kd_update[:,i] = interp1d((z_matrix[:,i]-z_matrix[0,i]), Kd[:,i],fill_value='extrapolate')(z_matrix_update[:,i]-z_matrix_update[0,i]) # 

    l         =  map(list, zip(*l))
    Kl_update = np.asarray(l[0])
    Kh_update = np.asarray(l[1])
    Kd_update = np.asarray(l[2])

    Kl_update = np.transpose(Kl_update)
    Kh_update = np.transpose(Kh_update)
    Kd_update = np.transpose(Kd_update)

    return Kl_update,Kh_update,Kd_update
        
        
def C_reslice(Cl,Ch,Cb, z_matrix,z_matrix_update,dh_ST,nrows, ncols):
    
    Cl_update = np.zeros((soil_layer_num,nrows*ncols))
    Ch_update = np.zeros((soil_layer_num,nrows*ncols))
    Cb_update = np.zeros((soil_layer_num,nrows*ncols))
    
    #deposition site Index
    p_dh_ind = np.where(dh_ST>=0.0)
    #erosion site Index
    n_dh_ind = np.where(dh_ST<0.0)
    
    #-----at the deposition site
    Cl_update[0,p_dh_ind[0]] =  Cl[0,p_dh_ind[0]]*1.0
    Ch_update[0,p_dh_ind[0]] =  Ch[0,p_dh_ind[0]]*1.0
    Cb_update[0,p_dh_ind[0]] =  Cb[0,p_dh_ind[0]]*1.0
    
    dz_update_p = (z_matrix_update[1:soil_layer_num,p_dh_ind[0]]-z_matrix_update[0:soil_layer_num-1,p_dh_ind[0]])
    zi_p        = z_matrix_update[1:soil_layer_num,p_dh_ind[0]] - (z_matrix[0:soil_layer_num-1,p_dh_ind[0]]+dh_ST[p_dh_ind[0]])
    zim1_p      = (z_matrix[0:soil_layer_num-1,p_dh_ind[0]]+dh_ST[p_dh_ind[0]])-z_matrix_update[0:soil_layer_num-1,p_dh_ind[0]]
    
    Cl_update[1:soil_layer_num,p_dh_ind[0]] = (Cl[1:soil_layer_num,p_dh_ind[0]]*zi_p + Cl[0:soil_layer_num-1,p_dh_ind[0]]*zim1_p )/dz_update_p
    Ch_update[1:soil_layer_num,p_dh_ind[0]] = (Ch[1:soil_layer_num,p_dh_ind[0]]*zi_p + Ch[0:soil_layer_num-1,p_dh_ind[0]]*zim1_p )/dz_update_p
    Cb_update[1:soil_layer_num,p_dh_ind[0]] = (Cb[1:soil_layer_num,p_dh_ind[0]]*zi_p + Cb[0:soil_layer_num-1,p_dh_ind[0]]*zim1_p )/dz_update_p
    
    #-----at the erosion site
    Cl_update[0,n_dh_ind[0]]                = (Cl[0,n_dh_ind[0]]*(z_matrix[0,n_dh_ind[0]]+dh_ST[n_dh_ind[0]]) + Cl[1,n_dh_ind[0]]*(z_matrix_update[0,n_dh_ind[0]]-(z_matrix[0,n_dh_ind[0]]+dh_ST[n_dh_ind[0]]) ) )/z_matrix_update[0,n_dh_ind[0]]
    Ch_update[0,n_dh_ind[0]]                = (Ch[0,n_dh_ind[0]]*(z_matrix[0,n_dh_ind[0]]+dh_ST[n_dh_ind[0]]) + Ch[1,n_dh_ind[0]]*(z_matrix_update[0,n_dh_ind[0]]-(z_matrix[0,n_dh_ind[0]]+dh_ST[n_dh_ind[0]]) ) )/z_matrix_update[0,n_dh_ind[0]]
    Cb_update[0,n_dh_ind[0]]                = (Cb[0,n_dh_ind[0]]*(z_matrix[0,n_dh_ind[0]]+dh_ST[n_dh_ind[0]]) + Cb[1,n_dh_ind[0]]*(z_matrix_update[0,n_dh_ind[0]]-(z_matrix[0,n_dh_ind[0]]+dh_ST[n_dh_ind[0]]) ) )/z_matrix_update[0,n_dh_ind[0]]

    Cl_update[soil_layer_num-1,n_dh_ind[0]] =  Cl[soil_layer_num-1,n_dh_ind[0]]*1.0
    Ch_update[soil_layer_num-1,n_dh_ind[0]] =  Ch[soil_layer_num-1,n_dh_ind[0]]*1.0
    Cb_update[soil_layer_num-1,n_dh_ind[0]] =  Cb[soil_layer_num-1,n_dh_ind[0]]*1.0
    
    dz_update_n  = (z_matrix_update[1:soil_layer_num-1,n_dh_ind[0]]-z_matrix_update[0:soil_layer_num-2,n_dh_ind[0]])
    zi_n         = (z_matrix[1:soil_layer_num-1,n_dh_ind[0]]+dh_ST[n_dh_ind[0]])-z_matrix_update[0:soil_layer_num-2,n_dh_ind[0]]
    zip1_n       = z_matrix_update[1:soil_layer_num-1,n_dh_ind[0]] - (z_matrix[1:soil_layer_num-1,n_dh_ind[0]]+dh_ST[n_dh_ind[0]])
    
    Cl_update[1:soil_layer_num-1,n_dh_ind[0]] = (Cl[1:soil_layer_num-1,n_dh_ind[0]]*zi_n + Cl[2:soil_layer_num,n_dh_ind[0]]*zip1_n)/dz_update_n
    Ch_update[1:soil_layer_num-1,n_dh_ind[0]] = (Ch[1:soil_layer_num-1,n_dh_ind[0]]*zi_n + Ch[2:soil_layer_num,n_dh_ind[0]]*zip1_n)/dz_update_n
    Cb_update[1:soil_layer_num-1,n_dh_ind[0]] = (Cb[1:soil_layer_num-1,n_dh_ind[0]]*zi_n + Cb[2:soil_layer_num,n_dh_ind[0]]*zip1_n)/dz_update_n


    return Cl_update ,Ch_update,Cb_update       
        
          
        
        