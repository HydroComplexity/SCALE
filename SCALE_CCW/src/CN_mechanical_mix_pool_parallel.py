# -*- coding: utf-8 -*-
"""
Created on Sat Nov 12 21:16:25 2016

@author: a
"""
from parameters import *
from scipy.interpolate import interp1d
import multiprocessing
import numpy as np

def f_mm_parallel(params):

    cell,z,dz,Cl_cell,Ch_cell,Cb_cell = params
    
    til_list = np.where(z<mix_soil_depth)
    til_list = np.asarray(til_list[0],dtype = int)
    layer_end = til_list[len(til_list)-1] 
    
    
    Cl_at_mix_soil_depth = interp1d((z), Cl_cell)(mix_soil_depth) #
    Ch_at_mix_soil_depth = interp1d((z), Ch_cell)(mix_soil_depth) #
    Cb_at_mix_soil_depth = interp1d((z), Cb_cell)(mix_soil_depth) #

    
    dh_extra = mix_soil_depth-z[layer_end]
    k = layer_end+1

    Cl_sum_toplayer = np.sum(dz[0:k]*Cl_cell[0:k])+ dh_extra*Cl_at_mix_soil_depth
    Ch_sum_toplayer = np.sum(dz[0:k]*Ch_cell[0:k])+ dh_extra*Ch_at_mix_soil_depth
    Cb_sum_toplayer = np.sum(dz[0:k]*Cb_cell[0:k])+ dh_extra*Cb_at_mix_soil_depth    
    

    Cl_top_ave =   Cl_sum_toplayer/mix_soil_depth
    Ch_top_ave =   Ch_sum_toplayer/mix_soil_depth
    Cb_top_ave =   Cb_sum_toplayer/mix_soil_depth
    
    return layer_end,Cl_top_ave,Ch_top_ave,Cb_top_ave

def C_mechamical_mix(Cl_all,Ch_all,Cb_all, z_matrix,dz_matrix,nrows, ncols):


    Cl_all_T    = np.transpose(Cl_all)
    Ch_all_T    = np.transpose(Ch_all)
    Cb_all_T    = np.transpose(Cb_all)
    z_matrix_T  = np.transpose(z_matrix)
    dz_matrix_T = np.transpose(dz_matrix)
                            
    params = zip(np.arange(nrows*ncols),z_matrix_T,dz_matrix_T,Cl_all_T,Ch_all_T,Cb_all_T) #z_matrix,Cv,

    cores = multiprocessing.cpu_count()
    pool  = multiprocessing.Pool(processes=cores)
    l     = pool.map(f_mm_parallel, params) 
    
    pool.close()
    pool.terminate()
    pool.join()
    
    l =  map(list, zip(*l))

    layer_end  = np.asarray(l[0])

    Cl_top_ave = np.asarray(l[1])
    Ch_top_ave = np.asarray(l[2])
    Cb_top_ave = np.asarray(l[3])
    
    Cl_top_ave = np.transpose(Cl_top_ave)
    Ch_top_ave = np.transpose(Ch_top_ave)
    Cb_top_ave = np.transpose(Cb_top_ave)
    
        
    for j in xrange(nrows*ncols):
        kk = layer_end[j]+1
        Cl_all[0:kk,j]= Cl_top_ave[j]
        Ch_all[0:kk,j]= Ch_top_ave[j]
        Cb_all[0:kk,j]= Cb_top_ave[j]
    return Cl_all,Ch_all,Cb_all


#np.where( (wd < hmin) | (wd_pN < hmin) | ((wd+wd_pN) < hmin) | (Sss < delta)  )