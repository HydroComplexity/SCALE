# -*- coding: utf-8 -*-
"""
Created on Sat Nov 12 21:16:25 2016

@author: a
"""
from parameters import *
from scipy.interpolate import interp1d
import numpy as np

def C_mechamical_mix(Cl_all,Ch_all,Cb_all, z_matrix,dz_matrix,nrows, ncols,level):
    #mix_soil_depth = 0.3 #m
    if level ==1:
        mix_soil_depth = mix_soil_depth_fast
    else:
        mix_soil_depth = mix_soil_depth_slow

    layer_end = np.zeros(nrows*ncols,dtype=int)
    Cl_top_ave = np.zeros(nrows*ncols)
    Ch_top_ave = np.zeros(nrows*ncols)
    Cb_top_ave = np.zeros(nrows*ncols)
    
    
    for i in xrange(soil_layer_num-1):
        for j in xrange(nrows*ncols):
            if z_matrix[i,j]<mix_soil_depth and z_matrix[i+1,j]>mix_soil_depth:
                layer_end[j] = i

    for j in xrange(nrows*ncols):
        k = layer_end[j]

        
        Cl_top_ave[j] = np.sum(dz_matrix[0:k+1,j]*Cl_all[0:k+1,j])/z_matrix[k,j]
        Ch_top_ave[j] = np.sum(dz_matrix[0:k+1,j]*Ch_all[0:k+1,j])/z_matrix[k,j]
        Cb_top_ave[j] = np.sum(dz_matrix[0:k+1,j]*Cb_all[0:k+1,j])/z_matrix[k,j]
    
    for j in xrange(nrows*ncols):
        k = layer_end[j]+1
        Cl_all[0:k,j]= Cl_top_ave[j]
        Ch_all[0:k,j]= Ch_top_ave[j]
        Cb_all[0:k,j]= Cb_top_ave[j]
    return Cl_all,Ch_all,Cb_all



def C_mechamical_mix1(Cl_all,Ch_all,Cb_all, z_matrix,dz_matrix,nrows, ncols):
    #mix_soil_depth = 0.3 #m
    print 'tillage'
    layer_end = np.zeros(nrows*ncols,dtype=int)
    Cl_sum_toplayer = np.zeros(nrows*ncols)
    Ch_sum_toplayer = np.zeros(nrows*ncols)
    Cb_sum_toplayer = np.zeros(nrows*ncols)
    
    
    for i in xrange(soil_layer_num-1):
        for j in xrange(nrows*ncols):
            if z_matrix[i,j]<mix_soil_depth and z_matrix[i+1,j]>mix_soil_depth:
                layer_end[j] = i
            
    for j in xrange(nrows*ncols):
        Cl_at_mix_soil_depth = interp1d((z_matrix[:,j]), Cl_all[:,j])(mix_soil_depth) #
        Ch_at_mix_soil_depth = interp1d((z_matrix[:,j]), Ch_all[:,j])(mix_soil_depth) #
        Cb_at_mix_soil_depth = interp1d((z_matrix[:,j]), Cb_all[:,j])(mix_soil_depth) #
        k = layer_end[j]

        dh_extra = mix_soil_depth-z_matrix[k,j]
    
        Cl_sum_toplayer[j] = np.sum(dz_matrix[0:k+1,j]*Cl_all[0:k+1,j])+ dh_extra*Cl_at_mix_soil_depth
        Ch_sum_toplayer[j] = np.sum(dz_matrix[0:k+1,j]*Ch_all[0:k+1,j])+ dh_extra*Ch_at_mix_soil_depth
        Cb_sum_toplayer[j] = np.sum(dz_matrix[0:k+1,j]*Cb_all[0:k+1,j])+ dh_extra*Cb_at_mix_soil_depth    
    
    Cl_top_ave =   Cl_sum_toplayer/mix_soil_depth
    Ch_top_ave =   Ch_sum_toplayer/mix_soil_depth
    Cb_top_ave =   Cb_sum_toplayer/mix_soil_depth
       
    for j in xrange(nrows*ncols):
        k = layer_end[j]+1
        Cl_all[0:k,j]= Cl_top_ave[j]
        Ch_all[0:k,j]= Ch_top_ave[j]
        Cb_all[0:k,j]= Cb_top_ave[j]
    return Cl_all,Ch_all,Cb_all