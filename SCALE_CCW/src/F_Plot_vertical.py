# -*- coding: utf-8 -*-
"""
Created on Wed Feb 10 06:23:10 2016

@author: a
"""
from parameters import * 
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
# from landlab.plot.imshow import imshow_cell_grid
# from landlab import RasterModelGrid
import matplotlib.cm as cm

def SOC_ver_profile(cell_ero_num, cell_depo_num, C_matrix, C_matrix_ini, z_matrix, z_matrix_ini):

    C_ini  = C_matrix_ini[:,cell_ero_num]
    z_ini  = z_matrix_ini[:,cell_ero_num]-z_matrix_ini[0,cell_ero_num]
    
    C_ero  = C_matrix[:,cell_ero_num]
    z_ero  = z_matrix[:,cell_ero_num]-z_matrix[0,cell_ero_num]
    C_depo = C_matrix[:,cell_depo_num]
    z_depo = z_matrix[:,cell_depo_num]-z_matrix[0,cell_depo_num]
    
    fig, ax = plt.subplots(figsize=(4,6)) 
    ax.plot(C_ini/1000.0,z_ini,'o--',label='ini')
    ax.plot(C_ero/1000.0,z_ero, 'o-',label='ero')
    ax.plot(C_depo/1000.0,z_depo,'o-', label='depo')
    ax.legend(loc='lower right')
    ax.set_xlabel(r'[$KgC m^{-3}$]')
    ax.set_ylabel(r'Soil Depth [m]')
    ax.invert_yaxis()
    
    return ()

    
def SOC_ver_profile2(cell_ero_num, cell_depo_num, C_matrix, C_matrix_ini, z_matrix, z_matrix_ini):

    C_ini  = C_matrix_ini[:,cell_ero_num]
    z_ini  = z_matrix_ini[:,cell_ero_num]#-z_matrix_ini[0,cell_ero_num]
    
    C_ero  = C_matrix[:,cell_ero_num]
    z_ero  = z_matrix[:,cell_ero_num]#-z_matrix[0,cell_ero_num]
    C_depo = C_matrix[:,cell_depo_num]
    z_depo = z_matrix[:,cell_depo_num]#-z_matrix[0,cell_depo_num]
    
    fig, ax = plt.subplots(figsize=(4,6)) 
    ax.plot(C_ini/1000.0,z_ini[soil_layer_num-1]-z_ini,'o--',label='ini')
    ax.plot(C_ero/1000.0,z_ero[soil_layer_num-1] - z_ero, 'o-',label='ero')
    ax.plot(C_depo/1000.0,z_depo[soil_layer_num-1] - z_depo,'o-', label='depo')
    ax.legend(loc='lower right')
    ax.set_xlabel(r'[$KgC m^{-3}$]')
    ax.set_ylabel(r'Soil Depth [m]')
#    ax.invert_yaxis()
    
    return ()
