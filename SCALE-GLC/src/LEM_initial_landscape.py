# -*- coding: utf-8 -*-
"""
Created on Mon Aug 17 22:20:18 2015

@author: a
"""
from parameters import *

from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import numpy as np
import gc
import os
import gdal

#---------------------------
# create initial elevation
#---------------------------
def ini_landscape():
    nrows = 20
    ncols = 50
# -------------Hill-----------------------------
    Y = np.arange(0, nrows*dy, dy)
    X = np.arange(0, ncols*dx, dx)
    X, Y = np.meshgrid(X, Y)
#    eta = (2.5-2*np.cos(4*np.pi*X/np.max(X))-np.cos(4*np.pi*Y/np.max(Y)))/20 + 100.0 #+0.1*mg.node_y
    eta = 12.0*(np.cos(4.0*np.pi*X/np.max(X))+2) #+0.1*mg.node_y
    eta = eta+(np.amax(Y)-Y)*0.075+ 0.01*np.random.rand(nrows,ncols) #+6.0
    

#    slope_vector = slope_direction(eta_vector)[0]
#    slope = slope_vector.reshape(nrows,ncols)
#    eta = eta+slope*dx*np.random.rand(nrows,ncols)

##------------valley----------------------------
#    Y = np.arange(0, nrows*dy, dy)
#    X = np.arange(0, ncols*dx, dx)
#    X_mid = ncols//2
#    
#
#    a = 10.0
#    b = 5.0
#    c = 1. # 300.0
#    X, Y = np.meshgrid(X, Y)
#    x1 = X[:,0:X_mid]
#    x2 = X[:,X_mid:ncols]
#    xx1 = x1/np.amax(x1)*a
#    xx2 = (x2-np.amax(x1))/(np.amax(x2)-np.amax(x1))*a
#    z1 = c/(1+np.exp(-xx1+b))
#    z2 = c/(1+np.exp(-xx2+b))
#    z1 = -z1 + np.amax(z1)+np.amin(z2)
#    
#    eta = np.zeros((nrows, ncols))
#    eta[:,0:X_mid] = z1
#    eta[:,X_mid:ncols] =z2
#    eta = eta+(np.amax(Y)-Y)*0.01 +0.001*np.random.rand(nrows,ncols) #(np.amax(Y)-Y)

    return X,Y,eta


def loadDEM(tre,unt):
    
    current_dir = os.getcwd()
    os.chdir(current_dir+"/"+'/Topo_SOCini_Landcover/')
    #
    if tre:
        gdata = gdal.Open("Treated_DEM_2m_2.tif") #CCW_2mDEM_deepStr, DEM_square_nofill2 saybrook_rec3_m  rec_test rec_test3_ditch rec_test3_stream
    if unt:
        gdata = gdal.Open("Untre_dem_2m.tif") #Untre_dem_2m
    geo = gdata.GetGeoTransform()
    Topography_rec = gdata.ReadAsArray()
#    Topography_rec[Topography_rec<0.0] = np.max(Topography_rec)
#    Topography_rec[115,118] =  Topography_rec[115,118]+0.99
    nrows, ncols = Topography_rec.shape
    Y = np.arange(0, nrows*dy, dy)
    X = np.arange(0, (ncols)*dx, dx) 
    X, Y = np.meshgrid(X, Y)
    
    if tre:
        Manndata = gdal.Open("landuse_tre_2m_clip.tif") # 
    if unt:
        Manndata = gdal.Open("landuse_untre_2m_clip.tif")# untreated_landuse_clip #landuse_untre_2m_clip #untreated_landuse_com_2grass 
    Mann_vg   = Manndata.ReadAsArray()
    
    mann_vg              = Mann_vg*1.0
    mann_vg[mann_vg==0] = -15
    Mann_largest = int(np.max(mann_vg))
    
    mannings_coeff = np.array([0.094, 0.125, 0.034, 0.188, 0.078,0.282,0.125])
    mannings_coeff = mannings_coeff/(9.81**0.5)/3600.0 # vegetation cover /24.0/3600.0 # unit is s/m^(1/3) -> hr/m^(1/3)

    for i in xrange(1,Mann_largest+1): 
        mann_vg[mann_vg==i] = mannings_coeff[i-1]
    
    mann_vg[mann_vg==-15] = 15
    
    os.chdir(current_dir)
    return X, Y, Topography_rec,mann_vg

#==============================================================================
# the main code (which does not go into fn)
#==============================================================================
if __name__ == "__main__":
    tre = False
    unt = True
    X,Y,Z,mann_vg = loadDEM(tre,unt) # ini_landscape  loadDEM
    #X,Y,Z = ini_landscape()
    
    nrows, ncols = Z.shape
    Z[Z<0.0] = np.nan
    np.save('eta_ini',Z)
    #Z_vector = np.array(Z).flatten()
    #direction = slope_direction(Z_vector)[1]
    #area_vector = drainge_area(direction)
    #area = area_vector.reshape(nrows,ncols)
    
#-------------------------------------------------------------
# find the locations of sampling point in the real landscape
#-------------------------------------------------------------   
#    current_dir = os.getcwd()
#    os.chdir(current_dir+"/"+'/Topo_SOCini_Landcover/')
#    #
#    gdata = gdal.Open("CCW_2mDEM.tif") #CCW_2mDEM_deepStr, DEM_square_nofill2 saybrook_rec3_m  rec_test rec_test3_ditch rec_test3_stream
#    geo = gdata.GetGeoTransform()
#    
#    x_easting  = np.array([589111.44,589125.07,589133.22,589156.32,589247.27])
#    y_northing = np.array([4620990.43,4621009.11,4621022.17,4621059.47,4621103.17])
#    xres = geo[1]
#    yres = geo[5]
#    xmin = geo[0] #+ xres * 0.5
#    ymax = geo[3] #- yres * 0.5
#    x_convt = (x_easting  - xmin)/xres
#    y_convt = (y_northing - ymax)/yres
#    col_samp = np.round(x_convt).astype(int)
#    row_samp = np.round(y_convt).astype(int)
#---------------------------------------------------
# plot 2D view of topogrpahy in the real landscape
#---------------------------------------------------
    fig, ax = plt.subplots(figsize=(12,6 ), dpi=120)     
    plt.rc("font", size=18)
    figplot  = ax.matshow(Z, extent=[0,ncols*dx,nrows*dy,0],cmap=cm.terrain) #,vmin=0,vmax=0.032 #pcolor , matshow, imshow
    cbar = fig.colorbar(figplot,shrink=0.6, pad = 0.02)  #if shrink, shrink=0.9,
    cbar.set_label(r'Elevation [m]')

    myLocator = mticker.MultipleLocator(200)
    ax.xaxis.set_major_locator(myLocator)
    ax.set_xlabel('x [m]', fontsize=18)
    ax.set_ylabel('y [m]', fontsize=18)
#    ax.scatter(X[row_samp,col_samp], Y[row_samp,col_samp], s=30,c ='orange' , lw=0.0,clip_on=True)
    
#---------------------------------------------------
# plot 2D view of topogrpahy in the sine wave
#---------------------------------------------------

#    fig = plt.figure(figsize=plt.figaspect(0.3))
##    fig, ax = plt.subplots(figsize=(8,4 ), dpi=300) 
#    ax  = fig.add_subplot(111, projection='3d')
#    im  = ax.plot_surface(Y,X,Z,cmap=cm.terrain,rstride=2,cstride=1) #,rstride=1,cstride=1,linewidth=0,antialiased=True
#    ax.set_xlabel('X (m)')
#    ax.set_ylabel('Y (m)')
#    ax.set_title('Initial Surface Topography')
#    ax.set_zlabel('Elevation (m)')
##    cb = fig.colorbar(im, shrink=0.5, pad = 0.2)
##    ax.view_init(30,340) # for CCW
#    ax.view_init(30,12) # for sine wave landform
##    fig = plt.figure()
##    ax = fig.gca(projection='3d')

#---------------------------
# plot stream line/river network
#---------------------------
#    river_lines = plt.matshow(area) #plt.matshow(np.random.rand(64, 64), fignum=100, cmap=plt.cm.gray)
#    plt.colorbar(river_lines)
#    plt.show()
    
#    fig, ax = plt.subplots(figsize=(8,4 ), dpi=90)  
#    ax.plot(np.arange(0, ncols*dx, dx),Z[19,:])
#    
#    fig, ax = plt.subplots(figsize=(8,4 ), dpi=90)  
#    ax.plot(np.arange(0, nrows*dy, dy),Z[:,12])
