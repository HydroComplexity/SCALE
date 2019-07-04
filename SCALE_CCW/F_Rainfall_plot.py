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

#from parameters import *
import numpy as np
#from scipy.sparse.linalg import spsolve
#import timeit
import matplotlib.pyplot as plt
#from matplotlib import cm

def accu_rain(rain_yr):
    
    rain_accu = np.zeros(len(rain_yr))
    rain_accu[0]=rain_yr[0]
    for j in xrange(1,len(rain_yr)):
        rain_accu[j] = rain_accu[j-1]+ rain_yr[j]
    return rain_accu

#...................Rainfall
rain = np.load('rainfall_100yr3_daily.npy') # mm/(day)
#rain = rain/1000.0
#    rain50 = rain62[0:50*365]/1000.0
#    rain = np.hstack([rain50,rain50[::-1]]) 

year       = 100
d          = 365
time_steps = year*d
months     = 12*year
month      = 12

rain_obsev = rain[365*90:365*100]
rain_matrix      = np.reshape(rain_obsev,(10,365))
rain_annual_mean = np.mean(rain_matrix,0)

###### bar plot for rainfall 100 yr in one row
#fig, ax = plt.subplots(figsize=(8,1.5 ), dpi=90)     
#plt.rc("font", size=10)
#plt.bar(np.arange(time_steps),1000*rain)
#ax.set_xlim([0,time_steps])
##    ax.set_ylim([0,50]) 
#ax.set_xlabel('DOY', fontsize=12)
#ax.set_ylabel('Rainfall [mm/day]', fontsize=12)

#==============bar plot of rainfall in 365 days==========
#fig, ax = plt.subplots(1,1,figsize=(8.4,2.6), dpi=300 )  
##plt.yscale('log')
#plt.rc("font", size=14)
#
#for y in xrange(0,year):
#    ax.bar(np.arange(d),rain[y*d:(y+1)*d],edgecolor ='grey',color="grey",alpha=0.05,label=r"rain fall intensity")
#
#ax.bar(np.arange(d),rain[98*d:(98+1)*d],edgecolor ='royalblue',color="royalblue",label=r"rain fall intensity")
#
#ax.set_xlabel('DOY', fontsize=14)
#ax.set_ylabel('Rainfall [mm/day]', fontsize=14)
###ax.set_ylim([0,80])
#ax.set_xlim([0,365])
#fig.tight_layout()
#fig.savefig("rainfall_paper.png", bbox_inches="tight")


#==============mean value of observation data==========
#ax1 = ax.twinx()
#
#ax1.bar(np.arange(d),rain_annual_mean,edgecolor ='k',color="k",label=r"rain fall intensity")
##
###ax1.set_ylim([0.0, 0.8])
#ax1.set_xlim([0,365])
#ax1.set_ylabel(r"Accumulated rain [mm]", fontsize=12, color="black")
#for label in ax1.get_yticklabels():
#    label.set_color("black")         

#==============accumulating data==========
#ax1 = ax.twinx()
#for y in xrange(0,year):
#
#    ax1.plot(np.arange(d),accu_rain(rain[y*d:(y+1)*d]),lw=2, color="grey",alpha=0.3,label=r"rain fall depth")
#    
#ax1.plot(np.arange(d),accu_rain(rain_annual_mean),lw=2, color="k",label=r"rain fall depth")
##
###ax1.set_ylim([0.0, 0.8])
#ax1.set_xlim([0,365])
#ax1.set_ylabel(r"Accumulated rain [mm]", fontsize=12, color="black")
#for label in ax1.get_yticklabels():
#    label.set_color("black")         

#fig.savefig("Rainfall_3.png" )

#==============bar plot of rainfall in 12 months==========
#DOY_at_month_end = np.array([31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334, 365])
#DOY_at_month_beg = np.array([0,  31, 59, 90,  120, 151, 181, 212, 243, 273, 304, 334])
#rain_month   = np.zeros(months)
#
##rain[d*i_year+DOY_at_month_beg[0]:DOY_at_month_end[0]+d*i_year]
#for i_year in xrange(year):
#    for i_month in xrange(month):
#        numofmonth = i_month+month*i_year
#        rain_month[numofmonth] = np.sum(rain[d*i_year+DOY_at_month_beg[i_month]:DOY_at_month_end[i_month]+d*i_year])
#
#fig, ax = plt.subplots(1,1,figsize=(8,3), dpi=300 )  
##plt.yscale('log')
#plt.rc("font", size=14)
#
#for y in xrange(0,year):
##    ax.bar(np.arange(month),rain_month[y*month:(y+1)*month],edgecolor ='grey',color="grey",alpha=0.1,label=r"rain fall intensity")
#    ax.plot(np.arange(month),rain_month[y*month:(y+1)*month],color = 'grey',alpha=0.3)
#    
##----plot the mean value--------
#rain_month_matrix = np.reshape(rain_month,(year,month))
#rain_month_mean = np.mean(rain_month_matrix,0)
#
##ax.plot(np.arange(month),rain_month_mean, 'k')
#    
##y_real = 2 # 0 is 2015, 1 is 2016, 2 is 2017
##ax.bar(np.arange(month),rain_month[y_real*month:(y_real+1)*month],edgecolor ='royalblue',color="royalblue",label=r"rain fall intensity")
##ax.plot(np.arange(month),rain_month[y*month:(y+1)*month], 'b')
#
##ax.set_xlabel('DOY', fontsize=14)
#ax.set_xticks(np.arange(month))
#ax.set_xticklabels(("Jan", "Feb", "Mar", "Apr", "May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"))
##ax.set_xticks(np.arange(month),
##           ["Jan", "Feb", "Mar", "Apr", "May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"])
#
#ax.set_ylabel('Rainfall [mm/mo]', fontsize=14)
###ax.set_ylim([0,80])
#ax.set_xlim([0,11])
#ax.set_ylim([0,430])
#fig.tight_layout()
#fig.savefig("rainfall_monthly.png", bbox_inches="tight")
