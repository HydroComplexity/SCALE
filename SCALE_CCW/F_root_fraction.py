# -*- coding: utf-8 -*-
"""
Created on Wed Oct 19 11:21:59 2016

@author: a
"""

import numpy as np
import matplotlib.pyplot as plt

from parameters import *

#zr = np.arange(0,1,0.1)

#dz = (initial_soil_depth-La)*dzz/np.sum(dzz) # delta z, soil layer thickness
#z = np.zeros(soil_layer_num) #soil_layer_depth
#z[0]=dz[0]+La
#for i in xrange(1,soil_layer_num):
#    z[i]=dz[i]+z[i-1]
        
def root_fr(z,dz):     
#    z50 = 0.1  # 0.1
#    z95 = 0.6 # 0.6
#    z50= 0.28
#    z95 = 0.92
    
    # MLCAN, SOY and maize
    z50 = 0.2  # 0.1
    z95 = 0.8 # 0.6
    
    crf = 1.27875 / (np.log10(z50) - np.log10(z95)); 
    Frooti = -dz*crf/z50*(z/z50)**(crf-1)*(1+(z/z50)**crf)**(-2)

    return Frooti

if __name__ == "__main__":
    dz = (initial_soil_depth-La)*dzz/np.sum(dzz) # delta z, soil layer thickness
    z = np.zeros(soil_layer_num) #soil_layer_depth
    z[0]=dz[0]
    for i in xrange(1,soil_layer_num):
        z[i]=dz[i]+z[i-1]
    
    Frooti = root_fr(z,dz)
    fig, ax = plt.subplots()
    
    ax.plot(Frooti, z, 'o', Frooti, z)
    plt.gca().invert_yaxis()

#ax.plot(x_prim, pN, label=r"pN(x)")
#ax.legend(loc=9) # upper left corner
#ax.set_xlabel(r'x')
#fig.savefig("Q6_a.png")