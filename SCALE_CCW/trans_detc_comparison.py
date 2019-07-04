#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 31 13:51:32 2018

@author: qina
"""

import numpy as np
import matplotlib.pyplot as plt

Kr    = 0.005       # s/m
tau_c = 5.6         # Pa
pf    = 2.5         # 2.5 or 1.5, try it
k     = 20          # [-] sediment transport coefficient
g     = 9.81        # m/s^2
rho_s = 1.34*1000.0 # kg/m^3
rho_w = 1000.0      # kg/m^3
D50   = 0.0024*1e-3 # m

R     = rho_s/rho_w-1.0
Kf    = k*(g*R*D50**3)**0.5/(rho_w*g*R*D50)**pf

tau   = np.arange(tau_c,10)

Dc    = Kr/rho_s*(tau-tau_c)
qs    = Kf*(tau-tau_c)**pf

qs_MP = 8*g**0.5/R*(tau/(rho_w*g)-tau_c/(rho_w*g))**1.5
fig, ax = plt.subplots(figsize=(5.3,2.5 ), dpi=90)    
ax.plot(tau,Dc   ,'r',alpha = 0.5,label = r'$D_c$' )
ax.plot(tau,qs   ,'g--',alpha = 0.5,label = r'$q_s$' )
#ax.plot(tau,qs_MP,'b--',alpha = 0.5,label = r'$q_sMP$' )

ax.set_xlabel(r'$\tau_b$ $[Pa]$ ',fontsize=16) #, linespacing=3.2
ax.set_ylabel(r'flux $[m/s]$',fontsize=16)
ax.legend(loc=2)
#ax.set_xlim() 

