# -*- coding: utf-8 -*-
"""
Copyright (C) 2017, Qina Corporation
All rights reserved.

Distributed Hydrologicc and Regional Analysis (DHARA) Model
DHARA model is made available as a restricted, non-exclusive, 
non-transferable license for education and research purpose only, 
and not for commercial use. See the LICENSE.txt for more details.

Author: Qina.yan@gmail.com (Qina Yan)
"""

from parameters import *
import numpy as np
from scipy.sparse import diags
import matplotlib.pyplot as plt
from scipy.sparse import dia_matrix

"""
{ item_description }
matrix_A:
matrix_BC:
"""

def linear_diffusion(nrows,ncols,level):
    if level ==1:
        Dx = Dx1
        Dy = Dy1
    else:
        Dx = Dx2
        Dy = Dy2

#------------------------------------------------------------------------
# THE KEY PARAMETERS OF MATRIX A
#------------------------------------------------------------------------
    Sx = dt*Dx/dx**2
    Sy = dt*Dy/dy**2
    Sxy = -2*Sx-2*Sy
    
    Sxe = dt*Dx/dx/2
    Sye = dt*Dy/dy/2
    
    Sxm = dt/dx
    Sym = dt/dy
    
    Scx = dt*Dx/2.0/dx**2
    Scy = dt*Dy/2.0/dy**2

    Nx = nrows
    Ny = ncols

#------------------------------------------------------------------------
# BUIL MATRIX
#------------------------------------------------------------------------
    idia_u = np.arange(0,Nx*Ny,Ny)
    idia_l = np.arange(Ny-2,Nx*Ny,Ny)
  
    MainD                     = Sxy*np.ones(Nx*Ny)
    Dia_top                   = Sy*np.ones(Ny*(Nx-1))
    Dia_top[0:Ny]             = 2.0*Sy  #*np.ones(Ny)
    Dia_up                    = Sx *np.ones(Nx*Ny-1)
    Dia_up[idia_u]            = 2.0*Sx  #*np.ones(Nx)
    Dia_low                   = Sx *np.ones(Nx*Ny-1)
    Dia_low[idia_l]           = 2.0*Sx  #*np.ones(Nx)
    Dia_bot                   = Sy *np.ones(Ny*(Nx-1))
    Dia_bot[Ny*(Nx-2):Ny*(Nx-1)]  = 2.0*Sy   #*np.ones(Ny)
    
    indZero            = np.arange(Ny-1, Nx*Ny-1, Ny)
    Dia_low[indZero] = 0.0
    Dia_up[indZero] = 0.0 
    
    matrix_A = diags([Dia_bot,Dia_low,MainD,Dia_up,Dia_top], [-Ny,-1,0,1,Ny],shape=(Nx*Ny, Nx*Ny))
    matrix_A = matrix_A.tocsr()
    

#------------------------------------------------------------------------
# BUIL MATRIX for carbon and soil diffusion part Ax and Ay
#------------------------------------------------------------------------
# this is for (i,j+1) - (i,j-1)
#    idia_t2                  = np.arange(Ny-1, Nx*Ny, Ny)
#    
#    DCx_main                  = np.zeros(Nx*Ny)
#    DCx_main[idia_u]          = -Scx
#    DCx_main[idia_t2]         = Scx           
#    DCx_up                    = Scx*np.ones(Nx*Ny-1)
#    DCx_low                   = -Scx*np.ones(Nx*Ny-1)
#        
#    DCx_low[indZero]          = 0.0
#    DCx_up[indZero]           = 0.0 
#
#    matrix_Ax = diags([DCx_low,DCx_main,DCx_up], [-1,0,1],shape=(Nx*Ny, Nx*Ny))
#    matrix_Ax = matrix_Ax.tocsr()
#    
# # this is for (i+1,j) - (i-1,j)   
#    DCy_main                  = np.zeros(Nx*Ny)
#    DCy_main[0:Ny]            = -Scy
#    DCy_main[Nx*Ny-Ny:]       = Scy
#    DCy_top                   = Scy*np.ones(Nx*Ny-Ny)
#    DCy_bot                   = -Scy*np.ones(Nx*Ny-Ny)
#
#    matrix_Ay = diags([DCy_bot,DCy_main,DCy_top], [-Ny,0,Ny],shape=(Nx*Ny, Nx*Ny))
#    matrix_Ay = matrix_Ay.tocsr()

# this is for (i, j+1)-(i,j), Axp1
    DCxp1_up            = Scx*np.ones(Nx*Ny-1)
    DCxp1_up[indZero]   = 0.0 
    indZero2            = np.arange(Ny-2, Nx*Ny-1, Ny)
    DCxp1_low           = np.zeros(Nx*Ny-1)
    DCxp1_low[indZero2] = Scx
    
    matrix_Axp1 = diags([DCxp1_low, -Scx, DCxp1_up], [-1,0,1], shape=(Nx*Ny, Nx*Ny))#.toarray()
    matrix_Axp1 = matrix_Axp1.tocsr()
       
# this is for (i, j)-(i,j-1), Axm1
    DCxm1_low           = -Scx*np.ones(Nx*Ny-1)
    DCxm1_low[indZero]  = 0.0 
    indZero3            = np.arange(0, Nx*Ny-1, Ny)
    DCxm1_up            = np.zeros(Nx*Ny-1)
    DCxm1_up[indZero3]  = -Scx 
    matrix_Axm1 = diags([DCxm1_low, Scx, DCxm1_up], [-1,0,1], shape=(Nx*Ny, Nx*Ny))#.toarray()
    matrix_Axm1 = matrix_Axm1.tocsr()
    
# this is for (i+1, j)-(i,j), Ayp1
    DCyp1_low           = np.zeros(Nx*Ny-Ny)
    indZero4            = np.arange(Ny*Nx-2*Ny, Nx*Ny-Ny)
    DCyp1_low[indZero4] = Scy
    matrix_Ayp1 = diags([DCyp1_low,-Scy,Scy], [-Ny,0,Ny],shape=(Nx*Ny, Nx*Ny))#.toarray()
    matrix_Ayp1 = matrix_Ayp1.tocsr()    

# this is for (i, j)-(i-1,j), Aym1
    DCym1_up            = np.zeros(Nx*Ny-Ny)
    indZero5            = np.arange(0, Ny)
    DCym1_up[indZero5]  = -Scy
    matrix_Aym1 = diags([-Scy,Scy,DCym1_up], [-Ny,0,Ny],shape=(Nx*Ny, Nx*Ny))#.toarray()
    matrix_Aym1 = matrix_Aym1.tocsr()    

    
    ############################################################
    ##### Build the boundary condition 
    ############################################################
    #------------------------------------------------------------------------
    # for landscape diffusion 
    #------------------------------------------------------------------------
    matrix_BC = np.zeros(Nx*Ny)
    
    matrix_BC[0:Ny]        =  2*dy*Sy*Bu
    matrix_BC[Nx*Ny-Ny:]   = -2*dy*Sy*Bd
    matrix_BC[0:Nx*Ny:Ny]=matrix_BC[0:Nx*Ny:Ny]+2*dx*Sx*Bl
    matrix_BC[Ny-1:Nx*Ny:Ny]=matrix_BC[Ny-1:Nx*Ny:Ny]-2*dx*Sx*Br
#    
#    matrix_BC = np.zeros(Nx*Ny)
#    
#    matrix_BC[0:Ny-1]=-2*dy*Sy*Cu
#    matrix_BC[0]=matrix_BC[0]-2*dx*Sx*Cl
#    matrix_BC[Ny-1]=matrix_BC[Ny-1]-2*dx*Sx*Cr
#    
#    matrix_BC[Ny]=-2*dx*Sx*Cl
#    matrix_BC[Nx*Ny-1-Ny]=-2*dx*Sx*Cr
#    
#    matrix_BC[Nx*Ny-Ny:Nx*Ny-1]=-2*dy*Sy*Cd
#    matrix_BC[Nx*Ny-Ny]=matrix_BC[Nx*Ny-Ny]-2*dx*Sx*Cl
#    matrix_BC[Nx*Ny-1]=matrix_BC[Nx*Ny-1]-2*dx*Sx*Cr
#    
    #------------------------------------------------------------------------
    # for SOC associated with landscape diffusion
    #------------------------------------------------------------------------
    BCx_p1                = np.zeros(Nx*Ny)
    BCx_p1[Ny-1:Nx*Ny:Ny] = BCx_p1[Ny-1:Nx*Ny:Ny]-2*dx*Scx*Br
    BCx_m1                = np.zeros(Nx*Ny)
    BCx_m1[0:Nx*Ny:Ny]    = BCx_m1[0:Nx*Ny:Ny]   -2*dx*Scx*Bl
    BCy_p1                = np.zeros(Nx*Ny)
    BCy_p1[Nx*Ny-Ny:]     = BCy_p1[Nx*Ny-Ny:]    -2*dy*Scy*Bd
    BCy_m1                = np.zeros(Nx*Ny)
    BCy_m1[0:Ny]          = BCy_m1[0:Ny]         -2*dy*Scy*Bu
    
    return matrix_A, matrix_BC, matrix_Axp1, matrix_Axm1, matrix_Ayp1, matrix_Aym1,BCx_p1, BCx_m1, BCy_p1,BCy_m1


def linear_diffusion1D(Zr,dZr):
    zh  = Zr[soil_layer_num-1]
    L   = np.arange(soil_layer_num)
    Lp1 = np.roll(L, -1)
    #Lm1 = np.roll(L, 1)
    
    zpf = (Zr+Zr[Lp1])/2.0
    #zmf = np.roll(zpf, 1) #(Zr[Lm1]+Zr)/2.0
    
    eZpf = dt*Do/dZr**2*np.exp(-zpf/zh)
    eZmf = dt*Do/dZr**2*np.roll(eZpf, 1)#np.exp(-zmf/zh)
    
    eZpf[0]                = 0.0
    eZpf[soil_layer_num-1] = 0.0
    eZmf[0]                = 0.0
    eZmf[soil_layer_num-1] = 0.0
    
    A = np.roll(eZpf,1) 
    B = -(eZpf+eZmf)
    C = np.roll(eZmf,-1) 

    ### build the matrix A, Ax = b
    A1d_C = dia_matrix( ([C,B,A],[-1,0,1]),shape=(nz, nz))
    A1d_C = A1d_C.tocsr()
    return A1d_C

def linear_diffusion1D_old(Zr):
    
    k1 = np.ones(soil_layer_num-1)
    k2 = np.ones(soil_layer_num)
    k3 = np.ones(soil_layer_num-1)
    for x in xrange(1,soil_layer_num-1):
        k1[x] = dt*Do*(2.0*np.exp(-Zr[x]/Zr[soil_layer_num-1])/(dZr[x+1]*(dZr[x+1]*dZr[x]))-np.exp(-Zr[x]/Zr[soil_layer_num-1])/Zr[soil_layer_num-1]/dZr[x+1])
        k2[x] = 1+dt*Do*(-2.0*np.exp(-Zr[x]/Zr[soil_layer_num-1])/(dZr[x+1]*(dZr[x+1]*dZr[x]))-2.0*np.exp(-Zr[x]/Zr[soil_layer_num-1])/(dZr[x]*(dZr[x+1]*dZr[x]))+np.exp(-Zr[x]/Zr[soil_layer_num-1])/Zr[soil_layer_num-1]/dZr[x+1])
        k3[x-1] = dt*Do*2.0*np.exp(-Zr[x-1]/Zr[soil_layer_num-1])/(dZr[x-1]*(dZr[x]*dZr[x-1]))
    k1[0]=0.0
    k2[0]=1.0
    k3[0]=0.0
    
    k1[soil_layer_num-2]=0.0 
    k2[soil_layer_num-1]=1.0
    k3[soil_layer_num-2]=0.0        
    matrix_1Diff_C = diags([k3,k2,k1], [-1,0,1], shape=(soil_layer_num,soil_layer_num))
    return matrix_1Diff_C


def Overland_SedimentTransport(wd_vector,slope,direction,Ca_vector,nrows,ncols,level,mask_ID):
#    print Ca_vector.shape
    if level ==1:
        Kw = Kw1
    else:
        Kw = Kw2
    
    wd = wd_vector.reshape(nrows,ncols)
    Cal = Ca_vector[0,:].reshape(nrows,ncols)
    Cah = Ca_vector[1,:].reshape(nrows,ncols)
    Cab = Ca_vector[2,:].reshape(nrows,ncols)
    
    xn = [-1,0,1,-1,1,-1,0,1,0]
    yn = [1,1,1,0,0,-1,-1,-1,0]
    dn = [2.**.5,1.,2.**.5,1.,1.,2.**.5,1.,2.**.5,0]
    Sw = Kw*dt/dx

#    alfa_e*dt*(area_vector**m3)*(slope_vector**m4)
#    eta_w = np.zeros((nrows,ncols))
    eta_w_vector = np.zeros(nrows*ncols)
    qs_upstream = np.zeros((nrows,ncols))
    qsCa_upstreaml = np.zeros((nrows,ncols))
    qsCa_upstreamh = np.zeros((nrows,ncols))
    qsCa_upstreamb = np.zeros((nrows,ncols))
    for i in xrange (0, nrows):
        for j in xrange (0,ncols):
            n=direction[i][j]
            if n == 8:
                tau = 0.0
            else:
                xw = i + xn[n]
                yw = j + yn[n]
                tau = (wd[i][j])*(slope[i][j]) 
            if tau > tau_c/(rho_w*g):
                AlSp = (tau -tau_c/(rho_w*g))**(1.6)
            else:
                AlSp = 0.0
            qs_upstream[xw][yw]= AlSp + qs_upstream[xw][yw] #eta = eta_vector.reshape(nrows,ncols)
            qsCa_upstreaml[xw][yw]= AlSp*Cal[i][j] + qsCa_upstreaml[xw][yw]
            qsCa_upstreamh[xw][yw]= AlSp*Cah[i][j] + qsCa_upstreamh[xw][yw]
            qsCa_upstreamb[xw][yw]= AlSp*Cab[i][j] + qsCa_upstreamb[xw][yw]

    
    qs_upstream_vector = np.array(qs_upstream).flatten()
    qsCa_upstream_vectorl = np.array(qsCa_upstreaml).flatten()
    qsCa_upstream_vectorh = np.array(qsCa_upstreamh).flatten()
    qsCa_upstream_vectorb = np.array(qsCa_upstreamb).flatten()
    qsCa_upstream_vector  = np.array([qsCa_upstream_vectorl,qsCa_upstream_vectorh, qsCa_upstream_vectorb])*1.0
    
    tau                               = wd*slope
    qs_downstream                     = np.zeros((nrows,ncols))
    qs_downstream[tau>(tau_c/(rho_w*g))] = (tau[tau>(tau_c/(rho_w*g))]-tau_c/(rho_w*g))**(1.6)
    qs_downstream_vector              = np.array(qs_downstream).flatten()
#    qs_upstream_vector[mask_ID]       = qs_upstream_vector[mask_ID]*0.0
#    qs_downstream_vector[mask_ID]     = qs_downstream_vector[mask_ID]*0.0
    qsCa_downstream_vector            = qs_downstream_vector*Ca_vector
    
#    # boundary condition
##    area_vector_bn = area_vector[(nrows-1)*ncols:nrows*ncols]
##    slope_vector_bn = slope_vector[(nrows-1)*ncols:nrows*ncols]
##    AS_bn = area_vector_bn**l*slope_vector_bn**p
#    q_out = 0.0
#    AS_bn = q_out*np.ones(ncols)
#    AS[(nrows-1)*ncols:nrows*ncols] = AS_bn
    eta_w_vector = Sw*(qs_upstream_vector -qs_downstream_vector)
    eta_w_Ca_vector = Sw*(qsCa_upstream_vector -qsCa_downstream_vector)

#    eta_w_vector_update = alfa_e*dt/dn[]*(area_vector**m3)*(slope_vector**m4)
#    eta_w_vector_update = Swx1*(area_vector)**(l-1)*matrix_Ax.dot(area_vector)*(slope_vector)**p+ \
#                          Swy1*(area_vector)**(l-1)*matrix_Ay.dot(area_vector)*(slope_vector)**p+ \
#                          Swx2*(area_vector)**(l)*(slope_vector)**(p-1)*(matrix_Bx.dot(eta_vector))+ \
#                          Swy2*(area_vector)**(l)*np.abs(slope_vector)**(p-1)*(matrix_By.dot(eta_vector))
#    eta_w_vector = np.array(eta_w).flatten()
    
    return eta_w_vector,eta_w_Ca_vector
    
#def Overland_SedimentTransport(wd,Sse,Ssw,Ssn,Sss,A_p1,A_m1,A_pN,A_mN,Ca_vector,nrows,ncols,level):
#    if level ==1:
#        Kw = Kw1
#    else:
#        Kw = Kw2
#    a = 17.0/6.0
#    b = 2.5
#    Swx = Kw*dt/dx
#    Swy = Kw*dt/dy
#    
#    wd_p1 = A_p1.dot(wd)
#    wd_m1 = A_m1.dot(wd)
#    wd_pN = A_pN.dot(wd)
#    wd_mN = A_mN.dot(wd)

##    eta_w = np.zeros((nrows,ncols))
#    eta_w_vector = Swy*(((wd_pN+wd)*0.5)**a*Sss**b-((wd_mN+wd)*0.5)**a*Ssn**b)
##Swx*(((wd_p1+wd)*0.5)**a*Sse**b-((wd_m1+wd)*0.5)**a*Ssw**b)#
#    return -eta_w_vector,np.zeros((nrows,ncols)) #eta_w_Ca_vector

def rill_erosion(wd_vector,slope,Ca_vector,nrows,ncols):

    
    wd = wd_vector.reshape(nrows,ncols)
    
    
    eta_r             = -dt*Kr/rho_b*(rho_w*g*wd*slope -tau_c)
    no_move_id        = np.where(rho_w*g*wd*slope<tau_c)
    eta_r[no_move_id] = 0.0
    eta_r_vector      = np.array(eta_r).flatten()
    eta_r_Ca_vector   = eta_r_vector*Ca_vector

    
    return eta_r_vector,eta_r_Ca_vector