import numpy as np
from parameters import *
from scipy.sparse import dia_matrix
import matplotlib.pyplot as plt
from scipy.sparse.linalg import spsolve
from scipy.sparse import csc_matrix
from numpy.linalg import matrix_rank
from scipy.sparse import diags
import timeit

import matplotlib.pyplot as plt

def overland_setup(eta_vector):

    wd = 0.0 * np.ones(eta_vector.size)
    we = eta_vector + wd
    #mann = manning_value * np.ones(eta.shape)
    return wd, we

def matrix_for_rolling(nrows,ncols):
    Nx = ncols
    My = nrows
# H_i+1,j, south, pN
    diag_pN = np.zeros(Nx*My,dtype = 'int')
    diag_pN[Nx*My-Nx :] = 1
    
    A_pN = diags([diag_pN, 1], [0,Nx],shape=(Nx*My, Nx*My))
    A_pN = A_pN.tocsr()
    
    # H_i-1,j, north value,mN
    diag_mN = np.zeros(Nx*My,dtype = 'int')
    diag_mN[0:Nx] = 1
    A_mN = diags([1,diag_mN], [-Nx,0],shape=(Nx*My, Nx*My))
    A_mN = A_mN.tocsr()
    
    # H_i,j+1, East value, p1
    diag_p1_u = np.ones(Nx*My-1,dtype = 'int')
    diag_p1_u[Nx-1:Nx*My-1:Nx] = 0
    
    diag_p1_l = np.zeros(Nx*My,dtype = 'int')
    diag_p1_l[Nx-1:Nx*My:Nx] = 1
    
    A_p1 = diags([diag_p1_l,diag_p1_u], [0,1],shape=(Nx*My, Nx*My))
    A_p1 = A_p1.tocsr()
    
    # H_i,j-1, west, m1
    diag_m1_u = np.zeros(Nx*My,dtype = 'int')
    diag_m1_u[0:Nx*My:Nx] = 1
    
    diag_m1_l = np.ones(Nx*My-1,dtype = 'int')
    diag_m1_l[Nx-1:Nx*My-1:Nx] = 0
    
    A_m1 = diags([diag_m1_u,diag_m1_l], [0,-1],shape=(Nx*My, Nx*My))
    A_m1 = A_m1.tocsr()
    
    # H_i+1,j+1, southeast, p1 and pN
    diag_p1N_u1 = np.ones(Nx*My-Nx-1,dtype = 'int')
    diag_p1N_u1[Nx-1:Nx*My-Nx-1:Nx] = 0
    
    diag_p1N_dia = np.zeros(Nx*My,dtype = 'int')
    diag_p1N_dia[Nx-1:Nx*My:Nx] = 1
    diag_p1N_dia[Nx*My-Nx :]=1
    A_p1N = diags([diag_p1N_dia, diag_p1N_u1], [0,Nx+1],shape=(Nx*My, Nx*My))
    A_p1N = A_p1N.tocsr()
    
    
    # H_i-1,j-1, northwest, pm1 and mN
    diag_m1N_dia = np.zeros(Nx*My,dtype = 'int')
    diag_m1N_dia[0:Nx] = 1
    diag_m1N_dia[0:Nx*My:Nx] = 1
    
    diag_m1N_l2 = np.ones(Nx*My-Nx-1,dtype = 'int')
    diag_m1N_l2[Nx-1:Nx*My-Nx-1:Nx] = 0
    A_m1N = diags([diag_m1N_l2,diag_m1N_dia], [-Nx-1,0],shape=(Nx*My, Nx*My))
    A_m1N = A_m1N.tocsr()
    
    # H_i-1,j+1, northeast, mNp1 
    diag_mNp1_l1 = np.ones(Nx*My-Nx+1,dtype = 'int')
    diag_mNp1_l1[0:Nx*My-Nx+1:Nx] = 0
    
    diag_mNp1_dia = np.zeros(Nx*My,dtype = 'int')
    diag_mNp1_dia[0:Nx] = 1
    diag_mNp1_dia[Nx-1:Nx*My:Nx] = 1
    A_mNp1 = diags([diag_mNp1_l1,diag_mNp1_dia], [-Nx+1,0],shape=(Nx*My, Nx*My))
    A_mNp1 = A_mNp1.tocsr()
    
    # H_i+1,j-1, southeast, m1pN 
    diag_pNm1_u2 = np.ones(Nx*My-Nx+1,dtype = 'int')
    diag_pNm1_u2[0:Nx*My-Nx+1:Nx] = 0
    
    diag_pNm1_dia = np.zeros(Nx*My,dtype = 'int')
    diag_pNm1_dia[Nx*My-Nx :] = 1
    diag_pNm1_dia[0:Nx*My:Nx] = 1
    A_pNm1 = diags([diag_pNm1_dia,diag_pNm1_u2], [0,Nx-1],shape=(Nx*My, Nx*My))
    A_pNm1 = A_pNm1.tocsr()
    
    ############################################################
    ##### Build the boundary condition for AH^n+1 = H^n+BC
    ############################################################
    BC_we                  = np.zeros(Nx*My)
    BC_we[0:Nx]            = -Qn_we*dts/dy
    BC_we[Nx*My-Nx :]      = -Qs_we*dts/dy
    BC_we[0:Nx*My:Nx]      =BC_we[0:Nx*My:Nx]- Qw_we*dts/dx
    BC_we[Nx-1:Nx*My:Nx]   =BC_we[Nx-1:Nx*My:Nx]- Qe_we*dts/dx 
    
    return A_p1,A_m1,A_pN,A_mN, A_p1N,A_m1N,A_pNm1,A_mNp1,BC_we
    
def EstimateKOverland(A_p1,A_m1,A_pN,A_mN, A_p1N,A_m1N,A_pNm1,A_mNp1,wd, we,mann):
    """
    @brief      This function estimate K for the overland flow equation
    @param      wd    water depth
    @param      we    water elevation
    @param      mann  Maning's coefficient
    
    @return     { description_of_the_return_value }
    """

#    mann = manning_value * np.ones(wd.shape)
    hmin = 1e-6
    delta = 1e-8
    delta_l = 1.0
#    eta   = we-wd
    
    #=========== calculate the slope we_p1, we_m1, we_pN, we_mN, we_pNp1,we_mNm1,we_mNp1,we_m1pN
    we_p1 =  A_p1.dot(we) 
    we_m1 =  A_m1.dot(we) 
    we_mN =  A_mN.dot(we) 
    we_pN =  A_pN.dot(we)
    
    we_pNp1 =  A_p1N.dot(we) 
    we_mNm1 =  A_m1N.dot(we)     
    we_mNp1 =  A_mNp1.dot(we) 
    we_pNm1 =  A_pNm1.dot(we) 
    
#    eta_p1 =  A_p1.dot(eta) 
#    eta_m1 =  A_m1.dot(eta) 
#    eta_mN =  A_mN.dot(eta) 
#    eta_pN =  A_pN.dot(eta)
#    
#    eta_pNp1 =  A_p1N.dot(eta) 
#    eta_mNm1 =  A_m1N.dot(eta)     
#    eta_mNp1 =  A_mNp1.dot(eta) 
#    eta_pNm1 =  A_pNm1.dot(eta) 

    
    wd_p1 = A_p1.dot(wd)
    wd_m1 = A_m1.dot(wd)
    wd_pN = A_pN.dot(wd)
    wd_mN = A_mN.dot(wd)
    
    mann_p1 = A_p1.dot(mann)
    mann_m1 = A_m1.dot(mann)
    mann_pN = A_pN.dot(mann)
    mann_mN = A_mN.dot(mann)

    #=========== calculate the slope S
    # inside the domain area
     
    Sse_i   = np.sqrt( ((we_p1 - we)/dx)**2. + ((we_pNp1+we_pN  -we_mNp1-we_mN)/(4*dy))**2. )
    Ssw_i   = np.sqrt( ((we - we_m1)/dx)**2. + ((we_pN  +we_pNm1-we_mN  -we_mNm1)/(4*dy))**2. )
    Ssn_i   = np.sqrt( ((we_mN - we)/dy)**2. + ((we_mNp1+we_p1  -we_mNm1-we_m1)/(4*dx))**2. )
    Sss_i   = np.sqrt( ((we - we_pN)/dy)**2. + ((we_p1  +we_pNp1-we_m1  -we_pNm1)/(4*dx))**2. )
    

    Sse_min = ((wd_p1+wd)*0.5)**(1.0/3.0)/((mann_p1+mann)*0.5)**2.0/g*((we_p1-we)/dx)**2.0
    Ssw_min = ((wd_m1+wd)*0.5)**(1.0/3.0)/((mann_m1+mann)*0.5)**2.0/g*((we-we_m1)/dx)**2.0
    Ssn_min = ((wd_mN+wd)*0.5)**(1.0/3.0)/((mann_mN+mann)*0.5)**2.0/g*((we-we_mN)/dy)**2.0
    Sss_min = ((wd_pN+wd)*0.5)**(1.0/3.0)/((mann_pN+mann)*0.5)**2.0/g*((we_pN-we)/dy)**2.0
    
    Sse     = np.maximum(Sse_i,Sse_min)
    Ssw     = np.maximum(Ssw_i,Ssw_min)
    Ssn     = np.maximum(Ssn_i,Ssn_min)
    Sss     = np.maximum(Sss_i,Sss_min)
#    
#    Sse = np.minimum(Sse1,Sse_eta)
#    Ssw = np.minimum(Ssw1,Ssw_eta)
#    Ssn = np.minimum(Ssn1,Ssn_eta)
#    Sss = np.minimum(Sss1,Sss_eta)


    ##### Calculate slope and diffusion coeff
    Ke = ( ((wd_p1+wd)/2.)**(5./3.) ) / ( 0.5*(mann_p1+mann) * np.sqrt(Sse) )    
    ind = np.where( (wd < hmin) | (wd_p1 < hmin) | ((wd+wd_p1) < hmin) | (Sse < delta)  )
    Ke[ind] = 0.0

#    Khe = np.zeros(np.size(Ke))
#    ind_pose = np.where(Ke>0.0)
#    Khe[ind_pose] = Ke[ind_pose]/(((wd_p1[ind_pose]+wd[ind_pose])/2.)**(11./10.))  

    Kw = ( ((wd_m1+wd)/2.)**(5./3.) ) / ( 0.5*(mann_m1+mann) * np.sqrt(Ssw) )
    ind = np.where( (wd < hmin) | (wd_m1 < hmin) | ((wd+wd_m1) < hmin) | (Ssw < delta)  ) # 4*wd_m1<wd
    Kw[ind] = 0.0
    
#    Khw = np.zeros(np.size(Kw))
#    ind_posw = np.where(Kw>0.0)
#    Khw[ind_posw] = Kw[ind_posw]/(((wd_m1[ind_posw]+wd[ind_posw])/2.)**(11./10.))  
    
    
    Kn = ( ((wd_mN+wd)/2.)**(5./3.) ) / ( 0.5*(mann_mN+mann) * np.sqrt(Ssn) )
    ind = np.where( (wd < hmin) | (wd_mN < hmin) | ((wd+wd_mN) < hmin) | (Ssn < delta) | (Ssn > delta_l) ) #4*wd_mN<wd 
    Kn[ind] = 0.0

#    Khn = np.zeros(np.size(Kn))
#    ind_posn = np.where(Kn>0.0)
#    Khn[ind_posn] = Kn[ind_posn]/(((wd_mN[ind_posn]+wd[ind_posn])/2.)**(11./10.))  
    
    Ks = ( ((wd_pN+wd)/2.)**(5./3.) ) / ( 0.5*(mann_pN+mann) * np.sqrt(Sss) )
    ind = np.where( (wd < hmin) | (wd_pN < hmin) | ((wd+wd_pN) < hmin) | (Sss < delta)  ) #4*wd<wd_pN
    Ks[ind] = 0.0
    

#    Ke_max = (g*((wd_p1+wd)/2.)**3./((we_p1-we)/dx)**2.0)**0.5
#    Ke_max[ind] = 0.0
#    Kw_max = (g*((wd_m1+wd)/2.)**3./((we-we_m1)/dx)**2.0)**0.5
#    Kw_max[ind] = 0.0
#    Kn_max = (g*((wd_mN+wd)/2.)**3./((we-we_mN)/dy)**2.0)**0.5
#    Kn_max[ind] = 0.0
#    Ks_max = (g*((wd_pN+wd)/2.)**3./((we_pN-we)/dy)**2.0)**0.5
#    Ks_max[ind] = 0.0
#    
#    Kee = np.minimum(Ke,Ke_max)
#    Kww = np.minimum(Kw,Kw_max)
#    Knn = np.minimum(Kn,Kn_max)
#    Kss = np.minimum(Ks,Ks_max)
    
    return Ke, Kw, Kn, Ks,Sse, Ssw, Ssn, Sss


def ImplicitSolver_1st_BC(we,eta_vector,ppt,ke,kw,kn,ks,nrows, ncols,BC_we):
    """
    Note that Qina uses parameter imported from another file: parameters.py
    Nx, My and some others parameters can be found there
    """
#    hmin = 1e-3
    dtodx2 = -dts/(dx**2)
    dtody2 = -dts/(dy**2)
    Nx = ncols
    My = nrows

    we_old = we 
    # Set up the LHS matrix A

#        
    a2d_top   = dtody2 * ks
    a2d_up    = dtodx2 * ke    
    a2d_dia   = 1.0 - dtodx2 * (ke+ kw) - dtody2 * (kn+ ks);
    a2d_low   = dtodx2 * kw
    a2d_bot   = dtody2 * kn
    
    # Index diagonal
    idialr = np.array([np.arange(Nx,Nx*(My-1),Nx), np.arange(2*Nx-1,Nx*(My-1),Nx)]).flatten()
    idia = np.concatenate((np.arange(0,Nx,1), idialr, np.arange(Nx*(My-1),Nx*My,1)),axis=0) 
    a2d_bot[idia] = 0.0
    a2d_low[idia] = 0.0
    a2d_dia[idia] = 1.0	
    a2d_up[idia] = 0.0
    a2d_top[idia] = 0.0
    #rhs2d[idia] = eta_vector[idia] + 0.02
    
    ### boundary condition
    BC = np.zeros(Nx*My)+ppt
    we_old[idia] = 0.001
    
    # Do not change after this line
    a2d_low_r = np.roll(a2d_low,-1)
    a2d_up_r  = np.roll(a2d_up,1)
    a2d_bot_r = np.roll(a2d_bot,-Nx)
    a2d_top_r = np.roll(a2d_top,Nx) 
    indlow            = np.arange(Nx-1, Nx*My, Nx)
    a2d_low_r[indlow] = 0.0
    indup             = np.arange(Nx, Nx*My, Nx)
    a2d_up_r[indup]   = 0.0
    
    ### build the matrix A, Ax = b
    a2d = dia_matrix( ([a2d_bot_r,a2d_low_r,a2d_dia,a2d_up_r,a2d_top_r], 
			   	   	   [-Nx,-1,0,1,Nx]),
					   shape=(Nx*My, Nx*My) )	# This is LHS matrix A of linear system Ax=b

    a2d = a2d.tocsr()
 
   ### Solve x, Ax = b
    #we_new = np.linalg.solve(a2d.toarray(),we_old+BC)
    #we_new = spsolve(a2d,we_old+BC)
    we_new = spsolve(a2d,we_old+BC)

    wd_new = we_new-eta_vector        
    
    ind_neg = np.where(wd_new<0.0)
    we_new[ind_neg] = eta_vector[ind_neg] 
    wd_new[ind_neg]=0.0
    
    return we_new,wd_new


def ImplicitSolver_2nd_BC(we,eta_vector,ppt,ke,kw,kn,ks,nrows, ncols,BC_we):
        
    """
    Note that Qina uses parameter imported from another file: parameters.py
    Nx, My and some others parameters can be found there
    
    """
#    hmin = 1e-3
    dtodx2 = -dts/(dx**2)
    dtody2 = -dts/(dy**2)
    Nx = ncols
    My = nrows

    ############################
    ### Set up the LHS matrix A
    ############################
    a2d_top   = dtodx2 * ks
    a2d_up    = dtody2 * ke    
    a2d_dia   = 1.0 - dtodx2 * (ke+ kw) - dtody2 * (kn+ ks)
    a2d_low   = dtody2 * kw
    a2d_bot   = dtodx2 * kn
    
    
    # Index diagonal
    idia_j1 = np.arange(0,Nx*My,Nx)
    idia_jnx = np.arange(Nx-1,Nx*My,Nx)
    
    ### Boundary conditions
    a2d_dia[0:Nx] = a2d_dia[0:Nx]+a2d_bot[0:Nx]
    a2d_dia[Nx*My-Nx:Nx*My] = a2d_dia[Nx*My-Nx:Nx*My] + a2d_top[Nx*My-Nx:Nx*My] 
    a2d_dia[idia_j1] = a2d_dia[idia_j1]+a2d_low[idia_j1]   
    a2d_dia[idia_jnx] = a2d_dia[idia_jnx]+a2d_up[idia_jnx]
    
#    # Index diagonal
#    idia_u = np.arange(0,Nx*My,Nx)
#    idia_l = np.arange(Nx-1,Nx*My,Nx)
#    
##    idia_top = np.arange(0,Nx)
##    idia_bot = np.arange(Nx*My-Nx,Nx*My)
#    
#    a2d_top[0:Nx] = dtodx2 * (ks[0:Nx]+kn[0:Nx])
#    a2d_up[idia_u] = dtody2 * (ke[idia_u]+kw[idia_u])    
#    a2d_low[idia_l] = dtody2 * (ke[idia_l]+kw[idia_l])  
#    a2d_bot[Nx*My-Nx:Nx*My] = dtodx2 * (kn[Nx*My-Nx:Nx*My]+ks[Nx*My-Nx:Nx*My])
#   

    BC = BC_we + ppt

    # Do not change after this line
    a2d_low_r = np.roll(a2d_low,-1)
    a2d_up_r  = np.roll(a2d_up,1)
    a2d_bot_r = np.roll(a2d_bot,-Nx)
    a2d_top_r = np.roll(a2d_top,Nx)

    ### Build Matirx A, Ax = b
    indlow            = np.arange(Nx-1, Nx*My, Nx)
    a2d_low_r[indlow] = 0.0
    indup             = np.arange(Nx, Nx*My, Nx)
    a2d_up_r[indup]   = 0.0   
    a2d = dia_matrix( ([a2d_bot_r,a2d_low_r,a2d_dia,a2d_up_r,a2d_top_r], 
			   	   	   [-Nx,-1,0,1,Nx]),
					   shape=(Nx*My, Nx*My) )	
    a2d = a2d.tocsr()
 
    ### Solve x, Ax = b
    #we_new = np.linalg.solve(a2d.toarray(),we_old+BC)
    #we_new = spsolve(a2d,we_old+BC)

    we_new = spsolve(a2d,we+BC)    
    
    wd_new = we_new-eta_vector
    
    # old type
#    ind_neg = np.where(wd_new<0.0)
#    we_new[ind_neg] = eta_vector[ind_neg] 
#    wd_new[ind_neg]=0.0
    
    # new stype
    ave_del_dept = np.mean(we_new - we)
    we_new = we_new - (ave_del_dept - ppt)
    wd_new = wd_new - (ave_del_dept - ppt)
     
    ind_neg = np.where(wd_new<0.0)
    ind_pos = np.where(wd_new>0.0)
    need_minu = np.sum(wd_new[ind_neg])
    wd_new_pos_sum = np.sum(wd_new[ind_pos])
    if wd_new_pos_sum>0:
        wd_new[ind_pos] = (need_minu/wd_new_pos_sum + 1.0)*wd_new[ind_pos]
    we_new = wd_new + eta_vector
    we_new[ind_neg] = eta_vector[ind_neg] 
    wd_new[ind_neg]=0.0
    
    return we_new, wd_new

def ExplicitSolver_2st_BC(we,eta_vector,ppt,ke,kw,kn,ks,nrows, ncols,BC_we):
    
    dtodx2 = dts/(dx**2)
    dtody2 = dts/(dy**2)
    Nx = ncols
    My = nrows
    
    ############################
    ### Set up the LHS matrix A
    ############################
    a2d_top   = dtodx2 * ks
    a2d_up    = dtody2 * ke    
    a2d_dia   = 1.0 - dtodx2 * (ke+ kw) - dtody2 * (kn+ ks)
    a2d_low   = dtody2 * kw
    a2d_bot   = dtodx2 * kn
    
    # Index diagonal
    idia_j1 = np.arange(0,Nx*My,Nx)
    idia_jnx = np.arange(Nx-1,Nx*My,Nx)
    
    ### Boundary conditions
    a2d_dia[0:Nx] = a2d_dia[0:Nx]+a2d_bot[0:Nx]
    a2d_dia[Nx*My-Nx:Nx*My] = a2d_dia[Nx*My-Nx:Nx*My] + a2d_top[Nx*My-Nx:Nx*My] 
    a2d_dia[idia_j1] = a2d_dia[idia_j1]+a2d_low[idia_j1]   
    a2d_dia[idia_jnx] = a2d_dia[idia_jnx]+a2d_up[idia_jnx]

    BC = -BC_we + ppt
        
    # Do not change after this line
    a2d_low_r = np.roll(a2d_low,-1)
    a2d_up_r  = np.roll(a2d_up,1)
    a2d_bot_r = np.roll(a2d_bot,-Nx)
    a2d_top_r = np.roll(a2d_top,Nx)

    ### Build Matirx A, Ax = b
    indlow            = np.arange(Nx-1, Nx*My, Nx)
    a2d_low_r[indlow] = 0.0
    indup             = np.arange(Nx, Nx*My, Nx)
    a2d_up_r[indup]   = 0.0   
    a2d = dia_matrix( ([a2d_bot_r,a2d_low_r,a2d_dia,a2d_up_r,a2d_top_r], 
			   	   	   [-Nx,-1,0,1,Nx]),
					   shape=(Nx*My, Nx*My) )	
    a2d = a2d.tocsr()
 
    ### Solve x, Ax = b
    #we_new = np.linalg.solve(a2d.toarray(),we_old+BC)
    #we_new = spsolve(a2d,we_old+BC)

    #matrix_Ad.dot(eta_vector) 
    we_new = a2d.dot(we)+BC
    
    wd_new = we_new-eta_vector
    ind_neg = np.where(wd_new<0.0)
    we_new[ind_neg] = eta_vector[ind_neg] 
    wd_new[ind_neg]=0.0
    return we_new, wd_new
    
    
def FlowVelocity(we,we_p1,we_m1,we_mN,we_pN,wd,wd_p1,wd_m1,wd_mN,wd_pN,ke,kw,kn,ks):
#    u = np.zeros(np.size(wd))
#    v = np.zeros(np.size(wd))
#    inde = np.where(wd>0.0)
#    u[inde] = -(ke[inde]+kw[inde])/2.0/(wd[inde])*(we_p1[inde]-we_m1[inde])/2.0/dx
#    v[inde] = -(kn[inde]+ks[inde])/2.0/wd[inde]*(we_pN[inde]-we_mN[inde])/2.0/dy
    
    ue = np.zeros(np.size(wd))
    uw = np.zeros(np.size(wd))
    vn = np.zeros(np.size(wd))
    vs = np.zeros(np.size(wd))


    ide = np.where((wd_p1+wd)>0.01)
    idw = np.where((wd_m1+wd)>0.01)
    idn = np.where((wd_mN+wd)>0.01)
    ids = np.where((wd_pN+wd)>0.01)
    
    ue[ide[0]] = -ke[ide[0]]/((wd_p1[ide[0]]+wd[ide[0]])*0.5)*(we_p1[ide[0]]-we[ide[0]])/dx
    uw[idw[0]] = -kw[idw[0]]/((wd_m1[idw[0]]+wd[idw[0]])*0.5)*(we[idw[0]]-we_m1[idw[0]])/dx
    vn[idn[0]] = -kn[idn[0]]/((wd_mN[idn[0]]+wd[idn[0]])*0.5)*(we[idn[0]]-we_mN[idn[0]])/dy
    vs[ids[0]] = -ks[ids[0]]/((wd_pN[ids[0]]+wd[ids[0]])*0.5)*(we_pN[ids[0]]-we[ids[0]])/dy
    
    
    u = 0.5*(ue+uw)
    v = 0.5*(vn+vs)
    
#    nrows = 20
#    ncols = 50
#    wd       = (we_pN-we).reshape((nrows,ncols))#qdC_vector.reshape((nrows,ncols))
#    fig1, ax1 = plt.subplots(figsize=(5.3,2.5 ), dpi=90)
#    plt.rc("font", size=14)
#    flip_plot = False
#    if flip_plot:
#        figplot  = ax1.matshow(np.fliplr(wd), extent=[0,ncols*dx,nrows*dy,0],origin="lower",vmin=0,vmax=1.0)  
#    else:
#        figplot  = ax1.matshow(wd, extent=[0,ncols*dx,nrows*dy,0]) #,vmin=0,vmax=0.032 #pcolor , matshow, imshow
#    cbar = fig1.colorbar(figplot,shrink=0.9, pad = 0.1)
#    cbar.set_label('Water Depth [m]')
##        ax1.set_title('Water Depth [m]', fontsize=20)
#    ax1.set_xlabel('x [m]', fontsize=18)
#    ax1.set_ylabel('y [m]', fontsize=18)
    return u,v
    

