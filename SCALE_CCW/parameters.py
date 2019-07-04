#Parameters
import numpy as np

#day = 365*50+3 #202 #2557 #18253
#input_file_rain = 'rainfall50yr4' #'rainfall5' #'rainfall50yr4', rainfall_200DOY_saved1

dx  = 2.0 # 300.0 1.2 
dy  = 2.0 # 300.0 1.2
dt  = 24.0 #365.0 #mark for landscape evolution model, in days
dts = 12.0
Dx1_save = 4.0*0.036/365.0/24. #1  #inputs.read_float('linear_diffusivity') m^2/yr --> m^2/day
Dy1_save = 4.0*0.036/365.0/24. #1  #inputs.read_float('linear_diffusivity') m^2/yr --> m^2/day
Dx1 = 1.8*0.036/365.0/24. #3, 0.8  #inputs.read_float('linear_diffusivity') m^2/yr --> m^2/day
Dy1 = 1.8*0.036/365.0/24. #3, 0.8  #inputs.read_float('linear_diffusivity') m^2/yr --> m^2/day

Dx2 = 1.0*0.036/365.0/24. #1  #inputs.read_float('linear_diffusivity') m^2/yr --> m^2/day
Dy2 = 1.0*0.036/365.0/24. #1  #inputs.read_float('linear_diffusivity') m^2/yr --> m^2/day

#cell_area = dx*dy

# soil flux for diffusion processes at boundaries
# positive value means facing downward or right
Bl = 0.0/365.0 # -0.01/365.0
Br = 0.0/365.0 # 0.01/365.0
Bu = 0.0/365.0
Bd = 0.0/365.0
uplift = 0.0#0.0003/365. # m/day, 0.0003/365.
lamdaP = 0.58

# the water flux on the 4 sides (boundaries)
Qn_we =  0.00
Qe_we =  0.00
Qw_we =  0.00
Qs_we =  0.00

g = 9.81*(3600.0)**2 # m/s^2 --> m/hr^2
manning_bs = 0.025/3600.0 # bare soil, winter time #/24.0/3600.0 # unit is s/m^(1/3) -> hr/m^(1/3)
manning_vg = 0.09/3600.0 # vegetation cover /24.0/3600.0 # unit is s/m^(1/3) -> hr/m^(1/3)
D50 = 0.002*1e-3 #m


# rill erosion
Kr    = 3e-05*0.005/3600# s/m -> hr/m #6e-06*
rho_b = 1.34*1000.0 # kg/m^3
rho_w = 1000.0 #kg/m^3
tau_c = 5.6*3600.0**2 # kg/m/s^2 -> kg/m/hr^2

R = rho_b/rho_w-1 #0.6 #[-]
# overland flow
#Kw = 2.6/dx/R*Cf**0.5/20.0 #8.0/dx/R*Cf**0.5/20.0# 4.0*10**(-3)/365.0 # Iff*8*Cf**0.5*i/dx/R  # # #
Kw1      = 0.06*1.0e-5*4.93*g**0.5/R #0.4, 0.02
Kw1_save = 0.10*1.0e-5*4.93*g**0.5/R 
Kw2      = 0.02*1.0e-5*4.93*g**0.5/R # 0.2*1e-5
#tau_star = 850.0
outlet_slope = 0.022 # 0.005

#Po = 0.0  #50.0/1000000.0/365.0 # m/day
#ho = 0.5 #m
initial_soil_depth = 1.0 #('initial_soil_depth')

La = 0.0#0.00015 # ('thinnest_soil_layer')

soil_layer_num= 7

# initial dz and z
mm = np.linspace(0.3,0.9,soil_layer_num)
dzz= -0.1/(mm-1)


#############
# soil moisture
##############
#n= 0.42 #inputs.read_float('n') #porosity
Infil_max = 2.0/1000.0  # m/day
total_soil_depth = initial_soil_depth # [m]
nz      = soil_layer_num # number of vertical grid size
poros   = 0.477       # soil porosity
Ksat    = 1.08*1e-3 #*24 # Sat. hydr.conductivity [m/hr]
alpha   = 0.0015             # parameter related to the inverse of the air entry suction [1/cm]
theta_S = poros             # Saturated water content [-]
theta_R = 0.19 #0.125       # Residual water content [-]
n       = 2.0               # Pore-size distributions [-]
m       = 1.0-1.0/n;    
Ss      = 5.0*1e-4 		# Storage
## soil-veg. parameter
s_star  = 0.31 #inputs.read_float('s_star') # point of incipient stress
sw      = theta_R/theta_S #inputs.read_float('sw') # permanent wilting point
sh      = 0.01 #inputs.read_float('sh') #hygroscopic point
sfc     = 0.3/poros #0.25/poros #inputs.read_float('sfc') # field capacity
Evap    = 3.2/10000.0         #Rainfall rate in unit [m/day]
T_max   = 9.1/10000.0 #inputs.read_float('T_max')/1000.0 #mm/day -> m/day

Ksoc = 1.0
ADD_s= 0.00/24.0 #1.5/24.0
Add_b= 0.00/24.0 # 0.8/24.0
## C/N ratio
CN_adds= 22.0 #inputs.read_float('CN_add')
CN_addb= 27.0 
CN_b_ini= 11.5 #inputs.read_float('CN_b')
CN_l_ini=30.0 #inputs.read_float('CN_l') # not provided
CN_h_ini= 22.0 #inputs.read_float('CN_h')
##others of the model estimated through calibration
a_plus=0.05 #inputs.read_float('a_plus')
a_minus= 1.0 #inputs.read_float('a_minus')
DEM_plus=0.2/24.0 #inputs.read_float('DEM_plus')
DEM_minus= 0.5/24.0 #inputs.read_float('DEM_minus')

b= 11.0 #inputs.read_float('b')
beta=2*b+4
dd= 3.0 #inputs.read_float('dd')
F= 0.1/24.0 #inputs.read_float('F')
ki_plus= 1.0/24.0 #inputs.read_float('ki_plus') # m^3/d g/C
ki_minus= 1.0/24.0 #inputs.read_float('ki_minus') # m^3/d g/C
kn= 0.5/24.0 #inputs.read_float('kn')
rr= 0.6 #inputs.read_float('rr')
rh_ini = 0.25
#rh= min(0.8,CN_h_ini/CN_l_ini)
D_top = 2.0#inputs.read_float('Dtop'), 
Do = D_top/100.0/100.0/365.0/24.0 #cm^2/yr --> m^2/day
#
mix_soil_depth_fast = 0.18
mix_soil_depth_slow = 0.05

# Kl = np.array([5.50656191e-07,   6.02273630e-08,   1.15977391e-07,         1.49302132e-07,   1.49563044e-07,   1.38723501e-07,         1.58610324e-07])
# Kh = np.array([5.05818905e-08,   4.74230242e-09,   7.33951743e-09,         6.89156369e-09,   4.70995757e-09,   7.58576366e-09,         1.03759598e-08])
# Kd = np.array([0.00119767,  0.00011956,  0.00021059,  0.00024506,  0.00021601,        0.00016623,  0.00013102])

Klc = 2.80432222e-07 #5.50656248e-07 #.12378826e-05
Khc = 1.14117216e-08 #5.05818905e-08 #1.36623120e-07
Kdc = 9.90434574e-05 #0.00119767 #0.00171096

#data save date
DOY_at_month = np.array([31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334, 365])





