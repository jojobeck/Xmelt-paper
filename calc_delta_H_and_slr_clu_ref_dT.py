
import numpy as np
import dimarray as da
from matplotlib import pylab as plt
from matplotlib import ticker, cm, colors

from matplotlib.colors import LogNorm
import matplotlib.gridspec as gridspec
import matplotlib.animation as animation
import os
import matplotlib.pyplot as pl

rm_fac = 1.25

files = [
    'proj_T_average_July_20_'+str(rm_fac)+'times10ymy_nodyn_but_mass_prech_dev_5.0_const_ref_0.25_temp_dT_par_275.15/',\
    'proj_T_average_July_10_'+str(rm_fac)+'times10ymy_nodyn_but_mass_prech_dev_5.0_const_ref_0.25_temp_dT_par_275.15/',\
    'proj_T_average_July_5_'+str(rm_fac)+'times10ymy_nodyn_but_mass_prech_dev_5.0_const_ref_0.25_temp_dT_par_275.15/',\
   ]
CLU = True
if CLU:
    clu_path = '/home/beckmann/exp_pism/pdd_runs/projection_runs/4500m/'
else:
    clu_path = output_proj
nn = np.load(clu_path +'proj_step6K_130a/' +'val_0.npz' )
m_gr=nn['m_gr']
m_thin_ground=nn['m_thin_ground']

H=nn['H']
H_ground = nn['H_ground']
smb_0_ground=nn['smb_0_ground']
H_float=nn['H_float']



A_ocean_m = 362.5*10**6*10**6
rho_ice = 910.0 
rho_fw =1000.
rho_sw =1028.0
dx = 4500.
dy = 4500.
m3_to_slr= rho_ice/(rho_fw*A_ocean_m)

to_month= 1/12.
smb_to_vol = to_month*dx*dy/rho_ice





def create_nc_monthly_variable_1d(exp,varibale,fname,dim_name):
    months = np.arange(len(exp))

    slr_po=da.DimArray(exp.tolist(),axes = [months],dims = ['months'])
    month=da.DimArray(months.tolist(),axes = [months],dims = ['months'])
#     month_co=da.DimArray(mons_cons.tolist(),axes = [mons_cons],dims = ['months'])
#     year_date=da.DimArray(years_mon.tolist(),axes = [mons_cons],dims = ['months'])
#     dim_name = my_vari_name(variable)
    dataset =da.Dataset({variable:slr_po})
#     print(clu+fname+dim_name)
    dataset.write_nc(clu_path+fname+dim_name)


for k, fname in enumerate(files):
    nn = np.load(clu_path + files[k] +'calc_delta_smb_grounded.npz')
    delta_smb_mioc_ground = nn['smb']
    delta_H = delta_smb_mioc_ground *smb_to_vol /(dx*dy)

    ddH =delta_H.cumsum(axis = 0)

    # change of volume above gflotation ,calc sea level rise
    H_above_fl=H_ground -H_float[m_gr]
    H_new = H_above_fl +ddH
    for i in range(H_new.shape[0]):
        h_melted = H_new[i,:] <0
        H_new[i,h_melted]=0


    slr_new_melted =m3_to_slr*dy*dx*H_new.sum(axis =1)

    #save calculation
    dim_name = 'ts_calc_slr_pot_no_neg_mass.nc'
    variable = 'sea_level_rise_potential'
    create_nc_monthly_variable_1d(slr_new_melted,variable,files[k],dim_name)
    
    # change of real ice thickness:
    H_new_all = H_ground +ddH
    for i in range(H_new_all.shape[0]):
        h_melted = H_new_all[i,:] <0
        H_new_all[i,h_melted]=0
    
    H_new_3d_all=np.zeros((H_new.shape[0],H.shape[0],H.shape[1]))
    for i in range(H_new.shape[0]):
        H_new_3d_all[i,m_gr]=H_new_all[i,]
        
    np.savez(clu_path+files[k] + 'calc_thickness_grounded_from_total_thickness.npz',thk_last = H_new_3d_all[-1,],thk_0 = H_new_3d_all[0,])
