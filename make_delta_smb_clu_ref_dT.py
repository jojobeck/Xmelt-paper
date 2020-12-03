
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
CLU = True
if CLU:
    clu_path = '/home/beckmann/exp_pism/pdd_runs/projection_runs/4500m/'
else:
    clu_path = output_proj
nn = np.load(clu_path +'proj_step6K_130a/' +'val_0.npz' )
m_gr=nn['m_gr']
m_thin_ground=nn['m_thin_ground']

H=nn['H']
smb_0_ground=nn['smb_0_ground']
H_float=nn['H_float']
files = [
    'proj_T_average_July_20_'+str(rm_fac)+'times10ymy_nodyn_but_mass_prech_dev_5.0_const_ref_0.25_temp_dT_par_275.15/',\
    'proj_T_average_July_10_'+str(rm_fac)+'times10ymy_nodyn_but_mass_prech_dev_5.0_const_ref_0.25_temp_dT_par_275.15/',\
    'proj_T_average_July_5_'+str(rm_fac)+'times10ymy_nodyn_but_mass_prech_dev_5.0_const_ref_0.25_temp_dT_par_275.15/',\
   ]

for k, fname in enumerate(files):
    smb00 = da.read_nc(clu_path + files[k]+ 'extra_y1gris_4500.nc', 'climatic_mass_balance')
    t,y,x = smb00.shape
    smb = np.zeros((t*130,y,x))
    
    for i in range(130):
        smb[i*12:(i+1)*12,:,:] = da.read_nc(clu_path + files[k]+ 'extra_y'+str(i+1)+'gris_4500.nc', 'climatic_mass_balance').values
    smb_miroc_ground = smb[:,m_gr]
    smb_miroc_ground[:,m_thin_ground]=0
    delta_smb_mioc_ground = np.zeros(np.shape(smb_miroc_ground))
    for i in range(130):
        for j in range(12):
            delta_smb_mioc_ground[i*12 +j,:] = smb_miroc_ground[i*12 +j,] - smb_0_ground[j,:]
            

    np.savez(clu_path + files[k] +'calc_delta_smb_grounded.npz',smb = delta_smb_mioc_ground)
