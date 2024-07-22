z_min = 0.4
z_max = 6

lim_deltaz = 2
Ks_cut = 24.5

# from myPackages import *
from astropy.io import fits
import numpy as np

cat_dir = "/Users/sinataamoli/Desktop/COSMOS/COSMOS2020_catalog/"
LSS_dir = "/Users/sinataamoli/Desktop/COSMOS/COSMOS2020_LSS/"
filament_dir = "/Users/sinataamoli/Desktop/COSMOS/COSMOS2020_FILAMENTS/"

# Importing The Catalog

hdul = fits.open(cat_dir + 'COSMOS2020_FARMER_R1_v2.1_p3.fits')
hdr = hdul[0].header

cosmos2020 = hdul[1].data
# hdul[1].columns or cols.info()

# # 'MODEL_FLAG': Flag (0: OK, 1: failed to converge, 2: drifted >0.6" from detection) 
# # 'FLAG_HSC':  Flag indicating quality of HSC imaging (0:clean, 1:masked)
# # 'FLAG_SUPCAM': Flag indicating quality of Suprime-Cam imaging (0:clean, 1:masked)
# # 'FLAG_UDEEP': Flag for the UltraVISTA ultra-deep regions (0:ultra-deep, 1:deep)
# # 'FLAG_UVISTA': Flag for the UltraVISTA region (0:inside, 1:outside)
# # 'FLAG_COMBINED': Comb. FLAG_UVISTA, FLAG_HSC, FLAG_SUPCAM (0:clean and inside UVISTA)
# # 'lp_zBEST': LePhare photo-z (=lp_zPDF if galaxy, NaN if X-ray source or masked)
# # 'lp_type': LePhare type (0: galaxy, 1: star, 2: Xray sour., -9: failure in fit)
# # 'lp_zp_2':  LePhare 2nd photo-z solution if a 2nd peak detected w. P>5% in PDF

print('COSMOS2020: %d sources \n' %len(cosmos2020))
print('Magnitude-Cut on Ks is %.2f \n' %Ks_cut)

# # Redshift between zmin and zmax
Data = cosmos2020[cosmos2020['lp_zPDF'] > z_min]
Data = Data[Data['lp_zPDF'] < z_max]

print('out of z_range (%.1f - %.1f): %d' %(z_min, z_max, (len(cosmos2020) - len(Data))))

# # "Conditions on redshift and uncertainty"
Data = Data[(Data['lp_zPDF_u68'][:] - Data['lp_zPDF_l68'][:]) < lim_deltaz]
Data = Data[(Data['lp_zPDF_u68'][:] - Data['lp_zPDF_l68'][:]) > 0]

print('z_uncertainty (>2): %d' %(len(cosmos2020[(cosmos2020['lp_zPDF_u68'] - cosmos2020['lp_zPDF_l68']) > lim_deltaz])))

# # No masked, No X-Ray source
Data = Data[Data['lp_zBEST']>0]

print('masked / X-ray sources: %d' %(len(cosmos2020) - len(cosmos2020[cosmos2020['lp_zBEST'] > 0])))

# # Magnitude 
Data = Data[Data['UVISTA_Ks_MAG'] > 0]
Data = Data[Data['UVISTA_H_MAG'] > 0]

Data = Data[Data['UVISTA_Ks_MAG'] < Ks_cut]
# # Final Sample

delta_z = np.array(Data['lp_zPDF_u68'] - Data['lp_zPDF_l68'])

print('\n Final Sample: \n')
print('** Data size: %d' %len(Data))
print('** z_uncertainty size: %d \n' %len(delta_z))

print('galaxies with second photo-z sol. (in Final Sample): %d' %len(Data[Data['lp_zp_2'] > 0]))
print('FLAG_UDEEP: 0 (ultra-deep): %d, 1: (deep): %d' %(len(Data[Data['FLAG_UDEEP'] == 0]), len(Data[Data['FLAG_UDEEP'] == 1])))

print('All other flags (HSC, SUPCAM, UDEEP, UVISTA, COMBINED) are checked / type: all are galaxies')
print('All sources have measured Ks and H magnitude')