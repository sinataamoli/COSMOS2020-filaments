from myPackages import *
import filament_functions as func

#########################################################################
################## Cosmological/Physical Parameters #####################
#########################################################################

c0 = 3e5
H_0 = 70
Omega_l = 0.7
Omega_m = 0.3

z_min = 0.4
z_max = 6

lim_deltaz = 2
Ks_cut = 24.5

cosmo = FlatLambdaCDM(H0 = H_0, Om0 = Omega_m)

#########################################################################
############################# Directories ###############################
#########################################################################

cat_dir = "/Users/sinataamoli/Desktop/COSMOS/COSMOS2020_catalog/"
LSS_dir = "/Users/sinataamoli/Desktop/COSMOS/COSMOS2020_LSS/"
filament_dir = "/Users/sinataamoli/Desktop/COSMOS/COSMOS2020_FILAMENTS/"

outputs_dir, plots_dir, eigenvalues_dir, signals_dir, inputs_dir = func.setup()
#########################################################################
############################# Import Data ###############################
#########################################################################

# %run import_COSMOS2020.py

# Data = Table(Data)
# Data.write(cat_dir + 'Data_magLimit245_Taamoli24.fits', overwrite=True)

Data = Table.read(cat_dir + 'Data_magLimit245_Taamoli24.fits')
delta_z = np.array(Data['lp_zPDF_u68'] - Data['lp_zPDF_l68'])


#########################################################################
############################# z slices ##################################
#########################################################################

physical_width = 35 # h^-1 Mpc

slice_centers, z_edges = func.redshift_bins(z_min, z_max, physical_width)

z_width = z_edges[:, 1] - z_edges[:, 0]

#########################################################################
############################# Covered Area ##############################
#########################################################################

alpha_min, alpha_max, delta_min, delta_max = 149.397, 150.786, 1.603, 2.816

height = delta_max - delta_min
width = alpha_max - alpha_min
field_area = height * width
holes_area = 0.211 * (height * width) 
corrected_area = field_area - holes_area

wh_ratio = width / height

# circles = np.load('inputs/circles.npy')
circles = np.load(LSS_dir+'inputs/circles_unmasked.npy')

#########################################################################
############################# Weights ###################################
#########################################################################
weight_threshold = 0.05

weights_block = np.load(LSS_dir + 'inputs/WeightsBlock_magcut245.npy')
W = np.load(LSS_dir + 'inputs/WeightsBlock-magcut245-W-th05-normalized.npy')
ind = np.array([np.where(weights_block > weight_threshold)[0], np.where(weights_block > weight_threshold)[1]]).T

gals_bin = []
for i in range(0, len(slice_centers)):
    gals_bin.append(np.unique(ind[ind[:, 1] == i, 0]))
    
inds_th = []
for i in range(0, len(slice_centers)):
    inds_th.append(np.where(weights_block[:, i] >= weight_threshold)[0])


print('\n length of Data:')
print(len(Data))
print('\n Covered Area:')
print('alpha_min-alpha_max: %.3f - %.3f | delta_min-delta_max: %.3f - %.3f \n' %(alpha_min, alpha_max, delta_min, delta_max))
print('height: %.3f - width: %.3f - masked regions: %d' %(height, width, len(circles)))