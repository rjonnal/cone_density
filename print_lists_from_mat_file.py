import numpy as np
import scipy.io as sio
import sys,os

# For better or worse the original data is in this horribly ill-conceived
# MATLAB file that I conceived.
# This program sorts out that mess and prints some nice, neat python
# code, which I then used to embed the data in __init__.py, to avoid
# having to keep a separage data file.
# Important conventions in the output: nasal and superior are negative
# and temporal and inferior are positive. All eccentricity values
# are in mm, even though the ConeDensityInterpolator.get_density_and_rowspacing
# function in __init__.py takes eccentricities in degrees.

mat_contents = sio.loadmat('./curcio_data.mat')
data = mat_contents['curcio_data']
sy,sx = data.shape

directions = {}
directions['superior'] = mat_contents['superior'][0][0]-1
directions['inferior'] = mat_contents['inferior'][0][0]-1
directions['nasal'] = mat_contents['nasal'][0][0]-1
directions['temporal'] = mat_contents['temporal'][0][0]-1
locations = ['central','peripheral']
positive_direction_names = ['temporal','inferior']
negative_direction_names = ['nasal','superior']
orientation_keys = ['horizontal','vertical']

central_index = 0
peripheral_index = 1

for pdn,ndn,ik in zip(positive_direction_names,negative_direction_names,orientation_keys):
    positive_direction_index = directions[pdn]
    negative_direction_index = directions[ndn]

    # pool the central and peripheral curves for interpolators
    # start with the positives:
    positive_central_ecc_mm = np.array(data[positive_direction_index,central_index][0],dtype=np.float)[:,0]
    positive_central_density = np.array(data[positive_direction_index,central_index][1],dtype=np.float)[:,0]
    # start indexing at 1 to skip duplicate 1 mm measurements
    positive_peripheral_ecc_mm = np.array(data[positive_direction_index,peripheral_index][0],dtype=np.float)[1:,0] 
    positive_peripheral_density = np.array(data[positive_direction_index,peripheral_index][1],dtype=np.float)[1:,0]
    # now move to the negatives, making sure to reverse signs for ecc:
    # start indexing at 1 to skip duplicate 1 mm measurements
    negative_central_ecc_mm = -np.array(data[negative_direction_index,central_index][0],dtype=np.float)[1:,0]
    negative_central_density = np.array(data[negative_direction_index,central_index][1],dtype=np.float)[1:,0]
    # start indexing at 1 to skip duplicate 1 mm measurements
    negative_peripheral_ecc_mm = -np.array(data[negative_direction_index,peripheral_index][0],dtype=np.float)[1:,0] 
    negative_peripheral_density = np.array(data[negative_direction_index,peripheral_index][1],dtype=np.float)[1:,0]

    ecc_mm = np.hstack((negative_peripheral_ecc_mm[::-1],negative_central_ecc_mm[::-1],positive_central_ecc_mm,positive_peripheral_ecc_mm))
    density = np.hstack((negative_peripheral_density[::-1],negative_central_density[::-1],positive_central_density,positive_peripheral_density))
    def print_vec(vec,vec_name):
        print '%s = ['%vec_name.upper(),
        for item in vec[:-1]:
            print '%0.1f,'%item,
        print '%0.1f ]'%vec[-1]
        print
    print_vec(ecc_mm,'%s_ecc_mm'%ik)
    print_vec(density,'%s_density'%ik)
