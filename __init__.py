import numpy as np
import scipy.io as sio
import sys,os
from scipy.optimize import curve_fit
from scipy.interpolate import interp1d
from matplotlib import pyplot as plt

debug = False

class ConeDensityInterpolator:


    def __init__(self):

        mat_contents = sio.loadmat('./curcio_data.mat')

        data = mat_contents['curcio_data']
        sy,sx = data.shape

        directions = {}
        directions['superior'] = mat_contents['superior'][0][0]-1
        directions['inferior'] = mat_contents['inferior'][0][0]-1
        directions['nasal'] = mat_contents['nasal'][0][0]-1
        directions['temporal'] = mat_contents['temporal'][0][0]-1


        locations = ['central','peripheral']
        direction_names = ['nasal','superior','temporal','inferior']
        positive_direction_names = ['temporal','inferior']
        negative_direction_names = ['nasal','superior']
        self.interpolator_keys = ['horizontal','vertical']
        self.interpolators = {}

        central_index = 0
        peripheral_index = 1


        self.max_ecc_mm = 18.0

        for pdn,ndn,ik in zip(positive_direction_names,negative_direction_names,self.interpolator_keys):
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

            interpolator = interp1d(ecc_mm,density)
            self.interpolators[ik] = interpolator

    def plot(self):
        for ik in self.interpolator_keys:
            plt.figure()
            ior = self.interpolators[ik]
            plt.plot(ior.x,ior.y,'ks')
            plt.xlabel('ecc (mm)')
            plt.ylabel('density (mm^-2)')
            plt.title(ik)
        plt.show()


    def get_density_and_rowspacing(self,x_deg,y_deg):
        x = x_deg*.3
        y = y_deg*.3

        r = np.sqrt(x**2+y**2)
        if r>self.max_ecc_mm:
            sys.exit('ConeDensityInterpolator error: request for ecc of %0.1f mm exceeds limit of %0.1f.'%(r,self.max_ecc_mm))

        hdensity = self.interpolators['horizontal'](r)
        vdensity = self.interpolators['vertical'](r)

        if x:
            vweight = np.arctan(np.abs(y/x))*2/np.pi
        else:
            vweight = 1.0
        hweight = 1.0 - vweight

        density = hdensity*hweight+vdensity*vweight

        row_spacing = 0.9306 / np.sqrt(density)*1e-3
        
        return density,row_spacing
        
if __name__=='__main__':

    cdi = ConeDensityInterpolator()
    N = 129
    lim = 40
    xs = np.linspace(-lim,lim,N)
    ys = np.linspace(-lim,lim,N)

    d = np.zeros((N,N))
    rs = np.zeros((N,N))

    for xidx,x in enumerate(xs):
        print xidx
        for yidx,y in enumerate(ys):
            d[yidx,xidx],rs[yidx,xidx] = cdi.get_density_and_rowspacing(x,y)

    plt.figure()
    plt.imshow(d,interpolation='none')
    plt.colorbar()
    plt.figure()
    plt.imshow(rs,interpolation='none')
    plt.colorbar()

    plt.show()


