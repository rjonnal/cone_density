import numpy as np
import sys,os
from scipy.interpolate import interp1d
from matplotlib import pyplot as plt

MM_PER_DEG = 0.3
MAX_ECC_MM = 18.0

HORIZONTAL_ECC_MM = [ -22.0, -21.0, -20.0, -19.0, -18.0, -17.0, -16.0, -15.0, -14.0, -13.0, -12.0, -11.0, -10.0, -9.0, -8.0, -7.0, -6.0, -5.0, -3.0, -2.0, -1.0, -0.9, -0.8, -0.7, -0.6, -0.5, -0.4, -0.3, -0.2, -0.2, -0.1, -0.1, 0.0, 0.1, 0.1, 0.2, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 17.0, 18.0 ]

HORIZONTAL_DENSITY = [ 5200.0, 5000.0, 5300.0, 4900.0, 4800.0, 4800.0, 4800.0, 4800.0, 4900.0, 5000.0, 5100.0, 5400.0, 5500.0, 5900.0, 6200.0, 6300.0, 6600.0, 7100.0, 8700.0, 12400.0, 21000.0, 22000.0, 23000.0, 26000.0, 31000.0, 36000.0, 42000.0, 53000.0, 73000.0, 93000.0, 115000.0, 154000.0, 196000.0, 161000.0, 121000.0, 98000.0, 79000.0, 58000.0, 45000.0, 38000.0, 34000.0, 29000.0, 24000.0, 22000.0, 21000.0, 11600.0, 9000.0, 6900.0, 5700.0, 5100.0, 4800.0, 4400.0, 4200.0, 3900.0, 3700.0, 3600.0, 3400.0, 3400.0, 3400.0, 3200.0, 3100.0, 3400.0 ]

VERTICAL_ECC_MM = [ -20.0, -19.0, -18.0, -17.0, -16.0, -15.0, -14.0, -13.0, -12.0, -11.0, -10.0, -9.0, -8.0, -7.0, -6.0, -5.0, -4.0, -3.0, -2.0, -1.0, -0.9, -0.8, -0.7, -0.6, -0.5, -0.4, -0.3, -0.2, -0.2, -0.1, -0.1, 0.0, 0.1, 0.1, 0.2, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 17.0, 18.0, 19.0 ]

VERTICAL_DENSITY = [ 3300.0, 4100.0, 3900.0, 3900.0, 3800.0, 3900.0, 3900.0, 4000.0, 4000.0, 4200.0, 4300.0, 4300.0, 4500.0, 4700.0, 4900.0, 5300.0, 6000.0, 7600.0, 10100.0, 16000.0, 17000.0, 19000.0, 21000.0, 24000.0, 28000.0, 36000.0, 47000.0, 66000.0, 84000.0, 106000.0, 141000.0, 198000.0, 147000.0, 108000.0, 86000.0, 69000.0, 51000.0, 40000.0, 33000.0, 28000.0, 24000.0, 21000.0, 19000.0, 17000.0, 10400.0, 7600.0, 6200.0, 5500.0, 5100.0, 4800.0, 4500.0, 4300.0, 4100.0, 3900.0, 3800.0, 3700.0, 3700.0, 3800.0, 3800.0, 3900.0, 3900.0, 4700.0 ]

class ConeDensityInterpolator:
    
    def __init__(self):

        self.interpolator_keys = ['horizontal','vertical']
        self.interpolators = {}
        self.interpolators['horizontal'] = interp1d(HORIZONTAL_ECC_MM,HORIZONTAL_DENSITY)
        self.interpolators['vertical'] = interp1d(VERTICAL_ECC_MM,VERTICAL_DENSITY)
        
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
        
        x = x_deg*MM_PER_DEG
        y = y_deg*MM_PER_DEG

        # magnitude of position vector:
        r = np.sqrt(x**2+y**2)

        # 
        if r>MAX_ECC_MM:
            sys.exit('ConeDensityInterpolator error: request for ecc of %0.1f mm exceeds limit of %0.1f.'%(r,MAX_ECC_MM))

        hdensity = self.interpolators['horizontal'](r)
        vdensity = self.interpolators['vertical'](r)

        # vweight is proportional to the inverse tangent of vertical/horizontal components
        # check to make sure the horizontal component x is non-zero; o.w. inverse tangent undefined
        if x:
            vweight = np.arctan(np.abs(y/x))*2/np.pi
        else:
            vweight = 1.0
        hweight = 1.0 - vweight

        # density is weighted sum of vertical and horizontal densities
        density = hdensity*hweight+vdensity*vweight

        # compute the row spacing 
        row_spacing = 0.9306 / np.sqrt(density)*1e-3
        return density,row_spacing
        
    def test(self):
        cdi = ConeDensityInterpolator()
        N = 129
        lim = 10 # compute over a 20x20 deg region
        xs = np.linspace(-lim,lim,N)
        ys = np.linspace(-lim,lim,N)

        d = np.zeros((N,N))
        rs_m = np.zeros((N,N))

        for xidx,x in enumerate(xs):
            for yidx,y in enumerate(ys):
                d[yidx,xidx],rs_m[yidx,xidx] = self.get_density_and_rowspacing(x,y)

        # convert row spacing in meters to microns for easier reading of colorbar
        rs_microns = rs_m*1e6

        plt.figure()
        plt.imshow(d,interpolation='none',extent=[-lim,lim,-lim,lim])
        plt.xlabel('eccentricity (deg)')
        plt.ylabel('eccentricity (deg)')
        plt.title('cone density ($mm^{-2}$)')
        plt.savefig('./maps/density.png')
        plt.colorbar()

        plt.figure()
        plt.imshow(rs_microns,interpolation='none',extent=[-lim,lim,-lim,lim])
        plt.xlabel('eccentricity (deg)')
        plt.ylabel('eccentricity (deg)')
        plt.title('cone row spacing ($\mu m$)')
        plt.colorbar()
        plt.savefig('./maps/row_spacing.png')
        plt.show()

if __name__=='__main__':
    cdi = ConeDensityInterpolator()
    cdi.test()
