# cone_density

A simple tool for interpolating cone photoreceptor density in the human retina. Using data from Curcio et al.'s 1990 paper *Human photoreceptor topography*.

To use, download and unzip the directory somewhere in your python path. Then, to estimate the cone density and row spacing at a given eccentricity, instantiate a `ConeDensityInterpolator` object and call its `get_density_and_rowspacing` function on the target x- and y-coordinates, specified in degrees of visual angle. The function returns a 2-tuple, `(density, row_spacing)`, where density is cone density in **mm<sup>-2</sup>** and row_spacing is in **m**.

Some conventions:

1. Nasal and superior eccentricities are specified with negative values, while temporal and inferior eccentricities are specified with positive ones.

2. The original paper gives densities in (**mm<sup>-2</sup>**) as a function of **mm**. The `get_density_and_rowspacing` function takes eccentricity in degrees and converts to **mm** using a factor of **300 &mu; m/deg**.

[Curcio, Christine A., et al. "Human photoreceptor topography." *Journal of Comparative Neurology* 292.4 (1990): 497-523.](https://www.ncbi.nlm.nih.gov/pubmed/2324310)
