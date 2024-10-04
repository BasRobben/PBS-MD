import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import make_interp_spline

# Load RDF data from the file
data = np.loadtxt("rdf_data.txt")
bin_centers = data[:, 0]
gij = data[:, 1]

# Create a smooth curve using spline interpolation
x_smooth = np.linspace(bin_centers.min(), bin_centers.max(), 500)
spl = make_interp_spline(bin_centers, gij, k=3)  # k=3 for cubic spline
gij_smooth = spl(x_smooth)

# Create the plot
plt.figure(figsize=(8, 6))
plt.plot(x_smooth, gij_smooth, linestyle='-', color='b')
plt.xlabel('Distance (r/σ)', size = 12)
plt.ylabel('Radial Distribution Function g(r)', size = 12)
plt.grid()
plt.show()