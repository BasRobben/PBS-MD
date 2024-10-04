import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import maxwell

# Load the data
data = pd.read_csv('velocity_data.txt', sep='\t')

# Separate the data by type
CH3_data = data[data['Type'] == 1]['Velocity']
CH4_data = data[data['Type'] == 0]['Velocity']

# Calculate the most probable speed vp for the Maxwell-Boltzmann distribution
kT = 2.01
m_ch3 = 0.94
m_ch4 = 1
vp_ch3 = np.sqrt(kT / m_ch3)
vp_ch4 = np.sqrt(kT / m_ch4)

# Plot the velocity distribution for CH3
plt.figure(figsize=(6, 8))

# Histogram for CH3 velocities
plt.subplot(2, 1, 1)
plt.hist(CH3_data, bins=30, density=True, alpha=0.6, color='b', label='CH3 Velocity Distribution')
plt.title('Distribution of Velocities (t=5000)', size = 12)

# Fit and plot the Maxwell-Boltzmann distribution for CH3
x = np.linspace(0, max(CH3_data), 100)
plt.plot(x, maxwell.pdf(x, scale=vp_ch3), 'r', label=f'Maxwell-Boltzmann (vp={vp_ch3:.2f})')
plt.ylabel('Probability Density', size = 12)
plt.legend()

# Histogram for CH4 velocities
plt.subplot(2, 1, 2)
plt.hist(CH4_data, bins=30, density=True, alpha=0.6, color='g', label='CH4 Velocity Distribution')

# Fit and plot the Maxwell-Boltzmann distribution for CH4
x = np.linspace(0, max(CH3_data), 100)
plt.plot(x, maxwell.pdf(x, scale=vp_ch4), 'r', label=f'Maxwell-Boltzmann (vp={vp_ch4:.2f})')
# plt.xlabel('Velocity', size = 12)
plt.xlabel('Velocity', size = 12)
plt.ylabel('Probability Density', size = 12)
plt.legend()

plt.tight_layout()
plt.show()