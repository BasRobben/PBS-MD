import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter

# Load data
data = np.loadtxt('data.txt', skiprows=1)
step = data[:, 0]
potential_energy = data[:, 2]
kinetic_energy = data[:, 3]
temperature = data[:, 5] * 298
pressure = data[:, 6]

# Create subplots
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(8, 8), sharex=True)

# # Plot potential energy on the first subplot
# ax1.plot(step, potential_energy, label='Potential Energy', color='blue')
# ax1.set_ylabel('Potential Energy [ε]', fontsize=12, color='blue')
# ax1.tick_params(axis='y', labelcolor='blue')
# ax1.legend(loc='upper left', fontsize=12)
# ax1.legend(loc='upper right', fontsize=12)

# # Plot kinetic energy on the second subplot
# ax2.plot(step, kinetic_energy, label='Kinetic Energy', color='orange')
# ax2.set_ylabel('Kinetic Energy [ε]', fontsize=12, color='orange')
# ax2.set_xlabel('Time Step [-]', fontsize=12)
# ax2.tick_params(axis='y', labelcolor='orange')
# ax2.legend(loc='upper left', fontsize=12)
# ax2.legend(loc='upper right', fontsize=12)

# Plot temperature on the first subplot
ax1.plot(step, temperature, label='Temperature', color='red')
ax1.set_ylabel('Temperature [K]', fontsize=12, color='red')
ax1.tick_params(axis='y', labelcolor='red')
ax1.legend(loc='upper left', fontsize=12)
ax1.legend(loc='upper right', fontsize=12)

# Plot pressure on the second subplot
ax2.plot(step, pressure, label='Pressure', color='green')
ax2.set_ylabel('Pressure [ε/σ³]', fontsize=12, color='green')
ax2.set_xlabel('Time Step [-]', fontsize=12)
ax2.tick_params(axis='y', labelcolor='green')
ax2.legend(loc='upper left', fontsize=12)
ax2.legend(loc='upper right', fontsize=12)

# Adjust layout
plt.tight_layout()
plt.show()