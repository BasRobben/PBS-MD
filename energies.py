import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter

# Load the data from energies.txt, skipping the header row
data = np.loadtxt('data.txt', skiprows=1)

# Assuming the columns are structured as step, time, potential_energy, kinetic_energy, total_energy
step = data[:, 0]
time = data[:, 1]
potential_energy = data[:, 2]
kinetic_energy = data[:, 3]
total_energy = data[:, 4]
temperature = data[:, 5] * 298

def calculate_rmsd(values):
    mean_value = np.mean(values)
    rmsd = np.sqrt(np.mean((values - mean_value) ** 2))
    return rmsd

def calculate_normalized_rmsd(values):
    rmsd = calculate_rmsd(values)
    mean_value = np.mean(values)
    normalized_rmsd = rmsd / mean_value
    return normalized_rmsd

def calculate_mean(values):
    return np.mean(values)

# Calculate RMSD and normalized RMSD for total energy
rmsd_total = calculate_rmsd(total_energy)
normalized_rmsd_total = calculate_normalized_rmsd(total_energy)
mean_kinetic_energy = calculate_mean(kinetic_energy[100:])

print(f"RMSD Total Energy: {rmsd_total}")
print(f"Normalized RMSD Total Energy: {normalized_rmsd_total}")
print(f"Mean Kinetic Energy: {mean_kinetic_energy}")

# Create subplots with shared x-axis
fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(8, 10), sharex=True)

# Scatter plot for potential energy
ax1.plot(step, potential_energy, label='Potential Energy', color='blue')
ax1.set_ylabel('Potential Energy [ε]', fontsize=12)
ax1.legend(fontsize=12)

# Scatter plot for kinetic energy
ax2.plot(step, kinetic_energy, label='Kinetic Energy', color='orange')
ax2.set_ylabel('Kinetic Energy [ε]', fontsize=12)
ax2.legend(fontsize=12)

# Scatter plot for total energy
ax3.plot(step, total_energy, label='Total Energy', color='green')
ax3.set_ylabel('Total Energy [ε]', fontsize=12)
ax3.set_xlabel('Time Step [-]', fontsize=12)
ax3.yaxis.set_major_formatter(ScalarFormatter(useMathText=True))
ax3.yaxis.get_major_formatter().set_useOffset(False)
ax3.legend(fontsize=12)

# # Scatter plot for temperature ratios
# ax3.plot(step, temperature, label='Temperature', color='green')
# ax3.set_ylabel('Temperature [K]', fontsize=12)
# ax3.set_xlabel('Time Step [-]', fontsize=12)
# ax3.legend(fontsize=12)


# Adjust layout
plt.tight_layout()
plt.show()