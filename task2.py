import numpy as np
import matplotlib.pyplot as plt

# Load the data from the files
data_with = np.loadtxt('data_t2with.txt', skiprows=1)
data_without = np.loadtxt('data_t2without.txt', skiprows=1)

# Extract columns for step, total energy, and temperature
step_with = data_with[:, 0]
etot_with = data_with[:, 4]
temp_with = data_with[:, 5] * 298

step_without = data_without[:, 0]
etot_without = data_without[:, 4]
temp_without = data_without[:, 5] * 298

# Create subplots
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8), sharex=True)

# Plot total energy comparison
ax1.plot(step_without, etot_without, label='Without Thermostat', color='orange')
ax1.plot(step_with, etot_with, label='With Thermostat', color='blue')
ax1.set_ylabel('Total Energy [ε]', fontsize=12)
ax1.legend(fontsize=12)
ax1.set_title('Total Energy Comparison', fontsize=14)

# Plot temperature comparison
ax2.plot(step_without, temp_without, label='Without Thermostat', color='orange')
ax2.plot(step_with, temp_with, label='With Thermostat', color='blue')
ax2.set_ylabel('Temperature [K]', fontsize=12)
ax2.set_xlabel('Time Step [-]', fontsize=12)
ax2.legend(fontsize=12)
ax2.set_title('Temperature Comparison', fontsize=14)

# Adjust layout
plt.tight_layout()
plt.show()