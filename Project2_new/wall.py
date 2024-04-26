import numpy as np
import matplotlib.pyplot as plt

# Clear all previous plots and close all figures
plt.close('all')

# Define the angle in radians
theta = 30 * np.pi / 180

# Components of vector t
tx = np.cos(theta)
ty = np.sin(theta)

# Plotting vector t
plt.plot([0, 6*tx], [0, 6*ty], '--k')
plt.grid(True)

# Define components of vector v1
u1 = 10 * np.cos(40 * np.pi / 180)
v1 = 10 * np.sin(40 * np.pi / 180)

# Plotting vector v1
plt.plot([0, u1], [0, v1], 'r')

# Calculate components of vector v0
u0 = u1 * (tx**2 - ty**2) + 2 * v1 * tx * ty
v0 = 2 * u1 * tx * ty + v1 * (ty**2 - tx**2)

# Plotting vector v0
plt.plot([0, u0], [0, v0], 'b')

# Set equal scaling
plt.axis('equal')

# Calculating magnitudes of vectors
V1 = np.sqrt(u1**2 + v1**2)
V0 = np.sqrt(u0**2 + v0**2)

# Display the plot
plt.show()
