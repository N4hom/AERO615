import matplotlib.pyplot as plt
import pandas as pd
import numpy as np


def plotForces(filename):

    data = pd.read_csv(filename)

    # Remove whitespace from column names
    data.columns = data.columns.str.strip()

    # Extract the time steps assuming it's based on the index values of data
    time_steps = data.index

    # Plotting Fx_bottom and Fy_bottom
    plt.figure(figsize=(10, 6))
    plt.plot(time_steps, data['Fx_bottom'], label='Fx_bottom', marker='o', linestyle='-')
    # plt.plot(time_steps, data['Fy_bottom'], label='Fy_bottom', marker='x', linestyle='-')

    plt.tick_params(axis='both', which='major', labelsize=13)

    plt.title('Bottom Forces Over Time')
    plt.xlabel('Time Step')
    plt.ylabel('Force')
    plt.legend()
    plt.grid(True)
    plt.show()


def plotResidual(filename , Mach, mesh):
    # Read the CSV file
    fileData = pd.read_csv(filename, header=None, delimiter=r"\s*,\s*", engine='python')

    # Assign column names (optional, if your file doesn't include headers)
    fileData.columns = ['Iterations', 'Mass', 'x-Momentum', 'y-Momentum', 'Energy']


    # Plotting the residuals with a logarithmic y-axis
    plt.figure(figsize=(10, 6))
    for column in fileData.columns[1:]:
        plt.plot(fileData['Iterations'], np.log10(fileData[column]), label=column)

    plt.tick_params(axis='both', which='major', labelsize=13)

    plt.xlabel('Number of Iterations' , fontsize = 17)
    plt.ylabel('Residuals (L-2 norm)' , fontsize = 17)
    # plt.yscale('log')  # Setting y-axis to logarithmic scale
    plt.title( 'Mach ' + str(Mach) + ' , ' +  mesh + ' mesh' , fontsize = 17)
    plt.legend(fontsize = 17)
    plt.grid(True)
    plt.show()

# Path to your CSV file
fileCoarseMach0_3 = 'results/residuals_coarse_0.300000.csv'
fileCoarseMach0_5 = 'results/residuals_coarse_0.500000.csv'
fileCoarseMach0_7 = 'results/residuals_coarse_0.700000.csv'
fileMediumMach0_3 = 'results/residuals_medium_0.300000.csv'
fileFineMach0_3 = 'results/residuals_fine_0.300000.csv'


gridSize = np.array([15*30 , 30*60 , 60*120])
forceBottom_x_03 = np.array([171.433 , 19.6695 ])
forceBottom_y_03 = np.array([-88360.8 , -91407.1 , -91356.6])

forceBottom_x_05 = np.array([343.898 , 42.968 , 9.4043 ])
forceBottom_y_05 = np.array([-70611.4 , -74976.3 , -74835])

forceBottom_x_07 = np.array([2809.21 , 2531.69 , 2498.29 ])
forceBottom_y_07 = np.array([-57298.2 , -55620 , -55728.3])


# MACH = 0.3
fig, axs = plt.subplots(2, 1)  # 2 rows, 1 column
# Plot on the first subplot
axs[0].plot(gridSize[:-1], forceBottom_x_03, 'b-o')  # 'b-o' is blue color with circle markers
axs[0].set_title('Mach = 0.3' , fontsize=17)
axs[0].set_ylabel('x-Force (N/m)', fontsize=17)
axs[0].grid(True)
axs[0].tick_params(axis='both', which='major', labelsize=13)

# Plot on the second subplot
axs[1].plot(gridSize, forceBottom_y_03, 'r-o')  # 'r-o' is red color with circle markers
axs[1].set_ylabel('y-Force (N/m)' , fontsize=17)
axs[1].set_xlabel('Number of cells' , fontsize=17)
axs[1].grid(True)
axs[1].tick_params(axis='both', which='major', labelsize=13)
# Improve the spacing between subplots
plt.tight_layout()
plt.show()


# MACH = 0.5
fig, axs = plt.subplots(2, 1)  # 2 rows, 1 column
# Plot on the first subplot
axs[0].plot(gridSize, forceBottom_x_05, 'b-o')  # 'b-o' is blue color with circle markers
axs[0].set_title('Mach = 0.5' , fontsize=17)
axs[0].set_ylabel('x-Force (N/m)', fontsize=17)
axs[0].grid(True)
axs[0].tick_params(axis='both', which='major', labelsize=13)

axs[1].plot(gridSize, forceBottom_y_05, 'r-o')  # 'r-o' is red color with circle markers
axs[1].set_ylabel('y-Force (N/m)' , fontsize=17)
axs[1].set_xlabel('Number of cells' , fontsize=17)
axs[1].grid(True)
axs[1].tick_params(axis='both', which='major', labelsize=13)
# Improve the spacing between subplots
plt.tight_layout()


# MACH = 0.7
fig, axs = plt.subplots(2, 1)  # 2 rows, 1 column
# Plot on the first subplot
axs[0].plot(gridSize, forceBottom_x_07, 'b-o')  # 'b-o' is blue color with circle markers
axs[0].set_title('Mach = 0.7' , fontsize=17)
axs[0].set_ylabel('x-Force (N/m)', fontsize=17)
axs[0].grid(True)
axs[0].tick_params(axis='both', which='major', labelsize=13)

axs[1].plot(gridSize, forceBottom_y_07, 'r-o')  # 'r-o' is red color with circle markers
axs[1].set_ylabel('y-Force (N/m)' , fontsize=17)
axs[1].set_xlabel('Number of cells' , fontsize=17)
axs[1].grid(True)
axs[1].tick_params(axis='both', which='major', labelsize=13)
# Improve the spacing between subplots
plt.tight_layout()

# Show the plot
plt.show()



# plt.plot(gridSize , forceBottom_y_03)

#   Mach = 0.3    force_x    force_y
#   coarse        171.433,  -88360.8,
#   medium        19.6695,  -91407.1,
#   fine

#   Mach = 0.5    force_x    force_y
#   coarse        343.898,  -70611.4, 
#   medium        42.968,   -74976.3
#   fine          9.4043    -74835

#   Mach = 0.7    force_x    force_y
#   coarse        2809.21,  -57298.2
#   medium        2531.69,  -55620
#   fine          2498.29   -55728.3


plotResidual("results/residuals_coarse_0.300000_.csv", 0.3 , 'coarse')
plotResidual("results/residuals_coarse_0.500000.csv", 0.5 , 'coarse')
plotResidual("results/residuals_coarse_0.700000.csv", 0.7 , 'coarse')
plotResidual("results/residuals_medium_01_008_0.300000.csv", 0.3 , 'medium')
plotResidual("results/residuals_medium_01_008_0.500000.csv", 0.5 , 'medium')
plotResidual("results/residuals_medium_01_008_0.700000.csv", 0.7 , 'medium')
plotResidual("results/residuals_medium_01_008_0.700000.csv", 0.7 , 'medium')
plotResidual("results/residuals_fine_1_005_0.300000.csv", 0.3 , 'fine')
plotResidual("results/residuals_fine_nu2_0.1_0.500000.csv", 0.5 , 'fine')
plotResidual("results/residuals_fine_0.3_0.005_0.700000.csv", 0.7 , 'fine')
plotForces("results/forces_fine_nu2_0.1_0.500000.csv")


# Path to your CSV file

# Read the CSV file
coarseMach0_3 = pd.read_csv(fileCoarseMach0_3, header=None, delimiter=r"\s*,\s*", engine='python')

# Assign column names (optional, if your file doesn't include headers)
coarseMach0_3.columns = ['Iterations', 'Mass', 'x-Momentum', 'y-Momentum', 'Energy']


# Plotting the residuals with a logarithmic y-axis
plt.figure(figsize=(10, 6))
for column in coarseMach0_3.columns[1:]:
    plt.plot(coarseMach0_3['Iterations'], coarseMach0_3[column], label=column)

plt.xlabel('Number of iterations' , fontsize = 17)
plt.ylabel('Residuals (L-2 norm)' , fontsize = 17)
plt.yscale('log')  # Setting y-axis to logarithmic scale
plt.title('Mach = 0.3, coarse mesh' , fontsize = 17)
plt.legend()
plt.grid(True)
plt.show()



# Read the CSV file
coarseMach0_5 = pd.read_csv(fileCoarseMach0_5, header=None, delimiter=r"\s*,\s*" ,engine = 'python')

# Assign column names (optional, if your file doesn't include headers)
coarseMach0_5.columns = ['Iterations', 'Mass', 'x-Momentum', 'y-Momentum', 'Energy']


# Plotting the residuals with a logarithmic y-axis
plt.figure(figsize=(10, 6))
for column in coarseMach0_5.columns[1:]:
    plt.plot(coarseMach0_5['Iterations'], coarseMach0_5[column], label=column)

plt.xlabel('Number of Iterations' , fontsize = 17)
plt.ylabel('Residual Values (L-2 norm)' , fontsize = 17)
plt.yscale('log')  # Setting y-axis to logarithmic scale
plt.title('Mach = 0.5, coarse mesh' , fontsize = 17)
plt.legend()
plt.grid(True)
plt.show()


# Read the CSV file
coarseMach0_7 = pd.read_csv(fileCoarseMach0_7, header=None, delimiter=r"\s*,\s*" ,engine = 'python')

# Assign column names (optional, if your file doesn't include headers)
coarseMach0_7.columns = ['Iterations', 'Mass', 'x-Momentum', 'y-Momentum', 'Energy']


# Plotting the residuals with a logarithmic y-axis
plt.figure(figsize=(10, 6))
for column in coarseMach0_7.columns[1:]:
    plt.plot(coarseMach0_7['Iterations'], coarseMach0_7[column], label=column)

plt.xlabel('Number of Iterations' , fontsize = 17)
plt.ylabel('Residual Values (L-2 norm)' , fontsize = 17)
plt.yscale('log')  # Setting y-axis to logarithmic scale
plt.title('Mach = 0.7, coarse mesh' , fontsize = 17)
plt.legend()
plt.grid(True)
plt.show()


# Read the CSV file
mediumMach0_3 = pd.read_csv(fileMediumMach0_3, header=None, delimiter=r"\s*,\s*", engine='python')

# Assign column names (optional, if your file doesn't include headers)
mediumMach0_3.columns = ['Iterations', 'Mass', 'x-Momentum', 'y-Momentum', 'Energy']


# Plotting the residuals with a logarithmic y-axis
plt.figure(figsize=(10, 6))
for column in mediumMach0_3.columns[1:]:
    plt.plot(mediumMach0_3['Iterations'], mediumMach0_3[column], label=column)

plt.xlabel('Number of Iterations' , fontsize = 17)
plt.ylabel('Residual Values (L-2 norm)' , fontsize = 17)
plt.yscale('log')  # Setting y-axis to logarithmic scale
plt.title('Mach = 0.3, medium mesh' , fontsize = 17)
plt.legend()
plt.grid(True)
plt.show()




# Read the CSV file
fineMach0_3 = pd.read_csv(fileFineMach0_3, header=None, delimiter=r"\s*,\s*", engine='python')

# Assign column names (optional, if your file doesn't include headers)
fineMach0_3.columns = ['Iterations', 'Mass', 'x-Momentum', 'y-Momentum', 'Energy']


# Plotting the residuals with a logarithmic y-axis
plt.figure(figsize=(10, 6))
for column in fineMach0_3.columns[1:]:
    plt.plot(fineMach0_3['Iterations'], fineMach0_3[column], label=column)

plt.xlabel('Number of Iterations' , fontsize = 17)
plt.ylabel('Residual Values (L-2 norm)' , fontsize = 17)
plt.yscale('log')  # Setting y-axis to logarithmic scale
plt.title('Mach = 0.3, medium mesh' , fontsize = 17)
plt.legend()
plt.grid(True)
plt.show()
