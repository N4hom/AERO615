import matplotlib.pyplot as plt
import pandas as pd


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

    plt.title('Bottom Forces Over Time')
    plt.xlabel('Time Step')
    plt.ylabel('Force')
    plt.legend()
    plt.grid(True)
    plt.show()


def plotResidual(filename):
    # Read the CSV file
    fileData = pd.read_csv(filename, header=None, delimiter=r"\s*,\s*", engine='python')

    # Assign column names (optional, if your file doesn't include headers)
    fileData.columns = ['Iterations', 'Mass', 'x-Momentum', 'y-Momentum', 'Energy']


    # Plotting the residuals with a logarithmic y-axis
    plt.figure(figsize=(10, 6))
    for column in fileData.columns[1:]:
        plt.plot(fileData['Iterations'], fileData[column], label=column)

    plt.xlabel('Number of Iterations' , fontsize = 17)
    plt.ylabel('Residual Values (L-2 norm)' , fontsize = 17)
    plt.yscale('log')  # Setting y-axis to logarithmic scale
    plt.title('Mach = 0.3, coarse mesh' , fontsize = 17)
    plt.legend()
    plt.grid(True)
    plt.show()

# Path to your CSV file
fileCoarseMach0_3 = 'results/residuals_coarse_0.300000.csv'
fileCoarseMach0_5 = 'results/residuals_coarse_0.500000.csv'
fileCoarseMach0_7 = 'results/residuals_coarse_0.700000.csv'
fileMediumMach0_3 = 'results/residuals_medium_0.300000.csv'
fileFineMach0_3 = 'results/residuals_fine_0.300000.csv'


plotForces("results/forces_medium_nu2_0.3_0.300000.csv")
plotResidual("results/residuals_medium_nu2_0.3_0.300000.csv")


# Path to your CSV file

# Read the CSV file
coarseMach0_3 = pd.read_csv(fileCoarseMach0_3, header=None, delimiter=r"\s*,\s*", engine='python')

# Assign column names (optional, if your file doesn't include headers)
coarseMach0_3.columns = ['Iterations', 'Mass', 'x-Momentum', 'y-Momentum', 'Energy']


# Plotting the residuals with a logarithmic y-axis
plt.figure(figsize=(10, 6))
for column in coarseMach0_3.columns[1:]:
    plt.plot(coarseMach0_3['Iterations'], coarseMach0_3[column], label=column)

plt.xlabel('Number of Iterations' , fontsize = 17)
plt.ylabel('Residual Values (L-2 norm)' , fontsize = 17)
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
