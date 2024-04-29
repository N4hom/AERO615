import matplotlib.pyplot as plt
import pandas as pd

# Path to your CSV file
fileCoarseMach0_3 = 'results/residuals_coarse_0.300000.csv'
fileCoarseMach0_5 = 'results/residuals_coarse_0.500000.csv'
fileCoarseMach0_7 = 'results/residuals_coarse_0.700000.csv'

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

fileCoarseMach0_3 = 'results/residuals_coarse_0.300000.csv'


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