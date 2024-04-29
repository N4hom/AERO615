import matplotlib.pyplot as plt
import pandas as pd

# Path to your CSV file
fileCoarseMach0_3 = 'residuals_coarse_0.3.csv'

# Read the CSV file
coarseMach0_3 = pd.read_csv(fileCoarseMach0_3, header=None, delimiter=r"\s*,\s*", engine='python')

# Assign column names (optional, if your file doesn't include headers)
data.columns = ['Iterations', 'Residual1', 'Residual2', 'Residual3', 'Residual4']


# Plotting the residuals with a logarithmic y-axis
plt.figure(figsize=(10, 6))
for column in data.columns[1:]:
    plt.plot(data['Iterations'], data[column], label=column)

plt.xlabel('Number of Iterations')
plt.ylabel('Residual Values')
plt.yscale('log')  # Setting y-axis to logarithmic scale
plt.title('Residuals vs. Number of Iterations')
plt.legend()
plt.grid(True)
plt.show()
