import numpy as np
import matplotlib.pyplot as plt

x = np.loadtxt('x', delimiter=',')
y = np.loadtxt('y', delimiter=',')

print(x)
print(y)

x = np.transpose(x)
y = np.transpose(y)

plt.scatter(x,y)
plt.grid()
plt.show()