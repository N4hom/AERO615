import numpy as np
import matplotlib.pyplot as plt

x = np.loadtxt('x', delimiter=',')
y = np.loadtxt('y', delimiter=',')

N = len(x[0]) - 2
M = len(x[:]) - 2

print("N " , N)
print("M " , M)

x = x/5
y = y/1

print(x)
print(y)


plt.plot(x,y , color='b')
x = np.transpose(x)
y = np.transpose(y)
plt.plot(x,y, color='b')
plt.title(f'Grid {M}x{N}')
plt.xlabel('x/L', fontsize=14)
plt.ylabel('y/H', fontsize=14)
plt.grid()
plt.show()