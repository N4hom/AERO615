import numpy as np
import matplotlib.pyplot as plt

xc = np.loadtxt('xc', delimiter=',')
yc = np.loadtxt('yc', delimiter=',')

x = np.loadtxt('x', delimiter=',')
y = np.loadtxt('y', delimiter=',')

N = len(x[0]) - 2
M = len(x[:]) - 2

print("N " , N)
print("M " , M)

x = x
y = y

print(xc)


plt.plot(xc,yc ,'.' ,color='b')
plt.title(f'Grid {M}x{N}')
plt.xlabel('x/L', fontsize=14)
plt.ylabel('y/H', fontsize=14)
plt.grid()
plt.show()


plt.plot(x,y , color='b')
x = np.transpose(x)
y = np.transpose(y)
plt.plot(x,y, color='b')
plt.plot(xc,yc , '.' , color='b')
xc = np.transpose(xc)
yc = np.transpose(yc)
plt.plot(xc,yc, '.' ,color='b')
plt.title(f'Grid {M}x{N}')
plt.xlabel('x/L', fontsize=14)
plt.ylabel('y/H', fontsize=14)
plt.grid()
plt.show()