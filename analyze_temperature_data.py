import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

grid_data_file = 'grid_data.txt'
nx, ny, Lx, Ly = np.genfromtxt(grid_data_file, unpack=True)

nx = int(nx)
ny = int(ny)

num_data_file = 'numerical_temperature_data' + '.txt'
X, Y, T = np.genfromtxt(num_data_file, unpack=True)

X_2d = np.zeros((nx, ny))
Y_2d = np.zeros((nx, ny))
T_2d = np.zeros((nx, ny))

for i in range(0, nx):
    for j in range(0, ny):
        num = j*nx+i
        X_2d[i][j] = X[num]
        Y_2d[i][j] = Y[num]
        T_2d[i][j] = T[num]

fig = plt.figure()
ax = fig.gca(projection='3d')
surf = ax.plot_surface(X_2d, Y_2d, T_2d)

ax.set_xlabel('X')
ax.set_xlim3d(0, 1)
ax.set_ylabel('Y')
ax.set_ylim3d(0, 1)
ax.set_zlabel('T')

plt.show()
