import numpy as np
import matplotlib.pyplot as plt
import matplotlib.pylab as plab
import scipy.stats
from numpy import *
import matplotlib.animation as manimation
import time
import sys
from mpl_toolkits.mplot3d import Axes3D

start_time = time.time()

grid_data_file = 'grid_data.txt'
nx, ny, Lx, Ly = np.genfromtxt(grid_data_file, unpack=True)

nx = int(nx)
ny = int(ny)

num_data_file = 'numerical_temperature_data' + '.txt'
X, Y, T = np.genfromtxt(num_data_file, unpack=True)

X_2d = np.zeros((nx, ny))
Y_2d = np.zeros((nx, ny))
T_2d = np.zeros((nx, ny))

num_data_file = 'analytical_temperature_data' + '.txt'
X_ana, Y_ana, T_ana = np.genfromtxt(num_data_file, unpack=True)

X_2d_ana = np.zeros((nx, ny))
Y_2d_ana = np.zeros((nx, ny))
T_2d_ana = np.zeros((nx, ny))

for i in range(0, nx):
    for j in range(0, ny):
        num = j*nx+i
        X_2d[i][j] = X[num]
        Y_2d[i][j] = Y[num]
        T_2d[i][j] = T[num]

for i in range(0, nx):
    for j in range(0, ny):
        num = j*nx+i
        X_2d_ana[i][j] = X_ana[num]
        Y_2d_ana[i][j] = Y_ana[num]
        T_2d_ana[i][j] = T_ana[num]
        
X_2d = np.array(X_2d)
Y_2d = np.array(Y_2d)
T_2d = np.array(T_2d)

fig = plt.figure()
ax = fig.gca(projection='3d')
surf = ax.plot_surface(X_2d, Y_2d, T_2d)
#surf = ax.plot_surface(X_2d_ana, Y_2d_ana, T_2d_ana)

ax.set_xlabel('X')
ax.set_xlim3d(0, 1)
ax.set_ylabel('Y')
ax.set_ylim3d(0, 1)
ax.set_zlabel('T')

plt.show()
