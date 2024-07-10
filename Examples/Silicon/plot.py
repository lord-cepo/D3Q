import numpy as np
import matplotlib.pyplot as plt


points = np.loadtxt('points.txt')
surf_c = np.loadtxt('surf_C.txt')

mask = (surf_c == -1.)

vectors_true = points[mask]
# Plotting
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

# Plot vectors with mask True
ax.scatter(vectors_true[:, 0], vectors_true[:, 1], vectors_true[:, 2], c='r', marker='o', label='True')
ax.set_xlabel('X axis')
ax.set_ylabel('Y axis')
ax.set_zlabel('Z axis')
ax.legend()

plt.savefig('out.png')