import numpy as np
import matplotlib.pyplot as plt
import os

for file in [file for file in os.listdir() if file.endswith('-test.dat')]:
    data = np.loadtxt(file)
    ref = np.loadtxt('./reference/' + file)
    if file.startswith('spf'):
        i = 7
    if file.startswith('self'):
        i = 9
    if file.startswith('dos'):
        i = 1
    if np.allclose(data, ref, rtol=1e-5):
        print(f"{file} matches reference.")
    else:
        print(f"---------------------> {file} does NOT match reference.")
        plt.plot(data[:, 0], data[:, i], label='test')
        plt.plot(ref[:, 0], ref[:, i], label='reference', linestyle='dashed')
        plt.title(file)
        plt.xlabel('Energy')
        plt.ylabel('Value')
        plt.legend()
        plt.savefig(file.replace('-test.dat', '.png'), dpi=600)