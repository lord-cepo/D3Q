import numpy as np
import matplotlib.pyplot as plt
import os

success = 0
for file in [file for file in os.listdir() if file.endswith('-test.dat')]:
    try:
        data = np.loadtxt(file)
    except:
        print(f"---------------------> Error reading {file}. Skipping.")
        continue
    ref = np.loadtxt('./reference/' + file)
    if file.startswith('spf'):
        i = 7
    if file.startswith('self'):
        i = 9
    if file.startswith('dos'):
        i = 1
    print(f"{file} read successfully.")
    if not np.allclose(data, ref, rtol=1e-5):
        print(f"---------------------> {file} does NOT match reference.")
        plt.plot(data[:, 0], data[:, i], label='test')
        plt.plot(ref[:, 0], ref[:, i], label='reference', linestyle='dashed')
        plt.title(file)
        plt.xlabel('Energy')
        plt.ylabel('Value')
        plt.legend()
        plt.savefig(file.replace('-test.dat', '.png'), dpi=600)
        plt.close()
    else:
        success += 1

print()
print("---------------------------")
if(success == 12):
    print("TEST SUCCESS")
else:
    print(f"TEST FAILED: {12 - success} files do not match reference.")
print("---------------------------")

