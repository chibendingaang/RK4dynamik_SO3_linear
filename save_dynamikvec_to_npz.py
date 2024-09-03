import numpy as np
import sys

filename = sys.argv[1]

# Read the temporary file
with open(f'./xpdrvn_NNN/L512/2emin3/J1only/{filename}.dat', 'rb') as f:
    # Read dimensions and data
    nx, ny, nz = np.fromfile(f, dtype=np.int32, count=3)
    data = np.fromfile(f, dtype=np.float128)

# Reshape the array
array_3d = data.reshape((nx, ny, nz), order='F')

print("Array shape, contents: ", array_3d.shape)
print(array_3d)

# Save as compressed npz
np.savez_compressed(f'{filename}.npz', array=array_3d, compression=4)
print(f'Array of shape {array_3d.shape} saved to output.npz with compression level 4')
