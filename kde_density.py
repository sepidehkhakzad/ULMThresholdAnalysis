
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde
from scipy.io import savemat
import scipy.io as sio

def calculate_kde_density(ground_truth, mapX, mapZ):
    """
    Calculate Gaussian KDE for object densities based on x and z values.

    Args:
        ground_truth (np.ndarray): The gt data where columns 3 and 5 contain z and x values respectively.
                                this is set to fit the way the gt data is stored in gt_scat_inside1.txt
        mapX (np.ndarray): Array for the X-axis values to evaluate KDE.
        mapZ (np.ndarray): Array for the Z-axis values to evaluate KDE.

    Returns:
        None
    """
    
    x_values = ground_truth[:, 2]  # z is in the 3rd column 
    z_values = ground_truth[:, 4]  # x is in the 5th column 

    # Stack x and z values for 2D KDE calculation
    values = np.vstack([x_values, z_values])

    kde = gaussian_kde(values)

    xx, zz = np.meshgrid(mapX, mapZ)
    grid_coords = np.vstack([xx.ravel(), zz.ravel()])
    density = kde(grid_coords).reshape(xx.shape)

    return density

# Load metadata
metadata = sio.loadmat('./metadata_simu1.mat')
mapX = metadata['PxSet']['mapX'][0][0][0]
mapZ = metadata['PxSet']['mapZ'][0][0][0]
dx = metadata['PxSet']['dx'][0][0][0][0]
dz = metadata['PxSet']['dz'][0][0][0][0]
cf = metadata['SimSet']['centre_frequency'][0][0][0][0]

lambdaa = 1540 / cf
tol = lambdaa/2

# Load ground truth data
gt_loc = np.loadtxt('./gt_scat_inside1.txt', delimiter=',')

density = calculate_kde_density(gt_loc, mapX, mapZ)
savemat("kde_density_simu1.mat", {
        'density': density,
        'mapX': mapX,
        'mapZ': mapZ
    })

# Plot the KDE density
plt.figure(figsize=(8, 6))
plt.imshow(density, origin='lower', aspect='auto',
           extent=[mapX.min(), mapX.max(), mapZ.min(), mapZ.max()],
           cmap='viridis')

plt.colorbar(label='Density')
plt.title('Gaussian KDE Density of Objects')
plt.xlabel('X')
plt.ylabel('Z')
plt.show()
