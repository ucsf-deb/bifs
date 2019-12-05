# A command line type example for using BIFS example

import numpy as np
from scipy import stats
import matplotlib.pyplot as plt
from bifs.bifs import Bifs
import bifs.bifs_util.util as bu

# Make 1D "image"
# Noise standard deviation in image space
noiseSD = 0.5
# Set y as a discrete set of values of even length
z = np.concatenate((np.zeros(14),np.arange(19)+1,np.zeros(7)+10,10-np.arange(10),np.zeros(14)))

y = z + stats.norm.rvs(size=len(z),loc=0.0,scale=noiseSD)

# Create mybifs BIFS object:
mybifs = Bifs()

# Can take a look at what functions and variables are available with, e.g.:
# dir(mybifs)

# Set a few things:
# Prior
mybifs.prior = "Gaussian" # Choices are currently: "Gaussian","Uniform"
# "Gaussian" is actually the default but for illustration...

# Lilelihood
mybifs.likelihood = "Gaussian" # Choices are currently: "Gaussian","Rician" 
# "Gaussian" is actually the default but again for illustration...

# Parameter Space Function
# mybifs.param_func_type = "Inverse Power Decay"
mybifs.param_func_type = "Linear Decay"
# Current choices are: "Inverse Power Decay","Banded Inverse Power Decay",
# "Linear Decay" with default "Inverse Power Decay",
# but again for illustration...

# Can check comments in bifs.py for description of other parametere to set

# Load the image - note, typically just start here re. loading a noisy image
mybifs.load_image(y)

# Run BIFS making sure that the initial image is loaded
if np.isscalar(mybifs.init_image):
  print("Error: Need to load an image into mybifs before running MAP")
else:
  print("Running BIFS_MAP() on image")
  mybifs.BIFS_MAP()

# Take a look at the current paramter function
bu.plot_param_func(mybifs)

# Look at the prior, liklelihood, and posterior at a voxel
bu.voxel_dist(mybifs,[mybifs.mod_image.shape[0]//2],do_plots=True)
  
# Plot the resulting "images", checking again that nothing went wrong...
if np.isscalar(mybifs.init_image):
  print("Error: Need to load an image into mybifs before running MAP")
else:
  # Initial noisy image
  plt.subplot(221)
  plt.axis('off')
  plt.title("Initial Image")
  plt.plot(mybifs.init_image)
  
  # Initial K-Space Image
  plt.subplot(223)
  plt.axis('off')
  plt.title("Initial K-Space Image")
  showim1k = np.roll(np.roll(mybifs.mod_image,mybifs.mod_image.shape[0]//2,0),1)
  plt.plot(np.log(showim1k))
             
  # Final K-Space Image after running BIFS
  plt.subplot(224)
  plt.axis('off')
  plt.title("Final K-Space Image")
  showim2k = np.roll(np.roll(mybifs.bifsk_image,mybifs.bifsk_image.shape[0]//2,0),1)
  plt.plot(np.log(showim2k))
             
  # Final Image after running BIFS
  plt.subplot(222)
  plt.axis('off')
  plt.title("Reconstructed Image")
  plt.plot(mybifs.final_image)

  plt.show()
