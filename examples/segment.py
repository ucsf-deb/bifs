
#! /usr/bin/python3
# Deal with segmentation of voxels by tissue type

# File: segement.py
# Author: Ross Boylan
# Created: 2020-05-13

import nibabel as nib
import matplotlib.pyplot as plt

SEG_FILE = r"C:\Users\rdboylan\Documents\Kornak\ExternalData\ycobigo\round3\TPM.nii"

s = nib.load(SEG_FILE)
print(s.header)
img = s.get_fdata()
print("Image dimensions are {}".format(img.shape))
# 6 (1-based) is outside head
fix, axes = plt.subplots(1, 6)
for i in range(6):
    axes[i].imshow(img[100, : , :, i])

plt.show()
