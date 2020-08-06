#! /usr/bin/env python3
# Build up empirical prior from ADNI images and save results (with masking)

# File: empirical_scan.py
# Author: Ross Boylan <ross.boylan@ucsf.edu>
# Created: 2020-08-05
#
# INPUTS
# All are in ..\ExternalData\ycobigo\round3\
#
# TPM.nii   segmentation map (-> mask)
# ADNI\
#      classified-subset.csv list of images and dx (excludes AD)
#      the image files  themselves
#
# OUTPUTS
# ep3.npz  The empirical prior
# vox3.npz sample  of voxels
#
## I think the code below means you can run from anywhere at or under \Kornak\bifs

ADNI = r"ExternalData\ycobigo\round3\ADNI"

import nibabel
import numpy as np
from pathlib import Path
import sys
import time

## Establish input and output paths

TOPDIR = Path.cwd()
while TOPDIR.name != "Kornak":
    BIFSDIR = TOPDIR
    TOPDIR = TOPDIR.parent
sys.path.insert(0, str(BIFSDIR))

ADNIDIR = TOPDIR / ADNI
SEGFILE = TOPDIR / "ExternalData/ycobigo/round3/TPM.nii"

EPNAME = "ep3.npz"
VOXNAME = "vox3.npz"

import BIFS
from BIFS.bifs_util.EmpiricalScanner import FeedScanner as Scnr


def mypaths():
    "Return a generator of the path names as <str> to operate on"
    with open(ADNIDIR / "classified-subset.csv", "rt") as fin:
        n=0
        for line in fin:
            n += 1
            if (n % 10) == 0:
                print("{} ".format(n),end='', flush=True)
            if n == 1:
                # skip initial header
                continue
            fn = line[:(line.find(","))]
            yield str(ADNIDIR / fn)


def get_mask():
    "Return a boolean mask for the images.  True means ignore this voxel"
    segs = nibabel.load(str(SEGFILE)).get_fdata()
    # True for voxels to set to 0
    mask = segs[..., 3:6].sum(3) > 0.80
    return mask

def go():
    mask = get_mask()
    samp_frac = 0.05
    unmask_frac = 1.0 - np.sum(mask)/mask.size
    print("Will sample {:.1%} of total voxels.  Since {:.1%} of total voxels are unmasked,".format(
        samp_frac, unmask_frac))
    start_t = time.perf_counter()
    print("will sample {:.1%} of unmasked voxels.".format(samp_frac/unmask_frac))
    s = Scnr(mypaths(), samp_frac, image_mask=mask)
    end_t = time.perf_counter()
    wall_t = end_t - start_t
    print("\n{:7.1f}s/image.  {} images in {:10.1f} seconds.".format(wall_t/s.nImages, s.nImages, wall_t))
    # and then output results to disk.
    np.savez(EPNAME, mean=s.mns, sd=s.sds)
    print("{} gets 2 {} matrices named mean and sd".format(EPNAME, s.mns.shape))
    np.savez(VOXNAME, samp=np.array(s.vox))
    print("{} gets {:,d} sample voxels named samp".format(VOXNAME, len(s.vox)))

if __name__ == "__main__":
    go()