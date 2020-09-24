#! /usr/bin/env python3
#  File: investigate07.py
#  Author: Ross Boylan
#  Created: 2020-09-23

# Do our local PETs and ADNI PETs have same distribution of focal intensities?
#
# INPUTS
#     vox3.npz   voxel intensities in focal area of ADNI PETS (large sample)
#     local PET images
#     TPM.nii segment map
#     subsamples.feather  file names and dx for input images
#         each row in subsample will produce a pdf
#     image files pointed to by subsample, MRI and PET

# OUTPUTS
#      investigate07.pdf  plots

import sys
import nibabel
import numpy as np
import pandas as pd
from pathlib import Path
from scipy import stats

import matplotlib
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.figure import Figure
from matplotlib.colors import NoNorm
from matplotlib import cm
import matplotlib.pyplot as plt

## Establish input and output paths

TOPDIR = Path.cwd()
while TOPDIR.name != "Kornak":
    BIFSDIR = TOPDIR
    TOPDIR = TOPDIR.parent
sys.path.insert(0, str(BIFSDIR))

VOXNAME = "vox3.npz"
SEGFILE = TOPDIR / "ExternalData/ycobigo/round3/TPM.nii"
SUBFILE = TOPDIR / "RKornak" / "subsamples.feather"
ODIR = TOPDIR 
if not ODIR.exists():
    ODIR.mkdir()

# ids are floating point without the assist
subsample = pd.read_feather(SUBFILE).astype({"id": int})
# to enhance blinding, reorder
subsample = subsample.sort_values(by=["id"])

segs = nibabel.load(str(SEGFILE)).get_fdata()
seg_names = ("GM",  "WM", "CSP", "SK", "BH", "OUT")
# segnames[i] is the probabilities given by segs[,,,i]

# True for voxels to set to 0
mask = segs[..., 3:6].sum(3) > 0.80
focus = np.logical_not(mask)

# if we parallelize, easier to load this once
refVoxels = np.load(VOXNAME)["samp"]

def print_summary(x, desc):
    ptiles = [1, 2.5, 5, 10, 25, 50, 75, 90, 95, 97.5, 99]
    print("Descriptive stats for {}\n {}\n Quantiles {} ({})".format(desc,
                                                                     stats.describe(x),
                                                                     stats.scoreatpercentile(x, ptiles).round(2),
                                                                     ptiles))

def do_one(i: int, pp):
    """plot subject i's PET distn vs ADNI
    pp is the backend"""
    global subsample, focus, refVoxels
    aPET = nibabel.load(subsample.at[i, "fn.pet"]).get_fdata()
    aPET[np.isnan(aPET)] = 0.0  # mirroring the processing of bifs.load_image
    # convert to list of focal voxels
    aPET = aPET[focus]
    dx = subsample.at[i, "DX"]
    subject = subsample.at[i, "id"]
    print_summary(aPET, f"{subject} ({dx})")

    fig = Figure(figsize=(10, 8), frameon = False)  # without this all the plots accumulate
    plt.title("PET intensities {} {} and reference, focal only".format(subject, dx))
    plt.hist((aPET, refVoxels, refVoxels[refVoxels>0.1]), bins=100, density=True, label=(subject, "reference", "reference > 0.1"))
    plt.legend()
    pp.savefig()
    plt.clf()

def go():
    global subsample
    ofile = ODIR / "investigate07.pdf"
    print_summary(refVoxels, "reference ADNI PETS (focal only)")
    pp = PdfPages(str(ofile))
    for i in range(0, subsample.shape[0]):
        do_one(i, pp)
    pp.close()

if __name__ == "__main__":
    go()
