#! /usr/bin/env python3
#  File: investigate08.py
#  Author: Ross Boylan
#  Created: 2020-10-02

# Explore the distribution of PET intensities in focal area
#
# INPUTS
#     vox3.npz   voxel intensities in focal area of ADNI PETS (large sample)
#     ep3.npz  The empirical prior
#     local PET images
#     TPM.nii segment map
#     subsamples.feather  file names and dx for input images
#         each row in subsample will produce a pdf
#     image files pointed to by subsample, MRI and PET

# OUTPUTS
#      investigate09.pdf  plots

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

import BIFS

EPNAME = "ep3.npz"
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
    print(" {} {}\n Quantiles {} ({})".format(desc, stats.describe(x), 
                                           stats.scoreatpercentile(x, ptiles).round(2),
                                           ptiles))

def referenceVoxels():
    """return  sorted array of voxels from images scanned"""
    return refVoxels

def adjustImage(img):
    """
    Adjust the distribution of pixels in the focus area of input image.
    Returns an input image of the same size with the intensity distribution adjusted 
    to match that in some scanned files.  Roughly the n'th percentile of image brightness in img
    will be assigned the n'th percentile of the reference images.

    All focal values < 5 are assigned the same value from the reference distribution.
    Rationale: All negative values are measurement errors.
    """
    ref = referenceVoxels()

    # directly updating the image does not work reliably
    img = np.copy(img, order='C')
    # argsort appears to break ties in a way that preserves the order it encounters elements
    cut = 5.0  # < cut in MRI space will be recoded to 0
    flat = img[focus]
    iFloor = flat < cut
    nFloor = np.sum(flat < cut)
    ix = flat.argsort()
    # The parentheses around ref.size/flat.size are essential.
    # Otherwise values wrap around
    subi = np.array(np.round(np.arange(flat.size)*(ref.size/flat.size)), dtype='int')
    flat[ix]  = ref[subi]
    vFloor = ref[subi][nFloor]
    flat[iFloor] = vFloor
    # reset min to 0
    flat = flat - vFloor
    img[focus] = flat
    return img

def do_dist(subj, dx, vs, desc, dists=None, labels=None):
    "handle all but final plot for a particular set of focal voxes vs"
    if dists is None:
        print(f"Results for {subj} ({dx}):")
        dists = []
    if labels is None:
        labels = []
    print_summary(vs, desc)
    dists.append(vs)
    labels.append(desc)
    return dists, labels


def do_one(i: int, pp):
    """plot subject i's PET distn vs ADNI
    pp is the backend"""
    global subsample, focus, refVoxels
    b = BIFS.bifs()
    b.load_image_file(subsample.at[i, "fn.cbf"])
    b._init_image[mask] = 0.0
    before = b.init_image()[focus]
    b.load_image(adjustImage(b._init_image))
    after = b.init_image()[focus]
    b.load_empirical(EPNAME)
    prior = b.prior_object()
    epmean = prior.mean()  # "origin" set to 0
    prior_as_image = np.real(b.bas.itxn(epmean*np.exp(1j*b.phase_image())))[focus]
    prior_no_phase = np.real(b.bas.itxn(epmean))[focus]



    aPET = nibabel.load(subsample.at[i, "fn.pet"]).get_fdata()
    aPET[np.isnan(aPET)] = 0.0  # mirroring the processing of bifs.load_image
    # convert to list of focal voxels
    aPET = aPET[focus]
    dx = subsample.at[i, "DX"]
    subject = subsample.at[i, "id"]

    rv = referenceVoxels()

    dists, labels = do_dist(subject, dx, aPET, "local PET")
    dists, labels = do_dist(subject, dx, rv[rv>0.1], "reference > 0.1", dists, labels )
    dists, labels = do_dist(subject, dx, prior_as_image, "prior w/image phase", dists, labels )
    dists, labels = do_dist(subject, dx, prior_no_phase, "prior, 0 phase", dists, labels )
    dists, labels = do_dist(subject, dx, after, "rescaled MRI", dists, labels)
    for scale in (1e-10, 4.0):
        prior.setScale(scale)
        b._invalidate_final()  # would be unnecessary in perfect world
        b.BIFS_MAP()  # unnecessary; call to final_image() triggers it anyway
        dists, labels = do_dist(subject, dx, after, f"posterior, scale={scale}", dists, labels)

    fig = Figure(figsize=(10, 8), frameon = False)  # without this all the plots accumulate
    i = 0
    while i<len(dists):
        j = min(len(dists), i+4)
        plt.subplot(2, 1, 1+i/4)
        plt.hist(dists[i:j], bins=100, density=True, label=labels[i:j])
        plt.legend()
        i = j
    plt.suptitle("Focal voxel intensities {} {}".format(subject, dx))

    pp.savefig()
    plt.clf()

def go():
    global subsample
    ofile = ODIR / "investigate08.pdf"
    pp = PdfPages(str(ofile))
    for i in range(0, 2): #subsample.shape[0]):
        do_one(i, pp)
    pp.close()

if __name__ == "__main__":
    go()

