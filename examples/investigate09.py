#! /usr/bin/env python3
#  File: investigate09.py
#  Author: Ross Boylan
#  Created: 2020-10-14

# Explore the distribution of PET intensities in the entire image
# Like investigate08, except no masking.  However, ep3 and vox3 were based on masks.
#
# INPUTS
#     vox3.npz   voxel intensities in focal area of ADNI PETS (large sample)
#     ep3.npz  The empirical prior
#     local PET images
#     subsamples.feather  file names and dx for input images
#         each row in subsample will produce a pdf
#     image files pointed to by subsample, MRI and PET

# OUTPUTS
#      investigate09.pdf  plots

import math
from math import pi
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


# if we parallelize, easier to load this once
refVoxels = np.load(VOXNAME)["samp"]

def print_summary(x, desc):
    ptiles = [1, 2.5, 5, 10, 25, 50, 75, 90, 95, 97.5, 99]
    print(" {} {}\n SS = {:10.5e}\n Quantiles {} ({})".format(desc, stats.describe(x), 
                                           np.sum(np.square(x)),
                                           stats.scoreatpercentile(x, ptiles).round(2),
                                           ptiles))

def referenceVoxels():
    """return  sorted array of voxels from images scanned"""
    return refVoxels

def adjustImage(img):
    """
    Returns an input image of the same size with the intensity distribution adjusted 
    to match that in some scanned files.  Roughly the n'th percentile of image brightness in img
    will be assigned the n'th percentile of the reference images.

    """
    ref = referenceVoxels()

    # directly updating the image does not work reliably
    img = np.copy(img, order='C')
    # argsort appears to break ties in a way that preserves the order it encounters elements
    flat = np.ravel(img)
    ix = flat.argsort()
    # The parentheses around ref.size/flat.size are essential.
    # Otherwise values wrap around
    subi = np.array(np.round(np.arange(flat.size)*(ref.size/flat.size)), dtype='int')
    flat[ix]  = ref[subi]
    img = np.reshape(flat, img.shape)
    return img

def do_dist(subj, dx, vs, desc, dists=None, labels=None):
    "handle all but final plot for a particular set of voxels vs"
    if dists is None:
        print(f"Whole image Results for {subj} ({dx}):")
        dists = []
    if labels is None:
        labels = []
    print_summary(vs, desc)
    dists.append(vs)
    labels.append(desc)
    return dists, labels


def do_one(i: int, pp, phase_generators):
    """plot subject i's PET distn vs ADNI
    pp is the backend
    phase_generators is a collection of functions that, when called with dimensions,
    return a random phase and a description"""
    global subsample, refVoxels
    b = BIFS.bifs()
    b.load_image_file(subsample.at[i, "fn.cbf"])
    dims = b.init_image().shape
    before = b.init_image().ravel()
    b.load_image(adjustImage(b._init_image))
    after = b.init_image().ravel()
    b.load_empirical(EPNAME)
    prior = b.prior_object()
    epmean = prior.mean()  # "origin" set to 0
    prior_as_image = np.real(b.bas.itxn(epmean*np.exp(1j*b.phase_image()))).ravel()
    prior_no_phase = np.real(b.bas.itxn(epmean)).ravel()



    aPET = nibabel.load(subsample.at[i, "fn.pet"]).get_fdata()
    aPET[np.isnan(aPET)] = 0.0  # mirroring the processing of bifs.load_image
    # convert to list of focal voxels
    aPET = aPET.ravel()
    dx = subsample.at[i, "DX"]
    subject = subsample.at[i, "id"]

    rv = referenceVoxels()

    dists, labels = do_dist(subject, dx, aPET, "local PET")
    dists, labels = do_dist(subject, dx, rv, "reference", dists, labels )
    dists, labels = do_dist(subject, dx, prior_as_image, "prior w/image phase", dists, labels )
    dists, labels = do_dist(subject, dx, prior_no_phase, "prior, 0 phase", dists, labels )
    dists, labels = do_dist(subject, dx, after, "rescaled MRI", dists, labels)
    for scale in (1e-10, 4.0):
        prior.setScale(scale)
        b._invalidate_final()  # would be unnecessary in perfect world
        b.BIFS_MAP()  # unnecessary; call to final_image() triggers it anyway
        dists, labels = do_dist(subject, dx, after, f"posterior, scale={scale}", dists, labels)
    dists, labels = do_dist(subject, dx, b.phase_image().ravel(), "phase", dists, labels)
    for f in phase_generators:
        m, des = f(dims)
        dists, labels = do_dist(subject, dx, m.ravel(), des, dists, labels)
        img = np.real(b.bas.itxn(epmean*np.exp(1j*m))).ravel()
        dists, labels = do_dist(subject, dx, img.ravel(), "img: "+des, dists, labels)

    fig = Figure(figsize=(10, 8), frameon = False)  # without this all the plots accumulate
    i = 0
    nDists = len(dists)
    nPanels = math.ceil(nDists/4)
    while i<len(dists):
        j = min(len(dists), i+4)
        plt.subplot(nPanels, 1, 1+i/4)
        plt.hist(dists[i:j], bins=100, density=True, label=labels[i:j])
        plt.legend()
        i = j
    plt.suptitle("Focal voxel intensities {} {}".format(subject, dx))

    pp.savefig()
    plt.clf()

## phase distributions
## each function takes a random seed as an argument and returns
## a function the takes dimensions as an input and returns
## a random phase array and a description
def gen_uniform(seed):
    dist = stats.uniform(loc= -pi, scale= 2*pi)
    dist.random_state =seed
    def f(dims):
        return (dist.rvs(size=dims), "uniform")
    return f

def gen_uniform_sym(seed):
    dist = stats.uniform(loc= -pi, scale= 2*pi)
    dist.random_state = seed
    def f(dims):
        m = dist.rvs(size=dims)
        # enforce Hermitian symmetry
        dims = np.array(dims) # so I can do math with it

        return (m, "uniform symtrc")
    return f

def go():
    global subsample
    ofile = ODIR / "investigate09.pdf"
    pp = PdfPages(str(ofile))
    phase_dists = [gen_uniform(123)]
    for i in range(0, 2): #subsample.shape[0]):
        do_one(i, pp, phase_dists)
    pp.close()

if __name__ == "__main__":
    go()

