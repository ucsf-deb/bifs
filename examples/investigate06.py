#! /usr/bin/env python3
# Investigate why there are so many zero voxels in the ADNI images
# File: investigate06.py
# Author: Ross Boylan
# Created: 2020-08-19

# First idea: orientation mismatch.  The mask and local images have one orientation,
# but at least some of the ADNI scans have a different one.

# INPUTS
# All are in ..\ExternalData\ycobigo\round3\
#
# TPM.nii   segmentation map (-> mask)
# ADNI\
#      classified-subset.csv list of images and dx (excludes AD)
#      the image files  themselves
#
# OUTPUTS
# investigation05.pdf plots
#
# HISTORY
# Code borrowed from empirical_scan.py and  presentation04-empirical.py.

ADNI = r"ExternalData\ycobigo\round3\ADNI"

import matplotlib
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.figure import Figure
from matplotlib.colors import NoNorm
from matplotlib import cm
import matplotlib.pyplot as plt

import nibabel
import numpy as np
from numpy.random import Generator, PCG64
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

import BIFS

def mypaths(seed=123, fraction=0.1):
    """Return a generator of the path names as <Path> to operate on
    randomly pick fraction of them to return"""
    rg = Generator(PCG64(seed))
    np.random.seed(seed)  # no longer best practice, but back compatible
    with open(ADNIDIR / "classified-subset.csv", "rt") as fin:
        n=0
        for line in fin:
            n += 1
            if n == 1:
                # skip initial header
                continue
            if rg.random() > fraction:
                continue
            fn = line[:(line.find(","))]
            yield ADNIDIR / fn


def get_mask():
    "Return a boolean mask for the images.  True means ignore this voxel"
    segs = nibabel.load(str(SEGFILE)).get_fdata()
    # True for voxels to set to 0
    mask = segs[..., 3:6].sum(3) > 0.80
    return mask

mask = get_mask()


def slice_fancy(image, ix=0, frac=0.5):
    "return rotated, masked 2D slice of 3D image."
    # rotate 90 degrees counter clockwise
    # ix and frac apply before rotation
    slice_index = np.int(np.round(frac*image.shape[ix]))
    delta = 0.5
    if ix == 0:
        im_slice = image[slice_index,:,:]
        im_mask = mask[slice_index,:,:]

    elif ix == 1:
        im_slice = image[:,slice_index,:]
        im_mask = mask[:,slice_index,:]
    elif ix == 2:
        im_slice = image[:,:,slice_index]
        im_mask= mask[:, :, slice_index]
    else:
        raise RuntimeError("Sorry slice index needs to be one of 0,1,2")
    vmin = im_slice.min()
    vmax = im_slice.max()
    adjust = np.full_like(im_slice, vmax)
    adjust[im_mask] = vmin
    im_slice = (im_slice + adjust)/2.0
    return np.rot90(im_slice)

def slice2(image, ix=0, frac=0.5):
    "return rotated, masked 2D slice of 3D image."
    # rotate 90 degrees counter clockwise
    # ix and frac apply before rotation
    slice_index = np.int(np.round(frac*image.shape[ix]))
    if ix == 0:
        im_slice = image[slice_index,:,:]
        im_mask = mask[slice_index,:,:]
    elif ix == 1:
        im_slice = image[:,slice_index,:]
        im_mask = mask[:,slice_index,:]
    elif ix == 2:
        im_slice = image[:,:,slice_index]
        im_mask= mask[:, :, slice_index]
    else:
        raise RuntimeError("Sorry slice index needs to be one of 0,1,2")
    m = np.zeros(im_mask.shape)
    m[im_mask] = 1.0
    return (np.rot90(im_slice), np.rot90(m))

def plot_prep(image, ix):
    "standard plot preparation for a 2D image"
    fig = Figure()
    plt.subplot(1, 3, ix+1)
    plt.rcParams["axes.grid"] = False # turn off grid lines for images
    plt.rcParams["xtick.color"] = (1,1,1,0)
    plt.rcParams["ytick.color"] = (1,1,1,0)
    plt.imshow(image, cmap = cm.Reds, alpha=0.6)


def plot_post(pp):
    "standard processing after all plotting of this page done"
    pp.savefig()
    # the text just the final plt.text persisted across figures without the next line
    plt.clf()

def go():
    b = BIFS.bifs()

    pp = PdfPages("investigate06.pdf")

    #info on slice to display
    ix = 2
    frac = 0.57

    for ix in range(3):
        init_image,  init_mask = slice2(np.zeros(mask.shape), ix = ix, frac = frac)
        plot_prep(init_mask, ix)
    plt.title("mask")
    plot_post(pp)

    n=0

    for fpath in mypaths():
        b.load_image_file(str(fpath))
        img = b.init_image()
        #img[mask] =0.0  leave it for now
        for ix in range(3):
            init_image, init_mask = slice2(img, ix = ix, frac = frac)
            plot_prep(init_image, ix)
            #plt.colorbar()
            plt.imshow(init_mask, cmap=cm.Blues, alpha=0.4)
        plt.title("ADNI for "+fpath.name)

        plot_post(pp)
        n += 1

    pp.close()

go()
