#! /usr/bin/env python3
#  File: presentation04-empirical.py
#  Author: Ross Boylan
#  Created: 2020-08-06
#  
# Prepare slides for evaluation by Howie
# uses analytic priors only on a selection of images
# One pdf per subject.
# pdf will include
#    MRI
#    PET
#    various transformed MRI
# All images are masked before processing to eliminate the skull, "big head" and outside voxels
#
# INPUTS
#     TPM.nii segment map
#     subsamples.feather  file names and dx for input images
#         each row in subsample will produce a pdf
#     image files pointed to by subsample, MRI and PET
#     ep3.npz  The empirical prior
#     vox3.npz sample  of voxels for rescaling
#
# OUTPUTS
#   presentation04/<id>-emp.pdf set of results
#   presentation04/guide-emp.csv  guide to files
# Python 3.5+ required by pathlib
#
# HISTORY
# 2020-08-06 copied from presentation04.py
#   borrow new path setup from empirical_scan.py
#
# 2020-09-02
#    reference voxels have ~1/3 that are near 0. eliminate them.
#
## I think the code below means you can run from anywhere at or under \Kornak\bifs

import concurrent.futures
import sys
import nibabel
import numpy as np
import pandas as pd
from pathlib import Path

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
ODIR = TOPDIR / "presentation04"
if not ODIR.exists():
    ODIR.mkdir()

# ids are floating point without the assist
subsample = pd.read_feather(SUBFILE).astype({"id": int})
# to enhance blinding, reorder
subsample = subsample.sort_values(by=["id"])
subsample.to_csv(str(ODIR / "guide-emp.csv"))
segs = nibabel.load(str(SEGFILE)).get_fdata()
seg_names = ("GM",  "WM", "CSP", "SK", "BH", "OUT")
# segnames[i] is the probabilities given by segs[,,,i]

# True for voxels to set to 0
mask = segs[..., 3:6].sum(3) > 0.80
focus = np.logical_not(mask)

# if we parallelize, easier to load this once
refVoxels = np.load(VOXNAME)["samp"]
# cut out the near 0, which we think is in error
refCut = 0.1
refVoxels = refVoxels[refVoxels>refCut]

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
    iFloor = np.logical_and(focus, img < 5.0)  # all these will be recoded to the floor value
    nFloor = np.sum(iFloor)
    flat = img[focus]
    ix = flat.argsort()
    # The parentheses around ref.size/flat.size are essential.
    # Otherwise values wrap around
    subi = np.array(np.round(np.arange(flat.size)*(ref.size/flat.size)), dtype='int')
    flat[ix]  = ref[subi]
    vFloor = ref[subi][nFloor]
    img[focus] = flat
    # reset all floor value
    img[iFloor] = vFloor
    return img


class OneImage:
    """data holder for a single image
    fname is a <pathlib.Path> for the file
    bifs  is a bifs object created for that file
    type is a string identifying the type of image this is
    """
    def __init__(self, fname, bifs, type):
        self.fname = fname
        self.bifs = bifs
        self.type = type

    def isLocal(self):
        "True if image is local"
        return self.type != "ADNI"

    def id(self):
        "subject id to match different images"
        # use the directory name in which the file resides
        return self.fname.parts[-2]

    def filename(self):
        return self.fname.name
   
    
def read_image(aPath : Path, type : str) -> OneImage:
    """return a <OneImage> created by reading from aPath,
    with type given by type.
    The image is masked"""
    b = BIFS.bifs()
    b.load_image_file(str(aPath))
    b._init_image[mask] = 0.0
    return OneImage(aPath, b, type)

def do_one(i : int):
    """Produce pdf file for i'th subject.
    This thin interface is designed to do parallel invocations.
    """
    global subsample
    aMRI = read_image(Path(subsample.at[i, "fn.cbf"]), "MRI")
    aPET = read_image(Path(subsample.at[i, "fn.pet"]), "PET")
    dx = subsample.at[i, "DX"]
    subject = subsample.at[i, "id"]
    ofile = ODIR / "{}-emp.pdf".format(subject )
    pp = PdfPages(str(ofile))

    #info on slice to display
    ix = 2
    frac = 0.57

    init_image = slice(aPET.bifs.init_image(), ix = ix, frac = frac)
    plot_prep(init_image)
    plt.title("FDG-PET")
    plot_post(pp)

    b = aMRI.bifs
    init_image = slice(b.init_image(), ix = ix, frac = frac)
    plot_prep(init_image)
    plt.title("Raw ASL perfusion CBF")
    plot_post(pp)

    before = b.init_image()[focus]
    b.load_image(adjustImage(b._init_image))
    after = b.init_image()[focus]
    init_image = slice(b.init_image(), ix = ix, frac = frac)
    img0 = np.copy(b.init_image())
    plot_prep(init_image)
    plt.title("MRI image with PET distribution values > {}".format(refCut))
    plot_post(pp)
    print("PET max={}, min={}".format(np.max(img0), np.min(img0)))

    plot_transform(before, after, pp)


    #im_slice = slice(MNI.bifs.init_image(), ix=ix, frac=frac)
    #plot_prep(im_slice)
    #plt.title("MNI Reference Image")
    #plot_post(pp)

    variant = 1
    b.load_empirical(EPNAME)
    prior = b.prior_object()
    for scale in (1e-10, 0.000001, 0.001, 0.01, 0.05, 0.4, 0.8, 1.0, 2.0, 4.0):
        prior.setScale(scale)
        b._invalidate_final()  # would be unnecessary in perfect world
        b.BIFS_MAP()  # unnecessary; call to final_image() triggers it anyway

        im_slice = slice(b.final_image(), ix=ix, frac=frac)
        deltatop = np.amax(img0-b.final_image())
        deltabot = np.amin(img0-b.final_image())
        print("With scale {} delta top = {}, bottom = {}.".format(scale, deltatop, deltabot))
        plot_prep(im_slice)
        myTitle = "BIFS reconstructed ASL perfusion CBF #{}".format(variant)
        plt.title(myTitle)
        plt.text(0, im_slice.shape[0], 
                 "Emprcl Prior (scale {:5.1e})\ndelta top = {:7.3G}, bottom = {:7.3G}.\nSlice {}% along axis {}".format(scale, deltatop, deltabot, round(frac*100), ix),
                 verticalalignment="top")
        plot_post(pp)
        variant += 1

    pp.close()  

def slice(image, ix=0, frac=0.5):
    "return rotated, masked 2D slice of 3D image."
    # rotate 90 degrees counter clockwise
    # ix and frac apply before rotation
    slice_index = np.int(np.round(frac*image.shape[ix]))
    if ix == 0:
        im_slice = image[slice_index,:,:]
        im_slice[mask[slice_index,:,:]] = 0.0
    elif ix == 1:
        im_slice = image[:,slice_index,:]
        im_slice[mask[:,slice_index,:]] = 0.0
    elif ix == 2:
        im_slice = image[:,:,slice_index]
        im_slice[mask[:, :, slice_index]] = 0.0
    else:
        raise RuntimeError("Sorry slice index needs to be one of 0,1,2")
    return np.rot90(im_slice)

def plot_prep(image):
    "standard plot preparation for a 2D image"
    fig = Figure()
    plt.rcParams["axes.grid"] = False # turn off grid lines for images
    plt.rcParams["xtick.color"] = (1,1,1,0)
    plt.rcParams["ytick.color"] = (1,1,1,0)
    plt.imshow(image, cmap = cm.Greys_r)

def plot_post(pp):
    "standard processing after all plotting of this page done"
    pp.savefig()
    # the text just the final plt.text persisted across figures without the next line
    plt.clf()

def plot_transform(before, after, pp):
    fig = Figure()  # RB not sure this is necessary
    plt.rcParams["axes.grid"] = True
    plt.rcParams["xtick.color"] = (0,0,0,0.5)
    plt.rcParams["ytick.color"] = (0,0,0,0.5)
    plt.minorticks_on()
    ix = np.argsort(before)
    ix = ix[np.linspace(0, ix.size-1, num=2000, endpoint = True, dtype = int)]
    plt.title("Transformation from MRI to PET in focal area")
    plt.plot(before[ix], after[ix], linestyle = "None", marker=".")
    #plt.hist(before) completely screws up previous graph
    pp.savefig()
    plt.clf()
    plt.title("Distribution of focal voxels")
    plt.minorticks_on()
    plt.hist(before, bins=500)
    pp.savefig()
    plt.clf()
    plt.title("Distribution of focal voxels in (0, 10]")
    plt.minorticks_on()
    # use & for and in next line didn't work
    plt.hist(before[ np.logical_and(before>0, before<=10)], bins=100)
    pp.savefig()
    plt.clf()


if __name__ == "__main__":
    with concurrent.futures.ProcessPoolExecutor(max_workers=7) as executor:
        executor.map(do_one, range(0, subsample.shape[0]))
    #do_one(0)