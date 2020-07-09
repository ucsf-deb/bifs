#! /usr/bin/env python3
#  File: presentation04.py
#  Author: Ross Boylan
#  Created: 2020-07-07
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
#
# OUTPUTS
#   presentation04/<id>.pdf set of results
#   presentation04/guide.csv  guide to files
# Python 3.5+ required by pathlib
# MUST be run from top level of project directory

import concurrent.futures
import sys
import numpy as np
from pathlib import Path

import matplotlib
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.figure import Figure
from matplotlib.colors import NoNorm
from matplotlib import cm
import matplotlib.pyplot as plt

sys.path.insert(0, ".")
import BIFS

import nibabel
import pandas as pd

## Establish input and output paths

TOPDIR = Path.cwd()
while TOPDIR.name != "Kornak":
    TOPDIR = TOPDIR.parent

SEGFILE = TOPDIR / "ExternalData/ycobigo/round3/TPM.nii"
SUBFILE = TOPDIR / "RKornak" / "subsamples.feather"
ODIR = TOPDIR / "presentation04"
if not ODIR.exists():
    ODIR.mkdir()

# ids are floating point without the assist
subsample = pd.read_feather(SUBFILE).astype({"id": int})
# to enhance blinding, reorder
subsample = subsample.sort_values(by=["id"])
subsample.to_csv(str(ODIR / "guide.csv"))
segs = nibabel.load(str(SEGFILE)).get_fdata()
seg_names = ("GM",  "WM", "CSP", "SK", "BH", "OUT")
# segnames[i] is the probabilities given by segs[,,,i]

# True for voxels to set to 0
mask = segs[..., 3:6].sum(3) > 0.80


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
    ofile = ODIR / "{}.pdf".format(subject )
    pp = PdfPages(str(ofile))

    #info on slice to display
    ix = 2
    frac = 0.57

    init_image = slice(aMRI.bifs.init_image(), ix = ix, frac = frac)
    plot_prep(init_image)
    plt.title("Initial MRI image")
    plot_post(pp)

    init_image = slice(aPET.bifs.init_image(), ix = ix, frac = frac)
    plot_prep(init_image)
    plt.title("Initial PET image")
    plot_post(pp)

    #im_slice = slice(MNI.bifs.init_image(), ix=ix, frac=frac)
    #plot_prep(im_slice)
    #plt.title("MNI Reference Image")
    #plot_post(pp)

    b = aMRI.bifs
    for pft in aMRI.bifs.param_func_choices:
        b.set_prior_func_type(pft)
        aprior = b.prior_object()
        scale = b.prior_object()._scale
        for sc in [scale/10, scale, scale*10]:
            aprior.setScale(sc)

            # ideally setScale would trigger this automatically
            b._invalidate_final()

            b.BIFS_MAP()  # unnecessary; call to final_image() triggers it anyway
            im_slice = slice(b.final_image(), ix=ix, frac=frac)
            plot_prep(im_slice)
            myTitle = aMRI.filename()
            if len(myTitle)>40:
                myTitle = myTitle[:40]+"...."
            plt.title(myTitle)
            plt.text(0, im_slice.shape[0]+10, "{} (scale {}). Slice {}% along axis {}".format(pft, sc, round(frac*100), ix))
            plot_post(pp)

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


if __name__ == "__main__":
    with concurrent.futures.ProcessPoolExecutor(max_workers=3) as executor:
        executor.map(do_one, range(0, subsample.shape[0]))
#    do_one(0)