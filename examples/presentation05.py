#! /usr/bin/env python3
#  File: presentation05.py
#  Author: Ross Boylan
#  Created: 2020-07-22
#  
# Show various enhancements to 2D breast image
# uses analytic priors only
# One pdf per subject.
# pdf will include
#    original image
#    various transformations
#
# INPUTS
#   sample-grayscale.png
#
# OUTPUTS
#   presentation05.pdf set of results
# Python 3.5+ required by pathlib
# MUST be run from top level of project directory

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

## Establish input and output paths

TOPDIR = Path.cwd()
while TOPDIR.name != "Kornak":
    TOPDIR = TOPDIR.parent

IMGFILE = TOPDIR / "ExternalData/BonnieJoeData/sample-grayscale.png"

OFILE = "presentation05.pdf"

def do_one(inf, outf):
    """Produce pdf outf from image inf
    """
    b = BIFS.bifs()
    b.load_image_file(str(inf))
    pp = PdfPages(str(outf))

    #info on slice to display
    ix = 2
    frac = 0.57

    init_image = b.init_image()
    plot_prep(init_image)
    plt.title("original image")
    plot_post(pp)

    variant = 1
    for pft in b.param_func_choices:
        if pft == "Banded Inverse Power Decay":
            # not interesting until we know know which frequencies to accentuate"
            continue
        b.set_prior_func_type(pft)
        aprior = b.prior_object()
        scale = b.prior_object()._scale
        scales = [scale/10, scale, scale*10]
        for sc in scales:
            aprior.setScale(sc)

            # ideally setScale would trigger this automatically
            b._invalidate_final()

            b.BIFS_MAP()  # unnecessary; call to final_image() triggers it anyway
            im_slice = b.final_image()
            plot_prep(im_slice)
            myTitle = "BIFS reconstructed #{}".format(variant)
            plt.title(myTitle)
            plt.text(0, im_slice.shape[0]+10, "{} (scale {:5.1e}).".format(pft, sc))
            plot_post(pp)
            variant += 1

    pp.close()  

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
    do_one(IMGFILE, OFILE)