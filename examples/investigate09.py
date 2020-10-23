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

from copy import deepcopy
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
    print(" {} {}\n    SS = {:10.5e}\n    Quantiles {} ({})".format(desc, stats.describe(x), 
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
        print(f"******Whole image Results for {subj} ({dx})******")
        dists = []
    if labels is None:
        labels = []
    print_summary(vs, desc)
    dists.append(vs)
    labels.append(desc)
    return dists, labels


class Enforcer:
    """Assure Hermitian symmetry
    It includes the following instance variables once setup:
    dims <ndarray> shape of arrays it will enforce Hermitian symmetry on. 
    sourcei tuple of indices from which the values will be drawn
    targeti corresponding indices that will get an appropriate conjugate value from the source
    nSource number of constraints imposed

    Anticipated useage pattern:
    Initialize once with the dimensions of arrays it will get.
    Call set_phase or set_magnitude repeatedly as appropriate on
    arrays of the same size as the initial dimensions.

    Note that set_* methods modify the array that is their argument to 
    enforce the desired symmetry.
    """
    def __init__(self, dims):
        "dims is a collection giving shape of arrays to be considered"
        dims = np.array(dims, dtype=np.int16)
        assert np.all(dims>0)
        self.dims = dims
        nDim = len(dims)
        freshIndices = []  # holds all indices not seen before
        touched = np.full(dims, False) # True if value is seen or set
        it = np.nditer(touched, flags = ['multi_index'])
        nTouched = 0
        maxTouched = touched.size
        while nTouched < maxTouched and not it.finished:
            ix = it.multi_index
            if touched[ix]:
                continue
            # without a tuple this returns a generator
            six = tuple( dims[i]-ix[i]-1 for i in range(nDim)) # symmetric index
            touched[ix] = True
            nTouched += 1
            # setting this first will prevent us from entering "0" distance elements in the list to recode
            if not touched[six]:
                touched[six] = True
                nTouched += 1
                freshIndices.append(ix)
            it.iternext()
        self.nSource = len(freshIndices)
        self.sourcei = tuple(np.zeros(self.nSource, dtype=np.int16) for i in range(nDim))
        self.targeti = deepcopy(self.sourcei)
        for ir in range(self.nSource):
            fi = freshIndices[ir]
            for id in range(nDim):
                self.sourcei[id][ir] = fi[id]
                self.targeti[id][ir] = dims[id] - fi[id] -1

    def set_phase(self, a):
        """Force a phase array a to obey Hermitian symmetry.
        updates a in place and returns it."""
        a[self.targeti] = -a[self.sourcei]
        return a

    def set_magnitude(self, a):
        """Force a magnitude/modulus array a to obey Hermitian symmetry.
        updates a in place and returns it."""
        a[self.targeti] = a[self.sourcei]
        return a

PHASE_GENERATOR = None  # will hold an instance of the class below

class PG:
    "produce random phase matrices with description"

    def __init__(self, dims, seed):
        """dims is the shape of the arrays to process
        seed is the random seed
        """
        self.dims = dims
        self.enf = Enforcer(dims)
        self.uniform = stats.uniform(loc= -pi, scale= 2*pi)
        self.uniform.random_state = seed

    def gen(self):
        p = self.uniform.rvs(size=self.dims)
        yield (p, "uniform")
        self.enf.set_phase(p)
        yield (p, "uniform symm")



def do_one(i: int, pp):
    """plot subject i's PET distn vs ADNI
    pp is the backend
    phase_generators is a collection of functions that, when called with dimensions,
    return a random phase and a description"""
    global PHASE_GENERATOR
    global subsample, refVoxels
    b = BIFS.bifs()
    b.load_image_file(subsample.at[i, "fn.cbf"])
    dims = b.init_image().shape
    if PHASE_GENERATOR is None:
        PHASE_GENERATOR = PG(dims, 12890)
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
    for  m, des in PHASE_GENERATOR.gen():
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


def go():
    global subsample
    ofile = ODIR / "investigate09.pdf"
    pp = PdfPages(str(ofile))
    for i in range(0, 2): #subsample.shape[0]):
        do_one(i, pp)
    pp.close()

if __name__ == "__main__":
    go()

