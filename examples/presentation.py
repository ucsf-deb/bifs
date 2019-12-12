# Prepare slides for a presentation
# First check if images are aligned
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

# hack the path to load bifs
sys.path.insert(0, Path.cwd() / "bifs")
import bifs

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
        
  
subjects = {}  # keys are ids, values are dictionaries of type, OneImage

def addFile(path : Path, type : str ):
    "populate the subjects dictionary. return the <OneImage>"
    b = bifs.bifs()
    b.load_image_file(str(path))
    b.image_file_loaded = True
    animage = OneImage(path, b, type)
    # I hope the {} is not shared...
    aSubject = subjects.setdefault(animage.id(), {})
    # assume no duplicates within type
    aSubject[type] = animage
    assert b.imdim == 3
    return animage


# Relative to top level project directory
DATATOP = r"..\ExternalData\ycobigo\registered"

## Gather images and do initial processing

MNI = None  # will hold a <OneImage>
for d in Path(DATATOP).iterdir():
    if d.is_dir() and not d.name.startswith("."):
        for f in d.glob("*.nii*"):
            if f.name.startswith("ICBM_ADNI"):
                addFile(f, "ADNI")
            elif f.name.startswith("mniwCBF_PVC_GM"):
                addFile(f, "MRI")
            elif f.name.startswith("mniwrsuvr_pons"):
                addFile(f, "PET")

MNI = addFile(Path(DATATOP) / "icbm152_GM.nii", "MNI")

def example01():
    """create example01.pdf with various slices of the original images.
    Intent is to check alignment and get a first cut at how PET and MRI
    images compare.
    """
    ## must be reworked to use subject to iterate
    pp = PdfPages('example01.pdf')
    for ix in range(3):
        for frac in [0.3, 0.5, 0.7]:
            # borrowed from bifs_gui.py show_initial_image
            for f, b in SUPER6:
                slice_index = np.int(np.round(frac*b.init_image.shape[ix]))
                if ix == 0:
                    init_im_slice = b.init_image[slice_index,:,:]
                elif ix == 1:
                    init_im_slice = b.init_image[:,slice_index,:]
                elif ix == 2:
                    init_im_slice = b.init_image[:,:,slice_index]
                else:
                    print("Sorry slice index needs to be one of 0,1,2")
                fig = Figure()
                plt.rcParams["axes.grid"] = False # turn off grid lines for images
                plt.rcParams["xtick.color"] = (1,1,1,0)
                plt.rcParams["ytick.color"] = (1,1,1,0)
                plt.imshow(init_im_slice, cmap = cm.Greys_r)
                myTitle = f.name
                if len(myTitle)>30:
                    myTitle = myTitle[:30]+"...."
                plt.title(myTitle)
                plt.text(0, init_im_slice.shape[0]+10, "Slice {}% along axis {}".format(frac*100, ix))
                pp.savefig()
                # the text just before savefig persisted across figures without the next line
                plt.clf()
    pp.close()

def slice(image, ix=0, frac=0.5):
    "return rotated 2D slice of 3D image"
    # rotate 90 degrees counter clockwise
    # ix and frac apply before rotation
    slice_index = np.int(np.round(frac*image.shape[ix]))
    if ix == 0:
        im_slice = image[slice_index,:,:]
    elif ix == 1:
        im_slice = image[:,slice_index,:]
    elif ix == 2:
        im_slice = image[:,:,slice_index]
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

def example02():
    """ graph various bifs transforms of the original images
    """
    pp = PdfPages('example02c.pdf')

    #info on slice to display
    ix = 0
    frac = 0.5

    for aSubject in subjects.values():
        if "MRI" in aSubject:
            break
    aMRI = aSubject["MRI"]
    aPET = aSubject["PET"]

    init_image = slice(aMRI.bifs.init_image, ix = ix, frac = frac)
    plot_prep(init_image)
    #plt.title("Initial MRI image")
    plot_post(pp)

    init_image = slice(aPET.bifs.init_image, ix = ix, frac = frac)
    plot_prep(init_image)
    #plt.title("Initial PET image")
    plot_post(pp)

    im_slice = slice(MNI.bifs.init_image, ix=ix, frac=frac)
    plot_prep(im_slice)
    #plt.title("MNI Reference Image")
    plot_post(pp)

    b = aMRI.bifs
    for pft in aMRI.bifs.param_func_choices:
        b.set_prior_func_type(pft)
        aprior = b.prior_object()
        scale = b.prior_object()._scale
        for sc in [scale/10, scale, scale*10]:
            aprior.setScale(sc)
            b.BIFS_MAP()
            im_slice = slice(b.final_image, ix=ix, frac=frac)
            plot_prep(im_slice)
            myTitle = aMRI.filename()
            if len(myTitle)>40:
                myTitle = myTitle[:40]+"...."
            #plt.title(myTitle)
            #plt.text(0, im_slice.shape[0]+10, "{} (scale {}). Slice {}% along axis {}".format(pft, sc, frac*100, ix))
            plot_post(pp)

    pp.close()

if __name__ == "__main__":

    example02()
