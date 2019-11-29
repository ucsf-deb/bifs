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

# Relative to top level project directory
DATATOP = r"..\ExternalData\ycobigo\registered"

## Gather images and do initial processing

# list of path-like objects of ADNI images, all PET
ADNI = []

# local PET images
PET = []

#local MRI images
MRI = []


for d in Path(DATATOP).iterdir():
    if d.is_dir() and not d.name.startswith("."):
        for f in d.glob("*.nii*"):
            if f.name.startswith("ICBM_ADNI"):
                ADNI.append(f)
            elif f.name.startswith("mniwCBF_PVC_GM"):
                MRI.append(f)
            elif f.name.startswith("mniwrsuvr_pons"):
                PET.append(f)

#print("ADNI: {}".format(ADNI))
#print("PET: {}".format(PET))
#print("MRI: {}".format(MRI))

ALL6 = ADNI+PET+MRI
SUPER6 = [] # holds pairs of f and bifs image
for f in ALL6:
    b = bifs.bifs()
    b.load_image_file(str(f))
    b.image_file_loaded = True
    assert b.imdim == 3
    SUPER6.append((f, b))

def example01():
    """create example01.pdf with various slices of the original images.
    Intent is to check alignment and get a first cut at how PET and MRI
    images compare.
    """
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

if __name__ == "__main__":
    example01()
