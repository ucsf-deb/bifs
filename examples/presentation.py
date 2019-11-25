# Prepare slides for a presentation
# First check if images are aligned
# Python 3.5+ required by pathlib
# MUST be run from top level of project directory

import sys
import numpy as np
from pathlib import Path

from PyQt5 import Qt, QtCore, QtGui, QtWidgets
from PyQt5.QtPrintSupport import QPrinter

sys.path.insert(0, Path.cwd() / "bifs")
import bifs

# Relative to top level project directory
DATATOP = r"..\ExternalData\ycobigo\registered"

# list of path-like objects of ADNI images, all PET
ADNI = []

# local PET images
PET = []

#local MRI images
MRI = []

app = QtWidgets.QApplication(sys.argv)
printer = QPrinter(QPrinter.HighResolution)
printer.setOutputFormat(QPrinter.PdfFormat)
# 1200 DPI resolution!
# setPaperSize always uses portrait orientation
printer.setPaperSize(QPrinter.Letter)
printer.setOrientation(QPrinter.Landscape)
printer.setOutputFileName("presentation01.pdf")
painter = QtGui.QPainter(printer)

# painter.scale() affects fonts and coordinates equally
# I want coordinates only.  So leave it alon.I hope this means I specify coordinates in inches
painter.setFont(QtGui.QFont("Times", pointSize = 20))
painter.drawText(1200, 1200, "Hi!")
painter.end()
sys.exit(0)

for d in Path(DATATOP).iterdir():
    if d.is_dir() and not d.name.startswith("."):
        for f in d.glob("*.nii*"):
            if f.name.startswith("ICBM_ADNI"):
                ADNI.append(f)
            elif f.name.startswith("mniwCBF_PVC_GM"):
                MRI.append(f)
            elif f.name.startswith("mniwrsuvr_pons"):
                PET.append(f)

print("ADNI: {}".format(ADNI))
print("PET: {}".format(PET))
print("MRI: {}".format(MRI))

ALL6 = ADNI+PET+MRI
for f in MRI:
    b = bifs.bifs()
    b.load_image_file(str(f))
    b.image_file_loaded = True
    assert b.imdim == 3
    # borrowed from bifs_gui.py show_initial_image
    slice_index = np.int(np.round(b.view3Dslice[1]*b.init_image.shape[b.view3Dslice[0]]))

    if b.view3Dslice[0] == 0:
        init_im_slice = b.init_image[slice_index,:,:]
    elif b.view3Dslice[0] == 1:
        init_im_slice = b.init_image[:,slice_index,:]
    elif b.view3Dslice[0] == 2:
        init_im_slice = b.init_image[:,:,slice_index]
    else:
        print("Sorry slice index needs to be one of 0,1,2")
    #self.ax1.imshow(init_im_slice, cmap = cm.Greys_r)