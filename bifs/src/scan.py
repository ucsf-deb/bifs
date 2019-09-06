import os
import numpy as np
import nibabel as nib
from math import isnan
import sys


IMAGETOP = r"C:\Users\rdboylan\Documents\Kornak\ExternalData\ycobigo\ASL"
#IMAGETOP = r"C:\Users\rdboylan\Documents\Kornak\ADNI\PET1"

def do_scan(theDir):
    mismatch = set()  # holds keys that had a mismatch
    benchmarkHdr = None
    for root, dirs, files in os.walk(theDir):
        if files:
            print("Dir {}:".format(root))
            for f in files:
                print("\t{}".format(f))
                if not (f.endswith(".nii") or f.endswith(".nii.gz")):
                    continue
                img = nib.load(os.path.join(root, f))
                hdr = img.header
                print("\t  Shape: {}.  Transform: {}.  Zooms: {}".format(img.shape, img.affine.shape,
                                                                         hdr.get_zooms()))
                if not benchmarkHdr:
                    # first header encountered
                    benchmarkHdr = hdr
                    # could not delete the following key
                    # it actually doesn't appear in the objects attributes
                    #del benchmarkHdr.__dict__['db_name']  # differences expected and no concern
                    print("\t  {}".format(hdr))
                else:
                    for key in benchmarkHdr:
                        if key == 'db_name':
                            continue

                        if key.startswith("scl_"):
                            # values were array(nan, dtype=float32) and I had no luck testing for them
                            # in various ways
                            continue
                        v1 = benchmarkHdr[key]
                        v2 = hdr[key]
                        if (v1 != v2).any():
                            mismatch.add(key)
                            print("\t  {}:{}".format(key, hdr[key]))
                print("\n")

    if mismatch:
        print("The following keys were not uniform in the files scanned: {}.".format(mismatch))
    else:
        print("Header keys are consistent across all images scanned.")


if __name__== "__main__":
    for dir in sys.argv[1:]:
        do_scan(dir)
else:
    do_scan(IMAGETOP)
