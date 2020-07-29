import itertools
import os
import re

# scans files to construct an empirical prior
from BIFS import bifs

# numpy >= 1.17
from numpy.random import Generator, PCG64
import numpy as np

class AbstractEmpiricalScanner:
    """ This class consumes a list of images to generate an empirical distribution by voxel.
    All images must have the same dimensions, and they should be aligned with each other
    for the results to be meaningful.

    Concrete classes provide particular ways to get images.  They then pass the images to _statsAccumulate and,
    optionally, _voxAccumulate (possibly different images for each) and call
    _post when done.  At that point, and only that point, are results available in self.mns and self.sds
    and, if requested, self.vox.

    mns and sds are arrays in the same shape as the images with the mean and sd at each pixel.
    vox is a 1-D array of voxel values, sorted in ascending order.  If sampleFraction<1 it will be a subset
    of all values seen.

    image_mask optionally indicates which areas of the image to ignore.
        It must be a boolean array with the same shape as image files.
        All voxels selected by image_mask are set to zero before doing BIFS processing.
        The mask applies to the original image NOT to the fourier space version, which will
        generally have non-0 values in the image_mask region.
        It is the subclass responsibility to implement these semantics.
        Note the "mask" here is not a mask in the numpy sense of a masked array, which
        concerns missing values.

    voxel sampling should be done only from the non-ignored regions, but the number sampled will be based on
        the total voxel count before masking.
    """
    def __init__(self, sampleFraction=0, seed=85792359, image_mask=None):
        """Setup for scan of images
        if sampleFraction is >0 (and it should be <=1) then that fraction of the image voxels will be retained.
        In that case, seed is used to set the random number generator.
        """
        self.nImages = 0
        self.sampleFraction = sampleFraction
        if sampleFraction>0:
            self._vox = []
            self._rg = Generator(PCG64(seed))
        self.masking = (image_mask is not None)
        if self.masking:
            self.image_mask = image_mask
            self.image_keep = np.logical_not(image_mask)


    def _statsAccumulate(self, m):
        """
        Accumulate running statistics on the values of m in different images
        self.nImages gives the current image number; it starts at 1.
        m is ordinarily the modulus, and must conform to numpy array protocols

        Updates self._mns and self._ss currently.
        """
        self.nImages += 1
        if self.nImages == 1:
            self._mns = m
            # ss will turn into a matrix later
            self._ss = 0.0
        else:
            lastdelta = m-self._mns
            self._mns += (lastdelta)/self.nImages
            # element by element multiplication in next line
            self._ss += lastdelta*(m-self._mns)

    def _voxAccumulate(self, m):
        """accumulate voxel values.
        In the most likely case, the voxels are from image space while the empirical prior
        is from k-space.  So we provide seperate functions for the 2 values.
        Calling this is pointless unless sampleFraction>0.
        """
        # always base number sampled on the complete image size
        nSamp = int(m.size*self.sampleFraction)

        if self.masking:
            self._vox.append(self._rg.choice(m[self.image_keep], nSamp))
        else:
            # m.ravel is not an acceptable first argument to choice
            # actually, it should have been np.ravel(m)
            # m.flatten and the mask selection above both create copies, unfortunately
            self._vox.append(self._rg.choice(m.flatten(), nSamp))


    def _statsPost(self):
        """
        Finalize computation of voxel by voxel statistics for all images.
        Call after all images have been seen.

        Results returned as arrays self.mns and self.sds.
        """
        self.mns = self._mns
        # element by element square root
        self.sds = np.sqrt(self._ss/(self.nImages-1))
        del self._mns
        del self._ss

    def _voxPost(self):
        """
        Finalize accumulated voxels.
        """
        if self.sampleFraction>0:
            self.vox = np.concatenate(self._vox)
            self.vox.sort()
            del self._vox

    def _post(self):
        "wrap up all processing"
        self._statsPost()
        self._voxPost()


class EmpiricalScanner(AbstractEmpiricalScanner):
    """Scan selected images on disk, ensuring they are alll compatible.

    topDir  path like object indicating where in the file system the scan should start
        all subdirectories will be scanned recursively unless they are excluded.
    matchFile  <String> regular expression for the file name of image files we want.
        Matching is on the file name only, not its full path.
    exclude     <String> optional regular expression.  Any directory matching this pattern is excluded.
        Any file that satisfies matchFile is excluded if it also matches exclude.

    See AbstractEmpiricalScanner for sampleFraction, seed and image_mask.

    The files are read in and converted to k-space.  We compute  the mean and sd of the k-space images,
    and optionally accumulate voxels from the original image.

    We also check that the headers are consistent.  This works for .nii files, and may or may not for others.
    """
    def __init__(self, sampleFraction=0, seed=85792359, topDir=".", matchFile="", exclude=None, image_mask=None):
        super().__init__(sampleFraction, seed)
        self._topDir = topDir
        self._matchRE = re.compile(matchFile, re.I)
        if exclude:
            self._excludeRE = re.compile(exclude, re.I)
        else:
            self._excludeRE = None
        self._benchmarkHdr = None  # checks for consistent headers
        self._mismatch = set()  # holds keys that had a mismatch
        self._bifs = bifs()
        self.go()

    def go(self):
        """Actually perform the scan.
        Note this is triggered by object initialization.
        Repeated calls may not work.
        """
        for root, dirs, files in os.walk(self._topDir):
            if self._excludeRE:
                # avoid directories with our target case for whom we are trying to predict
                iKill = [ i for i, d in zip(itertools.count(), dirs) if self._excludeRE.search(d)]
                if iKill:
                    nKill = 0
                    for i in iKill:
                        i -= nKill
                        print("Skipping {}".format(dirs[i]))
                        del dirs[i-nKill]
                        nKill += 1
            # look for files to import
            if files:
                for f in files:

                    if not self._matchRE.search(f):
                        continue
                    if self._excludeRE:
                        if self._excludeRE.search(f):
                            print("Skipping {}".format(f))
                            continue
                    self._bifs.load_image_file(os.path.join(root, f))
                    if self.masking:
                        # dirty trick.  But doesn't invalidate anything else in _bifs.
                        self._bifs._init_image[self.image_mask] = 0.0
                    self._statsAccumulate(self._bifs.mod_image())
                    if self.sampleFraction>0:
                        self._voxAccumulate(self._bifs.init_image())
                    hdr = self._bifs.read_imfile.header
                    if not self._benchmarkHdr:
                        # first header encountered
                        self._benchmarkHdr = hdr
                        # could not delete the following key
                        # it actually doesn't appear in the objects attributes
                        #del benchmarkHdr.__dict__['db_name']  # differences expected and no concern
                    else:
                        for key in self._benchmarkHdr:
                            if key == 'db_name':
                                continue

                            if key.startswith("scl_"):
                                # values were array(nan, dtype=float32) and I had no luck testing for them
                                # in various ways
                                continue
                            v1 = self._benchmarkHdr[key]
                            v2 = hdr[key]
                            if (v1 != v2).any():
                                mismatch.add(key)
        self._post()
