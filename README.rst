BIFS Module Repository
========================

This project is an implementation of Bayesian imaging in Fourier Space
(BIFS). For details see (ref.)

The basic idea is to perform something like a Bayesian image
restoration, that would ordinarily use an optimization technique
like Markov Chain Monte Carlo, but to try to do so in a more
efficient way by transforming to a space in which the degrees of
freedom can, at least approximately, be treated independently.
For a large class of images, transforming to a Fourier space
representation apears to accomplish this quite well.
And to emphasize a point that sometimes arises as a point of
confusion regarding this method, the goal is NOT to find a transform
that results in a parsimonious representation of the image but one
that reults in independence of the modes, so that the optimization
step can be performed independently on the modes, greatly increasing
effiency.

The package is seperated into a class representing the calculation
engine and a main GUI class. The files containted in the package
are:

bifs.py           - the main class containing the BIFS functions

bifs_gui.py       - the gui interface to bifs.py

bifs_cl_1D.py     - an example script showing how to perform
                    a BIFS analysis for 1D function in a
		    python or ipython shell (can also just be
		    run as a script)

bifs_cl_2D.py     - simlar to bifs_cl_1D.py but for a 2D image

bifs_cl_3D.py     - simlar to bifs_cl_1D.py but for a 3D data set

pset_dialogs.py   - set up functions for the BIFS Gui dialog boxes

bifs_util/util.py - some utility functions (used in the above scripts)

bases/fourier.py  - in principle a BIFS style analysis could be
                    performed in any basis (e.g. a wavelet basis) and
		    it is assumed that for various data sets, bases
		    may be found for which a better approximation to
		    independent components might be found. The
		    package is structured to allow for easy
		    implementation of other bases but the only basis
		    functions currently implemented are those for
		    Fourier bases. fourier.py contains the Fourier
		    specific functions.
		   
