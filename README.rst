BIFS Module Repository
========================

This project is an implementation of Bayesian imaging in Fourier Space
(BIFS). For details on the method see, e.g.:

* J. Kornak, K. Young, N. Schuff, A. Du, A. A. Maudsley, M. W. Weiner.
  K-Bayes Reconstruction for Perfusion MRI I: Concepts and Application. Journal of Digital Imaging. (2009) Feb 10
* J. Kornak, K. Young.
  K-Bayes Reconstruction for Perfusion MRI II: Modeling Motivation
  and Technical Development. Journal of Digital Imaging. (2009) Mar 10
* J. Kornak , K. Young, B.J. Soher, A.A. Maudsley.
  Bayesian k -space-time reconstruction of MR spectroscopic imaging for enhanced resolution. IEEE Trans Med Imaging. 2010 Jul;29(7):1333-50.

The basic idea is to perform something like a Bayesian image
restoration, that would ordinarily use an optimization technique
like Markov Chain Monte Carlo, but to try to do so in a more
efficient way by transforming to a space in which the degrees of
freedom can, at least approximately, be treated independently.
For a large class of images, transforming to a Fourier space
representation appears to accomplish this quite well.
And to emphasize a point that sometimes arises as a point of
confusion regarding this method, the goal is NOT to find a transform
that results in a parsimonious representation of the image, but one
that results in independence of the modes, so that the optimization
step can be performed independently on the modes, greatly increasing
efficiency.


Using the Package
-----------------

Once the dependencies are installed and the BIFS package has
been obtained from GitHub one can go to the directory containing
the class definition file, bifs.py, and the BIFS GUI, bifs_gui.py
and run the GUI via the command:

python bifs_gui.py

on MS-Windows, try
py bifs_gui.py
or
py -3 bifs_gui.py

This provides the easiest access to the BIFS package but provides
limited access to BIFS class variables.

To gain full access to the capabilities of BIFS one can: 1) run python,
2) load the BIFS class (import bifs), and 3) e.g. modify the
commands in one of the scripts (e.g. bifs_cl_2D.py) described below
to suit one's particular project.

Empirical Priors
~~~~~~~~~~~~~~~~

One can scan a set of images and use them to form a prior in Fourier space for
later analysis.  To do so, select "BIFS Operations" and then "Empirical Prior". 
Then select the top directory that holds all your images.  The code scans all subdirectories
and reads all files matching a regular expression, optionally skipping some of them.
It complains if the images don't all have the same dimensions.

You will probably want to do something different; you can edit the code in bifs_gui.py for
_scanImages.

The results are stored in a file ep1.npz in your current working directory.  This is defined
in EPFILE near the top of bifs_gui.py.  Saving the results means the scanning process
need only be done once; it it time and resource intensive.

To use the empirical prior select "Parameter Set" and then "Load Empirical Prior".  This reads
the file just mentioned and then uses it to form a prior.


Package Structure
-----------------

The package is separated into a class representing the calculation
engine and a main GUI class. The files contained in the package
are:

bifs/			- python package directory
	bifs.py           - the main class containing the BIFS functions

	bifs_gui.py       - the gui interface to bifs.py

	pset_dialogs.py   - set up functions for the BIFS Gui dialog boxes

	bifs.pyproj, bifs.sln  - MS Visual Studio 2019 project files

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

examples/		- sample code that uses BIFS package
	bifs_cl_1D.py     - an example script showing how to perform
						a BIFS analysis for 1D function in a
						python or ipython shell (can also just be
						run as a script)

	bifs_cl_2D.py     - similar to bifs_cl_1D.py but for a 2D image

	bifs_cl_3D.py     - similar to bifs_cl_1D.py but for a 3D data set

priors/	- definitions of different types of prior distributions
	FunctionalPrior.py  - priors based on smooth functions

tests/  -  test scripts.  For the pytest framework.
		   
Package Details
---------------

At the Python command line:

>>> import bifs
>>> help(bifs)

Or see the class documentation in BIFS/bifs.py.
