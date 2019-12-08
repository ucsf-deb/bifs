BIFS Module Repository
======================

This project is an implementation of Bayesian imaging in Fourier Space
(BIFS). For details on the method see, e.g.:

 - J. Kornak, K. Young, N. Schuff, A. Du, A. A. Maudsley, M. W. Weiner.
   K-Bayes Reconstruction for Perfusion MRI I: Concepts and Application. Journal of Digital Imaging. (2009) Feb 10
 - J. Kornak, K. Young.
   K-Bayes Reconstruction for Perfusion MRI II: Modeling Motivation
   and Technical Development. Journal of Digital Imaging. (2009) Mar 10
 - J. Kornak , K. Young, B.J. Soher, A.A. Maudsley.
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


Installation
------------

We are currently working on making the package fully PyPi compatible
so it will be easy to install, including installing all dependencies
In the meantime one can just get the package from GitHub (e.g. via
clone) and manually install the dependencies, e.g. via the pip command
(a good way to automatically install a number of the dependencies is to
start with an environment designed for scientific processing in
Python, e.g. Enthought's Canopy environment).

The currently required packages are (BIFS is known to work with
the following versions but will typically work with many earlier
versions as well):

::

 Python 3
 numpy 1.15.4-2
 scipy 1.2.1-1
 matplotlib 2.2.4-2
 imageio 2.4.1
 jsonpickle 1.2
 pyqt5 5.7.1-1
 (note: some version of PyQT 5 is required, i.e.
  the package will not work with PyQT 4)


Using the Package
-----------------

Once the dependencies are installed and the BIFS package has
been obtained from GitHub one can go to the directory containing
the class definition file, bifs.py, and the BIFS GUI, bifs_gui.py
(on typical install that will be ../bifs/src:) and run the GUI
via the command:

python bifs_gui.py

on MS-Windows, try
py bifs_gui.py
or
py -3 bifs_gui.py

This provides the easiest access to the BIFS package but provides
limited access to BIFS class variables.

To gain full access to the capabilities of BIFS one can:
 1) run python
 2) load the BIFS class (import bifs)
 3) e.g. modify the commands in one of the scripts (e.g. bifs_cl_2D.py) described below to suit one's particular project.

Empirical Priors
~~~~~~~~~~~~~~~~

One can scan a set of images and use them to form a prior in Fourier space for
later analysis.  To do so, select "BIFS Operations" and then "Empirical Prior". 
Then select the top directory that holds all your images.  The code scans all subdirectories
and reads all files named suvr_pons.nii.  It complains if they don't all have the same
dimensions.

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

bifs.py           - the main class containing the BIFS functions

bifs_gui.py       - the gui interface to bifs.py

bifs_cl_1D.py     - an example script showing how to perform
                    a BIFS analysis for 1D function in a
		    python or ipython shell (can also just be
		    run as a script)

bifs_cl_2D.py     - similar to bifs_cl_1D.py but for a 2D image

bifs_cl_3D.py     - similar to bifs_cl_1D.py but for a 3D data set

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
		   
Package Details
---------------

Much of the following information can be obtained at the Python
command line by importing the bifs package and typing help(bifs),
i.e. by issuing the following commands at the Python command line:

>>> import bifs
>>> help(bifs)

A number of class variables are available in the biffs class; many are
set to defaults by the class constructor and many are calculated and
set automatically when an initial image is loaded. To see how to set various
variables and perform an analysis using them, view the example scripts
provided with the package (and mentioned above).

Class variables available to constructor:

::

    init_image - initial loaded image
    k_image - initial k-space image
    mod_image - initial modulus image in k-space
    phase_image - initial phase image in k-space
    bifsk_image - final BIFS modulus image in k-space
    final_image - final reconstructed image in image space

    NOTE: All these images are currently stored in the BIFS object
    re. testing and experimentation - in future more
    parsimonious options may be provided re. production runs.

    image_file_loaded - whether an image is loaded (True,False)
    initial_image_file_name - file name of initial image

    imdim - int image dimension (1,2 or 3)
    imdim1 - int specifying size of 1st dimension of image
    imdim2 - int specifying size of 2nd dimension of image
    imdim3 - int specifying size of possible 3rd dimension of "image"

    kdist = distance function on the shifted k-space lattice

    view3Dslice - for 3D data this is a 2D array [a,b] where::
                   a = axis perpendicular to slice
                   b = fraction of maximum along that direction
                   for slice location
    
    prior - string specifying the prior distribution function to use
            current choices are:
              'Gaussian'

    prior_choices - list of current prior distribution
                    function choices (see above)
    prior_mean_init - prior mean before parameter space function
                      is set up (used for tests)
    prior_mean - the prior mean defined at each k-space point
                 by the k-space parameter function

    prior_std - the prior std defined at each k-space point
    prior_scale - the overall scale of the prior variance
    prior_scale_orig - prior scale at the origin - generally set huge
                       to allow the image data to determine overall scale

    likelihood - string specifying likelihood distribution function to use
                 current choices are:
                 'Gaussian'
                 'Rician'

    likelihood_choices - list of current choices (see above)
    likelihood_scale - the assumed (const) noise level in k-space

    bessel_approx_lims - limits for bessel approximation for rice
                         distribution - see paper referenced in code
    bessel_approx_array - array for bessel approximation for rice
                         distribution - see paper referenced in code

    rice_denom_cutoff - cutoff for the denominator of the closed form
                        of the posterior with a Gaussian prior and
                        Rician likelihood derived from bessel approximation
                        see paper referenced in code
    param_func_type - string specifying the k-space BIFS parameter
                      function to use
                      current choices are:
                      "Inverse Power Decay"
                      "Banded Inverse Power Decay"
                      "Linear Decay"
                      "Empirical"
		      
    param_func_choices - list of current choices (see above)
    decay - float decay exponent for the inverse power parameter function
    bvec - 2D float array specifying intercept and amplitude for parameter
           space functions 

    banded_cutoff - cutoff for banded, inverse power k-space parameter function
    basis - string specifying the basis to use - currently only choice
            is "Fourier"
    basis_choices - list of current choices (see above)

    bumps - dictionary containing set of "bump filters" to implement
            note: these "bump filters" are elements that are added
            to the parameter function to increase (or decrease if the
            amplitude is specified as negative) the sensitivity of the
            analysis to frequency ranges known in advance to be important
            (or missing) in the analyzed images. E.g. if there is a
            predominance of features of a give size, adding filters at
            wavelengths corresponding to that size could enhance the
            sensitivity of the analysis. The scipy.signal package
            provides a number of filters meant to applied in
            the time (image) domain to characterize properties in the
            Fourier domain. Providing these shapes for application in
            the Fourier domain for BIFS was straightforward and might
            be interesting to experiment with re. effective image
            feature enhancement.

    bump_types - set of choices for "bump" filter types to add to k-space
                 parameter function; uses scipy.signal window types
                 so consult that documentation for available types - 
                 currently only types that only require window type name
                 and size are used - current choices are: 
                 "boxcar"
                 "blackman"
                 "hann"
                 "bartlett"
                 "flattop"
                 "parzen"
                 "bohman"
                 "blackmanharris"
                 "nuttall"
                 "barthann"

    bump_default_type - the default window type used (currently "blackman")

    
