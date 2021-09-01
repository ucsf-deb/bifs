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


To Do
-----

These are for the package developers.  * as the first character in a headline means that is done.  This is different from * as an item mark, which
also happens.

  1. Settle overall package structure.
	  - invocation of test machinery
	  - capitalization of project name and class names, i.e., class names should start 
	    upper case.  Considerations:
		   * Python packge names should be lower case
		   * Class names should be capitalized
		   * github project names is 'bifs'
		   * having multiple things with the same name is confusing
		   * Old scheme was that the project name was 'bifs' in setup.py
		     but it was under the directory 'BIFS'.  It seems unlikely the package building machinery
		     would pick up on that.  import directives referred to the module as 'BIFS'.  And the main
		     class was 'bifs'. -> from BIFS import bifs
		   * New scheme switches it: from bifs import BIFS
	 - *Rename bifs.py -> bifscore.py to avoid import ambiguities.
  #. bifs_gui.py should be treated as a binary or installable script
      - Guido considers having executables in packages an anti-pattern
	  - hence relative imports, and scripts that might be run as either part of a package or separately
	    work awkwardly if at all.
	  - Concretely, this means bifs_gui should not use any relative imports.  
	    That seems to be the case already.
	  - How to test.
	  - How to package and install.
	      * should be in a different directory
		  * setup should list it as an entry point
		  * http://as.ynchrono.us/2007/12/filesystem-structure-of-python-project_21.html?showComment=1419123890864#c2535793882880782398 
		    argues entry point machinery is not great.  Note in particular the point about needing to install to run or test.
		  * I think there are entry points specific for GUIs
		  * How and where do I install it?
		  * Influential article says scripts in bin should be thin wrappers
		      + http://as.ynchrono.us/2007/12/filesystem-structure-of-python-project_21.html
			  + but no reason is ever given
  #. Other scripts?  Some of the Empirical Scanners might qualify.
      - bifs_util/EmpiricalScanner.py does not; it won't run by itself.
	  - bifs_util/util.py also does not
  #. Probably need to revive some of the scanner-related scripts
     In examples/
	   - empirical_scan.py
	   - scan.py
  #. __init__.py
      - *should be in every directory to import
	  - *was not true
	  - and yet the imports were working
	  - possibly related to relaxation of the __init_.py requirement in Python 3.3.
	  - but that relaxation was only for Namespace packages, which this is not
	  - What should go in __init__.py? opinions differ
	    https://www.reddit.com/r/Python/comments/1bbbwk/whats_your_opinion_on_what_to_include_in_init_py/
	     * Minimalist
		    + empty
			+ or just the version number
			+ maybe __all__, list of symbols to export
		* Maximalist
			+ key things is to regard it as a public interface
			+ this means using __all__ and maybe messing with nested modules to hoist some symbols up
			+ anything else that is appropriately global
  #. version number generally in 3 places, setup.py, PKG-INFO and __init__.py top level
      - must keep them in sync
	  - PKG-INFO appears to be generated automatically from other files.  It doesn't reflect the current version,
	  but I haven't gone through a packaging cycle with that version.
	  - there is a bump... package to do so, but it's not clear it worth the trouble
	       * since it requires naming the files to update
		   * though that info can go in a configuration file
  #. *convert all relative imports to absolute
      * generally recommended as more robust: The PEP 8 style guide recommends using absolute imports in general. 
  #. __main__.py, if present, executed automatically on package load
	Maybe only if started via `python -m`.
  #. Ensure existing tests work with current code.
	  - most seem to put tests directory top level
	  - http://as.ynchrono.us/2007/12/filesystem-structure-of-python-project_21.html?showComment=1419124050650#c5629056591479709945 argues against that
  #. Ensure existing examples work with current code.
  #. Consider which materials should be  distributed.
      - things to consider eliminating
	     * Visual Studio specific files and dirs
		 * some of the example code, i.e., Ross's
		 * git files
		 * other work files?
	  - decide how to eliminate them.  Could use a branch in git or a build script.
  #. Handling of multiple image loads in same session
  	- could be via GUI or via program, as in  presentation.py
  	- define desired functionality, which I think is and was to retain as many settings as possible.
  	- current behavior simply blows away previous settings
  #. Review and possibly remove bifs.copy_params
  #. Create  some general way to generate and use empirical priors.
  	- do not release with hard-coded path to store empirical prior.
	- literally I've done that since there is no path all! but name is hard-coded.
	- need a GUI to specify the empirical prior file, both construct and read
  #. Move much of the empirical prior logic out of bifs_gui.py.  Doneish?
  #. BIFS class should throw exceptions rather than print error messages
  #. bifs_gui should catch and display the exceptions.
  #. Continue delegation of model details to components of bifs.
  	- material currently stored in bifs should be obtained from
  		objects it holds, e.g., allowed function types, specific parameters.
  	- GUI will need to become more dynamic in response.
  #. Update to components should automatically invalidate the bifs object appropriately.
  	E.g., reset scale on prior causes bifs object to _invalidate_final().
  	The right way to do this is probably for the component, e.g., the prior, to hold a reference
  	to the parent, and send the parent _invalidate_final() as part of component's _mark_dirty().
  	For now we invalidate by default when handing out a prior object, and knowledgeable clients,
  	as in presentation.py that poke, plot, poke, replot must manually invalidate.
  #. @Karl Verify that the changes relating to isxcbanded are semantically correct; they are at least syntactically correct now.
  #. *If* more documentation outside the source code is desirable, e.g., this file used to have what was basically a copy
  	of the class comment for bifs, consider how to achieve that automatically.
  #. Review: NaNs in input file are now silently converted to 0 by bifs.load_image.  Is that desirable?
  	NaNs cause later processing to fail.
  #. Bug? After successfully running bifs_gui and selecting exit from its menu the console says
```
    File "J:\source\repos\bifs\bifs\bifs_gui.py", line 438, in close
    sys.exit(app.exec_())
	NameError: name 'app' is not defined
```
   In the debugger app was defined going into this step.

