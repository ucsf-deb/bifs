# BIFS
**Bayesian Image Analysis in Fourier Space**

BIFS is a framework and set of computational tools for implementation of Bayesian Image analysis in Fourier Space. The BIFS approach enables enhanced feature contrast and noise reduction, analogous to, but much more efficiently than Bayesian methods applied directly in image space, such as Markov Chain Monte Carlo methods. This efficiency is provided by the fact that for a large class of images the Fourier modes can reasonably be approximated as being independent, thus reducing the optimization problem from one involving a large number of coupled variables to one involving a product of independent one dimensional optimizations. BIFS is useful for a large class of images and in particular for medical imaging problems, including clinical imaging studies, screening, diagnosis, and treatment evaluation.


Installation
============
Prerequisites: 
   - Python3 installed on your computer.
   - git in order to follow these instructions.

A good way to automatically install a number of the dependencies is to
start with an environment designed for scientific processing in
Python, e.g. Enthought's Canopy environment.

For the impatient, here are steps which should install the package and
launch the GUI:

```shell
git clone https://github.com/ucsf-deb/bifs.git
python3 -m pip install --upgrade pip setuptools wheel
python3 -m pip install -e bifs
python3 bifs/bifs/bifs_gui.py   # use \ on MS-Windows.  launches program
```
bifs requires a number of other substantial software packages, including SciPy and PyQt5,
which the installation should take care of if they are not already present.

bifs is *not* currently available from the standard python package index.

See [this guide](https://packaging.python.org/tutorials/installing-packages) for much more
information about how to install python packages.  In particular, we endorse their
strong recommendation to *use virtual environments*, which are omitted from the steps above.

Fallback Installation of Dependencies
-------------------------------------

We are still working on the packaging, as well as the instructions, and it is possible 
you will encounter difficulties.  If the necessary dependencies are not installed automatically,
try changing to the top level of the repository and executing
```shell
python3 -m pip install -r requirements.txt
```
If even that doesn't work, you can manually tell it to install the packages listed in the 
[requirements file](requirements.txt), e.g.,
```shell
python3 -m pip install numpy
```

Fallback Package Setup
----------------------
If the regular setup fails it means that `import bifs` will not be able to find the package.
Since the bifs code itself does such imports, it won't run without it.

A workaround is to manually add an entry to the `PYTHONPATH` environment variable for the 
root directory of the repository you cloned from GitHub.

If you clone the project following the instructions above and you are in `/tmp`, it will create a directory `/tmp/bifs`.  It is this
directory, not `/tmp/bifs/bifs` that you will find under it, that should go in `PYTHONPATH`.
See how to do it graphically or from the command line on [MS-Windows](https://stackoverflow.com/questions/3701646/how-to-add-to-the-pythonpath-in-windows-so-it-finds-my-modules-packages)
or for [Mac or Linux](https://bic-berkeley.github.io/psych-214-fall-2016/using_pythonpath.html).  One can also manipulate `sys.path` inside a Python
program.

Known Compatible Software Versions
----------------------------------
BIFS is known to work with the following versions but will typically work with many earlier
versions as well:
 * Python 3
 * numpy 1.17 (no earlier)
 * scipy 1.2.1-1
 * matplotlib 2.2.4-2
 * imageio 2.4.1
 * jsonpickle 1.2
 * pyqt5 5.7.1-1 (note: some version of PyQT 5 is required, i.e.
   the package will not work with PyQT 4)
 * nibabel is listed in the requirements file, but the software will
   run without it.  It is only needed if you load an image type for which
   it is required, e.g., a .nii file.

Using BIFS
==========
See the more [detailed guide](README.rst) and the extensive comments in the code for how to use the package.