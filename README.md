# BIFS
**Bayesian Image Analysis in Fourier Space**

BIFS is a framework and set of computational tools for implementation of Bayesian Image analysis in Fourier Space. The BIFS approach enables enhanced feature contrast and noise reduction, analogous to, but much more efficiently than Bayesian methods applied directly in image space, such as Markov Chain Monte Carlo methods. This efficiency is provided by the fact that for a large class of images the Fourier modes can reasonably be approximated as being independent, thus reducing the optimization problem from one involving a large number of coupled variables to one involving a product of independent one dimensional optimizations. BIFS is useful for a large class of images and in particular for medical imaging problems, including clinical imaging studies, screening, diagnosis, and treatment evaluation.


Installation
============
Prerequisites: 
   - Python3 installed on your computer.

For the impatient, here are steps which should install the package and
launch the GUI:

```shell
#python may work instead of python3, depending on your environment
#py -3 should work on MS-Windows
python3 -m pip install --upgrade pip setuptools wheel
python3 -m pip install bifs
bifs_gui
```
bifs requires a number of other substantial software packages, including SciPy and PyQt5,
which the installation should take care of if they are not already present.

Virtual Environments
--------------------

See [this guide](https://packaging.python.org/tutorials/installing-packages) for much more
information about how to install python packages.  In particular, we endorse their
strong recommendation to *use virtual environments*, which are omitted from the steps above.

On Windows inside a virtual environment we have found `python` rather than `python3` is necessary.

We recommend against using the Windows-specific `py` command since it only respects the venv [sometimes, for example `py -3` will ignore it,](https://docs.python.org/3/using/windows.html#virtual-environments) under Python 3.9, and may have been even less reliable with earlier versions.

Qt Compatibility Problems
-------------------------

The graphical program `bifs_gui` depends on the `Qt` graphic framework.  This will be installed automatically if it is not present, 
but if it is already present the packages that `pip` installs may not be compatible with it.  We found this to be the case on Debian GNU/Linux 10 (buster).  If you get such a problem, you can fix it by specifying an explicit version for PyQt5 with `==`:

```bash
(testEnv) ross@barley:~/UCSF/Kornak$ python -m pip install PyQt5==5.11.2
Collecting PyQt5==5.11.2
  Downloading PyQt5-5.11.2-5.11.1-cp35.cp36.cp37.cp38-abi3-manylinux1_x86_64.whl (117.9 MB)
     |████████████████████████████████| 117.9 MB 84 kB/s 
Collecting PyQt5_sip<4.20,>=4.19.11
  Downloading PyQt5_sip-4.19.19-cp37-cp37m-manylinux1_x86_64.whl (67 kB)
     |████████████████████████████████| 67 kB 997 kB/s 
Installing collected packages: PyQt5-sip, PyQt5
  Attempting uninstall: PyQt5-sip
    Found existing installation: PyQt5-sip 12.9.0
    Uninstalling PyQt5-sip-12.9.0:
      Successfully uninstalled PyQt5-sip-12.9.0
  Attempting uninstall: PyQt5
    Found existing installation: PyQt5 5.15.4
    Uninstalling PyQt5-5.15.4:
      Successfully uninstalled PyQt5-5.15.4
Successfully installed PyQt5-5.11.2 PyQt5-sip-4.19.19
```

5.11.2 is the version of `Qt` that was already on the system, while 5.15.4 is the version that `pip` installed.

The current stable release of Debian is 11 (bullseye) and uses `Qt` 5.15.2, probably close enough to work with the default installation of `PyQt5`.

See https://github.com/ucsf-deb/bifs/issues/25 for more.

An alternative would be to manually install the required packages (see the `requirements.txt` or `setup.py` files) using your OS's package manager instead of `pip`, and then install `bifs` *without* using a  virtual environment.


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
See the more [detailed guide](https://github.com/ucsf-deb/bifs/blob/master/README.rst) and the extensive comments in the code for how to use the package.

Exploring the code
==================
https://github.com/ucsf-deb/bifs
