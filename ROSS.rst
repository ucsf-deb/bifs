Environment on ross1
====================

Overall Environment
-------------------

Debian 8, aka oldstable on 64 bit Intel

Debian Packages (name and version)
----------------------------------

- python3 3.3.2
- libopenimageio1.4 1.4.14
  there is no Debian package for python3-imageio in this release
- python3-jsonpickle 0.8.0
- python3-matplotlib 1.4.2
- python3-numpy	   1.8.2
- python3-scipy 0.14.0

did not install

- python3-pyqt4 4.11.2

Installed via pip3 install --user

- imageio 2.4.1
  Note the PKG-INFO says it runs on Python 3.4+, not what I have.
  It does not mention any dependency on a binary libimageio.


Results
=======

Non-Graphical Environment
-------------------------

Here's what happened when I tried to run something::

| ross@ross-node1:~/Kornak/bifs$ python3 bifs/src/bifs_cl_1D.py 
| Running BIFS_MAP() on image
| Traceback (most recent call last):
|   File "bifs/src/bifs_cl_1D.py", line 55, in <module>
|     bu.plot_param_func(mybifs)
|   File "/home/ross/Kornak/bifs/bifs/src/bifs_util/util.py", line 185, in plot_param_func
|     plt.title("K Space Parameter Function") 
|   File "/usr/lib/python3/dist-packages/matplotlib/pyplot.py", line 1349, in title
|     l =  gca().set_title(s, *args, **kwargs)
|   File "/usr/lib/python3/dist-packages/matplotlib/pyplot.py", line 828, in gca
|     ax =  gcf().gca(**kwargs)
|   File "/usr/lib/python3/dist-packages/matplotlib/pyplot.py", line 462, in gcf
|     return figure()
|   File "/usr/lib/python3/dist-packages/matplotlib/pyplot.py", line 435, in figure
|     **kwargs)
|   File "/usr/lib/python3/dist-packages/matplotlib/backends/backend_tkagg.py", line 81, in new_figure_manager
|     return new_figure_manager_given_figure(num, figure)
|   File "/usr/lib/python3/dist-packages/matplotlib/backends/backend_tkagg.py", line 89, in new_figure_manager_given_figure
|     window = Tk.Tk()
|   File "/usr/lib/python3.4/tkinter/__init__.py", line 1854, in __init__
|     self.tk = _tkinter.create(screenName, baseName, className, interactive, wantobjects, useTk, sync, use)
| _tkinter.TclError: no display name and no $DISPLAY environment variable


Graphical Environment
---------------------

With X enabled, it pops up a bunch of graphs and shows::

| ross@ross-node1:~/Kornak/bifs$ python3 bifs/src/bifs_cl_1D.py
| Running BIFS_MAP() on image
| Prior mean:  397.6
| Prior scale:  0.1988
| K-space data value: 2.88936330499
| Likelihood scale: 0.14446816525
| Maximum of Rice Likelihood at:  2.89296825122
| Maximum of Gaussian Likelihood at:  2.88936328963
| Conjugate Gaussian value:  139.297476848
| Gaussian Likelihood MAP estimate: 2.61803396034 at k-space co-ordinates: (32,)
| Rician Likelihood MAP estimate: 2.61803396034 at k-space co-ordinates: (32,)

i.e., it at least appears to work.
