# -*- coding: utf-8 -*-

from setuptools import setup, find_packages

with open('README.rst') as f:
    readme = f.read()

with open('LICENSE') as f:
    license = f.read()

packages_ = find_packages('src')

setup(
    name='bifs',
    version='0.1.0a5',
    description='Implementation of Bayesian Imaging in Fourier Space (BIFS)',
    long_description=readme,
    author='John Kornak, Karl Young, Ross Boylan',
    author_email='john.kornak@ucsf.edu',
    url='https://github.com/bifs/bifs',
    license='BSD',
    package_dir={'': 'src'},
    package_data={'bifs': ['images/*']},
    include_package_data=True,
    scripts=['src/bifs/gui/bifs_gui.py'],
    entry_points = {
        'console_scripts': ['bifs_gui=bifs.gui.bifs_gui:main'],
    },
    packages=packages_,
    # packages=find_packages(where='src') # does not work for me.

    install_requires=['imageio',
                      'jsonpickle',
                      'matplotlib',
                      'nibabel',
                      'numpy',
                      'PyQt5',
                      'scipy'],
    # Option items:
    # extras_require={...},
    # data_files=[(...)],

    classifiers=[
        # How mature is this project? Common values are
        #   3 - Alpha; 4 - Beta; 5 - Production/Stable
        'Development Status :: 3 - Alpha',

        # Indicate who your project is intended for
        'Intended Audience :: Developers',
        'Topic :: Scientific/Engineering',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'Topic :: Scientific/Engineering :: Image Recognition',

        # Pick your license as you wish (should match "license" above)
        'License :: OSI Approved :: BSD License',

        # Specify the Python versions you support here.
        # I assume we only support Python 3 and recent version.
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.5',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
        ],
    keywords='Bayesian Fourier Image development',
    python_requires='>=3.5',
)

