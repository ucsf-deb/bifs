# -*- coding: utf-8 -*-

from setuptools import setup, find_packages

with open('README.rst') as f:
    readme = f.read()

with open('LICENSE') as f:
    license = f.read()

setup(
    name='bifs',
    version='0.2.0',
    description='Implementation of Bayesian Imaging in Fourier Space (BIFS)',
    long_description=readme,
    author='John Kornak, Karl Young, Ross Boylan',
    author_email='john.kornak@ucsf.edu,kyoung21b2000@gmail.com',
    url='https://github.com/bifs',
    license=license,
    packages=['bifs'],
    install_requires=['imageio',
                      'jsonpickle',
                      'matplotlib',
                      'nibabel',
                      'numpy',
                      'PyQt5',
                      'scipy']
)

