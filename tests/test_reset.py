import bifs
import os
import numpy as np

def test_simple():
	b = bifs.bifs()
	print(os.getcwd())
	b.load_image_file("tests/images/lena512.bmp")
	b.BIFS_MAP()
	print(np.nonzero(b.prior_std2 - b.prior_std*b.prior_std))
	assert (b.prior_std2 == b.prior_std*b.prior_std).all()

def test_reset():
	b = bifs.bifs()
	b.load_image_file("tests/images/lena512.bmp")
	b.BIFS_MAP()

	b2 = bifs.bifs()
	b2.load_image_file("tests/images/lena512.bmp")
	b2.BIFS_MAP()
	b2.param_func_type = "Linear Decay"
	b2.BIFS_MAP()

	assert b.bas == b2.bas
	assert b.decay == b2.decay
	assert (b.prior_mean != b2.prior_mean).any()
	assert (b.bvec != b2.bvec).any()
	assert b.param_func_type != b2.param_func_type

def test_reset2():
	b = bifs.bifs()
	b.load_image_file("tests/images/lena512.bmp")
	b.BIFS_MAP()

	b2 = bifs.bifs()
	b2.param_func_type = "Linear Decay"
	b2.load_image_file("tests/images/lena512.bmp")
	b2.BIFS_MAP()


	assert b.bas == b2.bas
	assert b.decay == b2.decay
	assert (b.prior_mean != b2.prior_mean).any()
	assert (b.bvec != b2.bvec).any()
	assert b.param_func_type != b2.param_func_type

def test_inverse():
	b = bifs.bifs()
	b.param_func_type = "Banded Inverse Power Decay"
	b.load_image_file("tests/images/lena512.bmp")
	b.BIFS_MAP()

