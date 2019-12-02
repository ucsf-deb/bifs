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
test_simple()