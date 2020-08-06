# test problems mapping values between arrays
import matplotlib.pyplot as plt

import numpy as np
from numpy.random import Generator, PCG64

rg = Generator(PCG64(2398))

def quad(i, j):
    return i**2 + j**2

def step(i, j):
    return np.minimum(i, j)

dims = (2, 3)
a = np.fromfunction(step, dims)+rg.normal(0.0, 0.2, dims)
b = np.fromfunction(quad, dims)+rg.normal(0.0, 0.4, dims)
print(a)
print(b)
plt.subplot(221)
plt.title("a (step)")
plt.contourf(a)
  
plt.subplot(222)
plt.title("b (quadratic)")
plt.contourf(b)

a1 = a.reshape(-1)
a1 = np.sort(a1)

b1 = b.reshape(-1)
ix = b1.argsort()

c1 = np.empty_like(b1)
c1[ix]= a1
c = c1.reshape(b.shape)
plt.subplot(223)
plt.title("b with values from a")
plt.contourf(c)

plt.show()
pass
