# test timing of apparently slow  operation

## from fourier.py
import numpy as np
def kdist_calc(N):
    """
    kdist(N)

    kdist outputs centered k-space indicex array

    Inpute:
    
    N - Length of desired k-space index array

    Outputs:

    centered, length N k-space index array

    """
    kval = np.zeros((1,N))
    kval[0,0:N//2] = range(1,N//2+1)
    kval[0,0:N//2] = kval[0,0:N//2] - 1 # try this so you get an origin 
    kval[0,N//2:N] = range(-N//2,0)
        
    return kval
        
    
def kdist3D(N1,N2,N3):
    """
    kdist3D() generates array with distance from
    center of Fourier space but shifted so origin is at index (0,0,0) 

    Inpute:
    
    N1 - Size of desired array in kx
    N2 - Size of desired array in ky
    N3 - Size of desired array in kz 

    Outputs:

    shifted, size (N1,N2,N3) k-space distance array
    """
    xvec = kdist_calc(N1)
    yvec = kdist_calc(N2)
    zvec = kdist_calc(N3)
    kdmat = np.zeros((N1,N2,N3))
    for i in range(N1):
        for j in range(N2):
            for k in range(N3):
                kdmat[i,j,k] = np.sqrt(xvec[0,i]**2 + yvec[0,j]**2 + zvec[0,k]**2)
    return kdmat

def alt1(*ns):
    dist = [kdist_calc(n)[0, ]**2 for n in ns]
    a, b, c = np.ix_(*dist)
    kdmat = np.sqrt(a+b+c)
    return kdmat
    

## end fourier.py excerpt
import timeit
print(timeit.timeit("x=alt1(100, 256, 256)", number=3, globals=globals()))
#r=alt1(50, 25, 25)
#print("Result shape = {}".format(r.shape))
