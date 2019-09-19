
import numpy as np
from mayavi import mlab
from time import sleep

DATFILE=r"C:\Users\rdboylan\Documents\Kornak\bifs\bifs\src\ep1.npz"
x = np.load(DATFILE)
m = x["mean"]
sd = x["sd"]
x.close()
aSlice = mlab.volume_slice(m)
iMin = np.min(m.shape)
# we know it will be the x axis
iMin = m.shape[0]
m0 = m
m = np.zeros_like(m0)
m[m0>100] = 1.0
mlab.pipeline.volume(mlab.pipeline.scalar_field(m))
mlab.show_pipeline()
#mlab.show()

@mlab.animate(delay=200)
def anim():
    for i in range(iMin):
        aSlice.ipw.slice_index = i
        #mlab.process_ui_events()
        yield

mlab.outline()
mlab.orientation_axes()
anim()
mlab.show()
x=None