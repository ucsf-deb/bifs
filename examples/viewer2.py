
from mayavi import mlab
import numpy as np
from time import sleep
x, y, z = np.mgrid[-5:5:64j, -5:5:64j, -5:5:64j]
vs = mlab.volume_slice(x, y, z, x*x*0.5 + y*y + z*z*2)

@mlab.animate(delay=200)
def anim():
    # The attribute to modify is vs.ipw.slice_index
    for i in range(64):
        #print(i)
        vs.ipw.slice_index = i
        yield
        #mlab.process_ui_events()
mlab.outline()
mlab.orientation_axes()
anim()
mlab.show()