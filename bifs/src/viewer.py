
import numpy as np
from mayavi import mlab

DATFILE=r"C:\Users\rdboylan\Documents\Kornak\bifs\ep1.npz"
x = np.load(DATFILE)
m = x["mean"]
sd = x["sd"]
x.close()
aSlice = mlab.volume_slice(m)
mlab.pipeline.volume(mlab.pipeline.scalar_field(m))
mlab.show_pipeline()
mlab.show()
x=None