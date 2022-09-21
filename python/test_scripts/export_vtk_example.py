
from evtk import hl
import numpy as np
npoints = 10000
x = np.random.rand(npoints)
y = np.random.rand(npoints)
z = np.random.rand(npoints)
pressure = np.random.rand(npoints)
temp = np.random.rand(npoints)
hl.pointsToVTK("C:/Users/apoliti/Desktop/points.vu", x, y, z), data = {"temp" : temp, "pressure" : pressure})