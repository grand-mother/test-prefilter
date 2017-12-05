#!/usr/bin/env python
import sys
sys.path += ["deps/grand-tour/lib", "lib/python"]

import time
import numpy

from grand_tour import Topography
from default import SelectAntennas
from trigger import select_antennas

# Get a handle for the topography
topo = Topography(latitude=42.928056, longitude=86.741667,
                  path="share/topography", stack_size=25)

energy = 1E+09
direction = (1., 0., 0.)
position = (0., 0., topo.ground_altitude(0., 0.) + 1.)

t0 = time.time()
a = SelectAntennas(energy, direction, position, topo)
print "{:.1f} s".format(time.time() - t0)

t0 = time.time()
b = select_antennas(energy, direction, position, topo)
print "{:.1f} s".format(time.time() - t0)

for i in xrange(len(a)):
    print numpy.max(numpy.absolute(a[i] - b[i]))
