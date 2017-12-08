#!/usr/bin/env python
import sys
sys.path += ["deps/grand-tour/lib", "lib/python"]

import json
import time
import numpy

from grand_tour import Topography
from default import SelectAntennas
from preselector import Preselector

# Get a handle for the topography
topo = Topography(latitude=42.928056, longitude=86.741667,
                  path="share/topography", stack_size=25)

# Generate the antenna positions
def get_antenna():
    deltar = 200.
    antxmin = -3.E+04
    antxmax = 3.E+04
    antymin = -3.E+04
    antymax = 3.E+04

    x = numpy.arange(antxmin, antxmax, deltar)
    y = numpy.arange(antymin, antymax, deltar)
    r = numpy.zeros((len(x) * len(y), 3))
    for i, xi in enumerate(x):
        i *= len(y)
        for j, yj in enumerate(y):
            r[i + j, :] = xi, yj, topo.ground_altitude(xi, yj) + 3.
    return r

print "# Generating antenna positions ..."
t0 = time.time()
antenna = get_antenna()
with open("antenna-positions.json", "wb+") as f:
    json.dump(antenna.tolist(), f)
print "  --> Done in {:.1f} s".format(time.time() - t0)

# Shower settings
energy = 1E+09
azimuth, elevation = 180., -5.
theta, phi = map(lambda x: numpy.deg2rad(90. - x), (azimuth, elevation))
s = numpy.sin(theta)
direction = (s * numpy.cos(phi), s * numpy.sin(phi), numpy.cos(theta))
position = (0., 0., topo.ground_altitude(0., 0.) + 1000.)

# Call the default selection
print "# Running default selection ..."
t0 = time.time()
a = SelectAntennas(energy, direction, position, topo)
print "  --> Done in {:.1f} s".format(time.time() - t0)

# Call the preselector
print "# Running preselector ..."
t0 = time.time()
selector = Preselector(topo, { "position" : antenna })
b = selector(energy, position, direction)
print "  --> Done in {:.1f} s".format(time.time() - t0)

for i in xrange(len(a)):
    print numpy.max(numpy.absolute(a[i] - b[i]))
