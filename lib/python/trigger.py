import numpy

def select_antennas(energy, direction, position, topo):

    # Cone parameters
    gamma = (numpy.deg2rad(0.42) * numpy.log10(energy / 1E+08) +
             numpy.deg2rad(0.45))
    zcmin = 14000.
    zcmax = 27000. * numpy.log10(energy / 1E+08) + 22000.

    # Check if the shower crashes into a mountain early, before xmax
    deltar = 200
    s = numpy.arange(0., -zcmin - deltar, -deltar)
    xs, ys, zs = [position[i] + direction[i] * s for i in xrange(3)]
    zg = [topo.ground_altitude(xi, yi) for xi, yi in zip(xs, ys)]
    if (zs <= zg).any():
        return []

    def get_antennas():
        antxmin = -1.E+05
        antxmax = 1.E+05
        antymin = -1.E+05
        antymax = 1.E+05

        x = numpy.arange(antxmin, antxmax, deltar)
        y = numpy.arange(antymin, antymax, deltar)
        r = numpy.zeros((len(x) * len(y), 3))
        for i, xi in enumerate(x):
            i *= len(y)
            for j, yj in enumerate(y):
                r[i + j, :] = xi, yj, topo.ground_altitude(xi, yj)
        return r

    ra = get_antennas()

    # Select the antennas within the cone
    dr = ra - position
    zp = numpy.dot(dr, direction)
    rp2 = numpy.sum(dr**2, axis=1) - zp**2
    test_radius =  rp2 <= ((zp + zcmin) * numpy.tan(gamma))**2
    test_edgemin = zp <= -zcmin
    test_edgemax = zp >= -zcmax
    index= numpy.nonzero(test_radius & test_edgemin & test_edgemax)[0]
    if len(index) == 0:
        return index

    # Check for shadowing
    r0 = position - zcmin * direction

    def check_shadowing(i):
        u = ra[i, :] - r0
        d = numpy.linalg.norm(u)
        u /= d
        s = numpy.arange(0., d, deltar)
        xj, yj, zj = [r0[j] + u[j] * s for j in xrange(3)]
        zg = [topo.ground_altitude(x, y) for x, y in zip(xj, yj)]
        if (zj <= zg).any():
            return False
        return True

    return filter(check_shadowing, index)
