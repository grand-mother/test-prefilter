import numpy

def select_antenna(antenna, energy, direction, position, topo):

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

    # Select the antenna(s) within the cone
    dr = antenna - position
    zp = numpy.dot(dr, direction)
    rp2 = numpy.sum(dr**2, axis=1) - zp**2
    test_radius = rp2 <= ((zp + zcmin) * numpy.tan(gamma))**2
    test_edgemin = zp <= -zcmin
    test_edgemax = zp >= -zcmax
    index = numpy.nonzero(test_radius & test_edgemin & test_edgemax)[0]
    if len(index) == 0:
        return index

    # Check for shadowing
    r0 = position - zcmin * numpy.array(direction)

    def check_shadowing(i):
        u = antenna[i, :] - r0
        d = numpy.linalg.norm(u)
        u /= d
        s = numpy.arange(0., d, deltar)
        xj, yj, zj = [r0[j] + u[j] * s for j in xrange(3)]
        zg = [topo.ground_altitude(x, y) for x, y in zip(xj, yj)]
        if (zj <= zg).any():
            return False
        return True

    return filter(check_shadowing, index)
