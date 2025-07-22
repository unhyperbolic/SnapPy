from snappy import OrientableCuspedCensus, Manifold

import time
import sys

def surgery_description(M):
    G = M.fundamental_group(False)
    real_len, g = min((real_len, g) for g in G.generators() if (real_len := G.complex_length(g).real()) > 1e-6)
    if real_len > 0.4:
        return M

    try:
        N = M.drill_word(g, verified=True, bits_prec=1000)
    except:
        return M
    N.dehn_fill((1,0),-1)
    return N

def geometric_triangulation(M, use_surgery_description = True):
    if use_surgery_description:
        M = surgery_description(M)
    for i in range(1000):
        if M.solution_type() == 'all tetrahedra positively oriented':
            break
        M.randomize()
    return M

def print_systole_manifold(name, M, use_surgery_description = True):
    M_surgery = geometric_triangulation(M, use_surgery_description=use_surgery_description)
    isosig = M_surgery.triangulation_isosig(
        ignore_orientation=False,
        ignore_cusp_ordering=True,
        ignore_curves=True)

    if use_surgery_description:
        M = geometric_triangulation(M, use_surgery_description=False)
    else:
        M = M_surgery

    format_str = '"%s","%s",%s,%f,%f,%s'

    s = time.process_time()
    for bits_prec in [500]:
        try:
            L = M_surgery.length_spectrum_alt_gen(verified=True, bits_prec=bits_prec)
            systole = next(L)
            break
        except Exception as e:
            exception = e
        finally:
            time_new = time.process_time() - s
    else:
        print(
            format_str % (
                name,
                isosig,
                '', # systole
                time_new,
                0.0,
                '"%r"' % exception))

        sys.stdout.flush()
        return

    real_systole = systole.length.real()

    Mhigh = M.high_precision()

    if True:
        s = time.process_time()

        L = Mhigh.length_spectrum(real_systole.center() + 0.01)
        if not abs(L[0]['length'].real() - real_systole.center()) < 1e-6:
            raise Exception("Mismatch")

        time_old = time.process_time() - s

    if real_systole.is_NaN():
        raise Exception("NaN")

    if not real_systole > 0:
        raise Exception("Zero epsilon")

    if real_systole.diameter() > 1e-10:
        raise Exception("Big diameter")

    print(
        format_str % (
            name,
            isosig,
            repr(real_systole).replace('?',''),
            time_new,
            time_old,
            ''))

    sys.stdout.flush()

def print_systole_cusped_census(stride, i):
    for M in OrientableCuspedCensus[i * stride : (i + 1) *stride]:
        print_systole_manifold(str(M), M, use_surgery_description=True)

def print_systole_closed_census(stride, i):
    for i in range(i * stride, min((i + 1) * stride, 11031)):
        name = 'Closed%d.tri' % (i + 1)
        M = Manifold('ClosedManifolds/%s' % name)
        print_systole_manifold(name, M, use_surgery_description=False)

if __name__ == '__main__':
    import sys

    census = sys.argv[1]
    stride = int(sys.argv[2])
    i = int(sys.argv[3])

    if census == 'cusped':
        print_systole_cusped_census(stride, i)
    elif census == 'closed':
        print_systole_closed_census(stride, i)
    else:
        raise Exception("Unknonw census")
