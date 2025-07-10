from supremal_margulis_number import margulis_number

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

def print_margulis_manifold(name, M, use_surgery_description = True):
    M = geometric_triangulation(M, use_surgery_description=use_surgery_description)
    isosig = M.triangulation_isosig(
        ignore_orientation=False,
        ignore_cusp_ordering=True,
        ignore_curves=True)

    format_str = '"%s","%s",%s,"%s","%s",%f,%s'

    s = time.process_time()
    for bits_prec in [500, 1000, 2000, 5000]:
        try:
            epsilon, thin_part, collisions = margulis_number(
                M, bits_prec=bits_prec, verified=True, include_thin_part=True)
            break
        except Exception as e:
            exception = e
        finally:
            time_epsilon = time.process_time() - s
    else:
        print(
            format_str % (
                name,
                isosig,
                '', # epsilon
                '', # Thin part
                '', # Collisions
                time_epsilon,
                '"%r"' % exception))

        sys.stdout.flush()
        return

    if epsilon.is_NaN():
        raise Exception("NaN")

    if not epsilon > 0:
        raise Exception("Zero epsilon")

    if epsilon.diameter() > 1e-10:
        raise Exception("Big diameter")

    print(
        format_str % (
            name,
            isosig,
            repr(epsilon).replace('?',''),
            thin_part,
            collisions,
            time_epsilon,
            ''))

    sys.stdout.flush()

def print_margulis_cusped_census(stride, i):
    for M in OrientableCuspedCensus[i * stride : (i + 1) *stride]:
        print_margulis_manifold(str(M), M, use_surgery_description=True)

def print_margulis_closed_census(stride, i):
    for i in range(i * stride, min((i + 1) * stride, 11031)):
        name = 'Closed%d.tri' % (i + 1)
        M = Manifold('ClosedManifolds/%s' % name)
        print_margulis_manifold(name, M, use_surgery_description=False)

if __name__ == '__main__':
    import sys

    census = sys.argv[1]
    stride = int(sys.argv[2])
    i = int(sys.argv[3])
    
    if census == 'cusped':
        print_margulis_cusped_census(stride, i)
    elif census == 'closed':
        print_margulis_closed_census(stride, i)
    else:
        raise Exception("Unknonw census")
