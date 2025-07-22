"""


"""

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

def geometric_triangulation(M, use_surgery_description=True):
    if use_surgery_description:
        M = surgery_description(M)
    for i in range(1000):
        if M.solution_type() == 'all tetrahedra positively oriented':
            break
        M.randomize()
    return M

def print_systole_manifold(name, M, use_surgery_description=True, file=None):
    M = geometric_triangulation(M, use_surgery_description=use_surgery_description)
    isosig = M.triangulation_isosig(
        ignore_orientation=False,
        ignore_cusp_ordering=True,
        ignore_curves=True)

    format_str = '"%s","%s",%s,%s'

    for bits_prec in [500,1000,5000]:
        try:
            L = M.length_spectrum_alt_gen(verified=True, bits_prec=bits_prec)
            systole = next(L)
            break
        except Exception as e:
            exception = e
    else:
        print(
            format_str % (
                name,
                isosig,
                '', # systole
                '"%r"' % exception),
            file=file)

        file.flush()
        return

    real_systole = systole.length.real()

    print(
        format_str % (
            name,
            isosig,
            repr(real_systole).replace('?',''),
            ''),
        file=file)

    file.flush()

if __name__ == '__main__':
    import sys

    try:
        stride = int(sys.argv[1])
        i = int(sys.argv[2])
    except:
        print("Usage: sage print_systole_cusped_census.py [STRIDE] [i]")
        print()
        print("       Writes STRIDE manifolds to systole_cusped_[i].csv starting ")
        print("       with manifold STRIDE * i.")
        print("       Writes header to systole_cusped_0000.csv.")
        print()
        print("       Can be run with:")
        print("           parallel -j 4 sage print_systole_cusped_census.py 1000 -- `seq 0 61`")
        print()
        print("       parallel is a utility from moreutils package on Ubuntu Linux ")
        print("       or brew (Mac OS).")
        print("       -j is number of parallel processes")
        print()
        print("       After running, to get result that can be opened as a ")
        print("       spread sheet (for example with Numbers on Mac OS):")
        print("           cat systole_cusped_*.csv >systole_cusped.csv")

        sys.exit(1)

    with open('systole_cusped_%04d.csv' % i, 'w') as f:
        if i == 0:
            print('"Manifold","Triangulation","Systole","Exception"', file=f)

        for M in OrientableCuspedCensus[i * stride : (i + 1) *stride]:
            print_systole_manifold(str(M), M, use_surgery_description=True, file=f)
