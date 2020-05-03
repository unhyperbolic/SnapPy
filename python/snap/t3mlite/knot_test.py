"""
An example of the way in which the functions in knot and 
knot_simplification are supposed to be used.
"""

import snappy.snap.t3mlite.knot
import snappy.snap.t3mlite.knot_simplification

def attempt_trefoil_simplification(num_blowups, num_attempts):
    """
    Attempting to simplify the triangulation of K3a1 and its 
    longitude embedded as a PL arc in barycentric coordinates.
    Hopefully, it will result in a smaller triangulation, but
    with a more complicated PL arc.

    The end goal is to make the triangulation so simple it can be 
    embedded in R^3, and we can just connect up the arcs to see
    the knot as a PL knot in R^3, but I haven't figured out how
    to implement this part yet.
    """
    MC = snappy.snap.t3mlite.knot.knot_triangulation("K3a1")
    print("Starting with {} tetrahedra".format(len(MC)))
    snappy.snap.t3mlite.knot_simplification.two_three_simplify(MC, num_blowups, num_attempts)
    print("End with {} tetrahedra".format(len(MC)))
    print("All tetrahedra with arcs:")
    print([T.arcs for T in MC.Tetrahedra])
