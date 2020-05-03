from .knot import * 
"""
To do still: after simplifying a knot (embedded as a PL curve
in barycentric coordinates) to two tetrahdra, we need to still
embed this triangulation in R^3 to get a PL knot. Then, we can
recover the knot from the triangulation of the exterior.

I haven't been able to figure this out. There are three different
 two tetrahedron triangulations of S3 that result from the
 sequence of 2-3, 4-1, 3-2 and 1-4 moves that we can perform.
The isosigs are:
cMcabbgqs
cMcabbgdv
cMcabbgqv

There's a single 1 tetrahedron triangulation of S3, but I don't 
think you can get here without just the moves above, so I don't
think this will ever be the result of simplifying with the Pachner
moves.

Isosig:
bkaagj

To embed the triangulations above, the thought is to embed the two
remaining tetrahedra in some standard arrangement in R3, and 
connect the endpoints on the faces with affine maps taking each
face to the face it's glued to.

"""


def standard_tri_to_tri_matrix(a1, a2, a3):
    """
    Take the standard triangle (in z=1 hyperplane) with coordinates [0,0,0,1],
    [1,0,0,1], and [0,1,0,1] to the triangle with coordinates a1, a2, and a3 via
    a matrix fixing the z=1 hyperplane.
    """
    a1 = vector(a1)
    a2 = vector(a2)
    a3 = vector(a3)
    a1_3d = a1[:3]
    a2_3d = a2[:3]
    a3_3d = a3[:3]
    cross_3d = (a2_3d-a1_3d).cross_product(a3_3d-a1_3d)
    cross = list(cross_3d)
    cross.append(0)
    cross = vector(cross)
    return Matrix([a2-a1,a3-a1,cross,a1]).transpose()

def affine_triangle_map(a1,a2,a3,b1,b2,b3):
    """
    Given two triangles a1,a2,a3 and b1,b2,b3 in the z=1 hyperplane, compute
    a 4x4 map fixing the z=1 hyperplane and taking the first triangle to the 
    second. Takes a1 to b1, a2 to b2, and a3 to b3.
    """
    return standard_tri_to_tri_matrix(b1, b2, b3)*standard_tri_to_tri_matrix(a1, a2, a3).inverse()

        


