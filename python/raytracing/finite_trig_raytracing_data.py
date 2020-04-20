from snappy.snap import t3mlite as t3m
from snappy import Triangulation

from snappy.SnapPy import matrix, vector

from snappy.snap.mcomplex_base import *

from snappy.dev.vericlosed import compute_approx_hyperbolic_structure_orb

from .hyperboloid_utilities import *

from math import sqrt

__all__ = ['FiniteTrigRaytracingData']

class FiniteTrigRaytracingData(McomplexEngine):
    @staticmethod
    def from_triangulation(triangulation, areas = None, insphere_scale = 0.05):

        hyperbolicStructure = compute_approx_hyperbolic_structure_orb(triangulation)

        r = FiniteTrigRaytracingData(hyperbolicStructure.mcomplex)

        r.RF = hyperbolicStructure.edge_lengths[0].parent()

        r._compute_tet_vertices()

        return r

    def __init__(self, mcomplex, snappy_manifold = None):
        super(FiniteTrigRaytracingData, self).__init__(mcomplex)

    def _compute_tet_vertices(self):
        for tet in self.mcomplex.Tetrahedra:
            c = vector([1,0,0,0])

            tet.R13_vertices = {
                t3m.V0 : c,
                t3m.V1 : c,
                t3m.V2 : c,
                t3m.V3 : c }

    def get_uniform_bindings(self):

        R13Vertices = [
            tet.R13_vertices[V]
            for tet in self.mcomplex.Tetrahedra
            for V in t3m.ZeroSubsimplices ]

        return {
            'TetrahedraBasics.R13Vertices' :
                ('vec4[]', R13Vertices),

            'isNonGeometric' :
                ('bool', False),
            'nonGeometricTexture' :
                ('int', 0)}

    def get_compile_time_constants(self):
        return {
            b'##num_tets##' : len(self.mcomplex.Tetrahedra),
            b'##num_cusps##' : len(self.mcomplex.Vertices)
            }

    def initial_view_state(self):
        boost = matrix([[1.0,0.0,0.0,0.0],
                        [0.0,1.0,0.0,0.0],
                        [0.0,0.0,1.0,0.0],
                        [0.0,0.0,0.0,1.0]])
        tet_num = 0
        return (boost, tet_num)

    def update_view_state(self, boost_and_tet_num,
                          m = matrix([[1.0, 0.0, 0.0, 0.0], 
                                      [0.0, 1.0, 0.0, 0.0],
                                      [0.0, 0.0, 1.0, 0.0],
                                      [0.0, 0.0, 0.0, 1.0]])):
        boost, tet_num = boost_and_tet_num

        boost = matrix(boost, ring = self.RF)
        m = matrix(m, ring = self.RF)

        boost = O13_orthonormalize(boost * m)

        return boost, tet_num

        entry_F = -1

        for i in range(100):
            pos = boost.transpose()[0]
            tet = self.mcomplex.Tetrahedra[tet_num]

            amount, F = max(
                [ (R13_dot(pos, tet.R13_planes[F]), F)
                  for F in t3m.TwoSubsimplices ])

            if F == entry_F:
                break
            if amount < 0.0000001:
                break
            
            boost = O13_orthonormalize(tet.O13_matrices[F] * boost)
            tet_num = tet.Neighbor[F].Index
            entry_F = tet.Gluing[F].image(F)

        return boost, tet_num
