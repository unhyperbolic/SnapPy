from snappy.snap import t3mlite as t3m
from snappy import Triangulation

from snappy.SnapPy import matrix, vector

from snappy.snap.mcomplex_base import *

from snappy.dev.vericlosed import compute_approx_hyperbolic_structure_orb
from snappy.dev.vericlosed.truncatedComplex import *

from .hyperboloid_utilities import *

from math import sqrt

__all__ = ['FiniteTrigRaytracingData']

def _compute_path(letters, tet_num):
    path = []
    tet_and_perm = tet_num, t3m.Perm4([0,1,2,3])
    for letter in letters:
        edge = TruncatedComplex.Edge(letter, tet_and_perm)
        path.append(edge)
        tet_and_perm = edge.tet_and_perm_of_end()
    return path

class FiniteTrigRaytracingData(McomplexEngine):
    @staticmethod
    def from_triangulation(triangulation, areas = None, insphere_scale = 0.05):

        hyperbolic_structure = compute_approx_hyperbolic_structure_orb(triangulation)

        r = FiniteTrigRaytracingData(hyperbolic_structure)

        r.RF = hyperbolic_structure.edge_lengths[0].parent()

        r._compute_tet_vertices()
        r._compute_edge_ends()

        return r

    def __init__(self, hyperbolic_structure):
        super(FiniteTrigRaytracingData, self).__init__(
            hyperbolic_structure.mcomplex)
        self.hyperbolic_structure = hyperbolic_structure

    def _compute_tet_vertices(self):
        for tet in self.mcomplex.Tetrahedra:

            def _compute_vertex(path):

                c = vector([1,0,0,0])
                
                if not path:
                    return c
                
                m = self.hyperbolic_structure.pgl2_matrix_for_path(
                    _compute_path(path, tet.Index))
                
                return c * GL2C_to_O13(m)

            tet.R13_vertices = {
                t3m.V0 : _compute_vertex([]),
                t3m.V1 : _compute_vertex(['alpha']),
                t3m.V2 : _compute_vertex(['beta','gamma','alpha']),
                t3m.V3 : _compute_vertex(['beta','gamma','beta','gamma','alpha'])}

    def _compute_edge_ends(self):
        for tet in self.mcomplex.Tetrahedra:
            def _compute_edge_ends(path):
                cs = [ vector([1,  1, 0, 0]),
                       vector([1, -1, 0, 0]) ]
                if not path:
                    return cs

                m = self.hyperbolic_structure.pgl2_matrix_for_path(
                    _compute_path(path, tet.Index))
                
                return [ c * GL2C_to_O13(m) for c in cs ]

            tet.R13_edge_ends = {
                t3m.E01 : _compute_edge_ends([]),
                t3m.E02 : _compute_edge_ends(['beta']),
                t3m.E12 : _compute_edge_ends(['beta','alpha','beta']),
                t3m.E03 : _compute_edge_ends(['gamma','beta']),
                t3m.E13 : _compute_edge_ends(['alpha','gamma','beta']),
                t3m.E23 : _compute_edge_ends(['beta','alpha','gamma','beta']) }                

    def get_uniform_bindings(self):

        R13Vertices = [
            tet.R13_vertices[V]
            for tet in self.mcomplex.Tetrahedra
            for V in t3m.ZeroSubsimplices ]

        R13EdgeEnds = [
            edge_end
            for tet in self.mcomplex.Tetrahedra
            for E in t3m.OneSubsimplices
            for edge_end in tet.R13_edge_ends[E] ]
            
        print(R13EdgeEnds)

        return {
            'TetrahedraBasics.R13Vertices' :
                ('vec4[]', R13Vertices),
            'TetrahedraBasics.R13EdgeEnds' :
                ('vec4[]', R13EdgeEnds),
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
