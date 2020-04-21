from snappy.snap import t3mlite as t3m
from snappy import Triangulation

from snappy.SnapPy import matrix, vector

from snappy.snap.mcomplex_base import *

from snappy.dev.vericlosed import compute_approx_hyperbolic_structure_orb
from snappy.dev.vericlosed.truncatedComplex import *

from .hyperboloid_utilities import *

from snappy.verify.upper_halfspace import FinitePoint
from sage.all import RealField, ComplexField

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

def _path_finder(perm, tet_num):
    _already_seen = set()
    
    def _try(this_perm, path):
        if perm.tuple() == this_perm.tuple():
            return path
        for letter in ['alpha', 'beta', 'gamma']:
            edge = TruncatedComplex.Edge(letter, (tet_num, this_perm))
            dummy, new_perm = edge.tet_and_perm_of_end()

            if not new_perm.tuple() in _already_seen:
                _already_seen.add(new_perm.tuple())
                r = _try(new_perm, path + [ edge ])
                if not r is None:
                    return r

    return _try(t3m.Perm4([0,1,2,3]), [])

class FiniteTrigRaytracingData(McomplexEngine):
    @staticmethod
    def from_triangulation(triangulation, areas = None, insphere_scale = 0.05):

        hyperbolic_structure = compute_approx_hyperbolic_structure_orb(triangulation)

        r = FiniteTrigRaytracingData(hyperbolic_structure)

        r.RF = hyperbolic_structure.edge_lengths[0].parent()

        r._compute_tet_vertices()
        r._compute_edge_ends()
        r._compute_planes()
        r._compute_face_pairings()

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
                    _compute_path(path, tet.Index)).inverse()

                p = FinitePoint(
                    ComplexField()(0),
                    RealField()(1))

                #print(path)
                #print(p.translate_PGL(m))

                #print(path)
                #print(m)
                
                # m = matrix([[m[1,1],m[1,0]],[m[0,1],m[0,0]]])


                return GL2C_to_O13(m) * c, p.translate_PGL(m)

            #tet.R13_vertices = {
            #    t3m.V0 : _compute_vertex([]),
            #    t3m.V1 : _compute_vertex(['gamma', 'alpha']),
            #    t3m.V2 : _compute_vertex(['beta', 'alpha']),
            #    t3m.V3 : _compute_vertex(['beta','gamma','beta','alpha'])}

            tet.R13_vertices = {
                t3m.V0 : _compute_vertex(['gamma'])[0],
                t3m.V1 : _compute_vertex(['alpha'])[0],
                t3m.V2 : _compute_vertex(['beta', 'alpha', 'gamma'])[0],
                t3m.V3 : _compute_vertex(['beta', 'gamma', 'beta', 'alpha', 'gamma'])[0] }

            tet.A_vertices = {
                t3m.V0 : _compute_vertex(['gamma'])[1],
                t3m.V1 : _compute_vertex(['alpha'])[1],
                t3m.V2 : _compute_vertex(['beta', 'alpha', 'gamma'])[1],
                t3m.V3 : _compute_vertex(['beta', 'gamma', 'beta', 'alpha', 'gamma'])[1] }

            for i in t3m.ZeroSubsimplices:
                for j in t3m.ZeroSubsimplices:
                    if i != j:
                        print("edge", i, j)
                        edge_length = self.hyperbolic_structure.edge_lengths[tet.Class[i | j].Index]
                        print(edge_length)
                        print(R13_dot(tet.R13_vertices[i], tet.R13_vertices[j]))

            for v in t3m.ZeroSubsimplices:
                print("vertex", v, R13_time_vector_to_upper_halfspace(
                        tet.R13_vertices[v]))
                    

    def _compute_edge_ends(self):
        for tet in self.mcomplex.Tetrahedra:
            def _compute_edge_ends(path):
                cs = [ vector([1,  1, 0, 0]),
                       vector([1, -1, 0, 0]) ]
                if not path:
                    return cs

                m = self.hyperbolic_structure.pgl2_matrix_for_path(
                    _compute_path(path, tet.Index)).inverse()
                
                return [ GL2C_to_O13(m) * c for c in cs ]

            tet.R13_edge_ends = {
                t3m.E01 : _compute_edge_ends([]),
                t3m.E02 : _compute_edge_ends(['beta']),
                t3m.E12 : _compute_edge_ends(['beta','alpha','beta']),
                t3m.E03 : _compute_edge_ends(['gamma','beta']),
                t3m.E13 : _compute_edge_ends(['alpha','gamma','beta']),
                t3m.E23 : _compute_edge_ends(['beta','alpha','gamma','beta']) }                

    def _compute_planes(self):
        for tet in self.mcomplex.Tetrahedra:
            def _compute_plane(path):
                c = vector([0.0, 0.0, 0.0, -1.0])

                if not path:
                    return c

                m = self.hyperbolic_structure.pgl2_matrix_for_path(
                    _compute_path(path, tet.Index)).inverse()

                m = matrix([[ m[1,1],-m[0,1]],
                            [-m[1,0], m[0,0]]])

                v = c * GL2C_to_O13(m)
                
                return vector([-v[0], v[1], v[2], v[3]])

            tet.R13_planes = {
                t3m.F0 : _compute_plane(['beta', 'alpha', 'beta', 'gamma', 'beta']),
                t3m.F1 : _compute_plane(['beta', 'gamma', 'beta']),
                t3m.F2 : _compute_plane(['gamma']),
                t3m.F3 : _compute_plane(['beta']) }

    def _compute_face_pairing(self, tet, F):
        def first_perm():
            for p in t3m.Perm4.A4():
                if p.image(t3m.F3) == F:
                    return p

        tet0_perm = first_perm()

        #tet0_perm = t3m.Perm4(
        #    { v: k
        #      for k,v
        #      in first_perm().dict.items() })

        #print(tet0_perm)
        
        path0 = _path_finder(tet0_perm, tet.Index)

        # print(path0)

        if path0:
            m0 = self.hyperbolic_structure.pgl2_matrix_for_path(
                path0)
        else:
            m0 = matrix.identity(
                2, ComplexField())

        #print("face", F)

        #for V in t3m.ZeroSubsimplices:
        #    print("   V: ", tet.A_vertices[V].translate_PGL(m0.inverse()))

        return GL2C_to_O13(m0)
        

        tet1_perm = tet.Gluing[F] * tet0_perm
        tet1_perm = tet0_perm * tet.Gluing[F]

        path0 = _path_finder(tet0_perm, tet.Index)
        path1 = _path_finder(tet1_perm, tet.Neighbor[F].Index)

        if path0:
            m0 = self.hyperbolic_structure.pgl2_matrix_for_path(
                path0)
        else:
            m0 = matrix.identity(
                2, ComplexField())

        if path1:
            m1 = self.hyperbolic_structure.pgl2_matrix_for_path(
                path1)
        else:
            m1 = matrix.identity(
                2, ComplexField())

        m0Inv = matrix([[ m0[1,1],-m0[0,1]],
                        [-m0[1,0], m0[0,0]]])
            
        m1Inv = matrix([[ m1[1,1],-m1[0,1]],
                        [-m1[1,0], m1[0,0]]])

        #return GL2C_to_O13(m0 * m1Inv)
        #return GL2C_to_O13(m0Inv * m1)
        #return GL2C_to_O13(m1 * m0Inv)
        return GL2C_to_O13(m1Inv * m0)



    def _compute_face_pairings(self):
        for tet in self.mcomplex.Tetrahedra:
            tet.O13_matrices = {
                F : self._compute_face_pairing(tet, F)
                for F in t3m.TwoSubsimplices }

    def get_uniform_bindings(self):
        otherTetNums = [
            tet.Neighbor[F].Index
            for tet in self.mcomplex.Tetrahedra
            for F in t3m.TwoSubsimplices ]

        enteringFaceNums = [
            tet.Gluing[F][f]
            for tet in self.mcomplex.Tetrahedra
            for f, F in enumerate(t3m.TwoSubsimplices) ]

        SO13tsfms = [
            tet.O13_matrices[F]
            for tet in self.mcomplex.Tetrahedra
            for F in t3m.TwoSubsimplices ]

        R13Vertices = [
            tet.R13_vertices[V]
            for tet in self.mcomplex.Tetrahedra
            for V in t3m.ZeroSubsimplices ]

        R13EdgeEnds = [
            edge_end
            for tet in self.mcomplex.Tetrahedra
            for E in t3m.OneSubsimplices
            for edge_end in tet.R13_edge_ends[E] ]
            
        planes = [
            tet.R13_planes[F]
            for tet in self.mcomplex.Tetrahedra
            for F in t3m.TwoSubsimplices ]
        
        return {
            'otherTetNums' :
                ('int[]', otherTetNums),
            'enteringFaceNums' :
                ('int[]', enteringFaceNums),
            'TetrahedraBasics.SO13tsfms' :
                ('mat4[]', SO13tsfms),
            'TetrahedraBasics.planes' :
                ('vec4[]', planes),
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
