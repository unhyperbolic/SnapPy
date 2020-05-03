from snappy.snap.t3mlite.mcomplex import *
from sage.all import Rational, vector, Matrix
from snappy.snap.t3mlite.simplex import *
import snappy
import random
class BarycentricPoint(object):
    """
    A quadruple of Sage rational numbers whose sum is 1.
    """
    def __init__(self, c0, c1, c2, c3):
        c0 = Rational(c0)
        c1 = Rational(c1)
        c2 = Rational(c2)
        c3 = Rational(c3)
        self.vector = vector([c0,c1,c2,c3])

        if sum(self.vector) != 1:
            raise Exception("Barycentric point doesn't sum to 1")

    def __repr__(self):
        return self.vector.__repr__()

    def __eq__(self, other):
        return self.vector == other.vector

    def has_negative_coordinate(self):
        for l in self.vector:
            if l<0:
                return True
        return False

    def negative_coordinates(self):
        return [l for l in self.vector if l<0]
    
    def zero_coordinates(self):
        return [i for i in range(4) if self.vector[i] == 0]

    def is_on_boundary_face(self):
        return len(self.zero_coordinates()) == 1
    
    def convex_combination(self, other, t):
        c0,c1,c2,c3 = self.vector*(1-t)+other.vector*t
        return BarycentricPoint(c0,c1,c2,c3)

    def transform(self, matrix):
        v = self.vector
        new_v = matrix*v
        return BarycentricPoint(*new_v)

    def to_3d_point(self):
        return self.vector[0:3]

    def permute(self, perm):
        """
        Start with a permutation perm, which represents a map from the vertices of
        a tetrahedron T in which a point lies, to the vertices of another 
        tetrahedron S which is glued to T. Then, translate
        the barycentric coordinates of the point in T to the corresponding
        barycentric coordinates of the point in S.
        On the interior, this doesn't make much sense; it really should be 
        used for points on the common boundary triangle.
        """
        v = self.vector
        new_v = [0]*4
        for i in range(4):
            new_v[perm[i]] = v[i]
        return BarycentricPoint(*new_v)
    
P = BarycentricPoint

class BarycentricArc(object):
    """
    A line segment between two endpoints in barycentric coordinates.
    """
    def __init__(self, start, end):
        self.start = start
        self.end = end

    def reversed(self):
        return BarycentricArc(self.end, self.start)

    
    def trim_start(self):
        """
        If the start point has negative coordinates, intersect the
        arc with the 3-simplex. This is essentially moving the 
        start point along the arc until it intersects a face of
        the 3-simplex.
        """
        ts = []
        if not self.start.has_negative_coordinate():
            return self
        for i in range(4):
            start_i, end_i = self.start.vector[i], self.end.vector[i]
            if start_i == end_i:
                continue
            else:
                t = -start_i/(end_i-start_i)
                if 0 <= t <= 1:
                    ts.append(t)
                else:
                    continue
        # find the closest point to start which has
        # nonnegative coordinates along the line from start to end
        for t in sorted(ts):
            new_start = self.start.convex_combination(self.end,t)
            if new_start.has_negative_coordinate():
                continue
            else:
                return BarycentricArc(new_start, self.end)

        return None


    def trim_end(self):
        return self.reversed().trim_start().reversed()
        
    def __repr__(self):
        return '[{},{}]'.format(self.start, self.end)

    def __eq__(self, other):
        return (self.start == other.start) and (self.end == other.end)

    def direction(self):
        return self.end.vector-self.start.vector

    def transform(self, matrix):
        new_start = self.start.transform(matrix)
        new_end = self.end.transform(matrix)
        return BarycentricArc(new_start, new_end)

    def is_point(self):
        return self.start == self.end

    def to_3d_points(self):
        return (self.start.to_3d_point(), self.end.to_3d_point())



Arc = BarycentricArc

A = vector([Rational('1'),Rational('0'),Rational('0'),Rational('1')])
B = vector([Rational('-1'),Rational('1'),Rational('0'),Rational('1')])
C =vector([Rational('-1'),Rational('-1'),Rational('0'),Rational('1')])
N =vector([Rational('0'),Rational('0'),Rational('1'),Rational('1')])
S = vector([Rational('0'),Rational('0'),Rational('-1'),Rational('1')])


ABNS = Matrix([A,B,N,S]).T
BCNS = Matrix([B,C,N,S]).T
CANS = Matrix([C,A,N,S]).T

three = [ABNS, BCNS, CANS]


NABC = Matrix([N,A,B,C]).T
ABCS = Matrix([A,B,C,S]).T

two = [NABC, ABCS]

ABNSI = ABNS.inverse()
BCNSI = BCNS.inverse()
CANSI = CANS.inverse()

three_inv = [ABNSI, BCNSI, CANSI]

NABCI = NABC.inverse()
ABCSI = ABCS.inverse()

two_inv = [NABCI, ABCSI]

#trans_32_1 = [M*N for M in three for N in two_inv]
#trans_23_1 = [M*N for M in two for N in three_inv]
trans_32 = [N*M for M in three for N in two_inv]
trans_23 = [N*M for M in two for N in three_inv]


one_half = Rational('1/2')
one_third = Rational('1/3')
one_sixth = Rational('1/6')
c1 = Rational('1/5')
c2 = Rational('1/7')
c3 = Rational('23/35')
t3m_face_to_midpoint = {7:BarycentricPoint(one_third,one_third,one_third,0),
                        11:BarycentricPoint(one_third,one_third,0,one_third),
                        13:BarycentricPoint(one_third,0,one_third,one_third),
                        14:BarycentricPoint(0,one_third,one_third,one_third)}
t3m_face_to_point = {7:BarycentricPoint(c1,c2,c3,0),
                        11:BarycentricPoint(c1,c2,0,c3),
                        13:BarycentricPoint(c1,0,c2,c3),
                        14:BarycentricPoint(0,c1,c2,c3)}

def plot_arcs(arcs):
    import matplotlib as mpl
    from mpl_toolkits.mplot3d import Axes3D
    import matplotlib.pyplot as plt
    from itertools import combinations
    fig = plt.figure()
    ax = fig.gca(projection='3d')

    tet_points = [[0,0,0],[1,0,0],[0,1,0],[0,0,1]]
    for p1, p2 in combinations(tet_points,2):
        x = [p1[0],p2[0]]
        y = [p1[1],p2[1]]
        z = [p1[2],p2[2]]        
        ax.plot(x,y,z,color='black')
    for arc in arcs:
        p1, p2 = arc.to_3d_points()
        x = [p1[0],p2[0]]
        y = [p1[1],p2[1]]
        z = [p1[2],p2[2]]
        ax.plot(x,y,z)
    plt.show()

def smash(mc):    
    if len(mc.Vertices) > 1:
        for e in mc.Edges:
            if mc.smash_star(e):
                return 1
    return 0

def smash_all(mc):
    while len(mc.Vertices)>1:
        ans = smash(mc)


class BarycentricTetrahedronEmbedding(object):
    """
    A map from a tetrahedron with PL arcs in barycentric coordinates
    to R^3 embedded in R^4 as the z=1 hyperplane. The map is 
    described by choosing where the vertices go (giving three 3d 
    vectors, as the 4th coordinate is always 1).
    """
    def __init__(self, tetrahedron, vertex_images):
        self.tetrahedron = tetrahedron
        self.vertex_images = vertex_images
        assert len(vertex_images)==4
        for vi in vertex_images:
            assert len(vi)==3
            for i in range(3):
                vi[i] = Rational(vi[i])
            vi.append(Rational('1'))
        self.matrix = Matrix(vertex_images).T
        assert self.matrix.det() > 0
        
    def transfer_arcs_from(self, other_embedding):
        """
        Given two BarycentricTetrahedronEmbeddings, we get a way
        to map barycentric arcs in one tetrahedron to the other.
        To transfer from other to self, map the arcs in other to
        R^4 using other_embedding, and map those arcs back to
        barycentric coordinates in self using the inverse embedding
        of self.
        """
        transition_matrix = self.matrix.inverse()*other_embedding.matrix

        return [arc.transform(transition_matrix) for arc in other_embedding.tetrahedron.arcs]

    def info(self):
        self.tetrahedron.info()
        print(self.matrix)

    def plot(self):
        import matplotlib as mpl
        from mpl_toolkits.mplot3d import Axes3D
        import matplotlib.pyplot as plt
        from itertools import combinations
        fig = plt.figure()
        ax = fig.gca(projection='3d')

        
        tet_points = [v[:3] for v in self.vertex_images]
        for p1, p2 in combinations(tet_points,2):
            x = [p1[0],p2[0]]
            y = [p1[1],p2[1]]
            z = [p1[2],p2[2]]        
            ax.plot(x,y,z,color='black')
        for arc in self.tetrahedron.arcs:
            p1 = self.matrix * arc.start.vector
            p2 = self.matrix * arc.end.vector

            x = [p1[0],p2[0]]
            y = [p1[1],p2[1]]
            z = [p1[2],p2[2]]
            ax.plot(x,y,z)
        plt.show()



def barycentric_face_embedding(arrow):
    """
    The arrow here is a directed edge in a specified tetrahedron.
    It also specifies a face (the face which is disjoint from the
    arrow, except for the head).
    This helper function takes the arrow and embeds the face of the
    arrow in the xy-plane, with the two tetrahedra on either side
    of the face embedded in the upper and lower half-spaces.
    The specific coordinates are labeled below; A, B, and C are the
    vertices of the image of the face in the xy-plane, and N and S
    are are images of the vertices of the two tetrahedron not in
    the face. 
    """
    arrow = arrow.copy()
    next_arrow = arrow.glued()
    vertex_images = [None]*4
    vertex_images[VertexIndex[arrow.tail()]] = [0,0,1] #N
    vertex_images[VertexIndex[arrow.head()]] = [1,0,0] #A
    arrow.opposite()
    vertex_images[VertexIndex[arrow.tail()]] = [-1,-1,0] #C
    vertex_images[VertexIndex[arrow.head()]] = [-1,1,0] #B

    next_vertex_images = [None]*4
    next_vertex_images[VertexIndex[next_arrow.tail()]] = [1,0,0] #A
    next_vertex_images[VertexIndex[next_arrow.head()]] = [0,0,-1] #S
    next_arrow.opposite()
    next_vertex_images[VertexIndex[next_arrow.tail()]] = [-1,-1,0] #C
    next_vertex_images[VertexIndex[next_arrow.head()]] = [-1,1,0] #B

    return [BarycentricTetrahedronEmbedding(arrow.Tetrahedron,vertex_images),
            BarycentricTetrahedronEmbedding(next_arrow.Tetrahedron, next_vertex_images)]


def barycentric_edge_embedding(arrow):
    """
    Take the arrow corresponding to an edge of valence 3. This 
    function then creates an embedding of the three tetrahedra
    glued in pairs around the edge into R^3. The embedding is
    defined so that the edge goes from N to S, as labeled below, 
    The arrow goes from A to B; A, B, and C form a triangle in 
    the xy-plane.
    Note that this arrangement has the same image in R^3 as the 
    barycentric_face_embedding above -- that's by design, so that 
    we can use these two maps to transfer the arcs in barycentric
    coordinates under two-three and two-three moves.
    """
    assert len(arrow.linking_cycle())==3
    arrow = arrow.copy()
    next = arrow.glued()
    after_next = next.glued()

    vertex_images = [None]*4
    vertex_images[VertexIndex[arrow.tail()]] = [1,0,0] #A
    vertex_images[VertexIndex[arrow.head()]] = [-1,1,0] #B
    arrow.opposite()
    vertex_images[VertexIndex[arrow.tail()]] = [0,0,-1] #S
    vertex_images[VertexIndex[arrow.head()]] = [0,0,1] #N

    next_vertex_images = [None]*4
    next_vertex_images[VertexIndex[next.tail()]] = [-1,1,0] #B
    next_vertex_images[VertexIndex[next.head()]] = [-1,-1,0] #C
    next.opposite()
    next_vertex_images[VertexIndex[next.tail()]] = [0,0,-1] #S
    next_vertex_images[VertexIndex[next.head()]] = [0,0,1] #N

    after_next_vertex_images = [None]*4
    after_next_vertex_images[VertexIndex[after_next.tail()]] = [-1,-1,0] #C
    after_next_vertex_images[VertexIndex[after_next.head()]] = [1,0,0] #A
    after_next.opposite()
    after_next_vertex_images[VertexIndex[after_next.tail()]] = [0,0,-1] #S
    after_next_vertex_images[VertexIndex[after_next.head()]] = [0,0,1] #N

    return [BarycentricTetrahedronEmbedding(arrow.Tetrahedron,vertex_images),
            BarycentricTetrahedronEmbedding(next.Tetrahedron, next_vertex_images),
            BarycentricTetrahedronEmbedding(after_next.Tetrahedron,after_next_vertex_images)]
    
        


def knot_triangulation(name, add_arc = True):
    """
    Given the SnapPy manifold name of a knot exerior as a string,
    return an Mcomplex with the barycentric arcs representing
    the knot.
    """
    M = snappy.Triangulation(name+"(1,0)")
    MC = snappy.snap.t3mlite.Mcomplex(M._unsimplified_filled_triangulation())
    smash_all(MC)
    MC.rebuild()
    if add_arc:
        MC._add_arcs_around_valence_one_edge()
    return MC



    
