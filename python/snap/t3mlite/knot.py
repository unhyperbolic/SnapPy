from snappy.snap.t3mlite.mcomplex import *
from sage.all import Rational, vector, Matrix

class BarycentricPoint(object):
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
        return [l for l in self.vector if l == 0]

    def convex_combination(self, other, t):
        c0,c1,c2,c3 = self.vector*(1-t)+other.vector*t
        return BarycentricPoint(c0,c1,c2,c3)

    def transform(self, matrix):
        v = self.vector
        new_v = matrix*v
        return BarycentricPoint(*new_v)

    def to_3d_point(self):
        return self.vector[0:3]
    
P = BarycentricPoint

class BarycentricArc(object):
    def __init__(self, start, end):
        self.start = start
        self.end = end

    def reversed(self):
        return BarycentricArc(self.end, self.start)
    
    def trim_start(self):
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

        for t in sorted(ts): #find the closest point to start which has nonnegative coordinates along the line from start to end
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

def test():
    M = snappy.Triangulation('K3a1(1,0)')
    MC = snappy.snap.t3mlite.Mcomplex(M._unsimplified_filled_triangulation())
    MC._add_arcs_around_valence_one_edge()
    for i in range(100):
        MC._two_three_simplify(15)
        print([tet.arcs for tet in MC.Tetrahedra if tet.arcs])


    most_arcs = max(MC.Tetrahedra, key=lambda tet: len(tet.arcs)).arcs
    plot_arcs(most_arcs)
