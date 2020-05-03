from snappy import *
s = open('example1.orb').read()
T = Triangulation(s)
print(T.num_tetrahedra())
