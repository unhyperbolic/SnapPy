"""
Class hierarchy
---------------

>>> isinstance(Manifold("m004"), Triangulation)
True

>>> isinstance(ManifoldHP("m004"), TriangulationHP)
True

Low precision and high precision comparisons
--------------------------------------------
>>> M   = Manifold("m004")
>>> Mhp = Manifold("m004")
>>> N   = Manifold("m003")
>>> Nhp = Manifold("m003")
>>> M.is_isometric_to(M)
True
>>> M.is_isometric_to(Mhp)
True
>>> Mhp.is_isometric_to(M)
True
>>> Mhp.is_isometric_to(Mhp)
True
>>> M.is_isometric_to(N)
False
>>> M.is_isometric_to(Nhp)
False
>>> Mhp.is_isometric_to(N)
False
>>> Mhp.is_isometric_to(Nhp)
False

>>> O = Triangulation("mvvLALQQQhfghjjlilkjklaaaaaffffffff",
... remove_finite_vertices = False)
>>> O.has_finite_vertices()
True
>>> Ohp = TriangulationHP("mvvLALQQQhfghjjlilkjklaaaaaffffffff",
... remove_finite_vertices = False)
>>> Ohp.has_finite_vertices()
True

>>> len(O.isomorphisms_to(O))
8
>>> len(O.isomorphisms_to(Ohp))
8
>>> len(Ohp.isomorphisms_to(O))
8
>>> len(Ohp.isomorphisms_to(Ohp))
8
>>> len(M.isomorphisms_to(O))
0
>>> len(M.isomorphisms_to(Ohp))
0
>>> len(Mhp.isomorphisms_to(O))
0
>>> len(Mhp.isomorphisms_to(Ohp))
0
"""

if not __doc__:
    raise Exception("doc string with tests was not recognized.")

from .. import Manifold, ManifoldHP, Triangulation, TriangulationHP
