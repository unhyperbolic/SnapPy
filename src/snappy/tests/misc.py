"""

Canonical retriangulation
-------------------------

Some cases that should be rejected

>>> M = Manifold("m004(3,4)")
>>> M.canonical_retriangulation() # doctest: +ELLIPSIS
Traceback (most recent call last):
...
ValueError: Canonical retriangulation needs all cusps to be complete.

sage: M.canonical_retriangulation(verified=True) # doctest: +ELLIPSIS
Traceback (most recent call last):
...
ValueError: Canonical retriangulation needs all cusps to be complete.

Cusp areas
----------

>>> M = Manifold('o9_44210')
>>> M.cusp_areas(policy='greedy') # doctest: +NUMERIC9
[7.053940530873898, 3.2712450270, 2.7091590087]
>>> M.cusp_areas(policy='greedy', first_cusps=[]) # doctest: +NUMERIC9
[7.053940530873898, 3.2712450270, 2.7091590087]
>>> M.cusp_areas(policy='greedy', first_cusps=[0,]) # doctest: +NUMERIC9
[7.053940530873898, 3.2712450270, 2.7091590087]
>>> M.cusp_areas(policy='greedy', first_cusps=[0,1]) # doctest: +NUMERIC9
[7.053940530873898, 3.2712450270, 2.7091590087]
>>> M.cusp_areas(policy='greedy', first_cusps=[0,1,2]) # doctest: +NUMERIC9
[7.053940530873898, 3.2712450270, 2.7091590087]
>>> M.cusp_areas(policy='greedy', first_cusps=[0,2,1]) # doctest: +NUMERIC9
[7.053940530873898, 2.3513135103, 3.7690945490]
>>> M.cusp_areas(policy='greedy', first_cusps=[1,]) # doctest: +NUMERIC9
[2.30025338030798, 10.0315765558665, 0.883442685721903]

Cusp translations
-----------------

>>> M = Manifold("s776")
>>> M.cusp_translations(policy = 'greedy', first_cusps = [], bits_prec = 100) # doctest: +NUMERIC21
[(0.70710678118654752440084436210 + 1.8708286933869706927918743662*I, 2.8284271247461900976033774484), (0.35355339059327376220042218105 + 0.93541434669348534639593718308*I, 1.4142135623730950488016887242), (0.35355339059327376220042218105 + 0.93541434669348534639593718308*I, 1.4142135623730950488016887242)]


"""

if not __doc__:
    raise Exception("doc string with tests was not recognized.")

from .. import Manifold, ManifoldHP, Triangulation, TriangulationHP
