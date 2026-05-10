"""

>>> myOrbifold = Orbifold(os.path.join(test_files_paths[0], 'example.orb'))
>>> myOrbifold.volume() # doctest: +NUMERIC9
0.117838420347115

"""

if not __doc__:
    raise Exception("doc string with tests was not recognized.")

import os
from ..extensions.Orb import Orbifold
from .files import __path__ as test_files_paths
