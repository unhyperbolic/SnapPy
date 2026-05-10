"""

>>> M = Manifold(os.path.join(test_files_paths[0], 'example_m004.tri'))
>>> M.triangulation_isosig()
'cPcbbbiht_BaCB'

"""

if not __doc__:
    raise Exception("doc string with tests was not recognized.")

import os
from .. import Manifold
from .files import __path__ as test_files_paths
