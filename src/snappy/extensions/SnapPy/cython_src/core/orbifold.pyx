cdef class Orbifold(Triangulation):

    def __init__(self, spec=None, remove_finite_vertices=True):
        if self.c_triangulation != NULL:
            self.init_hyperbolic_structure()

    def init_hyperbolic_structure(self, force_recompute = False):
        if not self.c_triangulation:
            return
        if self.hyperbolic_structure_initialized and not force_recompute:
            return
        manual = False
        orb_find_hyperbolic_structure(self.c_triangulation, manual)
        self.hyperbolic_structure_initialized = True

    def _orb_cone_fill(self,
                       singular_order : Union[float, list[float]],
                       singular_index : Optional[SupportsIndex] = None) -> None:
        Triangulation._orb_cone_fill(self, singular_order, singular_index)
        manual = False
        orb_find_hyperbolic_structure(self.c_triangulation, manual)
        self._cache.clear(message='Manifold._orb_cone_fill')

    def volume(self):
        cdef Boolean ok

        return Real2Number(orb_volume(self.c_triangulation, &ok))
