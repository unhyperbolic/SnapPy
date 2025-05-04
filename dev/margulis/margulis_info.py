from typing import Optional

class MargulisInfo:
    pass

class MargulisTubeInfo(MargulisInfo):
    def __init__(self, word : str, radius, core_curve : Optional[int]):
        self.word = word
        self.radius = radius
        self.core_curve = core_curve
        self.type_ = 'geodesic'

    def __repr__(self):
        return "Tube about geodesic %s%s of radius %r" % (
            self.word, _format_core_curve(self.core_curve), self.radius)

class MargulisCuspNeighborhoodInfo(MargulisInfo):
    def __init__(self, cusp_index, cusp_area):
        self.cusp_index = cusp_index
        self.cusp_area = cusp_area
        self.type_ = 'cusp'

    def __repr__(self):
        return "Cusp Neighborhood for cusp %d of area %r" % (
            self.cusp_index, self.cusp_area)

def _format_core_curve(core_curve : Optional[int]):
    if core_curve is None:
        return ''
    return ' (Core curve of cusp %d)' % core_curve
