from . import hyperbolicity

__all__ = ['compute_hyperbolic_shapes']

from typing import Optional

def compute_hyperbolic_shapes_oriented(
        manifold, *,
        verified : bool, bits_prec : Optional[int] = None):

    # Get shapes, as intervals if requested
    shapes = manifold.tetrahedra_shapes('rect', intervals=verified,
                                        bits_prec=bits_prec)

    # Check it is a valid hyperbolic structure
    if verified:
        hyperbolicity.check_pos_imaginary_parts_and_lifts_to_log_gluing_equations(
            manifold, shapes)
    else:
        # If not verified, just ask SnapPea kernel for solution type
        sol_type = manifold.solution_type()
        if not sol_type == 'all tetrahedra positively oriented':
            raise RuntimeError(
                "Manifold has non-geometric solution type '%s'." % sol_type)

    return shapes

def compute_hyperbolic_shapes(
        manifold, *,
        verified : bool, bits_prec : Optional[int] = None):

    if verified or bits_prec is not None:
        if manifold.is_orientable():
            return compute_hyperbolic_shapes_oriented(
                manifold,
                verified = verified, bits_prec = bits_prec)
        else:
            return compute_hyperbolic_shapes_oriented(
                manifold.orientation_cover(),
                verified = verified, bits_prec = bits_prec)[::2]
    else:
        sol_type = manifold.solution_type()
        if not sol_type == 'all tetrahedra positively oriented':
            raise RuntimeError(
                "Manifold has non-geometric solution type '%s'." % sol_type)
        return manifold.tetrahedra_shapes('rect')
