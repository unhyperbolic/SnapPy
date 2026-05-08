from ..sage_helper import _within_sage, sage_method
from .. import snap
from . import exceptions

__all__ = [
    'check_pos_imaginary_parts_and_lifts_to_log_gluing_equations',
    'verify_hyperbolicity' ]

if _within_sage:
    from sage.symbolic.constants import pi
    from ..sage_helper import I, vector

class FalseTuple(tuple):
    def __nonzero__(self):
        return False

class NonIntegralFillingsError(RuntimeError):
    """
    Exception raised when Manifold has non-integral fillings, e.g.,
    for m004(1.1,1).
    """
    def __init__(self, manifold):
        self.manifold = manifold
        super.__init__('Manifold has non-integral Dehn-filings: %s' % manifold)

@sage_method
def check_pos_imaginary_parts_and_lifts_to_log_gluing_equations(
        manifold, shape_intervals) -> None:
    """
    Given a SnapPy manifold manifold and complex intervals for the shapes
    shape_intervals that are certified to contain a solution to the
    rectangular gluing equations, verify that the logarithmic gluing equations
    are also fulfilled and that all shapes have positive imaginary part.
    It will raise an exception if the verification fails.
    This is sufficient to prove that the manifold is indeed hyperbolic.

    Since the given interval are supposed to contain a true solution of
    the rectangular gluing equations, the logarithmic gluing equations
    are known to be fulfilled up to a multiple of 2 pi i. Thus it is enough
    to certify that the  absolute error of the logarithmic gluing
    equations is < 0.1. Using interval arithmetic, this function certifies
    this and positivity of the imaginary parts of the shapes::

        sage: from snappy import Manifold
        sage: M = Manifold("m019")
        sage: check_pos_imaginary_parts_and_lifts_to_log_gluing_equations(
        ...    M, M.tetrahedra_shapes('rect', intervals=True))


    The SnapPy triangulation of the hyperbolic manifold t02774 actually used
    to contain negatively oriented tetrahedra::

        sage: M = Manifold("iMzzzQcabcefghhhkxxjqoobo_abBa")
        sage: check_pos_imaginary_parts_and_lifts_to_log_gluing_equations(
        ...    M, M.tetrahedra_shapes('rect', intervals=True))    # doctest: +IGNORE_EXCEPTION_DETAIL
        Traceback (most recent call last):
        ...
        ShapePositiveImaginaryPartNumericalVerifyError: Numerical verification that shape has positive imaginary part has failed: Im(0.4800996900657? - 0.0019533695046?*I) > 0

    """

    for d in manifold.cusp_info():
        m, l = d['filling']
        if not (m.is_integer() and l.is_integer()):
            raise NonIntegralFillingsError(manifold)

    # A list
    #    z_0 z'_0 z''_0 z_1 z'_1 z''_1 ...
    all_shapes = [
        shape
        for z in shape_intervals
        for shape in [ z, 1 / (1 - z), ((z - 1) / z) ] ]

    # Check that the shapes have positive imaginary part.
    for shape in all_shapes:
        if not shape.imag() > 0:
            raise exceptions.ShapePositiveImaginaryPartNumericalVerifyError(
                shape)

    args = vector(shape.arg() for shape in all_shapes)

    # Compute the args of the LHS of the gluing
    # equations (in logarithmic form)
    #     a_0 * log(z_0) + b_0 * log(z'_0) + c_0 * log(z''_0) + ...
    # Also, see manifold.gluing_equations
    LHSs = manifold.gluing_equations() * args

    # Get the ComplexIntervalField of the shape intervals
    RIF = shape_intervals[0].real().parent()
    # 2 pi i in that field
    two_pi = RIF(2 * pi)
    epsilon = RIF(0.0625)

    num_edges = manifold.num_tetrahedra()

    # Compute the RHS of the gluing equations

    # Edge equations gives two_pi
    RHSs = num_edges * [ two_pi ]

    for complete in manifold.cusp_info('complete?'):
        if complete:
            # Complete cusp gives equation for meridian and longitude
            # that should each give 0.
            RHSs += 2 * [ 0 ]
        else:
            # Incomplete cusp gives equation for filling curve.
            RHSs += [ two_pi ]

    for i, (LHS, RHS) in enumerate(zip(LHSs, RHSs)):
        if not abs(LHS - RHS) < epsilon:
            if i < num_edges:
                raise exceptions.EdgeEquationLogLiftNumericalVerifyError(
                    LHS)
            else:
                raise exceptions.CuspEquationLogLiftNumericalVerifyError(
                    LHS, RHS)

@sage_method
def verify_hyperbolicity(manifold, verbose=False, bits_prec=None,
                         holonomy=False, fundamental_group_args=[], lift_to_SL=True):
    """
    Given an orientable SnapPy Manifold, verifies its hyperbolicity.

    Similar to HIKMOT's :py:meth:`verify_hyperbolicity`, the result is either
    ``(True, listOfShapeIntervals)`` or ``(False, [])`` if verification failed.
    ``listOfShapesIntervals`` is a list of complex intervals (elements in
    sage's ``ComplexIntervalField``) certified to contain the true shapes
    for the hyperbolic manifold.

    Higher precision intervals can be obtained by setting ``bits_prec``::

        sage: from snappy import Manifold
        sage: M = Manifold("m019")
        sage: M.verify_hyperbolicity() # doctest: +NUMERIC12
        (True, [0.780552527850? + 0.914473662967?*I, 0.780552527850? + 0.91447366296773?*I, 0.4600211755737? + 0.6326241936052?*I])

        sage: M = Manifold("t02333(3,4)")
        sage: M.verify_hyperbolicity() # doctest: +NUMERIC9
        (True, [2.152188153612? + 0.284940667895?*I, 1.92308491369? + 1.10360701507?*I, 0.014388591584? + 0.143084469681?*I, -2.5493670288? + 3.7453498408?*I, 0.142120333822? + 0.176540027036?*I, 0.504866865874? + 0.82829881681?*I, 0.50479249917? + 0.98036162786?*I, -0.589495705074? + 0.81267480427?*I])

    One can instead get a holonomy representation associated to the
    verified hyperbolic structure.  This representation takes values
    in 2x2 matrices with entries in the ``ComplexIntervalField``::

        sage: M = Manifold("m004(1,2)")
        sage: success, rho = M.verify_hyperbolicity(holonomy=True)
        sage: success
        True
        sage: trace = rho('aaB').trace(); trace # doctest: +NUMERIC9
        -0.1118628555? + 3.8536121048?*I
        sage: (trace - 2).contains_zero()
        False
        sage: (rho('aBAbaabAB').trace() - 2).contains_zero()
        True

    Here, there is **provably** a fixed holonomy representation rho0
    from the fundamental group G of M to SL(2, C) so that for each
    element g of G the matrix rho0(g) is contained in rho(g).  In
    particular, the above constitutes a proof that the word 'aaB' is
    non-trivial in G.  In contrast, the final computation is
    consistent with 'aBAbaabAB' being trivial in G, but *does not prove
    this*.

    A non-hyperbolic manifold (``False`` indicates that the manifold
    might not be hyperbolic but does **not** certify
    non-hyperbolicity. Sometimes, hyperbolicity can only be verified
    after increasing the precision)::

        sage: M = Manifold("4_1(1,0)")
        sage: M.verify_hyperbolicity()
        (False, [])

    Under the hood, the function will call the ``CertifiedShapesEngine`` to produce
    intervals certified to contain a solution to the rectangular gluing equations.
    It then calls ``check_pos_imaginary_parts_and_lifts_to_log_gluing_equations``
    to verify that the logarithmic gluing equations are fulfilled and that all
    tetrahedra are positively oriented.
    """

    try:
        shape_intervals = manifold.tetrahedra_shapes(
            'rect', bits_prec=bits_prec, intervals=True)
    except (ValueError, RuntimeError):
        if verbose:
            print("Could not certify solution to rectangular gluing equations")
        return FalseTuple((False, []))

    try:
        check_pos_imaginary_parts_and_lifts_to_log_gluing_equations(
            manifold, shape_intervals)
    except exceptions.NumericalVerifyError as e:
        if verbose:
            print(e)
        return FalseTuple((False, []))

    if holonomy:
        hol_rep = snap.interval_reps.holonomy_from_shape_intervals(
            manifold, shape_intervals, fundamental_group_args, lift_to_SL)
        hol_rep.shapes = shape_intervals
        return True, hol_rep
    else:
        return True, shape_intervals
