from geometric_object_neighborhood import GeometricObjectNeighborhood
from geodesic_neighborhood import GeodesicNeighborhood
from cusp_neighborhood_neighborhood import CuspNeighborhoodNeighborhood

from snappy.math_basics import correct_min, correct_max, upper

def epsilon_from_neighborhood_pair(
        neighborhood1 : GeometricObjectNeighborhood,
        neighborhood2 : GeometricObjectNeighborhood,
        object_distance, verified : bool):

    is_cusp1 = isinstance(neighborhood1, CuspNeighborhoodNeighborhood)
    is_cusp2 = isinstance(neighborhood2, CuspNeighborhoodNeighborhood)
    if is_cusp1 and is_cusp2:
        return epsilon_from_cusp_neighborhood_pair(
            neighborhood1, neighborhood2, object_distance)

    if neighborhood1.index == neighborhood2.index:
        return neighborhood1.epsilon_from_radius(object_distance / 2)

    return epsilon_from_neighborhood_pair_generic(
        neighborhood1, neighborhood2, object_distance, verified=verified)

def _my_max(values):
    n = len(values)
    if n == 0:
        raise RuntimeError("No values to take max")
    if n == 1:
        return values[0]
    return correct_max(values)

def epsilon_from_neighborhood_pair_generic(
        neighborhood1, neighborhood2, object_distance, verified):

    lower_bound_epsilon = _my_max([
        neighborhood.geodesic_info.length.real()
        for neighborhood in [ neighborhood1, neighborhood2 ]
        if isinstance(neighborhood, GeodesicNeighborhood) ])

    lower_bound_radius1 = neighborhood1.radius_from_epsilon(lower_bound_epsilon)
    lower_bound_radius2 = neighborhood2.radius_from_epsilon(lower_bound_epsilon)
    lower_bound_total_radius = lower_bound_radius1 + lower_bound_radius2

    if lower_bound_total_radius > object_distance:
        return lower_bound_epsilon

    if lower_bound_total_radius < object_distance:
        return epsilon_from_neighborhood_pair_implicit(
            neighborhood1, neighborhood2,
            object_distance,
            lower_bound_epsilon,
            verified)

    if verified:
        upper_bound_epsilon1 = (
            neighborhood1.epsilon_from_radius(object_distance))
        upper_bound_epsilon2 = (
            neighborhood2.epsilon_from_radius(object_distance))
        upper_bound_epsilon = correct_min([
            upper_bound_epsilon1, upper_bound_epsilon2])
        return lower_bound_epsilon.union(upper_bound_epsilon)
    else:
        return lower_bound_epsilon

def epsilon_from_neighborhood_pair_implicit(
        neighborhood1, neighborhood2,
        object_distance,
        lower_bound_epsilon,
        verified):

    initial_epsilon = correct_max([neighborhood1.epsilon,neighborhood2.epsilon])

    epsilon = epsilon_from_neighborhood_pair_newton(
        neighborhood1, neighborhood2,
        object_distance,
        initial_epsilon,
        lower_bound_epsilon,
        verified=verified)

    if verified:
        return epsilon_from_neighborhood_pair_newton_interval(
            neighborhood1, neighborhood2,
            object_distance,
            epsilon)
    else:
        return epsilon

_max_iterations_newton = 500

def epsilon_from_neighborhood_pair_newton(
        neighborhood1, neighborhood2,
        object_distance,
        initial_epsilon,
        lower_bound_epsilon,
        verified):

    RF = initial_epsilon.parent()
    epsilon = initial_epsilon

    pos_diff = None
    neg_diff = None

    abs_diff = None
    for i in range(_max_iterations_newton):
        radius1, d1 = neighborhood1.radius_and_derivative_from_epsilon(epsilon)
        radius2, d2 = neighborhood2.radius_and_derivative_from_epsilon(epsilon)

        radius = radius1 + radius2
        d = d1 + d2

        diff = object_distance - radius

        if diff >= 0:
            if pos_diff is not None:
                if not diff < pos_diff:
                    return epsilon
            pos_diff = diff
        else:
            if neg_diff is not None:
                if not -diff < -neg_diff:
                    return epsilon
            neg_diff = diff

        new_epsilon = epsilon + diff / d
        safe_epsilon = (epsilon - lower_bound_epsilon) / 8 + lower_bound_epsilon
        epsilon = correct_max([new_epsilon, safe_epsilon])

        if verified:
            epsilon = RF(epsilon.center())

    raise RuntimeError("Newton method did not converge")

def epsilon_from_neighborhood_pair_newton_interval(
        neighborhood1 : GeometricObjectNeighborhood,
        neighborhood2 : GeometricObjectNeighborhood,
        object_distance,
        epsilon_sample):

    radius1 = neighborhood1.radius_from_epsilon(epsilon_sample)
    radius2 = neighborhood2.radius_from_epsilon(epsilon_sample)
    radius = radius1 + radius2

    diff = object_distance - radius

    epsilon_interval = epsilon_sample

    for i in range(50):
        d1 = neighborhood1.radius_derivative_from_epsilon(epsilon_interval)
        d2 = neighborhood2.radius_derivative_from_epsilon(epsilon_interval)
        d = d1 + d2

        new_epsilon_interval = epsilon_sample + diff / d

        if new_epsilon_interval in epsilon_interval:
            return new_epsilon_interval

        epsilon_interval = epsilon_interval.union(new_epsilon_interval)

    raise RuntimeError(
        "Interval Newton method did not converge for neighborhood pair %r %r, object_distance = %r, epsilon_sample = %r, epsilon_interval = %r" % (
            neighborhood1, neighborhood2, object_distance, epsilon_sample, epsilon_interval))


def epsilon_from_cusp_neighborhood_pair(
        neighborhood1 : CuspNeighborhoodNeighborhood,
        neighborhood2 : CuspNeighborhoodNeighborhood,
        object_distance):
    # Radius1 + radius2 = minimum
    # log(h / _factor1) + log(h / _factor2) = minimum
    # log(h^2 / (_factor1 * _factor2)) = minimum
    # h^2 = exp(minimum) / (_factor1 * _factor2)

    len_prod = neighborhood1.euclidean_length * neighborhood2.euclidean_length
    h = (object_distance.exp() * len_prod).sqrt()
    return 2 * (h/2).arcsinh()
