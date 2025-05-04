from geometric_object_neighborhood import GeometricObjectNeighborhood
from geodesic_neighborhood import GeodesicNeighborhood
from cusp_neighborhood_neighborhood import CuspNeighborhoodNeighborhood
from epsilon_from_neighborhood_pair import epsilon_from_neighborhood_pair
from margulis_info import MargulisInfo

from snappy.verify.shapes import compute_hyperbolic_shapes
from snappy.geometric_structure.cusp_neighborhood.tiles_for_cusp_neighborhood import (
    mcomplex_for_tiling_cusp_neighborhoods)
from snappy.geometric_structure.cusp_neighborhood.complex_cusp_cross_section import (
    ComplexCuspCrossSection)
from snappy.geometric_structure.geodesic.add_core_curves import add_r13_core_curves
from snappy.geometric_structure import (
    add_r13_geometry, add_filling_information)
from snappy.len_spec.length_spectrum_geodesic_info import LengthSpectrumGeodesicInfo
from snappy.hyperboloid.distances import (
    distance_r13_horoballs, distance_r13_lines, distance_r13_horoball_line)
from snappy.hyperboloid.horoball import R13Horoball
from snappy.hyperboloid.line import R13Line
from snappy.math_basics import correct_min
from snappy.snap.t3mlite import Mcomplex
from snappy.sage_helper import Infinity
from snappy.tiling.tile import Tile

import heapq

from typing import Iterable, Union

def distance_r13_objects(object1 : Union[R13Line, R13Horoball],
                         object2 : Union[R13Line, R13Horoball]):
    is_horoball1 = isinstance(object1, R13Horoball)
    is_horoball2 = isinstance(object2, R13Horoball)
    if is_horoball1:
        if is_horoball2:
            return distance_r13_horoballs(
                object1.defining_vec, object2.defining_vec)
        else:
            return distance_r13_horoball_line(
                object1.defining_vec, object2)
    else:
        if is_horoball2:
            return distance_r13_horoball_line(
                object2.defining_vec, object1)
        else:
            return distance_r13_lines(
                object1, object2)

def distance_tiles(tile1 : Tile, tile2 : Tile):
    return distance_r13_objects(
        tile1.inverse_lifted_geometric_object,
        tile2.inverse_lifted_geometric_object)

class NeighborhoodPair:
    def __init__(self, infinity):
        self.distance_lifts = infinity
        self.finished : bool = False
        self.epsilon = None

class Neighborhoods:
    def __init__(self, mcomplex : Mcomplex):
        self.epsilon = mcomplex.infinity
        self.neighborhoods : list[GeometricObjectNeighborhood] = []
        self._indices_to_neighborhood_pair : list[list[NeighborhoodPair]] = []
        self.mcomplex : Mcomplex = mcomplex

    def add_neighborhood(
            self,
            neighborhood : GeometricObjectNeighborhood
            ) -> None:
        self.neighborhoods.append(neighborhood)
        self._indices_to_neighborhood_pair.append(
            [ NeighborhoodPair(self.mcomplex.infinity)
              for index in range(len(self.neighborhoods)) ])

    def get_neighborhood_pair(
            self,
            neighborhood1 : GeometricObjectNeighborhood,
            neighborhood2 : GeometricObjectNeighborhood
            ) -> NeighborhoodPair:
        i = neighborhood1.index
        j = neighborhood2.index
        if i >= j:
            return self._indices_to_neighborhood_pair[i][j]
        else:
            return self._indices_to_neighborhood_pair[j][i]

    def thin_part(self) -> list[MargulisInfo]:
        return [ neighborhood.info_for_epsilon(self.epsilon)
                 for neighborhood in self.neighborhoods ]

    def collisions(self) -> list[tuple[int, int]]:
        if self.mcomplex.verified:
            prec_epsilon = 0
        else:
            prec_epsilon = self.mcomplex.RF(1e-6)


        result = []
        for neighborhood1, neighborhood_pairs in zip(
                self.neighborhoods, self._indices_to_neighborhood_pair):
            for neighborhood2, neighborhood_pair in zip(
                    self.neighborhoods, neighborhood_pairs):
                if neighborhood_pair.epsilon is None:
                    continue
                if neighborhood_pair.epsilon > self.epsilon + prec_epsilon:
                    continue
                result.append((neighborhood2.index,neighborhood1.index))
        return result

def compute_cusp_shapes(M, *, bits_prec, verified):
    shapes = compute_hyperbolic_shapes(
        M, verified=verified, bits_prec=bits_prec)
    c = ComplexCuspCrossSection.fromManifoldAndShapes(M, shapes)
    return [
        ComplexCuspCrossSection.cusp_shape(v) if v.is_complete else None
        for v in c.mcomplex.Vertices ]

def mcomplex_for_margulis_number(M, bits_prec, verified):
    mcomplex = mcomplex_for_tiling_cusp_neighborhoods(
        M, bits_prec=bits_prec, verified=verified)
    add_filling_information(mcomplex, M)
    add_r13_core_curves(mcomplex, M)

    if verified:
        mcomplex.infinity = mcomplex.RF(Infinity)
    else:
        mcomplex.infinity = mcomplex.RF(1e20)

    return mcomplex

def add_cusp_to_queue_and_neighborhoods(
        neighborhood_queue : list[GeometricObjectNeighborhood],
        neighborhoods : Neighborhoods,
        vertex,
        cusp_shape
        ) -> None:
    index = len(neighborhoods.neighborhoods)
    neighborhood = CuspNeighborhoodNeighborhood(
        neighborhoods.mcomplex, index, vertex, cusp_shape)
    heapq.heappush(neighborhood_queue, neighborhood)
    neighborhoods.add_neighborhood(neighborhood)

def add_geodesic_to_queue(
        neighborhood_queue : list[GeometricObjectNeighborhood],
        neighborhoods : Neighborhoods,
        geodesic_info : LengthSpectrumGeodesicInfo
        ) -> None:
    index = len(neighborhoods.neighborhoods)
    neighborhood = GeodesicNeighborhood(
        neighborhoods.mcomplex, index, geodesic_info)
    heapq.heappush(neighborhood_queue, neighborhood)

def expand_next_neighborhood(
        neighborhood_queue : list[GeometricObjectNeighborhood],
        neighborhoods : Neighborhoods,
        len_spec : Iterable[LengthSpectrumGeodesicInfo]
    ) -> bool:

    neighborhood = heapq.heappop(neighborhood_queue)

    if neighborhood.epsilon > neighborhoods.epsilon:
        return False

    if isinstance(neighborhood, GeodesicNeighborhood):
        if not neighborhood.added:
            neighborhoods.add_neighborhood(neighborhood)
            neighborhood.added = True
            add_geodesic_to_queue(
                neighborhood_queue, neighborhoods, next(len_spec))

    next_tile : Tile = neighborhood.get_next_tile()
    neighborhood.update_radius_and_epsilon(next_tile)
    tet_index = next_tile.lifted_tetrahedron.tet.Index

    for other_neighborhood in neighborhoods.neighborhoods:
        neighborhood_pair = neighborhoods.get_neighborhood_pair(
            neighborhood, other_neighborhood)
        if neighborhood_pair.finished:
            continue

        for other_tile in other_neighborhood.tet_to_tiles[tet_index]:
            neighborhood_pair.distance_lifts = correct_min([
                neighborhood_pair.distance_lifts,
                distance_tiles(next_tile, other_tile)])

        total_radius = neighborhood.radius + other_neighborhood.radius

        if total_radius < neighborhood_pair.distance_lifts:
            continue

        neighborhood_pair.epsilon = epsilon_from_neighborhood_pair(
            neighborhood, other_neighborhood,
            neighborhood_pair.distance_lifts,
            verified=neighborhoods.mcomplex.verified)
        neighborhoods.epsilon = correct_min(
            [neighborhoods.epsilon, neighborhood_pair.epsilon])

        if neighborhood_pair.distance_lifts < total_radius:
            neighborhood_pair.finished = True

    neighborhood.add_next_tile(next_tile)
    heapq.heappush(neighborhood_queue, neighborhood)

    return True

def margulis_number(M, bits_prec=None, verified=False, include_thin_part=False):
    mcomplex = mcomplex_for_margulis_number(
        M, bits_prec=bits_prec, verified=verified)

    neighborhoods = Neighborhoods(mcomplex)
    neighborhood_queue : list[GeometricObjectNeighborhood] = []

    cusp_shapes = compute_cusp_shapes(
        M, bits_prec = bits_prec, verified = verified)

    for vertex, cusp_shape in zip(mcomplex.Vertices, cusp_shapes):
        if vertex.is_complete:
            add_cusp_to_queue_and_neighborhoods(
                neighborhood_queue, neighborhoods, vertex, cusp_shape)

    len_spec = M.length_spectrum_alt_gen(bits_prec=bits_prec, verified=verified)
    add_geodesic_to_queue(neighborhood_queue, neighborhoods, next(len_spec))

    while expand_next_neighborhood(
            neighborhood_queue, neighborhoods, len_spec):
        pass

    if include_thin_part:
        return (neighborhoods.epsilon,
                neighborhoods.thin_part(),
                neighborhoods.collisions())
    else:
        return neighborhoods.epsilon
