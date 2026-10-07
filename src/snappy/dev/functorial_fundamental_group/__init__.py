"""Functorial simplification of triangulations and fundamental groups."""

from collections import deque

from ...snap.t3mlite.simplex import FaceIndex, TwoSubsimplices
from .mcomplex import McomplexWithFundamentalGroup

__all__ = ['McomplexWithFundamentalGroup', 'simplify']


class IsomorphismImages(list):
    """A list of generator images with its domain and codomain attached."""

    def __init__(self, images, domain, codomain):
        super().__init__(images)
        self.domain = domain
        self.codomain = codomain


def _adjacency_from_info(info):
    return {
        (tet, TwoSubsimplices[face_index]):
            (entry['neighbors'][face_index],
             TwoSubsimplices[entry['gluings'][face_index][face_index]])
        for tet, entry in enumerate(info) for face_index in range(4)
    }


def _inverse_index_path(path, adjacency):
    return [(adjacency[step][0], adjacency[step][1])
            for step in reversed(path)]


def _generator_loops(info):
    """Give based dual-edge loops representing SnapPea's generators."""
    adjacency = _adjacency_from_info(info)
    roots = [entry['index'] for entry in info
             if entry['generator_path'] == -1]
    if len(roots) != 1:
        raise ValueError("SnapPea did not provide a unique generator base tet")
    root = roots[0]

    paths = {root: []}
    queue = deque([root])
    while queue:
        tet = queue.popleft()
        for face_index, face in enumerate(TwoSubsimplices):
            if info[tet]['generators'][face_index] != 0:
                continue
            neighbor, _ = adjacency[tet, face]
            if neighbor not in paths:
                paths[neighbor] = paths[tet] + [(tet, face)]
                queue.append(neighbor)
    if len(paths) != len(info):
        raise ValueError("The non-generator face pairings do not form a tree")

    count = max(abs(generator)
                for entry in info for generator in entry['generators'])
    loops = [None] * count
    for tet, entry in enumerate(info):
        for face_index, generator in enumerate(entry['generators']):
            if generator > 0:
                face = TwoSubsimplices[face_index]
                neighbor, _ = adjacency[tet, face]
                loops[generator - 1] = (
                    paths[tet] + [(tet, face)] +
                    _inverse_index_path(paths[neighbor], adjacency))
    if any(loop is None for loop in loops):
        raise ValueError("SnapPea generator data is incomplete")
    return loops


def _word(path, info, generator_names):
    generators = []
    for tet, face in path:
        generator = info[tet]['generators'][FaceIndex[face]]
        if generator:
            index = abs(generator) - 1
            signed_generator = index + 1 if generator > 0 else -index - 1
            if generators and generators[-1] == -signed_generator:
                generators.pop()
            else:
                generators.append(signed_generator)

    def render(signed_generator):
        name = generator_names[abs(signed_generator) - 1]
        return name if signed_generator > 0 else name.swapcase()

    return ''.join(render(generator) for generator in generators)


def _indexed_path(path, indices):
    return [(indices[tet], face) for tet, face in path]


def simplify(triangulation, jiggle_limit=None):
    """Return a simplified copy and the two induced isomorphisms.

    The return value is ``(simplified, forward, backward)``.  ``forward``
    maps generators of the input presentation to words in the simplified
    presentation; ``backward`` maps in the opposite direction.  The maps are
    list subclasses carrying ``domain`` and ``codomain`` attributes.

    The input triangulation is not mutated.
    """
    original = triangulation.copy().with_hyperbolic_structure()
    original._choose_generators(False, False)
    original_info = original._choose_generators_info()
    original_generator_names = original.fundamental_group(False).generators()

    mcomplex = McomplexWithFundamentalGroup(triangulation)
    mcomplex.simplify(jiggle_limit=jiggle_limit)
    mcomplex.rebuild()
    simplified = mcomplex.snappy_manifold()
    simplified._choose_generators(False, False)
    simplified_info = simplified._choose_generators_info()
    simplified_generator_names = (
        simplified.fundamental_group(False).generators())

    initial_indices = {
        tet: index for index, tet in enumerate(mcomplex._initial_tetrahedra)}
    current_indices = {
        tet: index for index, tet in enumerate(mcomplex.Tetrahedra)}

    forward_paths = mcomplex.initial_to_current['arrows']
    backward_paths = mcomplex.current_to_initial['arrows']

    forward_images = []
    for loop in _generator_loops(original_info):
        transported = [step
                       for tet, face in loop
                       for step in forward_paths[
                           (mcomplex._initial_tetrahedra[tet], face)]]
        forward_images.append(_word(
            _indexed_path(transported, current_indices), simplified_info,
            simplified_generator_names))

    backward_images = []
    for loop in _generator_loops(simplified_info):
        transported = [step
                       for tet, face in loop
                       for step in backward_paths[(mcomplex.Tetrahedra[tet], face)]]
        backward_images.append(_word(
            _indexed_path(transported, initial_indices), original_info,
            original_generator_names))

    forward = IsomorphismImages(forward_images, original, simplified)
    backward = IsomorphismImages(backward_images, simplified, original)
    return simplified, forward, backward
