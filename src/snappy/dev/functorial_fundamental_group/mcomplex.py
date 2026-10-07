"""Pachner moves which remember the induced map on the dual groupoid."""

from collections import deque

from ...exterior_to_link.mcomplex_with_expansion import McomplexWithExpansion
from ...snap.t3mlite.arrow import Arrow
from ...snap.t3mlite.mcomplex import edge_and_arrow
from ...snap.t3mlite.simplex import PickAnEdge, TwoSubsimplices


def _other_face(tet, face):
    return tet.Gluing[face].image(face)


def _inverse_path(path, adjacency):
    return [(adjacency[tet, face][0], adjacency[tet, face][1])
            for tet, face in reversed(path)]


def _path_inside(start, end, region, adjacency, boundary_faces=frozenset()):
    """Find a path through interior faces of an abstract Pachner ball."""
    if start is end:
        return []
    queue = deque([start])
    previous = {start: None}
    while queue:
        tet = queue.popleft()
        for face in TwoSubsimplices:
            if (tet, face) in boundary_faces:
                continue
            neighbor, _ = adjacency[tet, face]
            if neighbor in region and neighbor not in previous:
                previous[neighbor] = (tet, face)
                if neighbor is end:
                    ans = []
                    current = end
                    while previous[current] is not None:
                        source, source_face = previous[current]
                        ans.append((source, source_face))
                        current = source
                    return list(reversed(ans))
                queue.append(neighbor)
    raise ValueError("The tetrahedra in a Pachner ball are disconnected")


class McomplexWithFundamentalGroup(McomplexWithExpansion):
    """An ``McomplexWithExpansion`` carrying dual-groupoid equivalences.

    Paths are lists of ``(tetrahedron, outgoing_face)`` pairs.  The
    attributes :attr:`initial_to_current` and :attr:`current_to_initial`
    encode mutually inverse equivalences.  Their ``objects`` entries map
    objects, and their ``arrows`` entries map each oriented face pairing.

    The 2-0 moves of :class:`McomplexWithExpansion` are already expanded
    into 2-3 and 3-2 moves.  This class additionally transports across the
    4-4 moves used by ``simplify``.
    """

    def __init__(self, tetrahedron_list=None):
        super().__init__(tetrahedron_list)
        self._initial_tetrahedra = tuple(self.Tetrahedra)
        arrows = {(tet, face): [(tet, face)]
                  for tet in self.Tetrahedra for face in TwoSubsimplices}
        self._forward_objects = {tet: tet for tet in self.Tetrahedra}
        self._forward_arrows = dict(arrows)
        self._backward_objects = {tet: tet for tet in self.Tetrahedra}
        self._backward_arrows = dict(arrows)

    @staticmethod
    def _adjacency(tetrahedra):
        return {
            (tet, face): (tet.Neighbor[face], _other_face(tet, face))
            for tet in tetrahedra for face in TwoSubsimplices
        }

    @staticmethod
    def _local_equivalence(source, target, source_adj, target_adj,
                           boundary_map):
        """Return object and arrow maps between two Pachner balls.

        ``boundary_map`` matches occurrences in the boundary of the abstract
        Pachner ball.  It deliberately does not use the tetrahedron on the
        other side: that tetrahedron can itself belong to the move when
        boundary occurrences of the ball are identified in the ambient
        triangulation.
        """
        object_map = {}
        target_boundary = set(boundary_map.values())
        for tet in source:
            candidates = [target_face[0]
                          for source_face, target_face in boundary_map.items()
                          if source_face[0] is tet]
            object_map[tet] = candidates[0] if candidates else next(iter(target))

        arrow_map = {}
        for tet in source:
            for face in TwoSubsimplices:
                neighbor, _ = source_adj[tet, face]
                if (tet, face) not in boundary_map:
                    arrow_map[tet, face] = _path_inside(
                        object_map[tet], object_map[neighbor],
                        target, target_adj, target_boundary)
                else:
                    owner, boundary_face = boundary_map[tet, face]
                    path = (
                        _path_inside(object_map[tet], owner,
                                     target, target_adj, target_boundary) +
                        [(owner, boundary_face)])
                    target_neighbor, _ = target_adj[owner, boundary_face]
                    # If two boundary occurrences are identified, crossing
                    # the target occurrence remains inside the move ball.
                    # Whisker it to the chosen image of the actual neighbor.
                    if neighbor in source:
                        path += _path_inside(
                            target_neighbor, object_map[neighbor],
                            target, target_adj, target_boundary)
                    arrow_map[tet, face] = path
        return object_map, arrow_map

    def _finish_transport(self, before, before_adj, result, boundary_map):
        after = set(self.Tetrahedra)
        removed = before - after
        added = after - before
        after_adj = self._adjacency(after)
        if len(set(boundary_map.values())) != len(boundary_map):
            raise AssertionError(
                "A Pachner move did not match boundary occurrences bijectively")
        for tet in removed:
            for face in TwoSubsimplices:
                if ((tet, face) not in boundary_map and
                        before_adj[tet, face][0] not in removed):
                    raise AssertionError(
                        "A boundary occurrence of a Pachner move was not tracked")
        old_to_new_obj, old_to_new = self._local_equivalence(
            removed, added, before_adj, after_adj, boundary_map)
        new_to_old_obj, new_to_old = self._local_equivalence(
            added, removed, after_adj, before_adj,
            {new: old for old, new in boundary_map.items()})

        def forward_step(step):
            return old_to_new[step] if step[0] in removed else [step]

        self._forward_arrows = {
            arrow: [new_step
                    for step in path for new_step in forward_step(step)]
            for arrow, path in self._forward_arrows.items()
        }
        self._forward_objects = {
            obj: old_to_new_obj.get(image, image)
            for obj, image in self._forward_objects.items()
        }

        old_backward = self._backward_arrows
        self._backward_arrows = {
            (tet, face): [initial_step
                          for old_step in (new_to_old[tet, face]
                                           if tet in added else [(tet, face)])
                          for initial_step in old_backward[old_step]]
            for tet in after for face in TwoSubsimplices
        }
        self._backward_objects = {
            tet: self._backward_objects[new_to_old_obj[tet]]
                 if tet in added else self._backward_objects[tet]
            for tet in after
        }
        return result

    def two_to_three(self, face_or_arrow, tet=None, return_arrow=False,
                     must_succeed=False, unsafe_mode=False):
        if isinstance(face_or_arrow, Arrow):
            assert tet is None
            a = face_or_arrow.copy()
        else:
            a = Arrow(PickAnEdge[face_or_arrow], face_or_arrow, tet)
        b = a.copy().glued()
        if not unsafe_mode:
            possible, reason = self._face_permits_two_to_three(a, b)
            if not possible:
                if must_succeed:
                    raise ValueError(reason)
                return False

        before = set(self.Tetrahedra)
        before_adj = self._adjacency(before)
        a_orig = a.copy()
        new = self.new_arrows(3)
        for i in range(3):
            new[i].glue(new[(i + 1) % 3])
        boundary_map = {}
        a.reverse()
        for c in new:
            c.opposite()
            boundary_map[a.Tetrahedron, a.Face] = (c.Tetrahedron, c.Face)
            c.glue(a.glued())
            c.reverse()
            boundary_map[b.Tetrahedron, b.Face] = (c.Tetrahedron, c.Face)
            c.glue(b.glued())
            a.rotate(-1)
            b.rotate(1)
        for c in new:
            c.reverse().opposite()
        self._two_to_three_move_hook(a_orig, new)
        self.delete_tet(a.Tetrahedron)
        self.delete_tet(b.Tetrahedron)
        if not unsafe_mode:
            self.build_edge_classes()
        result = new[1].north_head().get_arrow() if return_arrow else True
        return self._finish_transport(
            before, before_adj, result, boundary_map)

    def three_to_two(self, edge_or_arrow, return_arrow=False,
                     must_succeed=False, unsafe_mode=False):
        edge, a = edge_and_arrow(edge_or_arrow)
        if not unsafe_mode:
            possible, reason = self._edge_permits_three_to_two(edge)
            if not possible:
                if must_succeed:
                    raise ValueError(reason)
                return False

        before = set(self.Tetrahedra)
        before_adj = self._adjacency(before)
        a_orig = a.copy()
        b = self.new_arrow()
        c = self.new_arrow()
        b.glue(c)
        b_orig = b.copy()
        b.reverse()
        b_to_return = b.copy()
        boundary_map = {}
        for _ in range(3):
            a.opposite()
            boundary_map[a.Tetrahedron, a.Face] = (b.Tetrahedron, b.Face)
            b.glue(a.glued())
            a.reverse()
            boundary_map[a.Tetrahedron, a.Face] = (c.Tetrahedron, c.Face)
            c.glue(a.glued())
            b.rotate(-1)
            c.rotate(1)
            a.reverse().opposite().next()

        self._three_to_two_move_hook(a_orig, (b_orig, b, c))
        if unsafe_mode:
            old_tets = [a_orig.Tetrahedron,
                        a_orig.copy().next().Tetrahedron,
                        a_orig.copy().next().next().Tetrahedron]
        else:
            old_tets = [corner.Tetrahedron for corner in edge.Corners]
        for old_tet in old_tets:
            if old_tet in self.Tetrahedra:
                self.delete_tet(old_tet)
        if not unsafe_mode:
            self.build_edge_classes()
        result = b_to_return if return_arrow else True
        return self._finish_transport(
            before, before_adj, result, boundary_map)

    def four_to_four(self, edge_or_arrow, must_succeed=False,
                     unsafe_mode=False):
        edge, a = edge_and_arrow(edge_or_arrow)
        possible, reason = self._edge_permits_four_to_four(edge)
        if not possible:
            if must_succeed:
                raise ValueError(reason)
            return False

        before = set(self.Tetrahedra)
        before_adj = self._adjacency(before)
        a_orig = a.copy()
        c = self.new_arrows(4)
        c_orig = [x.copy() for x in c]
        for i in range(4):
            c[i].glue(c[(i + 1) % 4])
        boundary_map = {}
        b = a.glued().reverse()

        a.rotate(1)
        boundary_map[a.Tetrahedron, a.Face] = (
            c[0].copy().opposite().Tetrahedron,
            c[0].copy().opposite().Face)
        c[0].opposite().glue(a.glued())
        b.rotate(-1)
        boundary_map[b.Tetrahedron, b.Face] = (
            c[1].copy().opposite().Tetrahedron,
            c[1].copy().opposite().Face)
        c[1].opposite().glue(b.glued())
        b.rotate(-1)
        boundary_map[b.Tetrahedron, b.Face] = (
            c[2].copy().opposite().Tetrahedron,
            c[2].copy().opposite().Face)
        c[2].opposite().glue(b.glued())
        a.rotate(1)
        boundary_map[a.Tetrahedron, a.Face] = (
            c[3].copy().opposite().Tetrahedron,
            c[3].copy().opposite().Face)
        c[3].opposite().glue(a.glued())

        a.rotate(1).reverse().next()
        b.rotate(-1).reverse().next()
        for new_arrow, old_arrow in ((c[0], a), (c[1], b)):
            new_arrow.reverse()
            old_arrow.rotate(-1 if new_arrow is c[0] else 1)
            boundary_map[old_arrow.Tetrahedron, old_arrow.Face] = (
                new_arrow.Tetrahedron, new_arrow.Face)
            new_arrow.glue(old_arrow.glued())
        for new_arrow, old_arrow in ((c[2], b), (c[3], a)):
            new_arrow.reverse()
            old_arrow.rotate(1 if new_arrow is c[2] else -1)
            boundary_map[old_arrow.Tetrahedron, old_arrow.Face] = (
                new_arrow.Tetrahedron, new_arrow.Face)
            new_arrow.glue(old_arrow.glued())

        self._four_to_four_move_hook(a_orig, c_orig)
        for old_tet in {corner.Tetrahedron for corner in edge.Corners}:
            self.delete_tet(old_tet)
        if not unsafe_mode:
            self.build_edge_classes()
        return self._finish_transport(
            before, before_adj, True, boundary_map)

    @property
    def initial_to_current(self):
        return {'objects': self._forward_objects,
                'arrows': self._forward_arrows}

    @property
    def current_to_initial(self):
        return {'objects': self._backward_objects,
                'arrows': self._backward_arrows}
