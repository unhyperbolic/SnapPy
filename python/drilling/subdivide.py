from .moves import one_four_move, two_three_move
from .tracing import GeodesicPiece, GeodesicPieceTracker

from . import debug

from ..snap.t3mlite import Mcomplex, Tetrahedron

from typing import Sequence, Dict

def traverse_geodesics_to_subdivide(
        mcomplex : Mcomplex,
        all_pieces : Sequence[Sequence[GeodesicPiece]]) -> Sequence[Tetrahedron]:

    """
    The function expects a triangulation and a list of list of geodesic pieces
    where each list was generated by calling trace_geodesic on the same
    triangulation. That is, the pieces of each list form a simple closed curve in
    the manifold that is geodesic except for one point. See trace_geodesic for
    further details on the input. Note that the function corrupts its input
    data.

    The function returns a list of tetrahedra that form a triangulation that is
    isometric to the given triangulation and such that the above simple closed
    curves all embed into the 1-skeleton. The triangulation will have ideal and
    finite vertices (that will be assigned light-like and (unnormalized)
    time-like vectors, respectively by tet.R13_vertices).

    An edge E of the new triangulation is part of one of the above simple
    closed curves if there is a pair of a tetrahedron tet and an edge e in that
    tetrahedron representing E and there is a GeodesicPiece in
    tet.geodesic_pieces with endpoints being the endpoints of e.
    """

    # We perform _traverse_geodesic_to_subdivide for each of the given simple
    # closed curves. A call to _traverse_geodesic_to_subdivide will use
    # 1-4 and 2-3 move to introduce finite vertices and edges connecting them
    # so that the given simple curve can be embedded into the 1-skeleton.
    #
    # While the details of this process are described in
    # _traverse_geodesic_to_subdivide, we describe here the nature of the
    # data structure to represent the intermediate triangulations with
    # line segments embedded in tetrahedra (that is GeodesicPiece's).
    #
    # We actually do not explicitly store the set or list of tetrahedra
    # and line segments making up the intermediate triangulations and
    # simple closed curves, respectively here. Instead, we just keep pointers to
    # some GeodesicPiece's and rely on the fact that the triangulation
    # and each simple closed curve is connected and we can find all
    # tetrahedra and pieces by traversal.
    #
    # The traversal uses that:
    # Each tetrahedron tet points to its four neighbors through tet.Neighbors.
    # Each tetrahedron tet points to all line segments embedded in it
    # through tet.geodesic_pieces.
    # Each GeodesicPiece points back to the tetrahedron where the line
    # segment it represents is embdedded in.
    # Each GeodesicPiece also has a pointer to the previous and next
    # GeodesicPiece in the simple closed curves. In other words, we have
    # a cyclically linked list of GeodesicPiece's.
    #
    # Note that the 1-4 and 2-3 moves will replace tetrahedra and
    # GeodesicPiece's with other tetrahedra and GeodesicPiece's.
    #
    # Thus, in general, it is not safe to keep a pointer to a tetrahedron
    # or GeodesicPiece around.
    #
    # The only exception are the "unprocessed" start pieces of the
    # simple closed curves. That is the first piece of a simple closed
    # curve before the simple closed curve became part of the 1-skeleton
    # through _traverse_geodesic_to_subdivide. Recall from trace_geodesic
    # that such a start piece starts in the interior of a tetrahedron and
    # ends at one of the faces of the tetrahedron. The code in the 1-4 and
    # 2-3 moves has special logic to detect pieces of this nature. For
    # such a piece, it actually modifies the GeodesicPiece in place
    # rather than replacing it with one or several other GeodesicPiece's
    # like it does for non-start GeodesicPiece's.
    #
    # Thus, we can store all start pieces here and iterate through them
    # (a start piece will be valid until we have given it to
    # _traverse_geodesic_to_subdivide but not after).

    # Add pointers to the geodesic pieces to the respective tetrahedra.
    for tet in mcomplex.Tetrahedra:
        tet.geodesic_pieces = []
    for pieces in all_pieces:
        for piece in pieces:
            piece.tet.geodesic_pieces.append(piece)

    # Store all start pieces.
    start_pieces = [ pieces[0] for pieces in all_pieces ]

    for start_piece in start_pieces:
        debug.check_consistency_segments(debug.flatten_link_list(start_piece))

    trackers = [ GeodesicPieceTracker(start_piece)
                 for start_piece in start_pieces ]

    # Iterate through start pieces.
    for tracker in trackers:
        # Make 1-skeleton contain the simple closed curve starting with
        # the start piece.
        last_piece : GeodesicPiece = _traverse_geodesic_to_subdivide(
            tracker.geodesic_piece, mcomplex.verified)

    # At this point, all start pieces have been processed and all elements
    # of start_pieces are invalid. Luckily, _traverse_geodesic_to_subdivide
    # gives us a valid GeodesicPiece we can use for traversal.

    return _find_and_index_all_tetrahedra(last_piece.tet)

def _traverse_geodesic_to_subdivide(
        start_piece : GeodesicPiece,
        verified : bool) -> GeodesicPiece:

    debug.check_consistency_2(start_piece)

    # We introduce the following notation for the current state
    # of the GeodesicPiece's that are cyclically linked to form
    # the simple closed curve that we are currently processing:
    #
    # end_piece
    #    v
    # F-F-T-F-...-F-
    #      ^
    #
    # X-Y denotes a GeodesicPiece where X and Y are V, F, or T to
    # indicate whether the start or end point, respectively, is
    # a vertex, on a face or in the interior of a tetrahedron.
    #
    # Since we have V = X for two consecutive pieces U-V X-Y, we
    # abbreviate to U-X-Y.
    #
    # Furthermore, since drawing circles in ASCII art is difficult,
    # X- indicates that we wrap back to the first letter in a line.
    #
    # "v" over a piece and "^" under a piece means that is currently
    # pointed to by a variable we are interested in where the name
    # the variable if it is not just "piece".


    # Following trace_geodesic, we start with
    #
    #    start_piece.prev
    #           v
    #        F-F-T-F-F-...-F-
    #             ^
    #         start_piece

    end_piece, piece = one_four_move(
        [start_piece],
        verified)

    # The first 1-4 move creates a vertex for the point where the
    # simple closed curve starts and ends:
    #
    #       end_piece
    #           v
    #        F-F-V-F-F-...-F-
    #             ^

    debug.check_consistency_2(piece)

    while True:
        # When entering the curve, we have
        #
        #
        #       V-...-V-F-...-F-
        #              ^
        # Note that there is only a single V if the curve
        # runs for the first time. And there might be only
        # a single F if we are about to finish.

        # Proceed to the next piece.
        piece = piece.next_

        # We now have a different picture depending on whether
        # there was only a single F.

        if piece.is_face_to_vertex():

            # We have
            #
            #       prev.piece
            #            v
            #     V-...-V-F-V-...-V-
            #              ^
            #
            # We can now achieve the goal by a single 2-3 move:
            # prev.piece and piece will turn into a single new piece
            # that is the edge created by the 2-3 move - the face
            # common to the tetrahedra containing prev.piece and piece
            # has disappeared.

            piece = two_three_move([piece.prev, piece], verified)
            debug.check_consistency_2(piece)

            # The entire simple closed curve consists of edges of
            # tetrahedra.
            #
            #      V-...-V-V-...-V-
            #             ^
            #
            # We are done.

            return piece

        # We have
        #
        #      V-...-V-F-F-...-F-
        #               ^
        #
        # We pick a point on the geodesic piece and use it as vertex
        # for a 1-4 move.

        piece, next_piece = one_four_move([piece], verified)

        debug.check_consistency_2(piece)

        # We have
        #
        #    piece.prev   next_piece
        #             v   v
        #      V-...-V-F-V-F-...-F-
        #               ^
        #
        # Remark: The current piece might not be the only geodesic
        # piece (of this simple closed curve) going through the tetrahedron.
        # The 1-4 move will split these other pieces into several piece.
        # But the pieces will be again F-F. And the process still finishes
        # in finite time. Similar for the 2-3 move.
        #
        # We perform a 2-3 move, see also above description.

        piece = two_three_move([piece.prev, piece], verified)
        debug.check_consistency_2(piece)

        # We have
        #
        #      V-...-V-V-F-...-F-
        #             ^
        #
        # Move to the next piece.

        piece = piece.next_

        # We have
        #
        #      V-...-V-V-F-...-F-
        #               ^
        #
        # Start next iteration.


def _find_and_index_all_tetrahedra(tet : Tetrahedron):
    """
    Recursively traverses neighbors of the given Tetrahedron
    to find all tetrahedra tet in the connected component.

    Assigns tet.Index to them.
    """
    result = []
    pending_tets = [tet]
    visited_tets = set()
    i = 0
    while pending_tets:
        tet = pending_tets.pop()
        if tet not in visited_tets:
            visited_tets.add(tet)
            tet.Index = i
            i += 1
            result.append(tet)
            for neighbor in tet.Neighbor.values():
                pending_tets.append(neighbor)

    return result
