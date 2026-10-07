"""Regression tests for functorial fundamental-group simplification."""

import snappy

from . import _word, simplify


def test_words_are_freely_reduced():
    info = [
        {'generators': (1, -1, 2, -2)},
    ]
    path = [(0, 14), (0, 13), (0, 11), (0, 13),
            (0, 14), (0, 7)]
    assert _word(path, info, ['x', 'y']) == ''

    path = [(0, 11), (0, 14), (0, 7)]
    assert _word(path, info, ['x', 'y']) == 'yxY'


def test_self_identified_pachner_boundary():
    """A move occupying all of m004 must not shortcut across its boundary."""
    manifold = snappy.Manifold('m004')
    manifold._two_to_three(0, 0)
    before = manifold._to_string()

    simplified, forward, backward = simplify(manifold)

    assert manifold._to_string() == before
    original_group = manifold.fundamental_group(False)
    simplified_group = simplified.fundamental_group(False)

    for index, image in enumerate(forward):
        generator = chr(ord('a') + index)
        assert abs(original_group.complex_length(generator) -
                   simplified_group.complex_length(image)) < 1e-6

    for index, image in enumerate(backward):
        generator = chr(ord('a') + index)
        assert abs(simplified_group.complex_length(generator) -
                   original_group.complex_length(image)) < 1e-6
