import xml.etree.ElementTree as ET

import pytest

from bacon import BaconError
from bacon.newick import ladderize, midpoint_root, parse, to_newick, to_svg


def leaves(node):
    return sorted(leaf.name for leaf in node.leaves())


def path_length(root, a, b):
    """Sum of branch lengths between two leaves."""
    def ancestors(name):
        out, stack = [], [(root, [])]
        while stack:
            node, path = stack.pop()
            if node.name == name and node.is_leaf():
                return path + [node]
            for c in node.children:
                stack.append((c, path + [node]))
        return out
    pa, pb = ancestors(a), ancestors(b)
    common = 0
    while common < min(len(pa), len(pb)) and pa[common] is pb[common]:
        common += 1
    return sum(n.length for n in pa[common:]) + sum(n.length for n in pb[common:])


def test_parse_lengths_supports_and_quotes():
    t = parse("((A:0.1,'B c':0.2)0.95:0.3,D:0.4);")
    assert leaves(t) == ["A", "B c", "D"]
    assert t.children[0].name == "0.95"
    assert t.children[0].length == pytest.approx(0.3)


@pytest.mark.parametrize("bad", ["(A,B)", "(A,B;", "(A,B));", "(A:x,B);"])
def test_parse_errors(bad):
    with pytest.raises(BaconError):
        parse(bad)


def test_roundtrip():
    text = "((A:0.1,B:0.2)0.9:0.3,C:0.4,D:0.5);"
    assert to_newick(parse(text)) == text


def test_midpoint_root_preserves_leaves_and_distances():
    text = "(A:1,B:2,(C:3,(D:4,E:10)0.8:1)0.7:1);"
    before = parse(text)
    rooted = midpoint_root(parse(text))
    assert leaves(rooted) == leaves(before)
    assert len(rooted.children) == 2
    for a, b in [("A", "E"), ("B", "D"), ("C", "E"), ("A", "B")]:
        assert path_length(rooted, a, b) == pytest.approx(path_length(before, a, b))
    # The longest path (E to B, C or D) is 14: the root sits 7 from E, on E's edge.
    e = next(ch for ch in rooted.children if ch.is_leaf())
    assert e.name == "E" and e.length == pytest.approx(7)


def test_midpoint_root_on_edge_halfway():
    rooted = midpoint_root(parse("(A:1,B:1,C:8);"))
    # Longest path is A-C or B-C = 9; midpoint 4.5 from C, on C's edge.
    c = next(ch for ch in rooted.children if ch.is_leaf() and ch.name == "C")
    assert c.length == pytest.approx(4.5)
    other = next(ch for ch in rooted.children if ch is not c)
    assert other.length == pytest.approx(3.5)


def test_midpoint_keeps_support_on_the_right_edge():
    rooted = midpoint_root(parse("((A:1,B:1)0.99:1,C:1,(D:1,E:20)0.5:1);"))
    labelled = {n.name: sorted(x.name for x in n.leaves()) for n in _internal(rooted) if n.name}
    # The clade (A,B) keeps its support whatever the rooting.
    assert labelled.get("0.99") == ["A", "B"]


def _internal(node):
    out = []
    for c in node.children:
        if c.children:
            out.append(c)
            out.extend(_internal(c))
    return out


def test_ladderize_and_svg_is_valid_xml():
    t = midpoint_root(parse("((A:0.1,B:0.2)0.95:0.3,(C:0.1,(D:0.1,E&F:0.2)1:0.1)0.9:0.2);"))
    ladderize(t)
    svg = to_svg(t, "title <x>")
    root = ET.fromstring(svg)
    texts = [el.text for el in root.iter("{http://www.w3.org/2000/svg}text")]
    assert "E&F" in texts and "title <x>" in texts
