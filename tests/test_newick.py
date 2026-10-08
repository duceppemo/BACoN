import xml.etree.ElementTree as ET

import pytest

from bacon import BaconError
from bacon.newick import depths, ladderize, midpoint_root, parse, to_newick, to_svg


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


def test_midpoint_keeps_the_support_of_the_split_edge():
    rooted = to_newick(midpoint_root(parse("((A:1,B:1)0.99:10,C:1,D:1);")))
    assert rooted == "((A:1,B:1)0.99:5,(C:1,D:1):5);"


def test_unterminated_quoted_label():
    import pytest

    from bacon import BaconError
    with pytest.raises(BaconError):
        parse("('A:1,B:1);")


def test_deep_tree_walks_are_iterative():
    # A caterpillar of 1,000 leaves is 999 levels deep: recursive walks would pass Python's recursion limit
    text = "".join("(" for _ in range(999)) + "L0:1" + "".join(f",L{i}:1):0.5" for i in range(1, 999)) + ",L999:1);"
    tree = parse(text)
    assert len(tree.leaves()) == 1000 and to_newick(tree) == text
    ladderize(tree)
    assert [n.name for n in tree.children] == ["L999", ""]  # The leaf first: fewer leaves
    rooted = midpoint_root(tree)
    assert leaves(rooted) == sorted(f"L{i}" for i in range(1000))
    assert to_svg(rooted).count("<text") >= 1000


def _caterpillar(prefix, m):
    """A caterpillar of `m` leaves whose internal nodes are all unnamed with length 0 (as FastTree writes the tree
    of identical sequences)."""
    return "(" * (m - 1) + f"{prefix}0:5" + "".join(f",{prefix}{i}:0.001):0" for i in range(1, m - 1)) \
        + f",{prefix}{m - 1}:0.001)"


def test_deep_tree_of_identical_internal_nodes():
    # Two 3,000-leaf caterpillars whose internal nodes have the same name and length: nodes compared by value
    # would recurse down the tree (RecursionError), and finding the path between the farthest leaves be quadratic
    text = f"({_caterpillar('A', 3000)}:0,{_caterpillar('B', 3000)}:0);"
    rooted = midpoint_root(parse(text))
    assert len(rooted.leaves()) == 6000 and len(rooted.children) == 2
    assert path_length(rooted, "A0", "B0") == pytest.approx(10)  # The longest path: the root halfway along it
    depth = depths(rooted)
    assert [depth[id(leaf)] for leaf in rooted.leaves() if leaf.name in ("A0", "B0")] == pytest.approx([5, 5])
    ladderize(rooted)
    assert to_svg(rooted).count("<text") >= 6000


def test_nodes_are_equal_only_to_themselves():
    a, b = parse("(A:1,B:1);"), parse("(A:1,B:1);")
    assert a != b and a == a and a.children[0] in a.children and a.children[0] not in b.children


@pytest.mark.parametrize("bad", ["(A:inf,B:1);", "(A:1,B:nan);", "(A:1,B:1e999);", "((A:1,B:1):-inf,C:1);",
                                 "((A:1,B:1)inf:1,C:1);", "((A:1,B:1)1e999:1,C:1);", "(A:1,B:1)nan;"])
def test_parse_rejects_non_finite_lengths_and_supports(bad):
    with pytest.raises(BaconError, match="not a finite number"):
        parse(bad)


def test_parse_keeps_named_internal_nodes():
    assert parse("((A:1,B:1)cladeX:1,C:1);").children[0].name == "cladeX"
