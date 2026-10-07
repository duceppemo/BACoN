"""Newick trees: parsing, midpoint rooting and drawing as SVG, with the standard library only. Every walk of a tree
is iterative: a deep tree (a caterpillar of thousands of leaves) does not reach Python's recursion limit."""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from xml.sax.saxutils import escape

from bacon import BaconError


@dataclass
class Node:
    name: str = ""
    length: float = 0.0
    children: list[Node] = field(default_factory=list)
    parent: Node | None = field(default=None, repr=False, compare=False)

    def is_leaf(self) -> bool:
        return not self.children

    def leaves(self) -> list[Node]:
        """The leaves, left to right."""
        return [n for n in preorder(self) if not n.children]

    def add(self, child: Node) -> Node:
        child.parent = self
        self.children.append(child)
        return child


def preorder(root: Node) -> list[Node]:
    """The nodes, each before its descendants, children left to right."""
    order, stack = [], [root]
    while stack:
        node = stack.pop()
        order.append(node)
        stack.extend(reversed(node.children))
    return order


def postorder(root: Node) -> list[Node]:
    """The nodes, each after its descendants."""
    return preorder(root)[::-1]


def parse(text: str) -> Node:
    """Parse one Newick tree. Internal node labels (support values) are kept as names."""
    text = "".join(text.split()) if "'" not in text else text.strip()
    if not text.endswith(";"):
        raise BaconError("Not a Newick tree (no final ';')")
    pos = 0

    def label() -> tuple[str, float]:
        nonlocal pos
        if pos < len(text) and text[pos] == "'":
            end = text.index("'", pos + 1)
            name = text[pos + 1: end]
            pos = end + 1
        else:
            start = pos
            while pos < len(text) and text[pos] not in ",():;":
                pos += 1
            name = text[start:pos]
        length = 0.0
        if pos < len(text) and text[pos] == ":":
            start = pos = pos + 1
            while pos < len(text) and text[pos] not in ",();":
                pos += 1
            try:
                length = float(text[start:pos])
            except ValueError:
                raise BaconError(f"Bad branch length in Newick tree: {text[start:pos]!r}") from None
        return name.strip(), length

    def tree() -> Node:
        """The tree starting at `pos`, read with a stack of the nodes whose children are being read."""
        nonlocal pos
        root = node = Node()
        opened: list[Node] = []
        while True:
            while text[pos] == "(":  # Down to the first child
                pos += 1
                opened.append(node)
                node = node.add(Node())
            node.name, node.length = label()
            while opened:  # Up past the closed parentheses, to the next sibling
                if text[pos] == ",":
                    pos += 1
                    node = opened[-1].add(Node())
                    break
                if text[pos] != ")":
                    raise BaconError(f"Unbalanced parentheses in Newick tree at position {pos}")
                pos += 1
                node = opened.pop()
                node.name, node.length = label()
            else:
                return root

    try:
        root = tree()
    except (IndexError, ValueError):  # ValueError: a quoted label without its closing quote
        raise BaconError("Truncated Newick tree") from None
    if text[pos:] != ";":
        raise BaconError(f"Unexpected text after the Newick tree: {text[pos:pos + 20]!r}")
    return root


def to_newick(node: Node) -> str:
    text: dict[int, str] = {}
    for n in postorder(node):
        name = n.name if all(c not in n.name for c in " ,():;'") else "'" + n.name.replace("'", "") + "'"
        inner = "(" + ",".join(text.pop(id(c)) for c in n.children) + ")" if n.children else ""
        text[id(n)] = f"{inner}{name}" + (f":{n.length:.8g}" if n.parent is not None else "")
    return text[id(node)] + ";"


def _distances_from(start: Node) -> dict[int, tuple[float, Node]]:
    """Path length from `start` to every node, walking the tree as an undirected graph."""
    dist = {id(start): (0.0, start)}
    stack = [start]
    while stack:
        node = stack.pop()
        d = dist[id(node)][0]
        neighbours = [(c, c.length) for c in node.children]
        if node.parent is not None:
            neighbours.append((node.parent, node.length))
        for other, length in neighbours:
            if id(other) not in dist:
                dist[id(other)] = (d + length, other)
                stack.append(other)
    return dist


def _path(a: Node, b: Node) -> list[Node]:
    ancestors = []
    n: Node | None = a
    while n is not None:
        ancestors.append(n)
        n = n.parent
    seen = {id(x) for x in ancestors}
    down = []
    n = b
    while id(n) not in seen:
        down.append(n)
        n = n.parent
    return ancestors[: ancestors.index(n) + 1] + down[::-1]


def _reroot_on_edge(child: Node, offset: float) -> Node:
    """A new tree rooted on the edge above `child`, `offset` away from `child`.

    Support values label edges: in the new tree each internal node carries the label of the edge above it.
    """
    top = child
    while top.parent is not None:
        top = top.parent
    adjacency: dict[int, list[tuple[Node, float, str]]] = {}
    original_leaves: set[int] = set()
    stack = [top]
    while stack:
        node = stack.pop()
        if node.is_leaf():
            original_leaves.add(id(node))
        for c in node.children:
            label = c.name if c.children else ""
            adjacency.setdefault(id(node), []).append((c, c.length, label))
            adjacency.setdefault(id(c), []).append((node, c.length, label))
            stack.append(c)

    def build(start: Node, came_from: Node, start_length: float, start_label: str) -> Node:
        """The subtree hanging from `start` away from `came_from`, built from its leaves up."""
        visits = []  # (node, came from, length, label), each node before the nodes beyond it
        stack = [(start, came_from, start_length, start_label)]
        while stack:
            visit = stack.pop()
            visits.append(visit)
            node, previous = visit[0], visit[1]
            if id(node) not in original_leaves:
                stack += [(nb, node, ln, lab) for nb, ln, lab in adjacency[id(node)] if nb is not previous]
        built: dict[int, Node] = {}
        for node, previous, length, label in reversed(visits):
            if id(node) in original_leaves:
                built[id(node)] = Node(node.name, length)
                continue
            kids = [built.pop(id(nb)) for nb, _, _ in adjacency[id(node)] if nb is not previous]
            if len(kids) == 1:  # The former root, now with one child: splice it out
                kids[0].length += length
                built[id(node)] = kids[0]
                continue
            new = Node(label, length)
            for k in kids:
                new.add(k)
            built[id(node)] = new
        return built[id(start)]

    parent = child.parent
    assert parent is not None
    root = Node()
    support = child.name if child.children else ""  # The split edge's support: shown on one of its halves
    root.add(build(child, parent, offset, support))
    root.add(build(parent, child, child.length - offset, ""))
    return root


def midpoint_root(root: Node) -> Node:
    """Root the tree halfway along the longest path between two leaves."""
    leaves = root.leaves()
    if len(leaves) < 3:
        return root
    far_a = max(_distances_from(leaves[0]).values(), key=lambda t: (t[1].is_leaf(), t[0]))[1]
    dist_a = _distances_from(far_a)
    far_b = max((t for t in dist_a.values() if t[1].is_leaf()), key=lambda t: t[0])[1]
    diameter = dist_a[id(far_b)][0]
    if diameter <= 0:
        return root
    half = diameter / 2
    path = _path(far_a, far_b)
    walked = 0.0
    for a, b in zip(path, path[1:]):
        edge_child, length = (a, a.length) if a.parent is b else (b, b.length)
        if walked + length >= half:
            offset_from_a = half - walked
            offset = offset_from_a if edge_child is a else length - offset_from_a
            return _reroot_on_edge(edge_child, offset)
        walked += length
    return root  # pragma: no cover


def ladderize(node: Node) -> None:
    """Sort the children of every node by their number of leaves, fewest first."""
    count: dict[int, int] = {}
    for n in postorder(node):
        count[id(n)] = sum(count[id(c)] for c in n.children) if n.children else 1
        n.children.sort(key=lambda c: count[id(c)])


def to_svg(root: Node, title: str = "") -> str:
    """Draw a rectangular phylogram with leaf names, support values and a scale bar."""
    leaves = root.leaves()
    row, left, top = 22, 36, 40 if title else 16
    label_width = 9 + 7.5 * max(len(leaf.name) for leaf in leaves)
    plot_width = 640.0
    depth = depths(root)
    max_depth = max(depth.values()) or 1.0
    scale = plot_width / max_depth
    y = heights(root, {id(leaf): top + i * row for i, leaf in enumerate(leaves)})
    width = left + plot_width + label_width + 20
    height = top + len(leaves) * row + 40
    parts = [f'<svg xmlns="http://www.w3.org/2000/svg" width="{width:.0f}" height="{height:.0f}" '
             f'viewBox="0 0 {width:.0f} {height:.0f}" font-family="Helvetica, Arial, sans-serif" font-size="13">',
             '<rect width="100%" height="100%" fill="white"/>']
    if title:
        parts.append(f'<text x="{left}" y="22" font-size="15" font-weight="bold">{escape(title)}</text>')
    lines, texts = [], []
    for n, entering in walk(root):
        x0 = left + depth[id(n)] * scale
        if entering and n.parent is not None:  # The branch to the node, before its subtree
            lines.append(f"M{left + depth[id(n.parent)] * scale:.1f},{y[id(n)]:.1f}H{x0:.1f}")
        if not entering:  # After its subtrees: an internal node's support
            if n.name and n.parent is not None:
                texts.append(f'<text x="{x0 - 3:.1f}" y="{y[id(n)] - 4:.1f}" text-anchor="end" font-size="10" '
                             f'fill="#666">{escape(n.name)}</text>')
        elif n.children:
            ys = [y[id(c)] for c in n.children]
            lines.append(f"M{x0:.1f},{min(ys):.1f}V{max(ys):.1f}")
        else:
            texts.append(f'<text x="{x0 + 5:.1f}" y="{y[id(n)] + 4.5:.1f}">{escape(n.name)}</text>')
    parts.append(f'<path d="{" ".join(lines)}" stroke="black" stroke-width="1.5" fill="none" '
                 f'stroke-linecap="square"/>')
    parts.extend(texts)
    bar = _nice_length(max_depth / 5)
    y_bar = top + len(leaves) * row + 12
    parts.append(f'<path d="M{left},{y_bar}H{left + bar * scale:.1f}" stroke="black" stroke-width="1.5"/>')
    parts.append(f'<text x="{left + bar * scale / 2:.1f}" y="{y_bar + 16}" text-anchor="middle" '
                 f'font-size="11">{bar:g} substitutions/site</text>')
    parts.append("</svg>")
    return "\n".join(parts) + "\n"


def depths(root: Node) -> dict[int, float]:
    """The distance of every node (by id) from the root, negative branch lengths counting as 0."""
    depth = {id(root): 0.0}
    for n in preorder(root):
        for c in n.children:
            depth[id(c)] = depth[id(n)] + max(c.length, 0.0)
    return depth


def heights(root: Node, leaf_y: dict[int, float]) -> dict[int, float]:
    """The vertical position of every node (by id): the leaves' as given, an internal node's halfway between its
    first and last children's extremes."""
    y = dict(leaf_y)
    for n in postorder(root):
        if n.children:
            ys = [y[id(c)] for c in n.children]
            y[id(n)] = (min(ys) + max(ys)) / 2
    return y


def walk(root: Node) -> list[tuple[Node, bool]]:
    """The nodes in the order of a recursive drawing: (node, True) on entering it, then its subtrees, then
    (node, False) on leaving it (internal nodes only)."""
    events, stack = [], [(root, True)]
    while stack:
        node, entering = stack.pop()
        events.append((node, entering))
        if entering and node.children:
            stack.append((node, False))
            stack.extend((c, True) for c in reversed(node.children))
    return events


def _nice_length(x: float) -> float:
    if x <= 0:
        return 1.0
    exponent = 10 ** math.floor(math.log10(x))
    for step in (1, 2, 5, 10):
        if step * exponent >= x:
            return step * exponent
    return 10 * exponent  # pragma: no cover
