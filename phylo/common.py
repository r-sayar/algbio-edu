"""Shared utilities: simple Tree node, Newick I/O, FASTA reader."""
from __future__ import annotations
from dataclasses import dataclass, field
from typing import Optional


@dataclass
class Node:
    """Phylogenetic tree node. children=[] for leaves."""
    name: str = ""
    branch_length: float = 0.0
    children: list = field(default_factory=list)
    parent: Optional["Node"] = None
    # Per-algorithm scratch fields:
    fitch_set: set = field(default_factory=set)        # Fitch up-pass character set
    fitch_state: str = ""                              # Fitch down-pass assigned state
    sankoff_costs: dict = field(default_factory=dict)  # Sankoff: state -> min cost
    likelihood: dict = field(default_factory=dict)     # Felsenstein: state -> conditional L

    def is_leaf(self) -> bool:
        return not self.children

    def leaves(self):
        if self.is_leaf():
            yield self
        else:
            for c in self.children:
                yield from c.leaves()

    def postorder(self):
        for c in self.children:
            yield from c.postorder()
        yield self

    def preorder(self):
        yield self
        for c in self.children:
            yield from c.preorder()


def read_fasta(path: str) -> dict[str, str]:
    seqs = {}
    name = None
    buf: list[str] = []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if name is not None:
                    seqs[name] = "".join(buf)
                name = line[1:].split()[0]
                buf = []
            elif line:
                buf.append(line)
        if name is not None:
            seqs[name] = "".join(buf)
    return seqs


def parse_newick(s: str) -> Node:
    """Parse a Newick string into a Node tree."""
    s = s.strip().rstrip(";")
    pos = [0]

    def parse() -> Node:
        node = Node()
        if pos[0] < len(s) and s[pos[0]] == "(":
            pos[0] += 1
            node.children.append(parse())
            while pos[0] < len(s) and s[pos[0]] == ",":
                pos[0] += 1
                node.children.append(parse())
            assert s[pos[0]] == ")"
            pos[0] += 1
        # Read name
        name = ""
        while pos[0] < len(s) and s[pos[0]] not in ":,)":
            name += s[pos[0]]; pos[0] += 1
        node.name = name
        # Read branch length
        if pos[0] < len(s) and s[pos[0]] == ":":
            pos[0] += 1
            num = ""
            while pos[0] < len(s) and s[pos[0]] not in ",)":
                num += s[pos[0]]; pos[0] += 1
            node.branch_length = float(num)
        for c in node.children:
            c.parent = node
        return node

    return parse()


def to_newick(node: Node, with_branch_lengths: bool = True) -> str:
    if node.is_leaf():
        s = node.name
    else:
        s = "(" + ",".join(to_newick(c, with_branch_lengths) for c in node.children) + ")"
        if node.name:
            s += node.name
    if with_branch_lengths and node.parent is not None:
        s += f":{node.branch_length:g}"
    return s


def name_internal_nodes(root: Node, prefix: str = "N") -> None:
    """Give each unnamed internal node a unique label N1, N2, …"""
    counter = [0]
    for n in root.preorder():
        if not n.is_leaf() and not n.name:
            counter[0] += 1
            n.name = f"{prefix}{counter[0]}"


def make_simple_tree(topology: str, seqs: dict[str, str]) -> Node:
    """Build a Node tree from a Newick topology and attach leaf sequences."""
    root = parse_newick(topology)
    name_internal_nodes(root)
    for leaf in root.leaves():
        if leaf.name in seqs:
            leaf.sequence = seqs[leaf.name]
    return root


def hamming_distance_matrix(seqs: dict[str, str]):
    """Return (names, n×n matrix of proportional Hamming distances)."""
    import numpy as np
    names = list(seqs)
    n = len(names)
    D = np.zeros((n, n))
    for i in range(n):
        for j in range(i+1, n):
            s1, s2 = seqs[names[i]], seqs[names[j]]
            d = sum(a != b for a, b in zip(s1, s2)) / len(s1)
            D[i, j] = D[j, i] = d
    return names, D
