"""Fitch's algorithm (1971): minimum number of substitutions needed to explain
character states at the leaves on a fixed tree topology.

Two passes:
  1. Bottom-up: at each internal node, take the intersection of children's sets if
     non-empty, else the union (and increment substitution count).
  2. Top-down: assign each internal node the parent's state if it's in its set,
     else any state from its set.

Time complexity: O(n) per alignment column.
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import Node, parse_newick, name_internal_nodes


def fitch_column(root: Node, leaf_states: dict[str, str]) -> int:
    """Run Fitch on one alignment column. Returns minimum number of substitutions.
    Annotates each node with .fitch_set (up-pass) and .fitch_state (down-pass)."""
    n_subs = 0

    # 1. Bottom-up (postorder)
    for node in root.postorder():
        if node.is_leaf():
            node.fitch_set = {leaf_states[node.name]}
        else:
            sets = [c.fitch_set for c in node.children]
            inter = set.intersection(*sets)
            if inter:
                node.fitch_set = inter
            else:
                node.fitch_set = set.union(*sets)
                n_subs += 1

    # 2. Top-down (preorder)
    for node in root.preorder():
        if node.parent is None:
            # Root: pick any element of its set (deterministic for reproducibility)
            node.fitch_state = sorted(node.fitch_set)[0]
        else:
            if node.parent.fitch_state in node.fitch_set:
                node.fitch_state = node.parent.fitch_state
            else:
                node.fitch_state = sorted(node.fitch_set)[0]

    return n_subs


def fitch_alignment(root: Node, seqs: dict[str, str]) -> tuple[int, list[int]]:
    """Run Fitch on every column; return (total tree length, per-column lengths)."""
    L = len(next(iter(seqs.values())))
    per_col = []
    for j in range(L):
        col = {name: s[j] for name, s in seqs.items()}
        per_col.append(fitch_column(root, col))
    return sum(per_col), per_col


# -------------------- Visualization --------------------

def visualize_fitch_column(root: Node, leaf_states: dict[str, str], outfile: str,
                           title: str = "Fitch algorithm — column trace"):
    """Run Fitch on a single column and draw the tree with up-pass sets and
    down-pass assignments shown at each node."""
    import matplotlib.pyplot as plt
    n_subs = fitch_column(root, leaf_states)

    # Layout
    name_internal_nodes(root)
    leaves = list(root.leaves())
    x_of = {l.name: i for i, l in enumerate(leaves)}

    def assign_x(node):
        if node.is_leaf():
            return x_of[node.name]
        xs = [assign_x(c) for c in node.children]
        x_of[node.name] = sum(xs) / len(xs)
        return x_of[node.name]
    assign_x(root)

    def depth(node):
        if node.is_leaf(): return 0
        return 1 + max(depth(c) for c in node.children)
    max_depth = depth(root)

    def y_of_node(node, d=0):
        if node.parent is None:
            return max_depth
        return max_depth - d

    # Recursively draw
    fig, ax = plt.subplots(figsize=(11, 7))

    def draw(node, d):
        x = x_of[node.name]
        y = max_depth - d
        if node.is_leaf():
            ax.text(x, y - 0.25, node.name, ha="center", fontsize=11, fontweight="bold")
            ax.text(x, y - 0.55, f"{leaf_states[node.name]}", ha="center", fontsize=12,
                    color="darkgreen", fontweight="bold")
        else:
            up_set = "{" + ",".join(sorted(node.fitch_set)) + "}"
            ax.text(x, y + 0.18, up_set, ha="center", fontsize=10, color="navy",
                    bbox=dict(boxstyle="round,pad=0.15", facecolor="lightyellow",
                              edgecolor="navy"))
            ax.text(x, y - 0.10, f"→{node.fitch_state}", ha="center", fontsize=11,
                    color="darkred", fontweight="bold")
        for c in node.children:
            cx = x_of[c.name]
            cy = max_depth - (d + 1)
            ax.plot([x, cx], [y, cy], "k-", lw=1.2)
            draw(c, d + 1)

    draw(root, 0)
    ax.set_title(f"{title}\nminimum substitutions = {n_subs}", fontsize=12)
    ax.text(0.02, 0.98,
            "Boxed sets = bottom-up (union if no intersection)\n"
            "Red letters = top-down assignment",
            transform=ax.transAxes, fontsize=9, va="top",
            bbox=dict(facecolor="white", alpha=0.8))
    ax.set_xlim(-0.6, len(leaves) - 0.4)
    ax.set_ylim(-1.0, max_depth + 0.7)
    ax.axis("off")
    plt.tight_layout()
    plt.savefig(outfile, dpi=140, bbox_inches="tight")
    plt.close(fig)
    return n_subs
