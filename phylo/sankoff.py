"""Sankoff's algorithm (1975): generalized parsimony with arbitrary substitution
costs. Dynamic programming bottom-up pass:

  l_c(v) = sum over children w of   min_{c'} [ l_{c'}(w) + cost(c', c) ]

The minimum entry of the root's cost map is the parsimony tree length for the column.
Sankoff supports arbitrary cost matrices — e.g., higher cost for transversions than
transitions — and trees with polytomies, both of which Fitch cannot handle.
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import Node


def sankoff_column(root: Node, leaf_states: dict[str, str],
                   states: list[str], cost: dict[tuple[str, str], float]) -> float:
    """Run Sankoff on one column. Returns minimum cost (= tree length for this column)."""
    INF = float("inf")
    for node in root.postorder():
        if node.is_leaf():
            s = leaf_states[node.name]
            node.sankoff_costs = {c: 0.0 if c == s else INF for c in states}
        else:
            costs = {}
            for c in states:
                total = 0.0
                for child in node.children:
                    best = min(child.sankoff_costs[cp] + cost[(cp, c)] for cp in states)
                    total += best
                costs[c] = total
            node.sankoff_costs = costs
    return min(root.sankoff_costs.values())


def jc_cost_matrix(states: list[str]) -> dict[tuple[str, str], float]:
    """Equal cost — recovers Fitch behaviour."""
    return {(a, b): (0.0 if a == b else 1.0) for a in states for b in states}


def kimura_cost_matrix(transition_cost: float = 1.0,
                       transversion_cost: float = 2.5) -> dict[tuple[str, str], float]:
    """Higher cost for transversions (A↔C, A↔T, G↔C, G↔T) than transitions (A↔G, C↔T)."""
    purines = {"A", "G"}
    pyrimidines = {"C", "T"}
    states = ["A", "C", "G", "T"]
    cost = {}
    for a in states:
        for b in states:
            if a == b:
                cost[(a, b)] = 0.0
            elif (a in purines) == (b in purines):
                cost[(a, b)] = transition_cost
            else:
                cost[(a, b)] = transversion_cost
    return cost


def sankoff_alignment(root: Node, seqs: dict[str, str], cost_matrix=None) -> tuple[float, list[float]]:
    states = ["A", "C", "G", "T"]
    if cost_matrix is None:
        cost_matrix = jc_cost_matrix(states)
    L = len(next(iter(seqs.values())))
    per_col = []
    for j in range(L):
        col = {name: s[j] for name, s in seqs.items()}
        per_col.append(sankoff_column(root, col, states, cost_matrix))
    return sum(per_col), per_col


# -------------------- Visualization --------------------

def visualize_sankoff_column(root: Node, leaf_states: dict[str, str], outfile: str,
                              cost_matrix=None,
                              title: str = "Sankoff algorithm — DP cost vectors"):
    import matplotlib.pyplot as plt
    states = ["A", "C", "G", "T"]
    if cost_matrix is None:
        cost_matrix = jc_cost_matrix(states)
    tree_len = sankoff_column(root, leaf_states, states, cost_matrix)

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

    fig, ax = plt.subplots(figsize=(12, 7.5))

    def fmt_costs(costs):
        return "  ".join(
            (f"{s}:∞" if v == float("inf") else f"{s}:{v:g}")
            for s, v in zip(states, [costs[s] for s in states])
        )

    def draw(node, d):
        x = x_of[node.name]
        y = max_depth - d
        if node.is_leaf():
            ax.text(x, y - 0.25, node.name, ha="center", fontsize=11, fontweight="bold")
            ax.text(x, y - 0.55, leaf_states[node.name], ha="center", fontsize=12,
                    color="darkgreen", fontweight="bold")
        else:
            ax.text(x, y + 0.20, fmt_costs(node.sankoff_costs), ha="center", fontsize=8.5,
                    family="monospace",
                    bbox=dict(boxstyle="round,pad=0.20", facecolor="lightyellow",
                              edgecolor="navy"))
            best = min(node.sankoff_costs.values())
            best_states = [s for s in states if node.sankoff_costs[s] == best]
            ax.text(x, y - 0.10, f"min={best:g} → {','.join(best_states)}",
                    ha="center", fontsize=9, color="darkred")
        for c in node.children:
            cx = x_of[c.name]
            cy = max_depth - (d + 1)
            ax.plot([x, cx], [y, cy], "k-", lw=1.2)
            draw(c, d + 1)

    draw(root, 0)
    ax.set_title(f"{title}\ntree length (column) = {tree_len:g}", fontsize=12)
    ax.set_xlim(-0.8, len(leaves) - 0.2)
    ax.set_ylim(-1.0, max_depth + 0.7)
    ax.axis("off")
    plt.tight_layout()
    plt.savefig(outfile, dpi=140, bbox_inches="tight")
    plt.close(fig)
    return tree_len
