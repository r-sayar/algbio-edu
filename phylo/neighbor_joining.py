"""Neighbor Joining (Saitou & Nei, 1987).

Does NOT assume a molecular clock — does assume the distances are (close to) additive.
For each iteration:
  r_i = sum_k d(i,k)
  Q_ij = (n-2) * d(i,j) - r_i - r_j      (equivalent to d_ij - (r_i + r_j)/(n-2))
  Pick (i,j) with min Q.
  d(u, k) = (d(i,k) + d(j,k) - d(i,j)) / 2
  d(u, i) = d(i,j)/2 + (r_i - r_j) / (2*(n-2))    [for visualization]

Time complexity: O(n^3).
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np


def neighbor_joining(names, D, record_steps=True):
    names = list(names); D = D.astype(float).copy()
    sizes = [1] * len(names)
    subtrees = list(names)
    steps = []

    while len(subtrees) > 2:
        n = D.shape[0]
        r = D.sum(axis=1)
        # Q matrix
        Q = (n - 2) * D - r[:, None] - r[None, :]
        np.fill_diagonal(Q, np.inf)
        # Find min off-diagonal
        idx = np.unravel_index(np.argmin(Q), Q.shape)
        i, j = idx
        if i > j: i, j = j, i

        d_ij = D[i, j]
        # Branch lengths for the join (step 4 in lecture)
        if n - 2 > 0:
            bli = d_ij / 2 + (r[i] - r[j]) / (2 * (n - 2))
        else:
            bli = d_ij / 2
        blj = d_ij - bli
        new_subtree = f"({subtrees[i]}:{bli:g},{subtrees[j]}:{blj:g})"

        if record_steps:
            steps.append({
                "matrix": D.copy(),
                "Q": Q.copy(),
                "names": list(subtrees),
                "joined": (i, j),
                "joined_names": (subtrees[i], subtrees[j]),
                "branch_lengths": (bli, blj),
                "new_subtree": new_subtree,
            })

        # Update distances
        new_row = np.zeros(n)
        for k in range(n):
            if k == i or k == j: continue
            new_row[k] = (D[i, k] + D[j, k] - d_ij) / 2

        keep = [k for k in range(n) if k != i and k != j]
        m = len(keep) + 1
        new_D = np.zeros((m, m))
        for a, ia in enumerate(keep):
            for b, ib in enumerate(keep):
                new_D[a, b] = D[ia, ib]
            new_D[a, -1] = new_D[-1, a] = new_row[ia]
        D = new_D
        subtrees = [subtrees[k] for k in keep] + [new_subtree]
        sizes = [sizes[k] for k in keep] + [sizes[i] + sizes[j]]

    # Final two: connect with one edge of length D[0,1]
    final_d = D[0, 1]
    newick = f"({subtrees[0]}:{final_d/2:g},{subtrees[1]}:{final_d/2:g});"
    return newick, steps


# -------------------- Visualization --------------------

def visualize_nj_steps(names, D, outfile_prefix: str, max_steps: int = 6):
    """One figure per merge step: distance matrix, Q matrix, current join."""
    import matplotlib.pyplot as plt
    import matplotlib.patches as mpatches
    newick, steps = neighbor_joining(names, D)
    out_files = []
    for k, step in enumerate(steps[:max_steps]):
        fig, axes = plt.subplots(1, 2, figsize=(13, 5.2))
        nms = step["names"]; M = step["matrix"]; Q = step["Q"]
        i, j = step["joined"]

        for ax, mat, ttl in [(axes[0], M, "Distance matrix D"),
                             (axes[1], Q, "Q matrix  Q_ij = (n−2)·d_ij − r_i − r_j")]:
            im = ax.imshow(np.where(np.isfinite(mat), mat, np.nan),
                           cmap="viridis_r")
            for a in range(mat.shape[0]):
                for b in range(mat.shape[1]):
                    val = mat[a, b]
                    if np.isfinite(val):
                        col = "white" if val > np.nanmax(np.where(np.isfinite(mat), mat, -np.inf))*0.5 else "black"
                        ax.text(b, a, f"{val:.2f}", ha="center", va="center",
                                fontsize=8, color=col)
            ax.add_patch(mpatches.Rectangle((j-0.5, i-0.5), 1, 1, fill=False,
                                            edgecolor="red", linewidth=3))
            ax.add_patch(mpatches.Rectangle((i-0.5, j-0.5), 1, 1, fill=False,
                                            edgecolor="red", linewidth=3))
            ax.set_xticks(range(len(nms))); ax.set_yticks(range(len(nms)))
            ax.set_xticklabels([nm[:6] for nm in nms], rotation=30, ha="right", fontsize=8)
            ax.set_yticklabels([nm[:6] for nm in nms], fontsize=8)
            ax.set_title(ttl, fontsize=10)
            plt.colorbar(im, ax=ax, fraction=0.04)

        bli, blj = step["branch_lengths"]
        fig.suptitle(f"NJ step {k+1}: pair with min Q is "
                     f"({step['joined_names'][0][:18]}, {step['joined_names'][1][:18]}) "
                     f"→ branch lengths {bli:.3f} / {blj:.3f}", fontsize=11)
        plt.tight_layout()
        f = f"{outfile_prefix}_step{k+1}.png"
        plt.savefig(f, dpi=130, bbox_inches="tight")
        plt.close(fig)
        out_files.append(f)
    return newick, out_files


def draw_nj_tree(newick: str, outfile: str, title="Neighbor Joining tree"):
    """Render a tree with branch lengths, parsed from Newick."""
    import matplotlib.pyplot as plt
    from common import parse_newick, name_internal_nodes
    root = parse_newick(newick); name_internal_nodes(root)

    leaves = list(root.leaves())
    y_of = {l.name: i for i, l in enumerate(leaves)}

    def assign_y(n):
        if n.is_leaf(): return y_of[n.name]
        ys = [assign_y(c) for c in n.children]
        y_of[n.name] = sum(ys) / len(ys)
        return y_of[n.name]
    assign_y(root)

    def x_of(n, x=0):
        # Use sum of branch lengths from root
        if n.parent is None:
            return 0.0
        return x_of(n.parent) + n.branch_length

    fig, ax = plt.subplots(figsize=(10, 6))

    def draw(n):
        x = x_of(n); y = y_of[n.name]
        if n.is_leaf():
            ax.text(x + 0.005, y, n.name, fontsize=10, va="center")
        for c in n.children:
            cx = x_of(c); cy = y_of[c.name]
            ax.plot([x, x], [y_of[c.name], y_of[n.name] if n.parent else cy], "k-", lw=1.2)
            ax.plot([x, cx], [cy, cy], "k-", lw=1.2)
            draw(c)

    # Draw using horizontal cladogram with branch lengths
    def render(n, x0):
        x = x0
        if n.is_leaf():
            ax.plot([x0 - n.branch_length, x0], [y_of[n.name], y_of[n.name]], "k-", lw=1.4)
            ax.text(x0 + 0.003, y_of[n.name], n.name, fontsize=10, va="center")
            return
        ymin = min(y_of[c.name] for c in n.children)
        ymax = max(y_of[c.name] for c in n.children)
        ax.plot([x0, x0], [ymin, ymax], "k-", lw=1.4)
        # Horizontal stub to parent
        if n.parent is not None:
            ax.plot([x0 - n.branch_length, x0], [y_of[n.name], y_of[n.name]], "k-", lw=1.4)
        for c in n.children:
            render(c, x0 + c.branch_length)

    render(root, 0.0)
    ax.set_title(title)
    ax.set_xlabel("Branch length (substitutions/site)")
    ax.set_yticks([]); ax.spines[["top","right","left"]].set_visible(False)
    plt.tight_layout(); plt.savefig(outfile, dpi=140, bbox_inches="tight"); plt.close(fig)
