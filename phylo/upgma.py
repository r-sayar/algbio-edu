"""UPGMA (Unweighted Pair Group Method with Arithmetic averages).

Agglomerative clustering. Repeatedly merge the closest cluster pair, taking the
weighted-by-size arithmetic mean of distances when forming new distances:
    d(k, u) = (n_i * d(k,i) + n_j * d(k,j)) / (n_i + n_j)
where u = i ∪ j. Height of new node = d(i,j) / 2.

Yields an ULTRAMETRIC tree (root-to-leaf distance constant) — assumes molecular clock.
Time complexity O(n^2). Returns a Newick string AND a step trace for visualization.
"""
import sys, os, copy
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np


def upgma(names, D, record_steps=True):
    """names: list of n leaf names; D: n×n distance matrix (numpy).
    Returns (newick_string, steps).
    Each step: dict(matrix=D, names=names, joined=(i,j), height=h)."""
    names = list(names)
    D = D.astype(float).copy()
    n = len(names)
    sizes = [1] * n
    heights = [0.0] * n  # height of each cluster's root
    # Each "name" entry is the Newick subtree string; we update as we merge.
    subtrees = list(names)
    steps = []

    while len(subtrees) > 1:
        # Find min pair
        m = D.shape[0]
        min_val = np.inf; mi, mj = -1, -1
        for i in range(m):
            for j in range(i+1, m):
                if D[i, j] < min_val:
                    min_val = D[i, j]; mi, mj = i, j
        h_new = D[mi, mj] / 2.0
        # Branch lengths for the children
        bli = h_new - heights[mi]
        blj = h_new - heights[mj]
        new_subtree = f"({subtrees[mi]}:{bli:g},{subtrees[mj]}:{blj:g})"

        if record_steps:
            steps.append({
                "matrix": D.copy(),
                "names": list(subtrees),
                "joined_indices": (mi, mj),
                "joined_names": (subtrees[mi], subtrees[mj]),
                "height": h_new,
                "new_subtree": new_subtree,
            })

        # Compute new row/col by weighted average
        ni, nj = sizes[mi], sizes[mj]
        new_row = (ni * D[mi] + nj * D[mj]) / (ni + nj)
        # Build the new D by deleting i,j and appending u
        keep = [k for k in range(m) if k != mi and k != mj]
        new_D = np.zeros((len(keep) + 1, len(keep) + 1))
        for a, ia in enumerate(keep):
            for b, ib in enumerate(keep):
                new_D[a, b] = D[ia, ib]
            new_D[a, -1] = new_D[-1, a] = new_row[ia]
        D = new_D
        subtrees = [subtrees[k] for k in keep] + [new_subtree]
        sizes = [sizes[k] for k in keep] + [ni + nj]
        heights = [heights[k] for k in keep] + [h_new]

    return subtrees[0] + ";", steps


# -------------------- Visualization --------------------

def visualize_upgma_steps(names, D, outfile_prefix: str, max_steps: int = 6):
    """Write one figure per merge step, showing the current distance matrix
    with the about-to-be-merged pair highlighted, and the running dendrogram so far."""
    import matplotlib.pyplot as plt
    import matplotlib.patches as mpatches
    newick, steps = upgma(names, D)
    n_steps = min(len(steps), max_steps)

    out_files = []
    for k, step in enumerate(steps[:n_steps]):
        fig, axes = plt.subplots(1, 2, figsize=(13, 5.5),
                                 gridspec_kw={"width_ratios": [1.1, 1]})

        # Left: distance matrix with min pair highlighted
        ax = axes[0]
        M = step["matrix"]; nms = step["names"]; mi, mj = step["joined_indices"]
        im = ax.imshow(M, cmap="viridis_r")
        for i in range(M.shape[0]):
            for j in range(M.shape[1]):
                color = "white" if M[i, j] > M.max() * 0.5 else "black"
                ax.text(j, i, f"{M[i,j]:.2f}", ha="center", va="center",
                        fontsize=9, color=color)
        # Highlight the chosen pair
        rect = mpatches.Rectangle((mj-0.5, mi-0.5), 1, 1, fill=False,
                                  edgecolor="red", linewidth=3)
        ax.add_patch(rect)
        rect2 = mpatches.Rectangle((mi-0.5, mj-0.5), 1, 1, fill=False,
                                   edgecolor="red", linewidth=3)
        ax.add_patch(rect2)
        ax.set_xticks(range(len(nms))); ax.set_yticks(range(len(nms)))
        ax.set_xticklabels([nm[:8] for nm in nms], rotation=30, ha="right", fontsize=9)
        ax.set_yticklabels([nm[:8] for nm in nms], fontsize=9)
        ax.set_title(f"Step {k+1}: merge {step['joined_names'][0][:18]} + "
                     f"{step['joined_names'][1][:18]}\nat height {step['height']:.3f}",
                     fontsize=11)
        plt.colorbar(im, ax=ax, fraction=0.04)

        # Right: cumulative dendrogram via scipy on the original matrix
        ax = axes[1]
        from scipy.cluster.hierarchy import linkage, dendrogram
        from scipy.spatial.distance import squareform
        Z = linkage(squareform(D), method="average")
        dendrogram(Z, labels=names, ax=ax, color_threshold=0)
        ax.set_title("Final UPGMA dendrogram (for reference)")
        ax.set_ylabel("Height")

        plt.tight_layout()
        f = f"{outfile_prefix}_step{k+1}.png"
        plt.savefig(f, dpi=130, bbox_inches="tight")
        plt.close(fig)
        out_files.append(f)

    return newick, out_files


def visualize_upgma_dendrogram(names, D, outfile: str, title="UPGMA dendrogram"):
    import matplotlib.pyplot as plt
    from scipy.cluster.hierarchy import linkage, dendrogram
    from scipy.spatial.distance import squareform
    fig, ax = plt.subplots(figsize=(10, 5.5))
    Z = linkage(squareform(D), method="average")
    dendrogram(Z, labels=names, ax=ax, leaf_rotation=30)
    ax.set_title(title); ax.set_ylabel("Height (= d/2)")
    plt.tight_layout(); plt.savefig(outfile, dpi=140, bbox_inches="tight"); plt.close(fig)
