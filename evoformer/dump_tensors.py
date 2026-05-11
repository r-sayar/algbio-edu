"""Run the toy Evoformer block and dump m / z snapshots as JSON
for the Three.js viewer to load.

Each snapshot is shape-flattened into a 1-D float array; the JSON also
records the shape so the renderer can rebuild the voxel grid.
"""

import os
import sys
import json
import torch

HERE = os.path.dirname(__file__)
sys.path.insert(0, os.path.abspath(os.path.join(HERE, "..", "src")))
from sublayers import (
    msa_row_attention_with_pair_bias,
    outer_product_mean,
    triangle_mult_outgoing,
    triangle_attention_starting_node,
)

# Slightly smaller than dataflow.png so the 3D voxel grids are legible.
torch.manual_seed(0)
N_seq, N_res, c = 3, 6, 4

m = torch.randn(N_seq, N_res, c)
z = torch.randn(N_res, N_res, c) * 0.3

snapshots = [{"name": "input",                  "m": m, "z": z}]

# 1. MSA row attention with pair bias  (pair -> MSA)
m1, _ = msa_row_attention_with_pair_bias(m, z)
m_running, z_running = m1, z
snapshots.append({"name": "after MSA row attn (pair->MSA bias)",
                  "m": m_running, "z": z_running})

# 4. Outer-product mean  (MSA -> pair)
opm = outer_product_mean(m_running)
z_running = z_running + opm
snapshots.append({"name": "after outer-product mean (MSA->pair)",
                  "m": m_running, "z": z_running})

# 5. Triangle multiplicative outgoing
z_running = triangle_mult_outgoing(z_running)
snapshots.append({"name": "after triangle mult outgoing",
                  "m": m_running, "z": z_running})

# 7. Triangle attention starting node
z_running, _ = triangle_attention_starting_node(z_running)
snapshots.append({"name": "after triangle attn starting node",
                  "m": m_running, "z": z_running})


def serialize(t: torch.Tensor):
    """JSON-friendly dump: list of floats + shape."""
    return {"shape": list(t.shape), "data": t.flatten().tolist()}


payload = {
    "N_seq": N_seq, "N_res": N_res, "c": c,
    "snapshots": [
        {"name": s["name"], "m": serialize(s["m"]), "z": serialize(s["z"])}
        for s in snapshots
    ],
}

out = os.path.join(HERE, "tensors.json")
with open(out, "w") as f:
    json.dump(payload, f)
print(f"wrote {out}  ({len(payload['snapshots'])} snapshots)")
