"""Blueprint gate-tree branch derivation for the interactive UMAP explorer."""
from __future__ import annotations


def branches_from_summary(summary):
    """Root->leaf paths of the blueprint gate tree, longest first.

    A non-root (gated) axis's parent is the axis whose positive module equals
    this axis's gate module. Singletons (length-1 paths) are skipped.
    Each branch: {"label", "leaf", "axes": [root..leaf]}.
    """
    axes = summary.get("blueprint", {}).get("axes", [])
    by_positive = {}
    for ax in axes:
        pos = ax.get("positive")
        if pos is not None:
            by_positive.setdefault(pos, ax["name"])
    parent_of = {}
    for ax in axes:
        gate = ax.get("gate")
        if not gate:
            continue
        parent = by_positive.get(gate)
        if parent and parent != ax["name"]:
            parent_of.setdefault(ax["name"], parent)
    parents = set(parent_of.values())
    leaves = [ax["name"] for ax in axes if ax["name"] not in parents]
    branches = []
    for leaf in leaves:
        path = [leaf]
        cur = leaf
        while cur in parent_of:
            cur = parent_of[cur]
            if cur in path:
                break
            path.append(cur)
        path = list(reversed(path))
        if len(path) < 2:
            continue
        branches.append({
            "label": f"{path[0]} → {path[-1]} ({len(path)} steps)",
            "leaf": leaf, "axes": path,
        })
    branches.sort(key=lambda b: (-len(b["axes"]), b["label"]))
    return branches
