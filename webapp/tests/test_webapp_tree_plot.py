import io

from matplotlib.figure import Figure

from webapp.pipeline import CELLTYPE_NOTE, TEMPORAL_NOTE
from webapp import tree_plot


def _bundle():
    """A 3-node cell-type tree: root (plain) -> child (gated) -> grandchild (gated)."""
    return {
        "axisOrder": ["imm_immuneRoot", "tcell_tCell", "tcell_cd8"],
        "axisMeta": {
            "imm_immuneRoot": {"kind": "plain", "positive": "immune",
                               "gate": None, "negatives": []},
            "tcell_tCell": {"kind": "gated", "positive": "tCell",
                            "gate": "immune", "negatives": ["bCell"]},
            "tcell_cd8": {"kind": "gated", "positive": "cd8Tcell",
                          "gate": "tCell", "negatives": ["cd4Tcell"]},
        },
    }


def test_build_edges_parentage():
    axes = tree_plot.axes_from_bundle(_bundle())
    edges, parent_of = tree_plot.build_edges(axes, stitch=True)
    assert parent_of == {"tcell_tCell": "imm_immuneRoot",
                         "tcell_cd8": "tcell_tCell"}
    # No ordinals in the names -> the stitch heuristic must not fire.
    assert all(k == "gate" for *_, k in edges)


def test_figure_from_bundle_shape_and_xlabel():
    b = _bundle()
    fig_ct = tree_plot.figure_from_bundle(b, "celltype")
    assert isinstance(fig_ct, Figure)
    # 3 nodes scattered (one PathCollection per axis).
    assert len(fig_ct.axes[0].collections) == 3
    xlabel_ct = fig_ct.axes[0].get_xlabel()
    xlabel_time = tree_plot.figure_from_bundle(b, "time").axes[0].get_xlabel()
    assert xlabel_ct != xlabel_time
    assert "broad class" in xlabel_ct and "earliest state" in xlabel_time
    # Unknown kind falls back to the cell-type label.
    assert tree_plot.figure_from_bundle(b, "bogus").axes[0].get_xlabel() == xlabel_ct


def test_detect_tree_kind():
    assert tree_plot.detect_tree_kind(TEMPORAL_NOTE) == "time"
    assert tree_plot.detect_tree_kind(CELLTYPE_NOTE) == "celltype"
    assert tree_plot.detect_tree_kind(None) == "celltype"
    assert tree_plot.detect_tree_kind("focus on myeloid") == "celltype"
    assert tree_plot.detect_tree_kind(f"custom\n\n{TEMPORAL_NOTE}") == "time"


def test_render_to_bytes_each_format():
    fig = tree_plot.figure_from_bundle(_bundle(), "celltype")
    magic = {"png": b"\x89PNG", "pdf": b"%PDF", "svg": b"<?xml"}
    for fmt, head in magic.items():
        buf = io.BytesIO()
        fig.savefig(buf, format=fmt, dpi=300, bbox_inches="tight")
        data = buf.getvalue()
        assert data, f"{fmt} produced no bytes"
        assert data[:len(head)] == head, f"{fmt} bad magic: {data[:8]!r}"
