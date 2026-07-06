"""Blueprint webapp design tokens + Plotly template."""
from __future__ import annotations

import plotly.graph_objects as go

ACCENT = "#2558d4"
BG = "#f6f7fb"
INK = "#172033"
MUTED = "#5f6b7a"
LINE = "#d9dee8"
GREEN = "#117a37"
RED = "#c4222c"
WHITE = "#FFFFFF"

# Okabe–Ito role colors for the blueprint sidebar gene chips.
ROLE_POSITIVE = "#0072B2"
ROLE_GATE = "#CC79A7"
ROLE_NEGATIVE = "#E69F00"

# Sequential: light bg -> blue -> green (positive-only quantities).
SEQUENTIAL = [
    [0.0, BG], [0.5, ROLE_POSITIVE], [1.0, GREEN],
]
# Diverging/symmetric for contrast/gated axes (neg=blue, mid=light, pos=red).
DIVERGING = [
    [0.0, ACCENT], [0.5, BG], [1.0, RED],
]

FONT_FAMILY = ("Inter, -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, "
               "Helvetica, Arial, sans-serif")


def plotly_template() -> go.layout.Template:
    return go.layout.Template(layout=dict(
        paper_bgcolor=WHITE,
        plot_bgcolor=WHITE,
        font=dict(family=FONT_FAMILY, color=INK),
        colorway=[ACCENT, GREEN, ROLE_GATE, ROLE_NEGATIVE, MUTED],
        margin=dict(l=40, r=20, t=40, b=40),
    ))
