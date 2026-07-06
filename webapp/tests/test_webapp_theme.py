from webapp import theme


def test_km_tokens_present():
    assert theme.ACCENT == "#2558d4"
    assert theme.BG == "#f6f7fb"
    assert theme.INK == "#172033"
    assert theme.MUTED == "#5f6b7a"
    assert theme.LINE == "#d9dee8"


def test_role_colors_present():
    assert theme.ROLE_POSITIVE == "#0072B2"
    assert theme.ROLE_GATE == "#CC79A7"
    assert theme.ROLE_NEGATIVE == "#E69F00"


def test_plotly_template_is_white():
    tpl = theme.plotly_template()
    assert tpl.layout.paper_bgcolor == "#FFFFFF"
    assert tpl.layout.plot_bgcolor == "#FFFFFF"
