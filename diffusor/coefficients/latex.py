r"""The equation of every law as LaTeX, built from the law's own parameter values.

The strings use the subset of LaTeX that matplotlib's mathtext renders (no TeX
installation is needed), so the interface can typeset them and the exports
can carry them as LaTeX source. Building them from ``params`` rather than
retyping them means the displayed equation is the one the code evaluates; the
test suite checks that every registered law and tracer family has one and that
each renders.

Conventions: D in m^2/s, T in K, R the gas constant, fO2 in the unit the law
states, X mole fractions.
"""
from __future__ import annotations

import math

R = r"R"


def num(v: float, digits: int = 6) -> str:
    """A number for LaTeX at full transcribed precision: 2.77\\times10^{-7}, 0.29, 262914."""
    v = float(v)
    if v == 0:
        return "0"
    a = abs(v)
    if 1e-2 <= a < 1e6:
        if a >= 1e4 and v.is_integer():
            return str(int(v))
        return f"{v:.{digits}g}"
    e = int(math.floor(math.log10(a)))
    m = v / 10 ** e
    ms = f"{m:.{digits}g}"
    if ms in ("1", "-1"):
        return ("-" if v < 0 else "") + rf"10^{{{e}}}"
    return rf"{ms}\times10^{{{e}}}"


def term(coef: float, var: str, digits: int = 6) -> str:
    """'+ 2.6\\,X' or '- 2.37\\,X' with the sign of ``coef``; a bare constant without ``var``."""
    sign = "-" if coef < 0 else "+"
    return rf" {sign} {num(abs(coef), digits)}" + (rf"\,{var}" if var else "")


def kj(q_kj: float) -> str:
    return rf"{num(q_kj)}\ \mathrm{{kJ\,mol^{{-1}}}}"


def arr_exp(q_kj: float, extra: str = "") -> str:
    r"""\exp\left(-\frac{Q}{RT}\right)"""
    top = kj(q_kj) + extra
    return rf"\exp\left(-\frac{{{top}}}{{{R}T}}\right)"


UNIT = r"\ \mathrm{m^2\,s^{-1}}"


def D_arrhenius(D0: float, q_kj: float, label: str = "D") -> str:
    return rf"{label} = {num(D0)}\,{arr_exp(q_kj)}{UNIT}"


def logD_arrhenius(logD0: float, q_kj: float, label: str = "D") -> str:
    return rf"\log_{{10}} {label} = {num(logD0)} - \frac{{{kj(q_kj)}}}{{2.303\,{R}T}}"


def _p(c, name):
    return c.params[name].value


# ---------------------------------------------------------------------------
def _dohmen_chakraborty(c):
    tamed = c.key.endswith("tamed")
    s = (rf"\log_{{10}} D_{{[001]}} = {num(_p(c, 'c0'))} - \frac{{{num(_p(c, 'Q') * 1e3, 6)} + (P - 10^5)\,"
         rf"{num(_p(c, 'dV'))}}}{{2.303\,{R}T}}")
    if tamed:
        s += rf" + \frac{{1}}{{6}}\log_{{10}}\frac{{f_{{O_2}}}}{{10^{{-7}}}}"
    s += rf" + {num(_p(c, 'm'))}\,(X_{{Fe}} - {num(_p(c, 'XFe_ref'))})"
    return s + r"\quad (P,\ f_{O_2}\ \mathrm{in\ Pa})"


def _oeser_axis(c):
    ax = c.reference_axis or "c"
    return (rf"\log_{{10}} D^*_{{\mathrm{{{c.species}}},{ax}}} = \log_{{10}} D_{{0,{ax}}} - "
            rf"\frac{{E_{{a,{ax}}}}}{{2.303\,{R}T}} + {num(_p(c, 'n'))}\,(X_{{Fe}} - {num(_p(c, 'XFe_ref'))})"
            + r";\ " + r",\ ".join(
                rf"{a}:\ {num(_p(c, f'logD0_{a}'))},\ {num(_p(c, f'Ea_{a}'))}\ \mathrm{{kJ\,mol^{{-1}}}}"
                for a in ("a", "b", "c")))


def _binary(c):
    a, b = c.derived_from
    from .registry import REGISTRY
    A, B = REGISTRY[a], REGISTRY[b]
    return (rf"D_{{\mathrm{{{c.species}}}}} = \frac{{D^*_{{\mathrm{{{A.species}}}}}\,D^*_{{\mathrm{{{B.species}}}}}}}"
            rf"{{X_{{\mathrm{{{A.species}}}}}\,D^*_{{\mathrm{{{A.species}}}}} + X_{{\mathrm{{{B.species}}}}}\,"
            rf"D^*_{{\mathrm{{{B.species}}}}}}}")


def _dohmen2016_fs9(c):
    return (rf"\log_{{10}} D_{{[001]}} = {num(_p(c, 'logD0'))} - \frac{{{kj(_p(c, 'Q'))}}}{{2.303\,{R}T}}"
            rf" + {num(_p(c, 'n'))}\,\log_{{10}} f_{{O_2}}[\mathrm{{Pa}}] + {num(_p(c, 'm'))}\,(X_{{Fe}} - 0.09)")


def _dias2025(c):
    return (rf"\mathrm{{for}}\ \log_{{10}} f_{{O_2}} > -10;\ D = {num(_p(c, 'D0_1'))}\left(\frac{{f_{{O_2}}}}{{10^{{-7}}}}\right)^{{{num(_p(c, 'n_1'))}}}"
            rf"{arr_exp(_p(c, 'Q_1'))}\,10^{{m_1(X_{{Fe}} - 0.1)}};\ \mathrm{{where}}\ m_1 = {num(_p(c, 'm1_slope'))}\,\frac{{10^4}}{{T}}"
            rf"{term(_p(c, 'm1_int'), '')};\ "
            rf"\mathrm{{for}}\ \log_{{10}} f_{{O_2}} \leq -10;\ D = {num(_p(c, 'D0_2'))}{arr_exp(_p(c, 'Q_2'))}\,10^{{m_2(X_{{Fe}} - 0.1)}};"
            rf"\ \mathrm{{where}}\ m_2 = {num(_p(c, 'm2_slope'))}\,\frac{{10^4}}{{T}}{term(_p(c, 'm2_int'), '')}\quad (f_{{O_2}}\ \mathrm{{in\ Pa}})")


def _dias2024(c):
    return (rf"D = {num(_p(c, 'D0'))}\,{arr_exp(_p(c, 'Q'))}\,10^{{m(X_{{Fe}} - 0.1)}};\ "
            rf"\mathrm{{where}}\ m = \frac{{{num(_p(c, 'm_slope'))}}}{{T}}{term(_p(c, 'm_int'), '')}")


def _dias_ree(c):
    s = rf"D_{{\mathrm{{{c.species}}}}} = {num(_p(c, 'D0'))}"
    if _p(c, "inv_n"):
        s += rf"\left(\frac{{f_{{O_2}}}}{{f_{{O_2}}^{{IW}}}}\right)^{{1/{round(1 / _p(c, 'inv_n'))}}}"
    return s + rf"\,{arr_exp(_p(c, 'Q'))}{UNIT}"


def _ganguly(c):
    s = rf"\log_{{10}} D = {num(_p(c, 'c0'))}{term(_p(c, 'a'), 'X_{Fe}')} - \frac{{{num(_p(c, 'b'))}}}{{T}}"
    if "n_fo2" in c.params:
        s += r" + \frac{1}{6}\log_{10}\frac{f_{O_2}}{f_{O_2}^{IW}}"
    return s


def _ln_xan(c):
    return (rf"\ln D = {num(_p(c, 'lnD0'))}{term(_p(c, 'a'), 'X_{An}')} - \frac{{{kj(_p(c, 'Q'))}}}{{{R}T}}")


def _log_xan(c):
    return (rf"\log_{{10}} D = {num(_p(c, 'a'))}\,X_{{An}}{term(_p(c, 'b'), '')} - "
            rf"\frac{{{num(_p(c, 'Q') * 1e3, 6)}}}{{2.303\,{R}T}}")


def _audetat(c):
    return _log_xan(c) + rf"{term(_p(c, 'c'), '(1 - a_{SiO_2})')}"


def _giletti(c):
    return rf"D = {num(_p(c, 'D0'))}\,{arr_exp(_p(c, 'Q'))}\,10^{{{num(_p(c, 'a'))}\,X_{{An}}}}{UNIT}"


def _feti(c):
    return rf"\ln D = {num(_p(c, 'c0'))}{term(_p(c, 'a'), 'x_{Ti}')} - \frac{{{num(_p(c, 'b'))}}}{{T}}"


def _table12(c):
    from .magnetite import TABLE12_PURE, TABLE12_XTI02
    v0, qv, i0, qi = TABLE12_PURE[c.species]
    s = (rf"D^* = {num(v0)}\,\exp\left(-\frac{{({num(qv)})\ \mathrm{{kJ\,mol^{{-1}}}}}}{{{R}T}}\right)a_{{O_2}}^{{2/3}}"
         rf" + {num(i0)}\,{arr_exp(qi)}\,a_{{O_2}}^{{-2/3}}")
    if c.species in TABLE12_XTI02:
        s += r"\ (x_{Ti}=0;\ \log D\ \mathrm{interpolated\ to}\ x_{Ti}=0.2)"
    return s + r",\ a_{O_2} = f_{O_2}/1\,\mathrm{atm}"


def _sievwright(c):
    from .magnetite import SIEVWRIGHT_TABLE5
    lv, li = SIEVWRIGHT_TABLE5[c.species][:2]
    s = rf"D = 10^{{{num(lv)}}}\,f_{{O_2}}^{{2/3}} + 10^{{{num(li)}}}\,f_{{O_2}}^{{-2/3}}\quad (f_{{O_2}}\ \mathrm{{in\ bar}},\ T = 1150^\circ\mathrm{{C}})"
    if c.fixed_temperature_K is None:
        s += r";\ \mathrm{each\ term\ scaled\ to\ }T\ \mathrm{with\ Table\ 12}\ Q_V, Q_I"
    return s


def _axis_law(c):
    fo2 = (rf" + {num(_p(c, 'n_fo2'))}\,\log_{{10}}\frac{{f_{{O_2}}}}{{10^{{-12}}\,\mathrm{{bar}}}}"
           if "n_fo2" in c.params else "")
    parts = []
    for a in ("a", "b", "c"):
        if f"logD0_{a}" in c.params:
            parts.append(rf"\log_{{10}} D_{a} = {num(_p(c, f'logD0_{a}'))} - "
                         rf"\frac{{{kj(_p(c, f'Q_{a}'))}}}{{2.303\,{R}T}}" + fo2)
    return r";\ ".join(parts)


def _borinski(c):
    return (rf"\log_{{10}} D^*_{{\mathrm{{{c.species}}}}} = \log_{{10}}({num(_p(c, 'D0'))}) - "
            rf"\frac{{{num(_p(c, 'Q') * 1e3, 6)} + {num(_p(c, 'dV'))}\,(P - 1)}}{{2.303\,{R}T}}\quad (P\ \mathrm{{in\ bar}})")


def _scalar(c):
    return logD_arrhenius(_p(c, "logD0"), _p(c, "Q"))


def _D0Q(c):
    return D_arrhenius(_p(c, "D0"), _p(c, "Q"))


def _logD0Q(c):
    return D_arrhenius(10 ** _p(c, "logD0"), _p(c, "Q"))


_BY_FUNC = {
    "_dohmen_chakraborty_tamed": _dohmen_chakraborty, "_dohmen_chakraborty_ped": _dohmen_chakraborty,
    "_chakraborty1997": _D0Q, "_plain_arrhenius": _D0Q, "_arrhenius_law": _D0Q, "_buffer_arrhenius": _D0Q,
    "_grove1984": _D0Q, "_schwandt": _logD0Q, "_dohmen2016_fs1": _logD0Q, "_pohl2024": _scalar,
    "_scalar": _scalar, "_dohmen2016_fs9": _dohmen2016_fs9, "_dias2025": _dias2025, "_dias2024": _dias2024,
    "_dias2025ree": _dias_ree, "_ganguly_tazzoli": _ganguly, "_ganguly_tazzoli_no_fo2": _ganguly,
    "_vanorman2014": _ln_xan, "_costa2003_mg": _ln_xan, "_grocolas_form": _log_xan, "_audetat2026_mg": _audetat,
    "_giletti_casserly1994": _giletti, "_feti_lnD": _feti, "_borinski": _borinski,
}


def coefficient_latex(c) -> str:
    """LaTeX (mathtext subset) for one registered law."""
    if c.derived_from:
        return _binary(c)
    q = c.func.__qualname__ if c.func is not None else ""
    if q.startswith("_oeser_axis"):
        return _oeser_axis(c)
    if q.startswith("_make_table12_func"):
        return _table12(c)
    if q.startswith("_make_sievwright_func"):
        return _sievwright(c)
    if q.startswith("_axis_law"):
        return _axis_law(c)
    name = q.split(".")[0]
    if name not in _BY_FUNC:
        raise KeyError(f"no LaTeX form for {c.key} ({q})")
    return _BY_FUNC[name](c)


def family_latex(f) -> str:
    """LaTeX for a garnet tracer family and the matrix it feeds."""
    if f.citation == "carlson2006":
        law = (r"\ln D^*_i = \ln D^*_{0,alm,i} + k_i\,(a_0 - 1.1525) - \frac{Q_i + P\,\Delta V_i}{RT}"
               r" + \frac{1}{6}\ln\frac{f_{O_2}}{f_{O_2}^{gr}},\ a_0 = \sum_j X_j\,a_{0,j}")
    else:
        law = r"\ln D^*_i = \ln D_{0,i} - \frac{Q_i}{RT} - \frac{(P - 1)\,\Delta V_i}{RT},\ D^*_{Ca} = 0.5\,D^*_{Fe}"
    matrix = (r"D_{ij} = D^*_i\,\delta_{ij} - \frac{D^*_i X_i}{\sum_k X_k D^*_k}\,(D^*_j - D^*_n)")
    return law + r";\quad " + matrix


# Typeset size: 15 pt at 100 dpi is a cap height of about 21 logical pixels, which keeps the
# subscripts and the fractions readable next to 13 px interface text.
BASE_DPI = 100
FONT_SIZE = 15
_TEXT_COLOUR = "#1f2a30"


def render_png(latex: str, dpi: float = BASE_DPI, fontsize: int = FONT_SIZE) -> bytes:
    """PNG bytes of the typeset equation (Computer Modern), for the interface."""
    import io
    import matplotlib
    from matplotlib.figure import Figure
    from matplotlib.backends.backend_agg import FigureCanvasAgg
    with matplotlib.rc_context({"mathtext.fontset": "cm", "text.color": _TEXT_COLOUR}):
        fig = Figure(figsize=(0.01, 0.01))
        FigureCanvasAgg(fig)
        fig.text(0, 0, f"${latex}$", fontsize=fontsize)
        buf = io.BytesIO()
        fig.savefig(buf, dpi=dpi, format="png", bbox_inches="tight", pad_inches=0.04, transparent=True)
    return buf.getvalue()


def equation_lines(latex: str) -> list:
    """Split a multi-branch equation at its top-level ';' separators (mathtext has no line breaks)."""
    out, depth, start, i = [], 0, 0, 0
    seps = (r";\quad ", r";\ ")
    while i < len(latex):
        ch = latex[i]
        if ch == "{":
            depth += 1
        elif ch == "}":
            depth -= 1
        elif depth == 0:
            for sep in seps:
                if latex.startswith(sep, i):
                    out.append(latex[start:i])
                    i += len(sep)
                    start = i
                    break
            else:
                i += 1
                continue
            continue
        i += 1
    out.append(latex[start:])
    return [s.strip() for s in out if s.strip()]


def png_size(png: bytes):
    """(width, height) in pixels from a PNG header."""
    import struct
    return struct.unpack(">II", png[16:24])


def _text_width(latex: str, dpi: float, fontsize: int = FONT_SIZE) -> float:
    """Width in pixels of the typeset equation at ``dpi``, measured without drawing it."""
    import matplotlib
    from matplotlib.backends.backend_agg import FigureCanvasAgg
    from matplotlib.figure import Figure
    from matplotlib.font_manager import FontProperties
    with matplotlib.rc_context({"mathtext.fontset": "cm"}):
        renderer = FigureCanvasAgg(Figure(dpi=dpi)).get_renderer()
        return renderer.get_text_width_height_descent(
            f"${latex}$", FontProperties(size=fontsize), ismath=True)[0]


def _break_points(latex: str) -> list:
    """Indices of the top-level ' + ' and ' - ' where a long line may be cut.

    Operators inside braces (exponents, subscripts, fractions), parentheses or between
    \\left and \\right are never cut, and neither is a sign that opens a term after '=',
    '(' or ','."""
    points, depth, fence, paren, i = [], 0, 0, 0, 0
    while i < len(latex):
        if latex.startswith(r"\left", i):
            fence += 1
            i += 5
            continue
        if latex.startswith(r"\right", i):
            fence -= 1
            i += 6
            continue
        ch = latex[i]
        if ch == "{":
            depth += 1
        elif ch == "}":
            depth -= 1
        elif ch == "(":
            paren += 1
        elif ch == ")":
            paren -= 1
        elif (depth == 0 and fence == 0 and paren == 0 and ch in "+-" and 0 < i < len(latex) - 1
              and latex[i - 1] == " " and latex[i + 1] == " "
              and latex[:i].rstrip()[-1:] not in "=(,:+-<>"
              and not latex[:i].rstrip().endswith(("\\", r"\leq", r"\geq"))):
            points.append(i)
        i += 1
    return points


def wrap_equation(latex: str, max_px: float, dpi: float = BASE_DPI) -> list:
    """One equation line cut into pieces that each fit in ``max_px`` pixels at ``dpi``.

    The cuts are at top-level operators, so every piece is valid on its own, and the
    continuation lines are indented. A line with no usable cut is returned whole."""
    if _text_width(latex, dpi) <= max_px:
        return [latex]
    points = _break_points(latex)
    pieces, start, prefix = [], 0, ""
    while True:
        rest = prefix + latex[start:]
        if _text_width(rest, dpi) <= max_px:
            pieces.append(rest.strip())
            return pieces
        later = [p for p in points if p > start]
        if not later:
            pieces.append(rest.strip())
            return pieces
        best = later[0]
        for p in later:
            if _text_width(prefix + latex[start:p], dpi) > max_px:
                break
            best = p
        pieces.append((prefix + latex[start:best]).strip())
        start, prefix = best, r"\quad "


def pixel_ratio() -> float:
    """Device pixels per logical pixel of the main screen (1 without Qt or a screen)."""
    try:
        from PySide6.QtGui import QGuiApplication
        screen = QGuiApplication.primaryScreen()
        return max(1.0, float(screen.devicePixelRatio())) if screen is not None else 1.0
    except Exception:
        return 1.0


_HTML_CACHE: dict = {}


_PIECE_CACHE: dict = {}


def equation_pieces(latex: str, max_width: int = 700, ratio: float = None) -> list:
    """One ``(png bytes, width, height)`` per drawn line of the equation, sizes in logical pixels.

    The pictures are drawn with one pixel for every pixel of the screen, so Qt shows
    them as they are. Pictures drawn at a fixed resolution and then stretched or
    shrunk by Qt look soft on a scaled display. Lines wider than ``max_width``
    logical pixels are cut at operators."""
    ratio = pixel_ratio() if ratio is None else float(ratio)
    key = (latex, max_width, ratio)
    if key in _PIECE_CACHE:
        return _PIECE_CACHE[key]
    dpi = BASE_DPI * ratio
    pieces = []
    for line in equation_lines(latex):
        for piece in wrap_equation(line, max_width * ratio, dpi):
            png = render_png(piece, dpi=dpi)
            w, h = png_size(png)
            if w > max_width * ratio:
                # nothing to cut at: draw it smaller instead of stretching the picture
                png = render_png(piece, dpi=dpi * max_width * ratio / w)
                w, h = png_size(png)
            pieces.append((png, round(w / ratio), round(h / ratio)))
    if len(_PIECE_CACHE) > 400:
        _PIECE_CACHE.clear()
    _PIECE_CACHE[key] = pieces
    return pieces


def equation_html(latex: str, max_width: int = 700, ratio: float = None) -> str:
    """<img> tags, one per line of the equation, as inline data for Qt rich text."""
    import base64
    ratio = pixel_ratio() if ratio is None else float(ratio)
    key = (latex, max_width, ratio)
    if key in _HTML_CACHE:
        return _HTML_CACHE[key]
    html = "<br>".join(
        f"<img src='data:image/png;base64,{base64.b64encode(png).decode()}' "
        f"width='{w}' height='{h}'>" for png, w, h in equation_pieces(latex, max_width, ratio))
    _HTML_CACHE[key] = html
    return html
