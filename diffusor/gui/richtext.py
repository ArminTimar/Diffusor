"""Formatted documents for the reading panes: methods, coefficients, log and help.

Qt renders a subset of HTML 4 and CSS, so the layout sticks to headings,
paragraphs, lists and tables, which it draws reliably.
"""
from __future__ import annotations

import html
import platform
from datetime import datetime, timezone
from pathlib import Path
from typing import Iterable, List, Optional, Sequence, Tuple

import numpy as np
from PySide6.QtCore import Qt
from PySide6.QtGui import QGuiApplication
from PySide6.QtWidgets import QDialog, QDialogButtonBox, QTextBrowser, QVBoxLayout

from .. import __version__
from ..references import cite, format_reference, get as get_reference
from ..thermo.units import human_time
from . import theme

CSS = f"""
body {{ font-family: 'Segoe UI', Arial, sans-serif; font-size: 13px; color: {theme.TEXT}; }}
h1 {{ font-size: 19px; font-weight: 600; margin: 0 0 2px 0; color: {theme.TEXT}; }}
h2 {{ font-size: 14px; font-weight: 600; margin: 18px 0 6px 0; color: {theme.ACCENT}; }}
h3 {{ font-size: 13px; font-weight: 600; margin: 10px 0 4px 0; color: {theme.TEXT}; }}
p {{ margin: 0 0 8px 0; line-height: 120%; }}
ul {{ margin: 2px 0 8px 0; }}
li {{ margin: 0 0 3px 0; line-height: 118%; }}
.sub {{ color: {theme.TEXT_MUTED}; font-size: 12.5px; }}
.muted {{ color: {theme.TEXT_MUTED}; }}
.faint {{ color: {theme.TEXT_FAINT}; font-size: 12px; }}
.mono {{ font-family: Consolas, 'Cascadia Mono', monospace; font-size: 12.5px; }}
.key {{ color: {theme.TEXT_MUTED}; padding-right: 14px; }}
.num {{ font-family: Consolas, 'Cascadia Mono', monospace; }}
th {{ color: {theme.TEXT_MUTED}; font-weight: 600; text-align: left; padding: 4px 12px 4px 0;
      border-bottom: 1px solid {theme.BORDER}; }}
td {{ padding: 4px 12px 4px 0; vertical-align: top; }}
a {{ color: {theme.ACCENT}; }}
"""

_BADGES = {
    "ok": (theme.OK, theme.OK_SOFT),
    "warn": (theme.WARN, theme.WARN_SOFT),
    "danger": (theme.DANGER, "#F8EBEB"),
    "accent": (theme.ACCENT, theme.ACCENT_SOFT),
    "muted": (theme.TEXT_MUTED, "#F1EEE9"),
}


def esc(text) -> str:
    return html.escape(str(text))


class Doc:
    """A small HTML document builder."""

    def __init__(self, title: str = "", subtitle: str = "", badges: Sequence[Tuple[str, str]] = ()):
        self.parts: List[str] = []
        if title:
            self.parts.append(f"<h1>{esc(title)}</h1>")
        b = "&nbsp;&nbsp;".join(badge(t, k) for t, k in badges)
        if subtitle or b:
            sep = "&nbsp;&nbsp;&nbsp;" if subtitle and b else ""
            self.parts.append(f"<p class='sub'>{subtitle}{sep}{b}</p>")

    def h(self, text: str):
        self.parts.append(f"<h2>{esc(text)}</h2>")
        return self

    def h3(self, text: str):
        self.parts.append(f"<h3>{esc(text)}</h3>")
        return self

    def p(self, text: str, cls: str = ""):
        c = f" class='{cls}'" if cls else ""
        self.parts.append(f"<p{c}>{text}</p>")
        return self

    def ul(self, items: Iterable[str]):
        items = list(items)
        if items:
            self.parts.append("<ul>" + "".join(f"<li>{i}</li>" for i in items) + "</ul>")
        return self

    def kv(self, rows: Iterable[Tuple[str, str]]):
        body = "".join(f"<tr><td class='key'>{esc(k)}</td><td>{v}</td></tr>" for k, v in rows)
        self.parts.append(f"<table cellspacing='0' cellpadding='0'>{body}</table>")
        return self

    def table(self, header: Sequence[str], rows: Iterable[Sequence[str]]):
        head = "".join(f"<th>{esc(h)}</th>" for h in header)
        body = "".join("<tr>" + "".join(f"<td>{c}</td>" for c in r) + "</tr>" for r in rows)
        self.parts.append(f"<table cellspacing='0' cellpadding='0' width='100%'>"
                          f"<tr>{head}</tr>{body}</table>")
        return self

    def box(self, text: str, kind: str = "accent"):
        fg, bg = _BADGES.get(kind, _BADGES["accent"])
        colour = theme.TEXT if kind in ("accent", "muted") else fg
        self.parts.append(
            f"<table width='100%' cellspacing='0' cellpadding='10' style='margin: 4px 0 10px 0;'>"
            f"<tr><td style='background-color:{bg}; color:{colour};'>{text}</td></tr></table>")
        return self

    def equation(self, text: str):
        self.parts.append(
            f"<table width='100%' cellspacing='0' cellpadding='10' style='margin: 2px 0 8px 0;'>"
            f"<tr><td style='background-color:{theme.SURFACE_ALT}; border: 1px solid {theme.BORDER};'>"
            f"<span class='mono'>{esc(text)}</span></td></tr></table>")
        return self

    def raw(self, fragment: str):
        self.parts.append(fragment)
        return self

    def html(self) -> str:
        return f"<html><head><style>{CSS}</style></head><body>{''.join(self.parts)}</body></html>"


def badge(text: str, kind: str = "muted") -> str:
    fg, bg = _BADGES.get(kind, _BADGES["muted"])
    return (f"<span style='background-color:{bg}; color:{fg}; font-size:11px; font-weight:600;'>"
            f"&nbsp;{esc(text)}&nbsp;</span>")


def show(parent, title: str, document: str, width: int = 820, height: int = 620,
         modal: bool = True) -> QDialog:
    """Open a document in a reading pane."""
    dlg = QDialog(parent)
    dlg.setWindowTitle(title)
    screen = QGuiApplication.primaryScreen().availableGeometry()
    dlg.resize(min(width, int(screen.width() * 0.9)), min(height, int(screen.height() * 0.85)))
    lay = QVBoxLayout(dlg)
    lay.setContentsMargins(16, 16, 16, 12)
    lay.setSpacing(10)
    view = QTextBrowser()
    view.setOpenExternalLinks(True)
    view.document().setDocumentMargin(22)
    view.setHtml(document)
    lay.addWidget(view)
    bb = QDialogButtonBox(QDialogButtonBox.Close)
    bb.rejected.connect(dlg.reject)
    lay.addWidget(bb)
    dlg.view = view
    if modal:
        dlg.exec()
    else:
        dlg.setAttribute(Qt.WA_DeleteOnClose)
        dlg.show()
    return dlg


def reference_items(keys: Iterable[str]) -> List[str]:
    out = []
    for k in keys:
        ref = get_reference(k)
        text = esc(format_reference(k))
        if ref.doi:
            url = f"https://doi.org/{ref.doi}"
            text = text.replace(esc(url), f"<a href='{url}'>{esc(url)}</a>")
        out.append(text)
    return out


def _range(r) -> str:
    return "not stated" if r.lo is None and r.hi is None else esc(str(r))


# ------------------------------------------------------------------ coefficient
def coefficient_html(c) -> str:
    badges = []
    if c.recommended:
        badges.append(("recommended", "ok"))
    badges.append(("verified", "accent") if c.verified else ("unverified", "warn"))
    if c.superseded_by or c.superseded_note:
        badges.append(("superseded", "danger"))
    d = Doc(c.label, esc(f"{c.species} in {c.mineral}. Source: {cite(c.citation)}"), badges)
    if c.superseded_by or c.superseded_note:
        msg = "A newer calibration replaces this one."
        if c.superseded_by:
            msg = f"Superseded by {esc(cite(c.superseded_by))}."
        if c.superseded_note:
            msg += " " + esc(c.superseded_note)
        d.box(msg, "danger")
    d.h("Equation")
    d.kv([("Transport kind", esc(c.kind)), ("Model family", esc(c.model_family)),
          ("Validation", esc(c.validation_level)), ("State variable", esc(c.transported_variable)),
          ("Reference state", esc(c.reference_state))])
    if c.calibration_notes:
        d.box("<br>".join(esc(s) for s in c.calibration_notes), "warn")
    d.equation(c.equation_text)
    d.p(f"<span class='muted'>Equation {esc(c.equation_number or 'unnumbered')} in the source. "
        f"D in m²/s.{' The fO2 term takes fO2 in ' + esc(c.fo2_unit) + '.' if c.needs_fo2 else ''}"
        f"</span>")
    if c.params:
        d.h("Parameters")
        rows = []
        for p in c.params.values():
            unc = f"± {p.sigma:g} ({'2σ' if p.sigma_level == '2s' else '1σ'})" if p.sigma else "not given"
            rows.append((esc(p.name), f"<span class='num'>{p.value:g}</span>", unc,
                         esc(p.unit), esc(p.description)))
        d.table(["Name", "Value", "Uncertainty", "Unit", "Meaning"], rows)
    d.h("Uncertainty in the Monte Carlo")
    if c.uncertainty_note:
        d.box(esc(c.uncertainty_note), "warn")
    items = []
    if c.covariance is not None:
        items.append("Correlated D0 and Q over " + esc(", ".join(c.cov_order)) + ".")
    else:
        items.append("No covariance of D0 and Q is published for this law.")
    if c.sigma_logD is not None:
        items.append(f"Scatter of log D at the working temperature: {c.sigma_logD:g} log units (1σ).")
    else:
        items.append("No scatter of log D is published. Only the individual parameters can be varied.")
    items.append("Best available: " + {"covariance": "correlated D0 and Q",
                                       "logD_at_T": "scatter of log D at T",
                                       "independent": "each parameter independently",
                                       "none": "held fixed; coefficient uncertainty not propagated"}[
        c.default_sampling_mode()] + ".")
    d.ul(items)
    d.h("Calibration range")
    d.kv([("Temperature", _range(c.T_range)), ("Pressure", _range(c.P_range)),
          ("Oxygen fugacity", _range(c.fo2_range)), ("Composition", _range(c.X_range))])
    if c.axis_factors:
        d.h("Anisotropy")
        d.ul(f"D along {esc(k)} = {v:g} × D along {esc(c.reference_axis)}"
             for k, v in c.axis_factors.items())
    if c.principal_funcs:
        d.h("Anisotropy")
        d.p("Independent Arrhenius laws for each principal direction; projected diffusivity varies with temperature.")
    if c.allowed_axes:
        d.p("Supported directions: " + esc(", ".join(c.allowed_axes)) + ". No unmeasured tensor components are inferred.")
    d.h("Verification")
    d.p(("Checked against " if c.verified else "Not yet checked against the original paper. ")
        + esc(c.verified_from or "") + ("." if c.verified_from and not c.verified_from.endswith(".") else ""))
    if c.notes:
        d.h("Notes")
        d.p(esc(c.notes))
    d.h("References")
    d.ul(reference_items([c.citation, *c.secondary_citations]))
    return d.html()


# ------------------------------------------------------------------ examples
def example_html(ds, coefficient=None) -> str:
    kind = ("measured", "ok") if ds.kind == "measured" else ("synthetic", "warn")
    d = Doc(ds.name.split(" (")[0], esc(f"{ds.mineral}, {ds.species}"), [kind])
    if ds.kind != "measured":
        d.box("Synthetic data made by Diffusor with a known answer. Not a measurement.", "warn")
    d.h("Where the data come from")
    d.p(esc(ds.provenance))
    s = ds.settings
    d.h("Conditions it sets")
    if "fo2_absolute" in s:
        fo2 = f"log fO2 {s['fo2_absolute']:g} ± {s.get('sigma_delta', 0):g} bar"
    else:
        fo2 = f"{s.get('buffer', 'NNO')} {s.get('delta_buffer', 0):+g} ± {s.get('sigma_delta', 0):g}"
    rows = [("Temperature", f"{s.get('T_C', 950):g} ± {s.get('sigma_T_K', 20):g} °C", "T"),
            ("Pressure", f"{s.get('P_MPa', 200):g} ± {s.get('sigma_P_MPa', 100):g} MPa", "P"),
            ("Oxygen fugacity", fo2, "fO2")]
    d.table(["", "Value (1σ)", "Source"],
            [(esc(k), f"<b>{esc(v)}</b>", esc(ds.sources.get(src, "Not recorded.")))
             for k, v, src in rows])
    if coefficient is not None:
        unused = unused_conditions(coefficient, s)
        if unused:
            d.p("<span class='muted'>" + esc(unused) + "</span>")
    if ds.expected:
        d.h("Expected answer")
        d.p(esc(ds.expected))
    if ds.notes:
        d.h("Notes")
        d.p(esc(ds.notes))
    if ds.citation:
        d.h("Reference")
        d.ul(reference_items([ds.citation]))
    return d.html()


def unused_condition_names(coef) -> list:
    """Which of "pressure" and "oxygen fugacity" the coefficient ignores."""
    from ..coefficients.base import Conditions
    missing = []

    def same(a, b):
        # diffusivities are around 1e-20, so compare them in log units: the default
        # absolute tolerance of np.isclose would call any two of them equal
        return abs(np.log10(a) - np.log10(b)) < 1e-9

    try:
        base = Conditions(T_K=1173.15, P_Pa=2.0e8, log_fo2_bar=-10.0,
                          X={k: 0.3 for k in coef.requires})
        d0 = float(np.mean(coef.D(base)))
        if same(float(np.mean(coef.D(base.replace(P_Pa=6.0e8)))), d0):
            missing.append("pressure")
        if not coef.needs_fo2 or same(
                float(np.mean(coef.D(base.replace(log_fo2_bar=-8.0)))), d0):
            missing.append("oxygen fugacity")
    except Exception:
        return []
    return missing


def unused_conditions(coef, settings=None) -> str:
    """Say which of pressure and fO2 the coefficient ignores."""
    missing = unused_condition_names(coef)
    if not missing:
        return ""
    who = cite(coef.citation)
    if len(missing) == 2:
        return (f"{who} has no pressure or oxygen fugacity term, so those two values do not "
                "change the result.")
    return f"{who} has no {missing[0]} term, so {missing[0]} does not change the result."


# ------------------------------------------------------------------ help
def format_html() -> str:
    from . import format_help as fh
    d = Doc("Input file format", esc(fh.SHORT))
    d.h("Example")
    lines = fh.EXAMPLE_TABLE.splitlines()
    header = lines[0].split(",")
    rows = [[f"<span class='num'>{esc(c)}</span>" for c in ln.split(",")] for ln in lines[1:]]
    d.table(header, rows)
    d.h("Rules")
    d.ul(f"<b>{esc(t)}.</b> {esc(b)}" for t, b in fh.RULES)
    d.h("Greyscale profiles")
    d.p(esc(fh.GREYSCALE))
    return d.html()


def boundaries_html() -> str:
    d = Doc("Choosing the boundaries",
            "What each end of the profile is in the crystal, and what that means for the model.")
    d.h("The four choices")
    d.table(["Choice", "Use it when", "Condition"], [
        ("<b>Plateau continues</b>",
         "The end of the traverse lies inside a zone that is still flat, far from the crystal "
         "edge and the centre. The zone boundary you model sits between two such plateaus.",
         "Concentration held at the plateau value. It stands in for a plateau of infinite "
         "length, the assumption of the closed-form solution (Crank 1975 eq. 2.14)."),
        ("<b>Crystal rim, held by the melt</b>",
         "The traverse starts at the crystal edge and the crystal sat in a large volume of "
         "melt that kept the rim composition constant.",
         "Fixed concentration, an open boundary (Costa et al. 2008, p. 555)."),
        ("<b>Crystal rim, closed</b>",
         "The rim touched a phase that takes up none of the element, or one in which it "
         "diffuses much more slowly, for example a crystal enclosed in another crystal.",
         "Zero flux. Times come out much shorter than with an open rim, by more than an "
         "order of magnitude in Costa et al. (2008) Fig. 6."),
        ("<b>Crystal centre</b>",
         "The traverse ends at the centre of a crystal that exchanged with the melt on "
         "both sides in the same way.",
         "Zero flux by symmetry. The gradient is zero at the centre of a symmetric profile, "
         "so half the crystal can be modelled (Crank 1975 section 4.3)."),
    ])
    d.h("Why fixed plateaus are so common")
    d.p("Most published profiles, including Ostorero et al. (2022) and Chamberlain et al. "
        "(2014), use the error-function solution for a step between two plateaus. That "
        "solution assumes both plateaus extend to infinity, which is the same as holding them "
        "fixed far away. Physically, the core composition is not fixed. The assumption "
        "only says that diffusion has not yet reached it.")
    d.p("While both plateaus survive in the data, the choice at the far ends makes no "
        "difference to the fitted time. When the diffusion front reaches the end of the "
        "traverse, Diffusor warns you. From then on the real geometry matters. Model the "
        "whole crystal with a rim and a centre, and pick the rim condition from the petrography.")
    d.h("Solver")
    d.p("Diffusor uses the closed-form solution when it is exact: two plateaus that continue, "
        "a sharp step, constant D and plane geometry. Anything else runs the Crank-Nicolson "
        "finite-difference solver, which handles composition-dependent D, crystal rims and "
        "centres, cylinders and spheres and cooling paths. Both give the same answer where "
        "both apply, and the test suite checks that.")
    d.h("References")
    d.ul(reference_items(["crank1975", "costa2008", "ostorero2022", "chamberlain2014"]))
    return d.html()


# ------------------------------------------------------------------ log
def log_html(entries: Sequence[Tuple[datetime, str, str]]) -> str:
    d = Doc("Log", f"{len(entries)} entries in this session")
    if not entries:
        d.p("Nothing has happened yet.", "muted")
        return d.html()
    rows = []
    for when, level, msg in entries:
        lines = str(msg).rstrip().splitlines() or [""]
        colour = {"warn": theme.WARN, "error": theme.DANGER}.get(level, theme.TEXT)
        head = f"<span style='color:{colour}'>{esc(lines[0])}</span>"
        if len(lines) > 1:
            if level == "error":
                rest = "<br>".join(esc(ln) for ln in lines[1:])
                head += f"<br><span class='mono faint'>{rest}</span>"
            else:
                head += "<ul>" + "".join(f"<li>{esc(ln.strip())}</li>" for ln in lines[1:]
                                         if ln.strip()) + "</ul>"
        rows.append((f"<span class='faint num'>{when.strftime('%H:%M:%S')}</span>", head))
    d.table(["Time", "Event"], rows)
    return d.html()


# ------------------------------------------------------------------ results
def mc_html(res, coefficient=None) -> str:
    d = Doc("Monte Carlo result", f"{res.n_draws} draws, seed {res.seed}")
    d.h("Time")
    d.kv([("Best fit", f"<b>{esc(human_time(res.t_best))}</b>"),
          ("Median", f"<b>{esc(human_time(res.median))}</b>"),
          ("68% interval", esc(f"{human_time(res.p16)} to {human_time(res.p84)}")),
          ("95% interval", esc(f"{human_time(res.p2_5)} to {human_time(res.p97_5)}")),
          ("Scatter", f"{res.sigma_log10:.3f} log10 units (1σ)")])
    d.p("Empirical percentiles are reported without assuming a distribution shape, instead of a mean and a "
        "symmetric error.", "faint")
    b = res.budget
    d.h("What was varied")
    names = {"temperature": f"Temperature, 1σ = {b.sigma_T_K:g} K",
             "fo2": (f"Oxygen fugacity, 1σ = {b.sigma_delta_buffer:g} log units on {b.buffer}, re-evaluated "
                     "at each sampled temperature" if b.fo2_mode == "buffer"
                     else f"fO2, 1σ = {b.sigma_log_fo2:g} log units"),
             "pressure": f"Pressure, 1σ = {b.sigma_P_Pa / 1e6:g} MPa",
             "diffusion_coefficient": "Diffusion coefficient",
             "measurement_noise": "Measurement noise, from the uncertainty columns",
             "distance_scale": f"Distance scale, 1σ = {b.sigma_distance_scale * 100:g}%",
             "boundary_compositions": "Plateau compositions"}
    d.ul(esc(names.get(s, s)) for s in b.active_sources())
    if res.contributions:
        d.h("What matters most")
        d.p("Spread of log10 t with one source varied at a time. The values are for ranking "
            "and do not add up to the total.", "muted")
        d.table(["Source", "σ log10 t"],
                [(esc(names.get(k, k).split(",")[0]), f"<span class='num'>{v:.3f}</span>")
                 for k, v in sorted(res.contributions.items(), key=lambda kv: -kv[1])])
    if res.n_failed:
        d.box(f"{res.n_failed} draws failed and were left out.", "warn")
    if res.warnings:
        d.h("Check before using the result")
        d.ul(esc(w) for w in dict.fromkeys(res.warnings))
    return d.html()


def compare_html(results) -> str:
    ok = {k: r for k, r in results.items() if not isinstance(r, Exception)}
    d = Doc("Coefficient comparison", f"The same profile fitted with {len(results)} coefficients")
    rows = []
    for k, r in results.items():
        if isinstance(r, Exception):
            rows.append((esc(k), "failed", "", "", esc(str(r))))
            continue
        c = r.model.coefficient
        flags = []
        if c.superseded_by or c.superseded_note:
            flags.append(badge("superseded", "danger"))
        if not c.verified:
            flags.append(badge("unverified", "warn"))
        if c.recommended:
            flags.append(badge("recommended", "ok"))
        rows.append((esc(k), f"<b>{esc(human_time(r.t_seconds))}</b>",
                     f"<span class='num'>{r.stats.reduced_chi2:.3g}</span>",
                     f"<span class='num'>{r.stats.r_squared:.4f}</span>", " ".join(flags)))
    d.table(["Coefficient", "Time", "Reduced χ²", "R²", ""], rows)
    if len(ok) > 1:
        ts = [r.t_seconds for r in ok.values()]
        d.box(f"The times span {esc(human_time(min(ts)))} to {esc(human_time(max(ts)))}, "
              f"a factor of {max(ts) / min(ts):.1f}.", "accent")
    return d.html()


def methods_html(fit_result, mc_result=None, profile=None) -> str:
    from ..dataio.export import collect_citations
    model = fit_result.model
    c = model.coefficient
    cond = model.conditions
    ok, why = model.can_use_analytical()
    d = Doc("Methods and references",
            esc(f"Diffusor {__version__}, Python {platform.python_version()}, "
                f"{datetime.now(timezone.utc).strftime('%Y-%m-%d %H:%M UTC')}"))

    d.h("Data")
    rows = []
    if profile is not None:
        rows.append(("File", esc(Path(profile.source).name if profile.source else "not saved")))
        if profile.spec:
            rows.append(("Columns", esc(profile.spec.describe())))
        rows.append(("Points", f"{len(profile)}, from {profile.x.min():.2f} to "
                               f"{profile.x.max():.2f} µm"))
    else:
        rows.append(("Points", f"{fit_result.x_data.size}"))
    d.kv(rows)
    if profile is not None and profile.notes:
        d.ul(esc(n) for n in profile.notes)

    d.h("Conditions")
    rows = [("Temperature", f"{cond.T_K - 273.15:.1f} °C ({cond.T_K:.2f} K)"),
            ("Pressure", f"{cond.P_Pa / 1e6:.1f} MPa")]
    if cond.log_fo2_bar is not None:
        f = f"log fO2 = {cond.log_fo2_bar:.3f} bar ({cond.log_fo2_Pa:.3f} Pa)"
        if mc_result is not None and mc_result.budget.fo2_mode == "buffer":
            f += f", set as {esc(mc_result.budget.buffer)} {mc_result.budget.delta_buffer:+.2f}"
        rows.append(("Oxygen fugacity", f))
    for k, v in cond.X.items():
        if np.ndim(v) == 0:
            rows.append((k, f"{float(v):.4f}"))
    if cond.axis:
        rows.append(("Traverse", f"along the {esc(cond.axis)} axis"))
    if cond.angles_deg:
        rows.append(("Traverse", f"at {esc(cond.angles_deg)} degrees to a, b and c, direction "
                                 "cosines after Costa & Chakraborty (2004)"))
    if model.history is not None and not model.history.is_isothermal:
        rows.append(("Thermal history", esc(model.history.label) + ". The integral of D(T(t)) dt "
                     "was evaluated numerically (Lasaga 1983)."))
    d.kv(rows)

    d.h("Diffusion coefficient")
    d.p(f"<b>{esc(c.label)}</b>, equation {esc(c.equation_number or 'unnumbered')}")
    d.equation(c.equation_text)
    if not c.verified:
        d.box("The numbers of this law have not been checked against the original paper.", "warn")

    d.h("Model")
    d.kv([("Geometry", esc(f"{model.geometry.kind}. {model.geometry.describe()}")),
          ("Initial profile", esc(model.initial.describe())),
          ("Left end", esc(model.bc_left.describe())),
          ("Right end", esc(model.bc_right.describe())),
          ("Solver", esc(("Closed form. " if ok else "Numerical Crank-Nicolson, "
                          f"{model.n_nodes} nodes, D at half-nodes (Dohmen et al. 2017). ") + why)),
          ("Beam", (f"model convolved with a Gaussian of σ = {model.beam_sigma_um:.2f} µm "
                    "(Ganguly et al. 1988)") if model.beam_sigma_um else "no correction")])

    d.h("Result")
    rows = [("Best-fit time", f"<b>{esc(human_time(fit_result.t_seconds))}</b> "
                              f"<span class='faint'>({fit_result.t_seconds:.6g} s)</span>"),
            ("Fit quality", esc(fit_result.stats.describe()))]
    for k, v in fit_result.free.items():
        if k != "t":
            rows.append((f"Fitted {k}", f"{v:.5g}"))
    if mc_result is not None:
        rows += [("Monte Carlo median", f"<b>{esc(human_time(mc_result.median))}</b>"),
                 ("68% interval", esc(f"{human_time(mc_result.p16)} to {human_time(mc_result.p84)}")),
                 ("95% interval", esc(f"{human_time(mc_result.p2_5)} to {human_time(mc_result.p97_5)}")),
                 ("Draws", f"{mc_result.n_draws}, seed {mc_result.seed}, varying "
                           + esc(", ".join(s.replace("_", " ").replace("fo2", "oxygen fugacity")
                                           for s in mc_result.budget.active_sources())))]
    d.kv(rows)

    warns = list(fit_result.warnings) + (list(mc_result.warnings) if mc_result else [])
    if warns:
        d.h("Caveats")
        d.ul(esc(w) for w in dict.fromkeys(warns))

    d.h("References")
    d.ul(reference_items(collect_citations(fit_result, mc_result)))
    return d.html()
