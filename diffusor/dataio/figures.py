"""Publication-friendly static reports, independent of the Qt canvas."""
import numpy as np
from matplotlib.figure import Figure
from matplotlib.backends.backend_agg import FigureCanvasAgg
from ..thermo.units import human_time
from ..references import get as reference

INK, TEAL, GOLD, MUTED = "#243746", "#007f82", "#dc994b", "#75848b"


def _style(ax):
    ax.set_facecolor("#ffffff")
    ax.spines[["top", "right"]].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_color("#ccd7d9")
    ax.tick_params(colors=INK, labelsize=9)
    ax.grid(axis="y", color="#edf1f2", linewidth=.7)
    ax.set_axisbelow(True)


def _profile(ax, residual_ax, fit, mc=None, name=None):
    order = np.argsort(fit.x_data)
    x, y = fit.x_data[order], fit.C_data[order]
    model = fit.C_model[order]
    sigma = None if fit.sigma is None else np.broadcast_to(fit.sigma, y.shape)[order]
    _style(ax)
    _style(residual_ax)
    if mc is not None and mc.profiles is not None:
        lo, hi = mc.envelope()
        ix = np.argsort(mc.x_profiles)
        ax.fill_between(np.asarray(mc.x_profiles)[ix], lo[ix], hi[ix],
                        color=TEAL, alpha=.16, linewidth=0, label="68% fitted-profile band")
    xi = np.linspace(x.min(), x.max(), 501)
    if fit.model.initial.kind == "step":
        x0 = fit.model.initial.params["x0"]
        if x.min() <= x0 <= x.max():
            xi = np.sort(np.r_[xi, np.nextafter(x0, -np.inf), x0, np.nextafter(x0, np.inf)])
    ax.plot(xi, fit.model.initial.evaluate(xi), color=MUTED, linestyle="--", linewidth=1.1, label="Initial state")
    ax.plot(x, model, color=TEAL, linewidth=2.6, label="Fitted model")
    ax.errorbar(x, y, yerr=sigma, fmt="o", markersize=4, color=INK,
                markerfacecolor="white", markeredgewidth=1., ecolor="#aab8bf",
                elinewidth=.8, capsize=1.5, label="Measurements")
    ax.set_title(name or f"{fit.model.coefficient.mineral.capitalize()} · {fit.model.coefficient.species}",
                 loc="left", fontsize=12, color=INK, fontweight="bold", pad=12)
    ax.set_ylabel("Concentration (input units)", color=INK)
    # outlined, so the legend's sample dot is not taken for a measured point
    ax.legend(loc="best", frameon=True, fancybox=False, framealpha=.95, facecolor="white",
              edgecolor="#ccd7d9", fontsize=8)
    if sigma is None:
        r, label = y-model, "Residual"
    else:
        r, label = (y-model)/sigma, "Residual / σ"
        residual_ax.axhspan(-2, 2, color="#edf5f5", zorder=0)
    residual_ax.axhline(0, color=MUTED, linewidth=.8)
    residual_ax.plot(x, r, "o", color=GOLD, markersize=3.6)
    residual_ax.set_ylabel(label, color=INK)
    residual_ax.set_xlabel("Distance (µm)", color=INK)
    ax.tick_params(labelbottom=False)


def report_figure(fit, mc=None):
    """Measured data, model, initial state, residuals and uncertainty summary."""
    fig = Figure(figsize=(10, 8), facecolor="#fafcfb", layout="constrained")
    FigureCanvasAgg(fig)
    grid = fig.add_gridspec(4, 1, height_ratios=(.5, 3.4, 1.1, .55))
    head = fig.add_subplot(grid[0]); head.axis("off")
    head.text(0, .8, "DIFFUSOR  /  PROFILE STUDY", color=TEAL, fontsize=11, weight="bold")
    summary = f"Best-fit duration  {human_time(fit.t_seconds)}"
    if mc is not None:
        summary += f"    |    MC median {human_time(mc.median)} · 68% {human_time(mc.p16)}–{human_time(mc.p84)}"
    head.text(0, .12, summary, fontsize=11, color=INK)
    ax = fig.add_subplot(grid[1])
    residual = fig.add_subplot(grid[2], sharex=ax)
    _profile(ax, residual, fit, mc)
    foot = fig.add_subplot(grid[3]); foot.axis("off")
    c = fit.model.coefficient
    foot.text(0, .8, f"{reference(c.citation).short()}  |  {c.key}", fontsize=8, color=MUTED)
    history = fit.model.history
    temperature = (f"{fit.model.conditions.T_K-273.15:g} °C" if history is None else
                   f"Prescribed history {np.min(history.temps_K)-273.15:g}–{np.max(history.temps_K)-273.15:g} °C")
    foot.text(0, .37, f"{temperature} · {fit.model.conditions.P_Pa/1e6:g} MPa · "
              f"{fit.model.geometry.kind} · RMSE {fit.stats.rmse:.3g}", fontsize=8, color=MUTED)
    foot.text(0, -.05, "Calibration limits, uncertainty assumptions and interpretation caveats accompany this figure in the report.", fontsize=8, color=MUTED)
    return fig


def joint_figure(result):
    """One pair of axes per clock, with an explicit shared-duration caption."""
    n = len(result.profiles)
    fig = Figure(figsize=(10, 2.9*n+1.1), facecolor="#fafcfb", layout="constrained")
    FigureCanvasAgg(fig)
    fig.suptitle(f"DIFFUSOR  /  SHARED DURATION\n{human_time(result.t_seconds)} · "
                 f"χ² / dof = {result.reduced_chi2:.3g}", x=.06, ha="left", color=INK, fontsize=14)
    grid = fig.add_gridspec(2*n, 1, height_ratios=[3, 1]*n)
    for i, (name, fit) in enumerate(zip(result.names, result.profiles)):
        ax = fig.add_subplot(grid[2*i])
        residual = fig.add_subplot(grid[2*i+1], sharex=ax)
        _profile(ax, residual, fit, name=name)
    return fig
