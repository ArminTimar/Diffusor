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


# validated categorical order (colour-vision-deficiency separation checked); a group
# also gets its own marker, so identity never rests on colour alone
GROUP_COLOURS = ("#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#4a3aa7", "#e34948")
GROUP_MARKERS = ("o", "s", "D", "^", "v", "P", "X", "h")
POPULATION_PANELS = ("ranked", "density", "cumulative")


def population_figure(records, panels=POPULATION_PANELS, unit=None):
    """Ranked times, summed densities and cumulative curves on one shared time axis."""
    from matplotlib.lines import Line2D
    from matplotlib.ticker import FuncFormatter
    from ..population import (TIME_UNITS, UNIT_NAMES, auto_unit, group_density,
                              group_name, groups_of, log_grid)
    if not records:
        raise ValueError("the population is empty")
    panels = [p for p in POPULATION_PANELS if p in panels] or ["ranked"]
    groups = groups_of(records)
    if len(groups) > len(GROUP_COLOURS):
        raise ValueError(f"{len(groups)} groups; a figure can tell at most {len(GROUP_COLOURS)} "
                         "apart. Merge groups or plot them in separate figures.")
    unit = unit or auto_unit(records)
    f = TIME_UNITS[unit]
    colour = {g: GROUP_COLOURS[i] for i, g in enumerate(groups)}
    marker = {g: GROUP_MARKERS[i] for i, g in enumerate(groups)}
    members = {g: [r for r in records if r.group == g] for g in groups}
    named = len(groups) > 1 or groups[0] != ""

    heights = {"ranked": max(2.6, min(6.5, 0.16 * len(records) + 1.4)), "density": 2.4, "cumulative": 2.2}
    fig = Figure(figsize=(9, sum(heights[p] for p in panels) + 1.0), facecolor="#fafcfb",
                 layout="constrained")
    FigureCanvasAgg(fig)
    t_all = np.array([r.t_best for r in records])
    fig.suptitle(f"DIFFUSOR  /  TIMESCALE POPULATION\n{len(records)} timescales · "
                 f"median {human_time(float(np.median(t_all)))}", x=.06, ha="left", color=INK, fontsize=14)
    grid = fig.add_gridspec(len(panels), 1, height_ratios=[heights[p] for p in panels])
    lg = log_grid(records)
    axes, first = {}, None
    for i, p in enumerate(panels):
        ax = fig.add_subplot(grid[i], sharex=first)
        first = first or ax
        _style(ax)
        ax.set_xscale("log")
        ax.grid(axis="x", which="major", color="#edf1f2", linewidth=.7)
        axes[p] = ax
    first.set_xlim(10 ** lg[0] / f, 10 ** lg[-1] / f)

    def title(ax, text):
        ax.set_title(text, loc="left", fontsize=11, color=INK, fontweight="bold", pad=8)

    if "ranked" in axes:
        ax = axes["ranked"]
        title(ax, "Ranked timescales")
        y, ticks, labels = 0, [], []
        for g in groups:
            for r in sorted(members[g], key=lambda r: r.t_best):
                c = colour[g]
                if r.p2_5 and r.p97_5:
                    ax.plot([r.p2_5 / f, r.p97_5 / f], [y, y], color=c, linewidth=.9, alpha=.55,
                            solid_capstyle="round")
                if r.has_interval:
                    ax.plot([r.p16 / f, r.p84 / f], [y, y], color=c, linewidth=2.4, solid_capstyle="round")
                ax.plot(r.t_best / f, y, marker[g], color=c, markersize=6.5,
                        markeredgecolor="white", markeredgewidth=.9, zorder=3)
                ticks.append(y)
                labels.append(r.label)
                y += 1
            y += 1 if len(groups) > 1 else 0      # a gap between groups
        ax.set_ylim(-1, y)
        if len(records) <= 30:
            ax.set_yticks(ticks, labels, fontsize=7)
        else:
            ax.set_yticks([])
            ax.set_ylabel("Crystals, shortest at the bottom", color=INK)
        ax.grid(axis="y", visible=False)

    if "density" in axes:
        ax = axes["density"]
        title(ax, "Summed probability density")
        x = 10 ** lg / f
        drawn = 0
        for g in groups:
            d, n = group_density(members[g], lg)
            if d is None:
                continue
            ax.fill_between(x, d, color=colour[g], alpha=.10, linewidth=0)
            ax.plot(x, d, color=colour[g], linewidth=2)
            drawn += 1
        if not drawn:
            ax.text(.5, .5, "No timescale has an uncertainty, so there is no density to draw.",
                    transform=ax.transAxes, ha="center", va="center", color=MUTED, fontsize=9)
        ax.set_ylim(bottom=0)
        ax.set_ylabel("Density per log$_{10}$ t", color=INK)

    if "cumulative" in axes:
        ax = axes["cumulative"]
        title(ax, "Cumulative distribution of best-fit times")
        for g in groups:
            t = np.sort([r.t_best for r in members[g]]) / f
            frac = np.arange(1, t.size + 1) / t.size
            ax.step(np.r_[t[0], t], np.r_[0, frac], where="post", color=colour[g], linewidth=2)
        ax.set_ylim(0, 1.03)
        ax.set_ylabel("Fraction of crystals", color=INK)

    fmt = FuncFormatter(lambda v, _: f"{v:g}")
    for p, ax in axes.items():
        ax.xaxis.set_major_formatter(fmt)
        if p != panels[-1]:
            ax.tick_params(labelbottom=False)
    axes[panels[-1]].set_xlabel(f"Time ({UNIT_NAMES[unit]})", color=INK)
    if named:
        handles = [Line2D([], [], color=colour[g], marker=marker[g], markersize=6.5, linewidth=2,
                          markeredgecolor="white", label=f"{group_name(g)} (n = {len(members[g])})")
                   for g in groups]
        axes[panels[0]].legend(handles=handles, loc="best", frameon=True, fancybox=False,
                               framealpha=.95, facecolor="white", edgecolor="#ccd7d9", fontsize=8)
    return fig
