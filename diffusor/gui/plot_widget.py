"""Matplotlib canvas: measured profile, initial condition, fit, Monte Carlo envelope."""
from __future__ import annotations

from pathlib import Path
from typing import Dict, Optional, Sequence

import numpy as np
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg
from matplotlib.backends.backend_qtagg import NavigationToolbar2QT
from matplotlib.collections import LineCollection
from matplotlib.figure import Figure
from PySide6.QtWidgets import QVBoxLayout, QWidget

from ..thermo.units import human_time
from . import theme


class SaveToolbar(NavigationToolbar2QT):
    """Matplotlib's toolbar, whose save dialog opens in the data's own folder.

    Matplotlib starts the dialog in ``savefig.directory``, which it overwrites
    with wherever the last figure went, and proposes a name taken from the
    window title. Set ``start_dir`` and ``default_name`` when data are loaded.
    """

    def __init__(self, canvas, parent=None):
        super().__init__(canvas, parent)
        self.start_dir: Optional[Path] = None
        self.default_name: Optional[str] = None

    def save_figure(self, *args):
        import matplotlib as mpl
        if self.start_dir is not None and Path(self.start_dir).is_dir():
            mpl.rcParams["savefig.directory"] = str(self.start_dir)
        if self.default_name:
            self.canvas.get_default_filename = lambda: self.default_name
        try:
            return super().save_figure(*args)
        finally:
            self.canvas.__dict__.pop("get_default_filename", None)


class ProfilePlot(QWidget):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.figure = Figure(figsize=(7, 5.5), layout="constrained",
                             facecolor=theme.SURFACE)
        self.canvas = FigureCanvasQTAgg(self.figure)
        self.toolbar = SaveToolbar(self.canvas, self)
        self.toolbar.setStyleSheet(
            f"QToolBar {{ background:{theme.SURFACE}; border:none; "
            f"border-bottom:1px solid {theme.BORDER}; padding:3px; }}")
        self.setStyleSheet(f"background:{theme.SURFACE};")
        lay = QVBoxLayout(self)
        lay.setContentsMargins(0, 0, 0, 0)
        lay.addWidget(self.toolbar)
        lay.addWidget(self.canvas)
        self.ax = None
        self.ax_res = None
        self._make_axes()

    def _make_axes(self):
        self._mc = None            # leaving the Monte Carlo view detaches its live panels
        self.figure.clear()
        gs = self.figure.add_gridspec(4, 1)
        self.ax = self.figure.add_subplot(gs[0:3, 0])
        self.ax_res = self.figure.add_subplot(gs[3, 0], sharex=self.ax)
        self.ax.set_ylabel("composition")
        self.ax_res.set_ylabel("residual")
        self.ax_res.set_xlabel("distance (um)")
        theme.apply_plot_style(self.figure, [self.ax, self.ax_res])
        self.ax_res.axhline(0.0, color=theme.BORDER_STRONG, lw=0.8)
        self.canvas.draw_idle()

    # ------------------------------------------------------------------
    def clear(self):
        self._make_axes()

    def show_data(self, x, C, sigma=None, label="measured", y_label="composition"):
        self._make_axes()
        self._plot_data(x, C, sigma, label)
        self.ax.set_ylabel(y_label)
        self.ax.legend(fontsize=9, frameon=False)
        theme.apply_plot_style(self.figure, [self.ax, self.ax_res])
        self.canvas.draw_idle()

    def _plot_data(self, x, C, sigma=None, label="measured"):
        if sigma is not None:
            self.ax.errorbar(x, C, yerr=sigma, fmt="o", ms=4, color=theme.PLOT_DATA,
                             ecolor=theme.PLOT_DATA_ERR, elinewidth=1, capsize=2, label=label, zorder=3)
        else:
            self.ax.plot(x, C, "o", ms=4.5, color=theme.PLOT_DATA, label=label, zorder=3)

    def show_fit(self, fit_result, mc_result=None, y_label="composition",
                 initial=True, title: str = ""):
        """Measured points, initial condition, best fit, Monte Carlo envelope, residuals."""
        self._make_axes()
        x, C, sig = fit_result.x_data, fit_result.C_data, fit_result.sigma
        self._plot_data(x, C, sig)

        if mc_result is not None:
            # The band is the spread of the profiles the Monte Carlo actually fitted.
            # It is narrow, because what the data fix is the diffusion length
            # sqrt(D t), and every draw is re-fitted to the same points. Temperature
            # and the diffusion coefficient barely move the curve. They convert that
            # length into a time, so their uncertainty belongs to t and is shown in
            # the histogram, not here. Drawing the profile at the ends of the time
            # interval with D held fixed would show curves that no draw ever fitted.
            self._plot_draw_envelope(mc_result)

        if initial:
            try:
                xf = np.linspace(x.min(), x.max(), 600)
                C0 = fit_result.model.initial.evaluate(xf)
                self.ax.plot(xf, C0, "--", color=theme.PLOT_INITIAL, lw=1.3,
                             label="initial condition", zorder=2)
            except Exception:
                pass

        xf = np.linspace(x.min(), x.max(), 400)
        try:
            Cf = fit_result.model.profile(fit_result.t_seconds, xf)
        except Exception:
            xf, Cf = x, fit_result.C_model
        label = f"fit: t = {human_time(fit_result.t_seconds)}"
        if mc_result is not None:
            label += (f"\n68% of t: {human_time(mc_result.p16)} to "
                      f"{human_time(mc_result.p84)}")
        self.ax.plot(xf, Cf, "-", color=theme.PLOT_MODEL, lw=2.2, label=label, zorder=4)

        self.ax.set_ylabel(y_label)
        self.ax.legend(fontsize=9, frameon=False, loc="best")
        if title:
            self.ax.set_title(title, fontsize=11, color=theme.TEXT, pad=10)

        res = C - fit_result.C_model
        if sig is not None:
            self.ax_res.errorbar(x, res, yerr=sig, fmt="o", ms=3.2, color=theme.PLOT_DATA,
                                 ecolor=theme.PLOT_DATA_ERR, elinewidth=0.8, capsize=1.5)
        else:
            self.ax_res.plot(x, res, "o", ms=3.2, color=theme.PLOT_DATA)
        self.ax_res.axhline(0.0, color=theme.PLOT_MODEL, lw=1)
        self.ax_res.set_xlabel("distance (um)")
        self.ax_res.set_ylabel("residual")
        theme.apply_plot_style(self.figure, [self.ax, self.ax_res])
        self.canvas.draw_idle()

    def _plot_draw_envelope(self, mc_result):
        xp = mc_result.x_profiles
        if xp is None:
            return
        order = np.argsort(np.asarray(xp, dtype=float))
        xs = np.asarray(xp, dtype=float)[order]
        for q_lo, q_hi, alpha, label in ((16, 84, 0.30, "68% of the re-fitted draws"),
                                         (2.5, 97.5, 0.15, "95% of the re-fitted draws")):
            env = mc_result.envelope(q_lo, q_hi)
            if env is None:
                return
            lo, hi = env[0][order], env[1][order]
            self.ax.fill_between(xs, lo, hi, color=theme.PLOT_BAND, alpha=alpha, lw=0,
                                 label=label, zorder=1)

    def show_comparison(self, results: Dict[str, object], y_label="composition"):
        """Overlay the fits from several diffusion coefficients."""
        self._make_axes()
        first = next((r for r in results.values() if not isinstance(r, Exception)), None)
        if first is None:
            self.canvas.draw_idle()
            return
        x, C, sig = first.x_data, first.C_data, first.sigma
        self._plot_data(x, C, sig)
        colors = theme.PLOT_SERIES
        xf = np.linspace(x.min(), x.max(), 400)
        for i, (key, r) in enumerate(results.items()):
            if isinstance(r, Exception):
                continue
            try:
                Cf = r.model.profile(r.t_seconds, xf)
            except Exception:
                continue
            self.ax.plot(xf, Cf, "-", lw=1.8, color=colors[i % len(colors)],
                         label=f"{key}: {human_time(r.t_seconds)}")
            self.ax_res.plot(x, C - r.C_model, "o", ms=3, color=colors[i % len(colors)])
        self.ax.set_ylabel(y_label)
        self.ax.legend(fontsize=8, frameon=False, loc="best")
        self.ax_res.axhline(0.0, color=theme.BORDER_STRONG, lw=1)
        self.ax_res.set_xlabel("distance (um)")
        self.ax_res.set_ylabel("residual")
        theme.apply_plot_style(self.figure, [self.ax, self.ax_res])
        self.canvas.draw_idle()

    def show_histogram(self, mc_result):
        """Histogram of the Monte Carlo times on a log axis."""
        self._make_axes()
        self.figure.clear()
        ax = self.figure.add_subplot(111)
        t = mc_result.times
        ax.hist(np.log10(t), bins=40, color=theme.PLOT_DATA, alpha=0.85)
        for q, c, lbl in ((mc_result.median, theme.PLOT_MODEL, "median"),
                          (mc_result.p16, theme.ACCENT, "16th / 84th"),
                          (mc_result.p84, theme.ACCENT, None)):
            ax.axvline(np.log10(q), color=c, lw=1.6, ls="--", label=lbl)
        ax.set_xlabel("log10 t (s)")
        ax.set_ylabel("draws")
        ax.set_title(f"Monte Carlo, {mc_result.n_draws} draws "
                     f"(median {human_time(mc_result.median)})", fontsize=11,
                     color=theme.TEXT, pad=10)
        ax.legend(fontsize=9, frameon=False)
        theme.apply_plot_style(self.figure, [ax])
        self.ax = ax
        self.ax_res = None
        self.canvas.draw_idle()

    # ------------------------------------------------------------------ live Monte Carlo
    def start_monte_carlo(self, x, C, sigma, y_label: str, n_total: int, sampled,
                          fixed: Dict[str, str]):
        """Lay out the live view: profile, time histogram and the sampled inputs.

        ``sampled`` names the sources being varied. ``fixed`` gives the text shown in
        a panel whose input is held constant, for example ``{"T": "850 °C"}``.
        """
        self.figure.clear()
        gs = self.figure.add_gridspec(3, 3, width_ratios=[1.0, 1.0, 0.95],
                                      height_ratios=[1.0, 1.0, 1.05])
        self.ax = self.figure.add_subplot(gs[0:2, 0:2])
        self.ax_res = None
        self._mc = {
            "x": np.asarray(x, float), "C": np.asarray(C, float),
            "sigma": None if sigma is None else np.asarray(sigma, float),
            "n_total": int(n_total), "sampled": set(sampled), "fixed": dict(fixed),
            "t": [], "T": [], "f": [], "D": [], "curves": 0, "clouds": 0,
            "ax_t": self.figure.add_subplot(gs[2, 0:2]),
            "ax_T": self.figure.add_subplot(gs[0, 2]),
            "ax_f": self.figure.add_subplot(gs[1, 2]),
            "ax_D": self.figure.add_subplot(gs[2, 2]),
            "y_label": y_label,
        }
        self._mc["order"] = np.argsort(self._mc["x"])
        ax = self.ax
        self._plot_data(x, C, sigma, label="measured")
        ax.set_ylabel(y_label)
        ax.set_xlabel("distance (um)")
        ax.plot([], [], "-", color=theme.PLOT_MODEL, alpha=0.5, lw=1, label="fit to each draw")
        if "measurement_noise" in self._mc["sampled"]:
            ax.plot([], [], "o", color=theme.PLOT_DATA_ERR, alpha=0.6, ms=3,
                    label="data as drawn with its noise")
        ax.legend(fontsize=8.5, frameon=False, loc="best")
        self._redraw_mc_panels()
        self.canvas.draw_idle()

    def add_monte_carlo_draws(self, draws, max_curves: int = 250):
        """Add a batch of draws (dicts from ``montecarlo.run(on_draw=...)``)."""
        mc = getattr(self, "_mc", None)
        if mc is None or not draws:
            return
        xs, ys, segments = [], [], []
        for dr in draws:
            mc["t"].append(dr["t"])
            mc["T"].append(dr["T_K"] - 273.15)
            mc["f"].append(np.nan if dr["log_fo2_bar"] is None else dr["log_fo2_bar"])
            mc["D"].append(dr["log10_D"])
            if mc["curves"] < max_curves:
                o = np.argsort(dr["x"])
                segments.append(np.column_stack([dr["x"][o], dr["C_model"][o]]))
                mc["curves"] += 1
            if mc["clouds"] < max_curves and "measurement_noise" in mc["sampled"]:
                xs.append(dr["x"])
                ys.append(dr["C"])
                mc["clouds"] += 1
        if segments:
            # one collection per batch draws far faster than a line per draw
            self.ax.add_collection(LineCollection(segments, colors=theme.PLOT_MODEL,
                                                  alpha=0.07, linewidths=1.0, zorder=1))
        if xs:
            self.ax.plot(np.concatenate(xs), np.concatenate(ys), "o", ms=2.2,
                         color=theme.PLOT_DATA_ERR, alpha=0.10, zorder=2, mew=0)
        self._redraw_mc_panels()
        self.canvas.draw_idle()

    def finish_monte_carlo(self, res, best_curve=None):
        mc = getattr(self, "_mc", None)
        if mc is None:
            return
        if best_curve is not None:
            xf, Cf = best_curve
            self.ax.plot(xf, Cf, "-", color=theme.PLOT_MODEL, lw=2.2, zorder=4,
                         label="best fit: " + human_time(res.t_best))
            self.ax.legend(fontsize=8.5, frameon=False, loc="best")
        mc["result"] = res
        self._redraw_mc_panels()
        self.canvas.draw_idle()

    def _panel_note(self, ax, text):
        ax.cla()
        for spine in ax.spines.values():
            spine.set_visible(False)
        ax.text(0.5, 0.5, text, ha="center", va="center", fontsize=9,
                color=theme.TEXT_MUTED, transform=ax.transAxes)
        ax.set_xticks([])
        ax.set_yticks([])

    def _redraw_mc_panels(self):
        mc = self._mc
        n = len(mc["t"])
        year = 365.25 * 86400.0
        t_yr = np.asarray(mc["t"], float) / year

        ax = mc["ax_t"]
        ax.cla()
        if n:
            lo, hi = np.log10(t_yr.min()), np.log10(t_yr.max())
            if hi - lo < 0.2:
                lo, hi = lo - 0.1, hi + 0.1
            ax.hist(t_yr, bins=np.logspace(lo, hi, 40), color=theme.PLOT_DATA, alpha=0.85)
            ax.set_xscale("log")
            p16, p50, p84 = np.percentile(t_yr, [16, 50, 84])
            ax.axvline(p50, color=theme.PLOT_MODEL, lw=1.8)
            for q in (p16, p84):
                ax.axvline(q, color=theme.PLOT_MODEL, lw=1.2, ls="--")
            ax.set_title("median " + human_time(p50 * year) + ", 68% "
                         + human_time(p16 * year) + " to " + human_time(p84 * year),
                         fontsize=9.5, color=theme.TEXT, loc="left")
        ax.set_xlabel("time (years)")
        ax.set_ylabel("draws")

        T = np.asarray(mc["T"], float)
        if "T" in mc["fixed"] or n == 0 or np.ptp(T) == 0:
            self._panel_note(mc["ax_T"], "temperature\nnot varied\n" + mc["fixed"].get("T", ""))
        else:
            mc["ax_T"].cla()
            mc["ax_T"].hist(T, bins=30, color=theme.ACCENT, alpha=0.8)
            mc["ax_T"].set_xlabel("temperature (°C)")

        f = np.asarray(mc["f"], float)
        good = np.isfinite(f)
        if "f" in mc["fixed"] or good.sum() < 2 or np.ptp(f[good]) == 0:
            self._panel_note(mc["ax_f"], "oxygen fugacity\nnot varied\n" + mc["fixed"].get("f", ""))
        else:
            ax = mc["ax_f"]
            ax.cla()
            ax.plot(T[good], f[good], "o", ms=2.5, color=theme.ACCENT, alpha=0.35, mew=0)
            ax.set_xlabel("temperature (°C)")
            ax.set_ylabel("log fO2 (bar)")

        D = np.asarray(mc["D"], float)
        D = D[np.isfinite(D)]
        if D.size < 2 or np.ptp(D) == 0:
            self._panel_note(mc["ax_D"], "diffusion coefficient\nnot varied")
        else:
            mc["ax_D"].cla()
            mc["ax_D"].hist(D, bins=30, color=theme.ACCENT, alpha=0.8)
            mc["ax_D"].set_xlabel("log10 D (m²/s)")

        state = "done" if mc.get("result") is not None else "running"
        self.ax.set_title(f"Monte Carlo, {n} of {mc['n_total']} draws ({state})",
                          fontsize=10.5, color=theme.TEXT, loc="left")
        theme.apply_plot_style(self.figure, [self.ax, mc["ax_t"]]
                               + [a for a in (mc["ax_T"], mc["ax_f"], mc["ax_D"]) if a.get_xticks().size])


class DataPreview(QWidget):
    """A small plot of the loaded profile, without the toolbar."""

    def __init__(self, parent=None):
        super().__init__(parent)
        self.figure = Figure(figsize=(5, 3), layout="constrained", facecolor=theme.SURFACE)
        self.canvas = FigureCanvasQTAgg(self.figure)
        lay = QVBoxLayout(self)
        lay.setContentsMargins(0, 0, 0, 0)
        lay.addWidget(self.canvas)
        self.clear()

    def clear(self, message: str = "Your profile will appear here."):
        self.figure.clear()
        ax = self.figure.add_subplot(111)
        ax.text(0.5, 0.5, message, ha="center", va="center", fontsize=10,
                color=theme.TEXT_FAINT, transform=ax.transAxes)
        ax.set_axis_off()
        self.canvas.draw_idle()

    def show_data(self, x, C, sigma=None, y_label="composition", title=""):
        self.figure.clear()
        ax = self.figure.add_subplot(111)
        if sigma is not None:
            ax.errorbar(x, C, yerr=sigma, fmt="o", ms=3.5, color=theme.PLOT_DATA,
                        ecolor=theme.PLOT_DATA_ERR, elinewidth=0.9, capsize=1.5)
        else:
            ax.plot(x, C, "o", ms=3.5, color=theme.PLOT_DATA)
        ax.set_xlabel("distance (um)")
        ax.set_ylabel(y_label)
        if title:
            ax.set_title(title, fontsize=10, color=theme.TEXT, loc="left")
        theme.apply_plot_style(self.figure, [ax])
        self.canvas.draw_idle()
