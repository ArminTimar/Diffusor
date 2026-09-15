"""Matplotlib canvas: measured profile, initial condition, fit, Monte Carlo envelope."""
from __future__ import annotations

from typing import Dict, Optional, Sequence

import numpy as np
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg
from matplotlib.backends.backend_qtagg import NavigationToolbar2QT
from matplotlib.figure import Figure
from PySide6.QtWidgets import QVBoxLayout, QWidget

from ..thermo.units import human_time


class ProfilePlot(QWidget):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.figure = Figure(figsize=(7, 5.5), layout="constrained")
        self.canvas = FigureCanvasQTAgg(self.figure)
        self.toolbar = NavigationToolbar2QT(self.canvas, self)
        lay = QVBoxLayout(self)
        lay.setContentsMargins(0, 0, 0, 0)
        lay.addWidget(self.toolbar)
        lay.addWidget(self.canvas)
        self.ax = None
        self.ax_res = None
        self._make_axes()

    def _make_axes(self):
        self.figure.clear()
        gs = self.figure.add_gridspec(4, 1)
        self.ax = self.figure.add_subplot(gs[0:3, 0])
        self.ax_res = self.figure.add_subplot(gs[3, 0], sharex=self.ax)
        self.ax.set_ylabel("composition")
        self.ax_res.set_ylabel("residual")
        self.ax_res.set_xlabel("distance (um)")
        self.ax_res.axhline(0.0, color="0.6", lw=0.8)
        self.canvas.draw_idle()

    # ------------------------------------------------------------------
    def clear(self):
        self._make_axes()

    def show_data(self, x, C, sigma=None, label="measured"):
        self._make_axes()
        self._plot_data(x, C, sigma, label)
        self.ax.legend(fontsize=8, frameon=False)
        self.canvas.draw_idle()

    def _plot_data(self, x, C, sigma=None, label="measured"):
        if sigma is not None:
            self.ax.errorbar(x, C, yerr=sigma, fmt="o", ms=4, color="#22456b",
                             ecolor="#8aa2bd", elinewidth=1, capsize=2, label=label, zorder=3)
        else:
            self.ax.plot(x, C, "o", ms=4, color="#22456b", label=label, zorder=3)

    def show_fit(self, fit_result, mc_result=None, y_label="composition",
                 initial=True, title: str = ""):
        """Measured points, initial condition, best fit, Monte Carlo envelope, residuals."""
        self._make_axes()
        x, C, sig = fit_result.x_data, fit_result.C_data, fit_result.sigma
        self._plot_data(x, C, sig)

        if mc_result is not None:
            # Band drawn at the 16th and 84th percentile *times*, under the
            # best-fit conditions. The envelope of the re-fitted Monte Carlo
            # profiles themselves is not informative: every draw is re-fitted,
            # so all of them pass through the data by construction. What the
            # reader wants to see is how different the profile would look at
            # the ends of the time interval.
            xf_b = np.linspace(x.min(), x.max(), 300)
            try:
                lo = fit_result.model.profile(mc_result.p16, xf_b)
                hi = fit_result.model.profile(mc_result.p84, xf_b)
                self.ax.fill_between(xf_b, lo, hi, color="#d1443f", alpha=0.20, lw=0,
                                     label="68% time interval", zorder=1)
                lo2 = fit_result.model.profile(mc_result.p2_5, xf_b)
                hi2 = fit_result.model.profile(mc_result.p97_5, xf_b)
                self.ax.fill_between(xf_b, lo2, hi2, color="#d1443f", alpha=0.10, lw=0,
                                     label="95% time interval", zorder=1)
            except Exception:
                pass

        if initial:
            try:
                xf = np.linspace(x.min(), x.max(), 600)
                C0 = fit_result.model.initial.evaluate(xf)
                self.ax.plot(xf, C0, "--", color="#6b6b6b", lw=1.2,
                             label="initial condition", zorder=2)
            except Exception:
                pass

        xf = np.linspace(x.min(), x.max(), 400)
        try:
            Cf = fit_result.model.profile(fit_result.t_seconds, xf)
        except Exception:
            xf, Cf = x, fit_result.C_model
        self.ax.plot(xf, Cf, "-", color="#d1443f", lw=2,
                     label=f"fit: t = {human_time(fit_result.t_seconds)}", zorder=4)

        self.ax.set_ylabel(y_label)
        self.ax.legend(fontsize=8, frameon=False, loc="best")
        if title:
            self.ax.set_title(title, fontsize=10)

        res = C - fit_result.C_model
        if sig is not None:
            self.ax_res.errorbar(x, res, yerr=sig, fmt="o", ms=3, color="#22456b",
                                 ecolor="#8aa2bd", elinewidth=0.8, capsize=1.5)
        else:
            self.ax_res.plot(x, res, "o", ms=3, color="#22456b")
        self.ax_res.axhline(0.0, color="#d1443f", lw=1)
        self.ax_res.set_xlabel("distance (um)")
        self.ax_res.set_ylabel("residual")
        self.canvas.draw_idle()

    def show_comparison(self, results: Dict[str, object], y_label="composition"):
        """Overlay the fits from several diffusion coefficients."""
        self._make_axes()
        first = next((r for r in results.values() if not isinstance(r, Exception)), None)
        if first is None:
            self.canvas.draw_idle()
            return
        x, C, sig = first.x_data, first.C_data, first.sigma
        self._plot_data(x, C, sig)
        colors = ["#d1443f", "#2f7d4f", "#7a4fa3", "#c47b16", "#1f7a8c", "#a3435f"]
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
        self.ax.legend(fontsize=7, frameon=False, loc="best")
        self.ax_res.axhline(0.0, color="0.5", lw=1)
        self.ax_res.set_xlabel("distance (um)")
        self.ax_res.set_ylabel("residual")
        self.canvas.draw_idle()

    def show_histogram(self, mc_result):
        """Histogram of the Monte Carlo times on a log axis."""
        self._make_axes()
        self.figure.clear()
        ax = self.figure.add_subplot(111)
        t = mc_result.times
        ax.hist(np.log10(t), bins=40, color="#22456b", alpha=0.8)
        for q, c, lbl in ((mc_result.median, "#d1443f", "median"),
                          (mc_result.p16, "#c47b16", "16th / 84th"),
                          (mc_result.p84, "#c47b16", None)):
            ax.axvline(np.log10(q), color=c, lw=1.6, ls="--", label=lbl)
        ax.set_xlabel("log10 t (s)")
        ax.set_ylabel("draws")
        ax.set_title(f"Monte Carlo, {mc_result.n_draws} draws "
                     f"(median {human_time(mc_result.median)})", fontsize=10)
        ax.legend(fontsize=8, frameon=False)
        self.ax = ax
        self.ax_res = None
        self.canvas.draw_idle()
