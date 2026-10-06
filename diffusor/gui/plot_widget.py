"""Matplotlib canvas: measured profile, initial condition, fit, Monte Carlo envelope."""
from __future__ import annotations

import math
from pathlib import Path
from typing import Dict, Optional, Sequence

import numpy as np
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg
from matplotlib.backends.backend_qtagg import NavigationToolbar2QT
from matplotlib.collections import LineCollection
from matplotlib.figure import Figure
from matplotlib.patches import Rectangle
from PySide6.QtCore import QPointF, QRectF, QSize, Qt, Signal
from PySide6.QtGui import QAction, QColor, QIcon, QPainter, QPen, QPixmap
from PySide6.QtWidgets import QVBoxLayout, QWidget

from ..dataio.vector import VECTOR_SUFFIXES, save_vector
from ..thermo.units import human_time
from . import theme

# The toolbar's own buttons are checkable (pan, zoom) and the cut tool is too. The
# checked one is filled and outlined, so the mode in use is always visible.
TOOLBAR_STYLE = (
    f"QToolBar {{ background:{theme.SURFACE}; border:none; "
    f"border-bottom:1px solid {theme.BORDER}; padding:3px; spacing:2px; }}"
    f"QToolButton {{ border:1px solid transparent; border-radius:5px; padding:3px; }}"
    f"QToolButton:hover {{ background:{theme.ACCENT_SOFT}; }}"
    f"QToolButton:checked, QToolButton:checked:hover {{ background:{theme.ACCENT_SELECTED}; "
    f"border:1px solid {theme.ACCENT}; }}")


def _tool_icon(kind: str, size: QSize, ratio: float) -> QIcon:
    """A toolbar icon in the style of matplotlib's own: dark line drawing, no fill.

    ``cut`` is a dashed selection box around a point with a cross; ``restore`` is an
    arrow running round counter-clockwise."""
    side = max(int(round(size.width() * ratio)), 16)
    pm = QPixmap(side, side)
    pm.fill(Qt.transparent)
    pm.setDevicePixelRatio(ratio)
    p = QPainter(pm)
    p.setRenderHint(QPainter.Antialiasing)
    k = side / ratio / 24.0
    p.scale(k, k)                       # the drawing below is on a 24 x 24 grid
    ink = QColor("#111111")
    pen = QPen(ink, 2.2)
    pen.setCapStyle(Qt.RoundCap)
    pen.setJoinStyle(Qt.RoundJoin)
    if kind == "cut":
        box = QPen(ink, 2.0, Qt.CustomDashLine)
        box.setDashPattern([2.0, 1.5])
        p.setPen(box)
        p.drawRect(QRectF(2.6, 4.6, 18.8, 14.8))
        p.setPen(pen)
        p.setBrush(ink)
        p.drawEllipse(QPointF(8.0, 12.0), 2.0, 2.0)
        p.setBrush(Qt.NoBrush)
        p.drawLine(QPointF(13.0, 8.8), QPointF(17.6, 15.2))
        p.drawLine(QPointF(17.6, 8.8), QPointF(13.0, 15.2))
    else:
        c, r = QPointF(12.0, 12.5), 7.5
        p.setPen(pen)
        # an arc from 60 degrees round to 330, counter-clockwise, ending in a filled arrow head
        p.drawArc(QRectF(c.x() - r, c.y() - r, 2 * r, 2 * r), 60 * 16, 270 * 16)
        end = math.radians(330.0)
        tx, ty = -math.sin(end), math.cos(end)          # direction of travel, y up
        nx, ny = ty, -tx
        ex, ey = c.x() + r * math.cos(end), c.y() - r * math.sin(end)    # y down on screen
        tip = QPointF(ex + 3.4 * tx, ey - 3.4 * ty)
        left = QPointF(ex - 2.6 * tx + 3.6 * nx, ey + 2.6 * ty - 3.6 * ny)
        right = QPointF(ex - 2.6 * tx - 3.6 * nx, ey + 2.6 * ty + 3.6 * ny)
        p.setPen(Qt.NoPen)
        p.setBrush(ink)
        p.drawPolygon([tip, left, right])
    p.end()
    return QIcon(pm)


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
        self.setStyleSheet(TOOLBAR_STYLE)

    def save_figure(self, *args):
        import matplotlib as mpl
        if self.start_dir is not None and Path(self.start_dir).is_dir():
            mpl.rcParams["savefig.directory"] = str(self.start_dir)
        if self.default_name:
            self.canvas.get_default_filename = lambda: self.default_name
        fig = self.canvas.figure

        def save(fname, *a, **k):
            # this stand-in goes first, so the figure saves itself through the plain method
            fig.__dict__.pop("savefig", None)
            # SVG, PDF and EPS are written for editing in Inkscape or CorelDRAW:
            # text stays text and every point is its own object
            if Path(str(fname)).suffix.lower() in VECTOR_SUFFIXES:
                return save_vector(fig, fname, **k)
            return fig.savefig(fname, *a, **k)
        fig.savefig = save
        try:
            return super().save_figure(*args)
        finally:
            self.canvas.__dict__.pop("get_default_filename", None)
            fig.__dict__.pop("savefig", None)


CUT_HINT = ("Cut points: drag a box around the points to leave out of the fit, "
            "Shift+drag to bring them back, click one point to toggle it.")
CLICK_PIXELS = 5            # a press and release closer than this is a click, not a box
PICK_PIXELS = 12            # how near a click must land to a point to pick it


def _legend(ax, fontsize: float = 9):
    """A legend on a white, outlined box.

    Without a frame the legend's sample dot looks like one more measured point."""
    leg = ax.legend(fontsize=fontsize, loc="best", frameon=True, fancybox=False,
                    framealpha=0.95, facecolor=theme.SURFACE, edgecolor=theme.BORDER_STRONG,
                    borderpad=0.6)
    leg.get_frame().set_linewidth(0.9)
    return leg


def keep_labels_in_view(canvas) -> None:
    """Draw once more when an axis label sticks out of the figure.

    Constrained layout makes room for the labels as they measure at the moment of
    drawing. If the figure's scale changes after that, which happens when the window
    meets another display scaling, the labels come out wider than the room made for
    them and the edge of the y label is cut off. A second draw lays the figure out
    again at the new scale. It is tried once per clipped state, so a label that is
    too long for the figure cannot start a loop of redraws.
    """
    state = {"again": False}

    def check(_event):
        try:
            renderer = canvas.get_renderer()
            out = False
            for ax in canvas.figure.axes:
                for label in (ax.yaxis.label, ax.xaxis.label):
                    if label.get_text():
                        box = label.get_window_extent(renderer)
                        out = out or box.x0 < -0.5 or box.y0 < -0.5
        except Exception:
            return
        if not out:
            state["again"] = False
        elif not state["again"]:
            state["again"] = True
            canvas.draw_idle()

    canvas.mpl_connect("draw_event", check)


class ProfilePlot(QWidget):
    # the new exclusion mask (one flag per point) after points were cut or restored
    points_cut = Signal(object)

    def __init__(self, parent=None):
        super().__init__(parent)
        self._pts = None           # (x, C, excluded) of every loaded point, for cutting
        self._drag = None          # (pixel x, pixel y, data x, data y) of a box being drawn
        self._band = None
        self.figure = Figure(figsize=(7, 5.5), layout="constrained",
                             facecolor=theme.SURFACE)
        self.canvas = FigureCanvasQTAgg(self.figure)
        keep_labels_in_view(self.canvas)
        self.toolbar = SaveToolbar(self.canvas, self)
        self.setStyleSheet(f"background:{theme.SURFACE};")
        lay = QVBoxLayout(self)
        lay.setContentsMargins(0, 0, 0, 0)
        lay.addWidget(self.toolbar)
        lay.addWidget(self.canvas)
        ratio = float(self.devicePixelRatioF())
        icon_size = self.toolbar.iconSize()
        self.act_cut = QAction(_tool_icon("cut", icon_size, ratio), "Cut points", self)
        self.act_cut.setCheckable(True)
        self.act_cut.setToolTip(CUT_HINT)
        self.act_cut.toggled.connect(self._cut_toggled)
        self.act_restore = QAction(_tool_icon("restore", icon_size, ratio), "Restore all", self)
        self.act_restore.setToolTip("Put every cut point back into the fit")
        self.act_restore.setEnabled(False)
        self.act_restore.triggered.connect(self._restore_all)
        # before the coordinate readout, which stretches and would push these into the overflow
        readout = next((a for a in self.toolbar.actions()
                        if self.toolbar.widgetForAction(a) is getattr(self.toolbar, "locLabel", None)),
                       None)
        self.toolbar.insertSeparator(readout)
        self.toolbar.insertAction(readout, self.act_cut)
        self.toolbar.insertAction(readout, self.act_restore)
        # one tool at a time: choosing pan or zoom ends the cut mode, and the other way round
        for key in ("pan", "zoom"):
            act = getattr(self.toolbar, "_actions", {}).get(key)
            if act is not None:
                act.toggled.connect(lambda on: on and self.act_cut.setChecked(False))
        self.canvas.mpl_connect("button_press_event", self._cut_press)
        self.canvas.mpl_connect("motion_notify_event", self._cut_motion)
        self.canvas.mpl_connect("button_release_event", self._cut_release)
        self.ax = None
        self.ax_res = None
        self._make_axes()

    def _make_axes(self):
        self._mc = None            # leaving the Monte Carlo view detaches its live panels
        self._drag = self._band = None
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
        _legend(self.ax, 9)
        theme.apply_plot_style(self.figure, [self.ax, self.ax_res])
        self.canvas.draw_idle()

    def _plot_data(self, x, C, sigma=None, label="measured"):
        if sigma is not None:
            self.ax.errorbar(x, C, yerr=sigma, fmt="o", ms=4, color=theme.PLOT_DATA,
                             ecolor=theme.PLOT_DATA_ERR, elinewidth=1, capsize=2, label=label, zorder=3)
        else:
            self.ax.plot(x, C, "o", ms=4.5, color=theme.PLOT_DATA, label=label, zorder=3)
        self._plot_cut()

    def _plot_cut(self):
        """Mark the points the user has cut, so it is clear what the fit leaves out."""
        if self._pts is None or not self._pts[2].any():
            return
        x, C, cut = self._pts
        self.ax.plot(x[cut], C[cut], "x", ms=6.5, mew=1.5, color=theme.TEXT_FAINT,
                     label="cut from the fit", zorder=2)

    # ------------------------------------------------------------------ cutting points
    def set_points(self, x, C, excluded=None):
        """Tell the plot every loaded point, so a box or click can cut them."""
        x = np.asarray(x, dtype=float)
        cut = (np.zeros(x.shape, bool) if excluded is None else np.asarray(excluded, bool))
        self._pts = (x, np.asarray(C, dtype=float), cut)
        n = int(cut.sum())
        self.act_restore.setEnabled(n > 0)
        self.act_cut.setToolTip(f"{n} point{'s' if n != 1 else ''} cut. {CUT_HINT}" if n else CUT_HINT)

    def _nav_mode(self) -> str:
        return str(getattr(self.toolbar.mode, "value", self.toolbar.mode))

    def _cut_toggled(self, on: bool):
        self.canvas.setCursor(Qt.CrossCursor if on else Qt.ArrowCursor)
        if on:
            # pan and zoom use the same mouse gestures
            mode = self._nav_mode()
            if mode == "pan/zoom":
                self.toolbar.pan()
            elif mode == "zoom rect":
                self.toolbar.zoom()
            self.toolbar.set_message(CUT_HINT)
        else:
            self._drop_band()

    def _cut_active(self) -> bool:
        return (self.act_cut.isChecked() and not self._nav_mode() and self._mc is None
                and self.ax_res is not None and self._pts is not None)

    def _drop_band(self):
        self._drag = None
        if self._band is not None:
            try:
                self._band.remove()
            except (ValueError, NotImplementedError):
                pass
            self._band = None
            self.canvas.draw_idle()

    def _cut_press(self, ev):
        if self._cut_active() and ev.button == 1 and ev.inaxes is self.ax:
            self._drag = (ev.x, ev.y, ev.xdata, ev.ydata)

    def _cut_motion(self, ev):
        if self._drag is None or ev.inaxes is not self.ax:
            return
        _, _, x0, y0 = self._drag
        if self._band is None:
            self._band = Rectangle((x0, y0), 0, 0, fill=True, alpha=0.18, lw=1.2,
                                   facecolor=theme.PLOT_MODEL, edgecolor=theme.PLOT_MODEL,
                                   zorder=6)
            self.ax.add_patch(self._band)
        self._band.set_bounds(min(x0, ev.xdata), min(y0, ev.ydata),
                              abs(ev.xdata - x0), abs(ev.ydata - y0))
        self.canvas.draw_idle()

    def _cut_release(self, ev):
        if self._drag is None:
            return
        px, py, x0, y0 = self._drag
        self._drop_band()
        x, C, cut = self._pts
        new = cut.copy()
        if np.hypot(ev.x - px, ev.y - py) < CLICK_PIXELS:
            xy = self.ax.transData.transform(np.column_stack([x, C]))
            d = np.hypot(xy[:, 0] - px, xy[:, 1] - py)
            i = int(np.argmin(d)) if d.size else -1
            if i < 0 or d[i] > PICK_PIXELS:
                return
            new[i] = not new[i]
        else:
            x1, y1 = self.ax.transData.inverted().transform((ev.x, ev.y))
            inside = ((x >= min(x0, x1)) & (x <= max(x0, x1))
                      & (C >= min(y0, y1)) & (C <= max(y0, y1)))
            new[inside] = ev.key != "shift"
        if not np.array_equal(new, cut):
            self.points_cut.emit(new)

    def _restore_all(self):
        if self._pts is not None and self._pts[2].any():
            self.points_cut.emit(np.zeros(self._pts[0].shape, bool))

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
        _legend(self.ax, 9)
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
        _legend(self.ax, 8)
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
        _legend(ax, 9)
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
        _legend(ax, 8.5)
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
            _legend(self.ax, 8.5)
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
        keep_labels_in_view(self.canvas)
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

    def show_data(self, x, C, sigma=None, y_label="composition", title="", cut=None):
        """``cut`` is an (x, C) pair of points the user left out of the fit."""
        self.figure.clear()
        ax = self.figure.add_subplot(111)
        if sigma is not None:
            ax.errorbar(x, C, yerr=sigma, fmt="o", ms=3.5, color=theme.PLOT_DATA,
                        ecolor=theme.PLOT_DATA_ERR, elinewidth=0.9, capsize=1.5)
        else:
            ax.plot(x, C, "o", ms=3.5, color=theme.PLOT_DATA)
        if cut is not None and len(cut[0]):
            ax.plot(cut[0], cut[1], "x", ms=6, mew=1.4, color=theme.TEXT_FAINT)
        ax.set_xlabel("distance (um)")
        ax.set_ylabel(y_label)
        if title:
            ax.set_title(title, fontsize=10, color=theme.TEXT, loc="left")
        theme.apply_plot_style(self.figure, [ax])
        self.canvas.draw_idle()
