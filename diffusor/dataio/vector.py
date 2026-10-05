"""Vector figures that open and edit cleanly in Inkscape, CorelDRAW and Illustrator.

A figure saved with matplotlib's defaults is awkward to edit afterwards. Text is
turned into outlines, so a label cannot be retyped, and every marker of a series
is a ``<use>`` copy of one shared shape, so one point cannot be recoloured or moved
without first unlinking it. The figures saved here differ in three ways.

* Text stays text (SVG ``svg.fonttype = none``, PDF and PostScript TrueType
  fonts), set in a font that a Windows machine has.
* Every data point is its own group with its own id, for example
  ``measured_point_007``, holding the marker and its error bar. One point can be
  selected, recoloured, enlarged or deleted, and the id can be searched for
  (Inkscape: Edit, Find; the XML editor shows the structure). The fit, the bands,
  the legend and each axis are named too.
* A marker is drawn as a plain path, not as a copy of a shared shape.

The figure on screen is not touched: the points are split on a copy of it.
"""
from __future__ import annotations

import pickle
import re
import xml.etree.ElementTree as ET
from pathlib import Path
from typing import Dict

import matplotlib as mpl
import numpy as np
from matplotlib.container import ErrorbarContainer
from matplotlib.lines import Line2D
from matplotlib.markers import MarkerStyle
from matplotlib.patches import PathPatch
from matplotlib.transforms import Affine2D, ScaledTranslation

VECTOR_SUFFIXES = (".svg", ".pdf", ".eps")

# Fonts a Windows machine has, then the one matplotlib ships. A text object keeps
# its font name in the file, so a viewer without that font replaces it.
VECTOR_RC: Dict[str, object] = {
    "svg.fonttype": "none",
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
    "font.sans-serif": ["Arial", "Liberation Sans", "Helvetica", "DejaVu Sans"],
}


def _slug(text: str) -> str:
    return re.sub(r"[^A-Za-z0-9]+", "-", str(text)).strip("-").lower()


def _named(label: str, fallback: str) -> str:
    return _slug(label) if label and not label.startswith("_") and _slug(label) else fallback


class _Ids:
    """Hands out ids that are unique in one figure, which an SVG needs."""

    def __init__(self):
        self.used = set()

    def get(self, base: str) -> str:
        base = base or "item"
        name, n = base, 1
        while name in self.used:
            n += 1
            name = f"{base}-{n}"
        self.used.add(name)
        return name


def _markers_only(line) -> bool:
    return (line.get_linestyle() in ("None", "none", "", " ")
            and line.get_marker() not in (None, "None", "none", "", " "))


def _color(value, fallback):
    return fallback if value is None or (isinstance(value, str) and value == "auto") else value


def _point_patches(ax, line, base: str, ids: _Ids, role: str = "marker") -> None:
    """Replace a markers-only line by one path per point, drawn the same size."""
    ms = MarkerStyle(line.get_marker(), fillstyle=line.get_fillstyle())
    shape = ms.get_path().transformed(ms.get_transform())
    filled = ms.is_filled() and line.get_fillstyle() != "none"
    colour = line.get_color()
    face = _color(line.get_markerfacecolor(), colour) if filled else "none"
    edge = _color(line.get_markeredgecolor(), colour)
    xs = np.asarray(line.get_xdata(orig=False), dtype=float)
    ys = np.asarray(line.get_ydata(orig=False), dtype=float)
    size_in = line.get_markersize() / 72.0
    digits = max(3, len(str(len(xs))))
    for i, (x, y) in enumerate(zip(xs, ys)):
        if not (np.isfinite(x) and np.isfinite(y)):
            continue
        # the marker's size is in points, so it is scaled to inches and from there
        # to pixels, and put at the point's place by a translation measured in data
        trans = (Affine2D().scale(size_in) + ax.figure.dpi_scale_trans
                 + ScaledTranslation(x, y, ax.transData))
        patch = PathPatch(shape, transform=trans, facecolor=face, edgecolor=edge,
                          linewidth=line.get_markeredgewidth(), alpha=line.get_alpha(),
                          zorder=line.get_zorder(), joinstyle="miter" if filled else "round")
        patch.set_gid(ids.get(f"{base}_point_{i + 1:0{digits}d}_{role}"))
        ax.add_artist(patch)
    line.remove()


def _bar_lines(ax, coll, base: str, ids: _Ids) -> None:
    """Replace the error bars of a series, one collection, by one line per bar."""
    segments = coll.get_segments()
    colours = coll.get_colors()
    widths = coll.get_linewidth()
    digits = max(3, len(str(len(segments))))
    for i, seg in enumerate(segments):
        line = Line2D(seg[:, 0], seg[:, 1], color=colours[i % len(colours)],
                      linewidth=widths[i % len(widths)], alpha=coll.get_alpha(),
                      zorder=coll.get_zorder(), solid_capstyle="butt")
        line.set_transform(ax.transData)
        line.set_gid(ids.get(f"{base}_point_{i + 1:0{digits}d}_errorbar"))
        ax.add_artist(line)
    coll.remove()


def _split_axes(ax, ids: _Ids) -> None:
    # (series name, part of a point) for each piece of an error-bar series
    parts: Dict[object, tuple] = {}
    for n, cont in enumerate(ax.containers):
        if not isinstance(cont, ErrorbarContainer):
            continue
        base = _named(cont.get_label(), f"series{n + 1}")
        data, caps, bars = cont.lines
        if data is not None:
            parts[data] = (base, "marker")
        for j, cap in enumerate(caps):
            parts[cap] = (base, f"cap{j + 1}")
        for bar in bars:
            parts[bar] = (base, "errorbar")
    for n, line in enumerate(list(ax.lines)):
        base, role = parts.get(line) or (_named(line.get_label(), f"line{n + 1}"), "marker")
        if _markers_only(line):
            _point_patches(ax, line, base, ids, role)
        else:
            line.set_gid(ids.get(base))
    for n, coll in enumerate(list(ax.collections)):
        if coll in parts and parts[coll][1] == "errorbar":
            _bar_lines(ax, coll, parts[coll][0], ids)
        else:
            coll.set_gid(ids.get(_named(coll.get_label(), f"collection{n + 1}")))
    legend = ax.get_legend()
    if legend is not None:
        legend.set_gid(ids.get("legend"))


def editable_copy(fig):
    """A copy of ``fig`` with each marker split into its own named path.

    A figure that cannot be copied is returned as it is, with its markers
    unsplit, rather than stopping the export."""
    try:
        clone = pickle.loads(pickle.dumps(fig))
    except Exception:
        return fig
    ids = _Ids()
    for k, ax in enumerate(clone.axes):
        ax.set_gid(ids.get(("profile-axes", "residual-axes")[k] if k < 2 else f"axes-{k + 1}"))
    for ax in clone.axes:
        _split_axes(ax, ids)
    return clone


_POINT_PART = re.compile(r"^(.+_point_\d+)_([a-z0-9]+)$")
_SVG = "http://www.w3.org/2000/svg"


def _group_points(svg_path: Path) -> None:
    """Gather the pieces of one point (marker, caps, error bar) into one SVG group.

    Matplotlib writes each piece as a group of its own, in a flat list. Selecting
    a point in an editor should take its error bar with it, so the pieces whose ids
    share the prefix ``<series>_point_<n>`` go into a group with that id."""
    for prefix, uri in (("", _SVG), ("xlink", "http://www.w3.org/1999/xlink"),
                        ("dc", "http://purl.org/dc/elements/1.1/"),
                        ("cc", "http://creativecommons.org/ns#"),
                        ("rdf", "http://www.w3.org/1999/02/22-rdf-syntax-ns#")):
        ET.register_namespace(prefix, uri)
    tree = ET.parse(svg_path)
    g_tag = "{%s}g" % _SVG
    for parent in list(tree.getroot().iter()):
        groups: Dict[str, ET.Element] = {}
        for child in list(parent):
            m = _POINT_PART.match(child.get("id", "")) if child.tag == g_tag else None
            if m is None:
                continue
            group = groups.get(m.group(1))
            if group is None:
                group = ET.Element(g_tag, {"id": m.group(1)})
                groups[m.group(1)] = group
                parent.insert(list(parent).index(child), group)
            parent.remove(child)
            group.append(child)
    tree.write(svg_path, encoding="utf-8", xml_declaration=True)


def save_vector(fig, path, individual_points: bool = True, **kwargs) -> str:
    """Write ``fig`` as an SVG, PDF or EPS file, chosen by the suffix of ``path``."""
    path = Path(path)
    if path.suffix.lower() not in VECTOR_SUFFIXES:
        raise ValueError(f"{path.suffix or 'this file type'} is not a vector format; "
                         f"use one of {', '.join(VECTOR_SUFFIXES)}")
    with mpl.rc_context(VECTOR_RC):
        work = editable_copy(fig) if individual_points else fig
        work.savefig(path, **kwargs)
    if individual_points and path.suffix.lower() == ".svg":
        _group_points(path)
    return str(path)
