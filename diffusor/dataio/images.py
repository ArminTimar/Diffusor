"""Reading micrographs and element maps, and turning pixels into numbers.

Profiles drawn on images start here. Three jobs:

* **Read the file.** BSE images and element maps arrive as 8-bit PNG or JPEG,
  16-bit or 32-bit float TIFF (often multi-page), ENVI or other raw binary
  dumps, NumPy arrays, or text grids of counts exported by microprobe
  software. :func:`load_image` reads all of these into one float array and
  never rescales the values, so a 16-bit BSE image keeps its 0-65535 range and a
  quantitative map keeps its wt% values.
* **Find the pixel size** where the instrument wrote it: ImageJ, Zeiss SmartSEM,
  Thermo Fisher (FEI) and Tescan TIFF tags, and JEOL or Hitachi text sidecars.
  It is only ever a suggestion that the user confirms.
* **Turn a pixel into one value.** A grey image is its own value. A colour
  image is reduced by luminance, by the mean of the channels, or by a single
  channel. An element map drawn in a rainbow or other colour scale is inverted
  through its legend with :class:`ColorScale`, and pixels whose colour is not on
  the legend (black cracks, white labels, grey epoxy) are flagged instead of
  being given a value.
"""
from __future__ import annotations

import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np

# Formats Pillow reads natively; TIFF goes through tifffile first when installed.
PILLOW_SUFFIXES = (".png", ".jpg", ".jpeg", ".jpe", ".jfif", ".bmp", ".dib", ".gif", ".tif",
                   ".tiff", ".webp", ".ppm", ".pgm", ".pbm", ".pnm", ".tga", ".pcx", ".jp2",
                   ".j2k", ".sgi", ".im", ".msp", ".dds", ".icns", ".ico", ".psd")
TEXT_SUFFIXES = (".txt", ".csv", ".tsv", ".asc", ".dat", ".xyz")
RAW_SUFFIXES = (".raw", ".img", ".bin", ".bsq", ".bil", ".bip")
ALL_SUFFIXES = PILLOW_SUFFIXES + TEXT_SUFFIXES + RAW_SUFFIXES + (".hdr", ".npy")

FILE_FILTER = ("Images and maps (" + " ".join("*" + s for s in ALL_SUFFIXES) + ");;"
               "TIFF (*.tif *.tiff);;Text grids (*.txt *.csv *.tsv *.asc *.dat);;"
               "ENVI or raw binary (*.hdr *.raw *.img *.bin);;All files (*)")

VALUE_MODES = {
    "luminance": "Grey value (luminance)",
    "mean": "Mean of the colour channels",
    "channel": "One channel",
    "colour_scale": "Colour scale (legend)",
}

LENGTH_UNITS_UM = {"m": 1e6, "mm": 1e3, "um": 1.0, "µm": 1.0, "μm": 1.0, "micron": 1.0,
                   "microns": 1.0, "micrometer": 1.0, "micrometre": 1.0, "nm": 1e-3,
                   "pm": 1e-6, "cm": 1e4}


@dataclass
class LoadedImage:
    """An image as a float array, with what is known about where it came from."""
    data: np.ndarray                         # (H, W) or (H, W, C), values as stored
    source: str = ""
    format: str = ""
    channel_names: List[str] = field(default_factory=list)
    value_range: Tuple[float, float] = (0.0, 255.0)   # nominal full scale of the storage type
    pixel_size_um: Optional[float] = None
    pixel_size_source: str = ""
    notes: List[str] = field(default_factory=list)

    @property
    def height(self) -> int:
        return int(self.data.shape[0])

    @property
    def width(self) -> int:
        return int(self.data.shape[1])

    @property
    def n_channels(self) -> int:
        return 1 if self.data.ndim == 2 else int(self.data.shape[2])

    @property
    def is_rgb(self) -> bool:
        return self.n_channels == 3 and self.channel_names[:3] == ["red", "green", "blue"]

    def channel(self, i: int) -> np.ndarray:
        return self.data if self.data.ndim == 2 else self.data[:, :, i]

    def rgb01(self) -> np.ndarray:
        """The colour image scaled to 0-1, for display and colour-scale lookups."""
        if not self.is_rgb:
            raise ValueError("this image has no red, green and blue channels")
        lo, hi = self.value_range
        return np.clip((self.data - lo) / (hi - lo if hi > lo else 1.0), 0.0, 1.0)

    def describe(self) -> str:
        ch = "1 channel" if self.n_channels == 1 else f"{self.n_channels} channels ({', '.join(self.channel_names)})"
        s = f"{self.width} x {self.height} px, {ch}, {self.format}"
        if self.pixel_size_um:
            s += f", {self.pixel_size_um:.5g} um/px from {self.pixel_size_source}"
        return s


# ----------------------------------------------------------------------- loading
def load_image(path, raw: Optional[dict] = None) -> LoadedImage:
    """Read any supported image or map into a :class:`LoadedImage`.

    ``raw`` describes a headerless binary file (keys ``width``, ``height`` and
    optionally ``dtype``, ``byte_order``, ``offset``, ``bands``, ``interleave``)
    and forces the raw reader. Without it, a file no reader recognises raises a
    ``ValueError`` saying so, and the caller can ask for those numbers.
    """
    p = Path(path)
    if raw is not None:
        return load_raw(p, **raw)
    suf = p.suffix.lower()
    if suf == ".npy":
        return _from_array(np.load(p, allow_pickle=False), p, "NumPy array")
    if suf == ".hdr" or (suf in RAW_SUFFIXES + ("",) and _envi_header_for(p) is not None):
        return load_envi(p)
    if suf in TEXT_SUFFIXES:
        return load_text_grid(p)
    errors = []
    if suf in (".tif", ".tiff"):
        try:
            return _load_tifffile(p)
        except ImportError:
            pass
        except Exception as exc:              # compressed TIFF without imagecodecs, etc.
            errors.append(f"tifffile: {exc}")
    try:
        return _load_pillow(p)
    except Exception as exc:
        errors.append(f"Pillow: {exc}")
    if _envi_header_for(p) is not None:
        return load_envi(p)
    raise ValueError(f"{p.name} is not an image format Diffusor recognises ("
                     + "; ".join(errors) + "). If it is a headerless binary dump, give its "
                     "width, height and data type to read it as raw data.")


def _from_array(arr: np.ndarray, p: Path, fmt: str, channel_names=None,
                value_range=None, notes=None) -> LoadedImage:
    arr = np.asarray(arr)
    notes = list(notes or [])
    if arr.ndim == 3 and arr.shape[0] <= 16 and arr.shape[2] > 16:
        arr = np.moveaxis(arr, 0, -1)          # (pages, H, W) -> (H, W, pages)
    if arr.ndim > 3:
        arr = arr.reshape(arr.shape[0], arr.shape[1], -1)
        notes.append("extra dimensions were flattened into channels")
    if arr.ndim not in (2, 3):
        raise ValueError(f"{p.name} holds a {arr.ndim}-D array, not an image")
    if arr.ndim == 3 and arr.shape[2] == 1:
        arr = arr[:, :, 0]
    if value_range is None:
        value_range = _nominal_range(arr)
    data = arr.astype(float)
    if channel_names is None:
        n = 1 if data.ndim == 2 else data.shape[2]
        if n == 1:
            channel_names = ["value"]
        elif n in (3, 4) and arr.dtype in (np.uint8, np.uint16):
            channel_names = ["red", "green", "blue", "alpha"][:n]
        else:
            channel_names = [f"channel {i + 1}" for i in range(n)]
    if data.ndim == 3 and channel_names[-1] == "alpha":
        data = data[:, :, :-1]
        channel_names = channel_names[:-1]
        notes.append("the alpha (transparency) channel was ignored")
    bits = {np.dtype(np.uint8): "8-bit", np.dtype(np.uint16): "16-bit", np.dtype(np.int16): "16-bit signed",
            np.dtype(np.uint32): "32-bit", np.dtype(np.int32): "32-bit signed",
            np.dtype(np.float32): "32-bit float", np.dtype(np.float64): "64-bit float"}
    fmt = f"{fmt} ({bits.get(arr.dtype, str(arr.dtype))})"
    return LoadedImage(data, str(p), fmt, list(channel_names), value_range, notes=notes)


def _nominal_range(arr: np.ndarray) -> Tuple[float, float]:
    if arr.dtype == np.uint8:
        return (0.0, 255.0)
    if arr.dtype == np.uint16:
        return (0.0, 65535.0)
    if arr.dtype == bool:
        return (0.0, 1.0)
    finite = np.asarray(arr, dtype=float)
    finite = finite[np.isfinite(finite)]
    if finite.size == 0:
        return (0.0, 1.0)
    return (float(finite.min()), float(finite.max()))


def _load_pillow(p: Path) -> LoadedImage:
    from PIL import Image
    with Image.open(p) as im:
        fmt = im.format or p.suffix.upper().lstrip(".")
        pages = []
        n_frames = getattr(im, "n_frames", 1)
        notes = []
        for i in range(n_frames):
            im.seek(i)
            pages.append(_pillow_frame_to_array(im))
        tags = dict(getattr(im, "tag_v2", {}) or {})
        info = dict(im.info)
    names = None
    if len(pages) > 1 and all(pg.shape == pages[0].shape for pg in pages) and pages[0].ndim == 2:
        if fmt == "GIF":                   # an animation, not a stack of maps
            arr = pages[0]
            notes.append(f"only the first of {len(pages)} GIF frames was read")
        else:
            arr = np.stack(pages, axis=-1)
            names = [f"page {i + 1}" for i in range(len(pages))]
            notes.append(f"{len(pages)} pages read as channels")
    else:
        arr = pages[0]
        if len(pages) > 1:
            notes.append(f"only the first of {len(pages)} pages was read (pages differ in size)")
    img = _from_array(arr, p, fmt, channel_names=names, notes=notes)
    size, where = pixel_size_from_tags(tags, info)
    if size is None:
        size, where = pixel_size_from_sidecar(p)
    img.pixel_size_um, img.pixel_size_source = size, where
    return img


def _pillow_frame_to_array(im) -> np.ndarray:
    mode = im.mode
    if mode in ("P", "PA"):
        im = im.convert("RGBA" if "transparency" in im.info or mode == "PA" else "RGB")
    elif mode in ("CMYK", "YCbCr", "LAB", "HSV"):
        im = im.convert("RGB")
    elif mode == "1":
        im = im.convert("L")
    elif mode == "LA":
        im = im.convert("L")
    elif mode.startswith("I;16"):
        return np.array(im, dtype=np.uint16)
    return np.array(im)


def _load_tifffile(p: Path) -> LoadedImage:
    import tifffile                       # optional; better for scientific TIFFs
    with tifffile.TiffFile(p) as tf:
        arr = tf.asarray()
        page = tf.pages[0]
        n_pages = len(tf.pages)
        tags = {t.code: t.value for t in page.tags.values()}
        photometric = getattr(page, "photometric", None)
    names, notes = None, []
    if arr.ndim == 3 and arr.shape[-1] in (3, 4) and str(photometric).endswith("RGB"):
        names = ["red", "green", "blue", "alpha"][:arr.shape[-1]]
    elif arr.ndim == 3 and n_pages > 1:
        names = [f"page {i + 1}" for i in range(n_pages)]
        notes.append(f"{n_pages} pages read as channels")
    elif arr.ndim == 3:
        names = [f"channel {i + 1}" for i in range(min(arr.shape[0], arr.shape[-1]))]
    img = _from_array(arr, p, "TIFF", channel_names=names, notes=notes)
    try:                                  # Pillow keeps the vendor tags as raw text
        from PIL import Image
        with Image.open(p) as im:
            tags = {**tags, **dict(getattr(im, "tag_v2", {}) or {})}
    except Exception:
        pass
    size, where = pixel_size_from_tags(tags, {})
    if size is None:
        size, where = pixel_size_from_sidecar(p)
    img.pixel_size_um, img.pixel_size_source = size, where
    return img


def load_text_grid(p) -> LoadedImage:
    """A matrix of numbers in a text file, one image row per line.

    Microprobe software (JEOL, Cameca, Probe for EPMA) exports maps like this.
    Separators may be tabs, spaces, commas or semicolons (with decimal commas).
    Header and footer lines, label columns and a leading row-number column are
    dropped when they are not part of the rectangular block of numbers.
    """
    p = Path(p)
    rows = []
    for line in p.read_text(errors="ignore").splitlines():
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        if ";" in line:
            toks = [t.strip().replace(",", ".") for t in line.split(";")]
        else:
            toks = [t for t in re.split(r"[,	 ]+", line) if t]
        vals = []
        for t in toks:
            try:
                vals.append(float(t))
            except ValueError:
                vals.append(np.nan)
        if vals and np.mean(np.isfinite(vals)) > 0.5:
            rows.append(vals)
    if not rows:
        raise ValueError(f"{p.name} does not hold a 2-D grid of numbers")
    lengths = np.array([len(r) for r in rows])
    modal = int(np.bincount(lengths).argmax())
    arr = np.array([r for r in rows if len(r) == modal], dtype=float)
    notes = []
    if len(arr) < len(rows):
        notes.append(f"{len(rows) - len(arr)} lines of a different length were skipped")
    arr = arr[:, ~np.all(np.isnan(arr), axis=0)]          # label columns
    if arr.shape[1] > 2:
        first = arr[:, 0]
        if np.allclose(first, np.arange(first.size)) or np.allclose(first, np.arange(1, first.size + 1)):
            arr = arr[:, 1:]
            notes.append("a leading row-number column was dropped")
    if arr.ndim != 2 or min(arr.shape) < 2:
        raise ValueError(f"{p.name} does not hold a 2-D grid of numbers")
    img = _from_array(arr, p, "text grid", notes=notes)
    img.format = f"text grid ({arr.shape[1]} x {arr.shape[0]} values)"
    size, where = pixel_size_from_sidecar(p)
    img.pixel_size_um, img.pixel_size_source = size, where
    return img


ENVI_TYPES = {1: np.uint8, 2: np.int16, 3: np.int32, 4: np.float32, 5: np.float64,
              12: np.uint16, 13: np.uint32, 14: np.int64, 15: np.uint64}


def _envi_header_for(p: Path) -> Optional[Path]:
    for cand in (p.with_suffix(".hdr"), Path(str(p) + ".hdr")):
        if cand.exists() and cand != p:
            try:
                if cand.read_text(errors="ignore").lstrip().upper().startswith("ENVI"):
                    return cand
            except OSError:
                pass
    return None


def load_envi(p) -> LoadedImage:
    """An ENVI header (.hdr) and its binary data file."""
    p = Path(p)
    hdr = p if p.suffix.lower() == ".hdr" else _envi_header_for(p)
    if hdr is None:
        raise ValueError(f"no ENVI header found for {p.name}")
    text = hdr.read_text(errors="ignore")
    fields = {}
    for m in re.finditer(r"^\s*([\w ]+?)\s*=\s*(\{[^}]*\}|[^\n]*)", text, flags=re.M):
        fields[m.group(1).strip().lower()] = m.group(2).strip()
    data = p
    if p.suffix.lower() == ".hdr":
        stem = p.with_suffix("")
        cands = [stem] + [stem.with_suffix(s) for s in RAW_SUFFIXES + (".dat",)]
        if fields.get("data file"):
            cands.insert(0, p.parent / fields["data file"])
        data = next((c for c in cands if c.exists() and c.is_file()), None)
        if data is None:
            raise ValueError(f"the ENVI header {p.name} has no data file beside it")
    dtype = ENVI_TYPES.get(int(fields.get("data type", 4)))
    if dtype is None:
        raise ValueError(f"unsupported ENVI data type {fields.get('data type')}")
    img = load_raw(data, width=int(fields["samples"]), height=int(fields["lines"]),
                   bands=int(fields.get("bands", 1)), dtype=np.dtype(dtype).name,
                   byte_order="big" if fields.get("byte order", "0").strip() == "1" else "little",
                   offset=int(fields.get("header offset", 0)),
                   interleave=fields.get("interleave", "bsq").lower())
    img.format = "ENVI " + img.format
    names = fields.get("band names")
    if names and names.startswith("{"):
        bn = [s.strip() for s in names.strip("{}").split(",")]
        if len(bn) == img.n_channels:
            img.channel_names = bn
    size = fields.get("pixel size")
    if size:
        m = re.search(r"([\d.eE+-]+)", size)
        unit = re.search(r"units\s*=\s*(\w+)", size)
        if m and unit and unit.group(1).lower() in ("meters", "metres", "m"):
            img.pixel_size_um, img.pixel_size_source = float(m.group(1)) * 1e6, "ENVI header"
    return img


def load_raw(p, width: int, height: int, dtype: str = "uint16", byte_order: str = "little",
             offset: int = 0, bands: int = 1, interleave: str = "bsq") -> LoadedImage:
    """A headerless binary image of known size and type."""
    p = Path(p)
    dt = np.dtype(dtype).newbyteorder("<" if byte_order == "little" else ">")
    n = int(width) * int(height) * int(bands)
    buf = np.fromfile(p, dtype=dt, count=n, offset=int(offset))
    if buf.size < n:
        raise ValueError(f"{p.name} holds {buf.size} values but {width} x {height} x {bands} "
                         f"= {n} were expected")
    if bands == 1:
        arr = buf.reshape(height, width)
    elif interleave == "bil":
        arr = buf.reshape(height, bands, width).transpose(0, 2, 1)
    elif interleave == "bip":
        arr = buf.reshape(height, width, bands)
    else:
        arr = buf.reshape(bands, height, width).transpose(1, 2, 0)
    img = _from_array(arr.astype(dt.newbyteorder("=")), p, "raw binary",
                      channel_names=None if bands > 1 else ["value"])
    if bands > 1:
        img.channel_names = [f"band {i + 1}" for i in range(bands)]
    size, where = pixel_size_from_sidecar(p)
    img.pixel_size_um, img.pixel_size_source = size, where
    return img


# ------------------------------------------------------------------- pixel size
def _to_um(value: float, unit: str) -> Optional[float]:
    f = LENGTH_UNITS_UM.get(unit.strip().lower().rstrip("."))
    return None if f is None else float(value) * f


def _text(v) -> str:
    if isinstance(v, bytes):
        return v.decode("latin-1", errors="ignore")
    if isinstance(v, (tuple, list)):
        return "".join(_text(x) for x in v)
    return str(v)


def pixel_size_from_tags(tags: Dict, info: Dict) -> Tuple[Optional[float], str]:
    """Pixel size in micrometres from TIFF tags written by SEM and imaging software."""
    # Zeiss SmartSEM, tag 34118: "Image Pixel Size = 24.39 nm"
    if 34118 in tags:
        m = re.search(r"Image Pixel Size'?\s*[=:,]\s*([\d.]+)'?,?\s*'?([a-zA-Zµ]+)", _text(tags[34118]))
        if m:
            um = _to_um(float(m.group(1)), m.group(2))
            if um:
                return um, "Zeiss SmartSEM metadata"
    # Thermo Fisher / FEI, tag 34682: INI text with PixelWidth in metres
    if 34682 in tags:
        m = re.search(r"PixelWidth'?\s*[=:]\s*([\d.eE+-]+)", _text(tags[34682]))
        if m and float(m.group(1)) > 0:
            return float(m.group(1)) * 1e6, "Thermo Fisher (FEI) metadata"
    # Tescan, tag 50431: PixelSizeX in metres
    if 50431 in tags:
        m = re.search(r"PixelSizeX'?\s*[=:]\s*([\d.eE+-]+)", _text(tags[50431]))
        if m and float(m.group(1)) > 0:
            return float(m.group(1)) * 1e6, "Tescan metadata"
    # ImageJ: calibrated unit in the description, resolution = pixels per unit
    desc = _text(tags.get(270, info.get("description", "")))
    if desc.startswith("ImageJ"):
        m = re.search(r"unit=(\S+)", desc)
        xres = tags.get(282)
        if m and xres:
            try:
                ppu = float(xres[0]) / float(xres[1]) if isinstance(xres, tuple) else float(xres)
            except (TypeError, ZeroDivisionError, ValueError):
                ppu = 0.0
            unit = m.group(1).replace("\\u00B5", "µ")
            if ppu > 0:
                um = _to_um(1.0 / ppu, unit)
                if um:
                    return um, "ImageJ calibration"
    return None, ""


def pixel_size_from_sidecar(p: Path) -> Tuple[Optional[float], str]:
    """Pixel size from the text file some instruments write beside each image."""
    for cand in (p.with_suffix(".txt"), p.with_suffix(".TXT")):
        if cand == p or not cand.exists():
            continue
        try:
            text = cand.read_text(errors="ignore")
        except OSError:
            continue
        # JEOL SEM: a scale bar of $$SM_MICRON_BAR pixels labelled $$SM_MICRON_MARKER
        bar = re.search(r"\$+SM_MICRON_BAR\s+([\d.]+)", text)
        mark = re.search(r"\$+SM_MICRON_MARKER\s+([\d.]+)\s*([a-zA-Zµ]+)", text)
        if bar and mark and float(bar.group(1)) > 0:
            um = _to_um(float(mark.group(1)), mark.group(2))
            if um:
                return um / float(bar.group(1)), f"JEOL sidecar {cand.name}"
        # Hitachi: PixelSize in nanometres
        m = re.search(r"^PixelSize\s*=\s*([\d.eE+-]+)", text, flags=re.M)
        if m and float(m.group(1)) > 0:
            return float(m.group(1)) * 1e-3, f"Hitachi sidecar {cand.name}"
    return None, ""


# --------------------------------------------------------------- pixel -> value
def srgb_to_lab(rgb) -> np.ndarray:
    """CIELAB (D65) from sRGB in 0-1, so colour distances match what the eye sees."""
    c = np.asarray(rgb, dtype=float)
    lin = np.where(c <= 0.04045, c / 12.92, ((c + 0.055) / 1.055) ** 2.4)
    m = np.array([[0.4124564, 0.3575761, 0.1804375],
                  [0.2126729, 0.7151522, 0.0721750],
                  [0.0193339, 0.1191920, 0.9503041]])
    xyz = lin @ m.T / np.array([0.95047, 1.0, 1.08883])
    f = np.where(xyz > (6 / 29) ** 3, np.cbrt(xyz), xyz / (3 * (6 / 29) ** 2) + 4 / 29)
    return np.stack([116 * f[..., 1] - 16, 500 * (f[..., 0] - f[..., 1]),
                     200 * (f[..., 1] - f[..., 2])], axis=-1)


@dataclass
class ColorScale:
    """A colour legend: an ordered run of colours, each standing for a value.

    A pixel is given the value of the nearest point on the legend, measured in
    CIELAB so that "nearest" means "looks most alike", and interpolated between
    the two legend colours it falls between. Pixels farther than
    ``max_distance`` (Delta E) from every legend colour are not on the scale and
    get no value: cracks, holes, epoxy, labels and the scale bar itself.

    A legend with discrete colour steps gives stepped values: the resolution is
    one step, whatever the number of pixels averaged later.
    """
    colors: np.ndarray                 # (n, 3) sRGB in 0-1, in legend order
    values: np.ndarray                 # (n,) value of each colour
    source: str = ""
    units: str = ""
    max_distance: float = 20.0         # Delta E (CIE76)

    def __post_init__(self):
        self.colors = np.clip(np.asarray(self.colors, dtype=float), 0.0, 1.0)
        self.values = np.asarray(self.values, dtype=float)
        if self.colors.ndim != 2 or self.colors.shape[1] != 3 or len(self.colors) < 2:
            raise ValueError("a colour scale needs at least two RGB colours")
        if self.values.shape != (len(self.colors),):
            raise ValueError("a colour scale needs one value per colour")

    @classmethod
    def from_colormap(cls, name: str, vmin: float, vmax: float, n: int = 256,
                      units: str = "", log: bool = False) -> "ColorScale":
        from matplotlib import colormaps
        cmap = colormaps[name]
        t = np.linspace(0.0, 1.0, n)
        vals = np.geomspace(vmin, vmax, n) if log else vmin + t * (vmax - vmin)
        return cls(cmap(t)[:, :3], vals, f"matplotlib colour map '{name}'", units)

    @classmethod
    def from_legend(cls, image: LoadedImage, start, end, v_start: float, v_end: float,
                    width_px: int = 3, n: Optional[int] = None, units: str = "",
                    log: bool = False) -> "ColorScale":
        """Sample a legend bar drawn in the image from ``start`` to ``end`` (x, y).

        The colours are averaged over ``width_px`` pixels across the bar, which
        smooths JPEG noise. Values run linearly (or logarithmically) from
        ``v_start`` at the first point to ``v_end`` at the second.
        """
        from scipy.ndimage import map_coordinates
        rgb = image.rgb01()
        p0, p1 = np.asarray(start, float), np.asarray(end, float)
        length = float(np.hypot(*(p1 - p0)))
        if length < 2:
            raise ValueError("the legend line is shorter than two pixels")
        n = int(n or max(2, round(length) + 1))
        t = np.linspace(0.0, 1.0, n)
        d = (p1 - p0) / length
        normal = np.array([-d[1], d[0]])
        offsets = np.arange(width_px) - (width_px - 1) / 2.0
        cols = np.zeros((n, 3))
        for o in offsets:
            pts = p0[None, :] + t[:, None] * (p1 - p0)[None, :] + o * normal[None, :]
            for k in range(3):
                cols[:, k] += map_coordinates(rgb[:, :, k], [pts[:, 1], pts[:, 0]], order=1,
                                              mode="nearest")
        cols /= len(offsets)
        if log:
            if v_start <= 0 or v_end <= 0:
                raise ValueError("a logarithmic legend needs positive end values")
            vals = np.geomspace(v_start, v_end, n)
        else:
            vals = v_start + t * (v_end - v_start)
        return cls(cols, vals, f"legend drawn in {Path(image.source).name} from "
                   f"({p0[0]:.0f}, {p0[1]:.0f}) to ({p1[0]:.0f}, {p1[1]:.0f})", units)

    @classmethod
    def from_table(cls, df, units: str = "") -> "ColorScale":
        """A table with columns value, R, G, B (0-1 or 0-255)."""
        low = {str(c).strip().lower(): c for c in df.columns}
        try:
            v = np.asarray(df[low["value"]], float)
            rgb = np.column_stack([np.asarray(df[low[k]], float) for k in ("r", "g", "b")])
        except KeyError:
            raise ValueError("a colour-scale table needs columns value, R, G and B")
        if rgb.max() > 1.0:
            rgb = rgb / 255.0
        return cls(rgb, v, "colour table", units)

    def to_values(self, rgb) -> Tuple[np.ndarray, np.ndarray]:
        """Values and colour distances (Delta E) for an array of RGB colours in 0-1.

        Values of colours farther than ``max_distance`` from the legend are NaN.
        Unique colours are looked up once, so a full element map is quick.
        """
        from scipy.spatial import cKDTree
        rgb = np.asarray(rgb, dtype=float)
        shape = rgb.shape[:-1]
        flat = np.clip(rgb.reshape(-1, 3), 0.0, 1.0)
        q = np.round(flat * 255).astype(np.int32)
        packed = (q[:, 0] << 16) | (q[:, 1] << 8) | q[:, 2]
        uniq, inverse = np.unique(packed, return_inverse=True)
        u_rgb = np.column_stack([(uniq >> 16) & 255, (uniq >> 8) & 255, uniq & 255]) / 255.0
        lab = srgb_to_lab(u_rgb)
        ref = srgb_to_lab(self.colors)
        _, idx = cKDTree(ref).query(lab)
        # refine between the nearest legend colour and whichever neighbour is closer
        best_v = self.values[idx].astype(float)
        best_d = np.linalg.norm(lab - ref[idx], axis=1)
        for step in (-1, 1):
            j = np.clip(idx + step, 0, len(ref) - 1)
            seg = ref[j] - ref[idx]
            L2 = np.einsum("ij,ij->i", seg, seg)
            with np.errstate(invalid="ignore", divide="ignore"):
                t = np.clip(np.einsum("ij,ij->i", lab - ref[idx], seg) / L2, 0.0, 1.0)
            t = np.where(L2 > 0, t, 0.0)
            proj = ref[idx] + t[:, None] * seg
            d = np.linalg.norm(lab - proj, axis=1)
            better = d < best_d
            best_d = np.where(better, d, best_d)
            best_v = np.where(better, self.values[idx] + t * (self.values[j] - self.values[idx]), best_v)
        best_v = np.where(best_d <= self.max_distance, best_v, np.nan)
        return best_v[inverse].reshape(shape), best_d[inverse].reshape(shape)

    def table(self):
        import pandas as pd
        return pd.DataFrame({"Position": np.arange(len(self.values)), "Value": self.values,
                             "R": self.colors[:, 0], "G": self.colors[:, 1], "B": self.colors[:, 2]})


def value_map(image: LoadedImage, mode: str = "luminance", channel: int = 0,
              color_scale: Optional[ColorScale] = None) -> Tuple[np.ndarray, np.ndarray, str]:
    """One number per pixel, a mask of pixels with no value, and a label for the values.

    ``luminance`` uses the Rec. 601 weights (0.299 R + 0.587 G + 0.114 B), the
    same as Pillow's greyscale conversion. ``mean`` is the unweighted mean that
    ImageJ uses by default. ``channel`` picks one channel or TIFF page (NIDIS
    reads the first channel only). ``colour_scale`` inverts a legend.
    """
    if mode == "colour_scale":
        if color_scale is None:
            raise ValueError("pick or draw a colour legend first")
        vals, _ = color_scale.to_values(image.rgb01())
        label = "legend value" + (f" ({color_scale.units})" if color_scale.units else "")
        return vals, ~np.isfinite(vals), label
    if image.n_channels == 1:
        vals, label = image.data.copy(), "grey value"
    elif mode == "channel":
        vals, label = image.channel(channel).copy(), f"{image.channel_names[channel]} value"
    elif mode == "mean":
        vals, label = image.data[:, :, :3].mean(axis=2), "mean channel value"
    elif image.n_channels >= 3:
        r, g, b = (image.data[:, :, i] for i in range(3))
        vals, label = 0.299 * r + 0.587 * g + 0.114 * b, "grey value"
    else:
        vals, label = image.data.mean(axis=2), "mean channel value"
    return vals, ~np.isfinite(vals), label
