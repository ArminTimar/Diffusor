"""Generate the synthetic images in examples/images/ for the profile extractor.

Both images are SYNTHETIC. They exist to try the extractor on, and to check
that what comes out of it is what went in:

* ``cpx_bse_zoned.tif``: a 16-bit BSE-like image of a zoned clinopyroxene with a
  curved Fe-rich rim, 0.05 um per pixel written as an ImageJ calibration. Grey
  values follow X_Fe = (grey - 60)/400 on the 8-bit scale (x 257 in 16 bits), so
  X_Fe 0.16 in the core and 0.26 in the rim are grey 124 and 164, across an
  error-function boundary 6 um inside the crystal face with erf half-width
  sqrt(4 D t) = 0.8 um. A dark crack, a bright oxide inclusion and an
  exsolution lamella cross the rim so the cleaning has something to remove.
* ``opx_mg_map_jet.png``: a 300 x 220 px MgO element map drawn with the 'jet'
  colour scale, with its legend (16 to 30 wt% MgO) in the lower left, a black
  crack and white label text. 0.5 um per pixel, stated here and on the image's
  scale bar (20 um = 40 px) only, as a microprobe map export would.

Run from the project root:  python scripts/make_example_images.py
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
from matplotlib import colormaps
from PIL import Image, ImageDraw
from PIL.TiffImagePlugin import ImageFileDirectory_v2
from scipy.ndimage import gaussian_filter
from scipy.special import erf

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "examples" / "images"
rng = np.random.default_rng(20260929)

PX_UM = 0.05
HALF_WIDTH_UM = 0.8
XFE_CORE, XFE_RIM = 0.16, 0.26


def bse_image():
    h, w = 520, 700
    yy, xx = np.mgrid[0:h, 0:w].astype(float)
    # a crystal face bowed outwards: distance from a large circle, in micrometres
    r = np.hypot(xx - 350, yy - 1450)
    inward_um = (1150 - r) * PX_UM                   # positive inside the crystal
    xfe = XFE_CORE + (XFE_RIM - XFE_CORE) * 0.5 * (1 - erf((inward_um - 6.0) / HALF_WIDTH_UM))
    grey8 = 60 + 400 * xfe
    grey8[inward_um < 0] = 25                        # epoxy outside the crystal
    grey8 = gaussian_filter(grey8, 0.8) + rng.normal(0, 2.5, grey8.shape)
    grey8[np.abs(yy - 0.35 * xx - 120) < 1.6] = 8    # crack
    grey8[(xx - 420) ** 2 + (yy - 330) ** 2 < 11 ** 2] = 238   # oxide inclusion
    lam = np.abs((xx - 250) * 0.97 + (yy - 250) * 0.24) < 2.2   # lamella across the rim
    grey8[lam & (inward_um > 0)] -= 30
    img16 = np.clip(grey8 * 257, 0, 65535).astype(np.uint16)
    ifd = ImageFileDirectory_v2()
    ifd[270] = "ImageJ=1.54f\nunit=micron\n"
    ifd[282] = 1.0 / PX_UM
    ifd[283] = 1.0 / PX_UM
    Image.fromarray(img16).save(OUT / "cpx_bse_zoned.tif", tiffinfo=ifd,
                                compression="tiff_adobe_deflate")


def mg_map():
    h, w = 220, 300
    yy, xx = np.mgrid[0:h, 0:w].astype(float)
    d = (xx * 0.8 - yy * 0.6 - 60) * 0.5             # um from a straight zone boundary
    mgo = 22 + 6 * erf(d / 3.0) + rng.normal(0, 0.5, (h, w))
    rgb = colormaps["jet"](np.clip((mgo - 16) / 14, 0, 1))[:, :, :3]
    rgb[np.abs(yy - 0.2 * xx - 40) < 1.2] = 0.0      # crack
    im = Image.fromarray((rgb * 255).round().astype(np.uint8))
    dr = ImageDraw.Draw(im)
    legend = (colormaps["jet"](np.linspace(0, 1, 120))[:, :3] * 255).astype(np.uint8)
    for i, c in enumerate(legend):
        dr.line([(10 + i, 196), (10 + i, 206)], fill=tuple(int(v) for v in c))
    dr.text((8, 184), "16", fill=(255, 255, 255))
    dr.text((118, 184), "30 wt% MgO", fill=(255, 255, 255))
    dr.line([(240, 204), (280, 204)], fill=(255, 255, 255), width=3)
    dr.text((244, 190), "20 um", fill=(255, 255, 255))
    im.save(OUT / "opx_mg_map_jet.png")


if __name__ == "__main__":
    OUT.mkdir(parents=True, exist_ok=True)
    bse_image()
    mg_map()
    print("written:", *sorted(p.name for p in OUT.iterdir()))
