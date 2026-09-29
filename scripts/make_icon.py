"""Draw Diffusor's application icon into diffusor/gui/icons/.

A teal tile in the interface's accent colour holding a zoned crystal whose
two zones blend across a gold zone boundary, and over it the diffusion
profile across that boundary: the error-function curve every Diffusor fit
comes down to. The crystal's shading follows the same error function.
Drawn at 1024 px and reduced, so the small sizes stay smooth; the curve is
thick enough to read at 16 px.

Writes ``diffusor.png`` (256 px, the window icon) and ``diffusor.ico``
(16 to 256 px, for Windows shortcuts).
Run from the project root:  python scripts/make_icon.py
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
from PIL import Image, ImageDraw
from scipy.special import erf

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "diffusor" / "gui" / "icons"
N = 1024

TEAL_TOP, TEAL_BOTTOM = (36, 132, 127), (20, 92, 89)     # around the accent #1B6E6B
CRYSTAL_LOW, CRYSTAL_HIGH = (62, 150, 143), (140, 204, 194)
GOLD = (233, 184, 84)
WHITE = (255, 255, 255)
ICO_SIZES = [16, 20, 24, 32, 40, 48, 64, 128, 256]


def tile() -> Image.Image:
    """The rounded square with a gentle top-to-bottom gradient."""
    t = np.linspace(0, 1, N)[:, None, None]
    grad = (np.array(TEAL_TOP) * (1 - t) + np.array(TEAL_BOTTOM) * t) * np.ones((N, N, 3))
    img = Image.fromarray(grad.astype(np.uint8)).convert("RGBA")
    mask = Image.new("L", (N, N), 0)
    ImageDraw.Draw(mask).rounded_rectangle([0, 0, N - 1, N - 1], radius=int(N * 0.22), fill=255)
    img.putalpha(mask)
    return img


def hexagon(cx, cy, r):
    ang = np.radians([0, 60, 120, 180, 240, 300])
    return [(cx + r * np.cos(a), cy + r * 0.86 * np.sin(a)) for a in ang]


def draw() -> Image.Image:
    img = tile()
    cx, cy = N * 0.5, N * 0.53
    bx, width = cx, N * 0.10                  # zone boundary and diffusion half-width
    # the crystal, shaded across the boundary by the same error function as the curve
    shape = Image.new("L", (N, N), 0)
    ImageDraw.Draw(shape).polygon(hexagon(cx, cy, N * 0.40), fill=255)
    f = 0.5 * (1 + erf((np.arange(N) - bx) / width))[None, :, None]
    shade = np.array(CRYSTAL_LOW) * (1 - f) + np.array(CRYSTAL_HIGH) * f
    crystal = Image.fromarray((shade * np.ones((N, 1, 1))).astype(np.uint8)).convert("RGBA")
    crystal.putalpha(shape)
    img.alpha_composite(crystal)
    d = ImageDraw.Draw(img)
    d.line([(bx, N * 0.17), (bx, N * 0.87)], fill=GOLD, width=int(N * 0.028))
    # the diffusion profile, stamped with a round brush so the stroke stays smooth
    x = np.linspace(N * 0.13, N * 0.87, 1500)
    y = cy - N * 0.25 * erf((x - bx) / width)
    r = N * 0.037
    for px, py in zip(x, y):
        d.ellipse([px - r, py - r, px + r, py + r], fill=WHITE)
    return img


if __name__ == "__main__":
    OUT.mkdir(parents=True, exist_ok=True)
    big = draw()
    big.resize((256, 256), Image.LANCZOS).save(OUT / "diffusor.png")
    big.save(OUT / "diffusor.ico", sizes=[(s, s) for s in ICO_SIZES])
    print("written:", OUT / "diffusor.png", OUT / "diffusor.ico")
