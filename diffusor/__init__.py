"""Diffusor -- unified diffusion chronometry toolkit.

Sub-packages
------------
constants     physical constants (CODATA 2018, IUPAC 2021)
references    citation registry; every equation/coefficient points here
thermo        oxygen-fugacity buffers and composition/unit conversions
minerals      mineral descriptions (species, axes, composition variables)
coefficients  literature diffusion coefficients (Arrhenius laws with fO2, P, X, axis terms)
solvers       analytical (Crank 1975) and numerical (Crank-Nicolson) 1-D solvers
fitting       weighted least squares for time and Monte Carlo error propagation
io            profile loading, greyscale calibration, result export
gui           PySide6 desktop application
"""

__version__ = "0.4.0"
