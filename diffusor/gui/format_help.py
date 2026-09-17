"""What an input file has to look like, in one place."""
from __future__ import annotations

SHORT = "CSV, TSV or Excel. A distance column and one or two composition columns."

FILE_TYPES = ".csv, .txt, .tsv or .xlsx"

EXAMPLE_TABLE = """Distance_um,FeO_wt,MgO_wt,FeO_err,MgO_err
0.0,20.9,22.4,0.15,0.20
1.7,20.8,22.5,0.15,0.20
3.3,21.0,22.3,0.15,0.20"""

RULES = [
    ("Distance",
     "Any unit. It can start anywhere and run in either direction. Spacing can be uneven."),
    ("Composition",
     "Two columns are modelled as the molar ratio A/(A+B), for example FeO and MgO giving "
     "X_Fe. One column is modelled as it stands, for example Sr in ppm or Fo in mol%."),
    ("Oxides",
     "Columns named after an oxide (FeO, MgO_wt, CaO (wt%)) are recognised and converted "
     "to cation moles before the ratio is taken."),
    ("Uncertainties",
     "Optional. Columns named like FeO_err or Sr_2s are matched to their values. Without "
     "them every point has the same weight."),
    ("Anything else",
     "Extra columns are ignored. Keep an anorthite column for plagioclase, it can drive "
     "the equilibrium initial profile."),
    ("Blank cells",
     "Rows with missing values are skipped and the log says how many."),
]

GREYSCALE = (
    "BSE grey values need at least two microprobe anchor points to convert grey value "
    "to composition.")


def as_plain_text() -> str:
    lines = ["Example", "", EXAMPLE_TABLE, ""]
    for title, body in RULES:
        lines += [title, f"    {body}", ""]
    lines += ["Greyscale profiles", f"    {GREYSCALE}", ""]
    return "\n".join(lines)
