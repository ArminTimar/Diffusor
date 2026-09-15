"""Visual theme: a warm off-white palette shared by the interface and the plots."""
from __future__ import annotations

# --- palette ------------------------------------------------------------------
BG = "#F7F5F2"            # page background, warm off-white
SURFACE = "#FFFFFF"       # cards and inputs
SURFACE_ALT = "#FBFAF8"   # subtle alternate fill
BORDER = "#E4E0D9"
BORDER_STRONG = "#D2CCC2"
TEXT = "#22201D"
TEXT_MUTED = "#6E6960"
TEXT_FAINT = "#96907F"

ACCENT = "#1B6E6B"        # deep teal
ACCENT_HOVER = "#175D5B"
ACCENT_SOFT = "#E7F0EF"
ACCENT_TEXT = "#FFFFFF"

WARN = "#A65423"
WARN_SOFT = "#FBF0E7"
DANGER = "#9B3535"
OK = "#2F6B45"

# --- plot colours (kept in step with the interface) ---------------------------
PLOT_DATA = "#2E4B66"
PLOT_DATA_ERR = "#9FB2C6"
PLOT_MODEL = "#C2562F"
PLOT_BAND = "#C2562F"
PLOT_INITIAL = "#8A8478"
PLOT_GRID = "#E4E0D9"
PLOT_SERIES = ["#C2562F", "#1B6E6B", "#6B4FA0", "#B0891C", "#2E4B66", "#9B3535"]

FONT_STACK = '"Segoe UI", "Inter", "Helvetica Neue", Arial, sans-serif'
MONO_STACK = '"Cascadia Mono", Consolas, "SF Mono", Menlo, monospace'


def stylesheet() -> str:
    return f"""
* {{
    font-family: {FONT_STACK};
    font-size: 13px;
    color: {TEXT};
}}
QMainWindow, QDialog {{ background: {BG}; }}
QWidget#Page {{ background: {BG}; }}

/* ---- cards -------------------------------------------------------------- */
QFrame#Card {{
    background: {SURFACE};
    border: 1px solid {BORDER};
    border-radius: 10px;
}}
QFrame#Divider {{ background: {BORDER}; max-height: 1px; border: none; }}

QLabel#H1 {{ font-size: 21px; font-weight: 600; color: {TEXT}; }}
QLabel#H2 {{ font-size: 15px; font-weight: 600; color: {TEXT}; }}
QLabel#Sub {{ font-size: 12.5px; color: {TEXT_MUTED}; }}
QLabel#Hint {{ font-size: 12px; color: {TEXT_FAINT}; }}
QLabel#FieldLabel {{ font-size: 12.5px; color: {TEXT_MUTED}; }}
QLabel#Mono {{ font-family: {MONO_STACK}; font-size: 12px; color: {TEXT_MUTED}; }}
QLabel#Warn {{
    background: {WARN_SOFT}; color: {WARN}; border: 1px solid #F0DCC9;
    border-radius: 7px; padding: 9px 11px; font-size: 12.5px;
}}
QLabel#Good {{ color: {OK}; font-size: 12.5px; }}
QLabel#Danger {{ color: {DANGER}; font-size: 12.5px; }}

/* ---- buttons ------------------------------------------------------------- */
QPushButton {{
    background: {SURFACE}; border: 1px solid {BORDER_STRONG};
    border-radius: 7px; padding: 7px 15px; color: {TEXT};
}}
QPushButton:hover {{ background: {SURFACE_ALT}; border-color: {TEXT_FAINT}; }}
QPushButton:pressed {{ background: {BORDER}; }}
QPushButton:disabled {{ color: {TEXT_FAINT}; background: {SURFACE_ALT}; border-color: {BORDER}; }}
QPushButton#Primary {{
    background: {ACCENT}; border: 1px solid {ACCENT}; color: {ACCENT_TEXT}; font-weight: 600;
    padding: 8px 20px;
}}
QPushButton#Primary:hover {{ background: {ACCENT_HOVER}; border-color: {ACCENT_HOVER}; }}
QPushButton#Primary:disabled {{ background: #B9CBCA; border-color: #B9CBCA; color: #EEF4F3; }}
QPushButton#Ghost {{ background: transparent; border: none; color: {ACCENT}; padding: 5px 9px; }}
QPushButton#Ghost:hover {{ background: {ACCENT_SOFT}; border-radius: 6px; }}
QPushButton#Chip {{
    background: {SURFACE}; border: 1px solid {BORDER}; border-radius: 16px;
    padding: 5px 13px; font-size: 12px; color: {TEXT_MUTED};
}}
QPushButton#Chip:hover {{ border-color: {ACCENT}; color: {ACCENT}; }}

/* ---- inputs --------------------------------------------------------------- */
QLineEdit, QDoubleSpinBox, QSpinBox, QComboBox, QTextEdit, QPlainTextEdit {{
    background: {SURFACE}; border: 1px solid {BORDER_STRONG}; border-radius: 7px;
    padding: 6px 9px; selection-background-color: {ACCENT_SOFT}; selection-color: {TEXT};
}}
QLineEdit:focus, QDoubleSpinBox:focus, QSpinBox:focus, QComboBox:focus {{
    border-color: {ACCENT};
}}
QComboBox::drop-down {{ border: none; width: 20px; }}
QComboBox::down-arrow {{
    image: none; border-left: 4px solid transparent; border-right: 4px solid transparent;
    border-top: 5px solid {TEXT_MUTED}; margin-right: 8px;
}}
QComboBox QAbstractItemView {{
    background: {SURFACE}; border: 1px solid {BORDER_STRONG};
    selection-background-color: {ACCENT_SOFT}; selection-color: {TEXT}; outline: none;
    padding: 4px;
}}
QDoubleSpinBox::up-button, QDoubleSpinBox::down-button,
QSpinBox::up-button, QSpinBox::down-button {{ width: 16px; border: none; background: transparent; }}

QCheckBox, QRadioButton {{ spacing: 8px; color: {TEXT}; padding: 2px 0; }}
QCheckBox::indicator, QRadioButton::indicator {{
    width: 16px; height: 16px; border: 1px solid {BORDER_STRONG};
    background: {SURFACE}; border-radius: 4px;
}}
QRadioButton::indicator {{ border-radius: 9px; }}
QCheckBox::indicator:checked, QRadioButton::indicator:checked {{
    background: {ACCENT}; border-color: {ACCENT};
}}
QCheckBox::indicator:hover, QRadioButton::indicator:hover {{ border-color: {ACCENT}; }}

/* ---- lists and tables ------------------------------------------------------ */
QListWidget {{
    background: {SURFACE}; border: 1px solid {BORDER}; border-radius: 8px;
    outline: none; padding: 4px;
}}
QListWidget::item {{ padding: 8px 8px; border-radius: 6px; }}
QListWidget::item:selected {{ background: {ACCENT_SOFT}; color: {TEXT}; }}
QListWidget::item:hover {{ background: {SURFACE_ALT}; }}

/* ---- scroll areas ---------------------------------------------------------- */
QScrollArea {{ background: transparent; border: none; }}
QScrollArea > QWidget > QWidget {{ background: transparent; }}
QScrollBar:vertical {{ background: transparent; width: 11px; margin: 2px; }}
QScrollBar::handle:vertical {{ background: {BORDER_STRONG}; border-radius: 5px; min-height: 30px; }}
QScrollBar::handle:vertical:hover {{ background: {TEXT_FAINT}; }}
QScrollBar::add-line, QScrollBar::sub-line {{ height: 0; width: 0; }}
QScrollBar::add-page, QScrollBar::sub-page {{ background: transparent; }}
QScrollBar:horizontal {{ background: transparent; height: 11px; margin: 2px; }}
QScrollBar::handle:horizontal {{ background: {BORDER_STRONG}; border-radius: 5px; min-width: 30px; }}

/* ---- step rail -------------------------------------------------------------- */
QLabel#StepDot {{
    font-size: 12px; color: {TEXT_FAINT}; padding: 5px 11px;
    border-radius: 14px; background: transparent;
}}
QLabel#StepDotActive {{
    font-size: 12px; font-weight: 600; color: {ACCENT_TEXT}; padding: 5px 11px;
    border-radius: 14px; background: {ACCENT};
}}
QLabel#StepDotDone {{
    font-size: 12px; color: {ACCENT}; padding: 5px 11px;
    border-radius: 14px; background: {ACCENT_SOFT};
}}

/* ---- summary sidebar --------------------------------------------------------- */
QWidget#Summary {{ background: {SURFACE}; }}
QLabel#SummaryKey {{ font-size: 11.5px; color: {TEXT_FAINT}; }}
QLabel#SummaryVal {{ font-size: 12.5px; color: {TEXT}; }}
QLabel#SummaryHead {{
    font-size: 11px; font-weight: 600; color: {TEXT_MUTED};
    letter-spacing: 1px; padding-top: 4px;
}}

QProgressBar {{
    background: {BORDER}; border: none; border-radius: 4px; height: 7px; text-align: center;
    color: transparent;
}}
QProgressBar::chunk {{ background: {ACCENT}; border-radius: 4px; }}

QStatusBar {{ background: {BG}; color: {TEXT_MUTED}; border-top: 1px solid {BORDER}; }}
QToolTip {{
    background: {TEXT}; color: {BG}; border: none; padding: 6px 9px; border-radius: 6px;
}}
QMenuBar {{ background: {BG}; border-bottom: 1px solid {BORDER}; }}
QMenuBar::item {{ padding: 6px 11px; background: transparent; border-radius: 6px; }}
QMenuBar::item:selected {{ background: {ACCENT_SOFT}; }}
QMenu {{ background: {SURFACE}; border: 1px solid {BORDER_STRONG}; border-radius: 8px; padding: 5px; }}
QMenu::item {{ padding: 7px 22px; border-radius: 6px; }}
QMenu::item:selected {{ background: {ACCENT_SOFT}; }}
QSplitter::handle {{ background: {BORDER}; width: 1px; }}
"""


def apply_plot_style(fig, axes):
    """Match a matplotlib figure to the interface palette."""
    fig.patch.set_facecolor(SURFACE)
    for ax in axes:
        if ax is None:
            continue
        ax.set_facecolor(SURFACE)
        ax.grid(True, color=PLOT_GRID, lw=0.8, alpha=0.9)
        ax.set_axisbelow(True)
        for side in ("top", "right"):
            ax.spines[side].set_visible(False)
        for side in ("left", "bottom"):
            ax.spines[side].set_color(BORDER_STRONG)
        ax.tick_params(colors=TEXT_MUTED, labelsize=9, length=3, width=0.8)
        ax.xaxis.label.set_color(TEXT_MUTED)
        ax.yaxis.label.set_color(TEXT_MUTED)
        ax.xaxis.label.set_fontsize(10)
        ax.yaxis.label.set_fontsize(10)
