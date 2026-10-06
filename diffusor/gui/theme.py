"""Visual theme: a warm off-white palette shared by the interface and the plots."""
from __future__ import annotations

from pathlib import Path

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
ACCENT_SOFT = "#E7F0EF"   # hover fill
ACCENT_SELECTED = "#CFE4E1"  # selected fill, one step stronger than hover
ACCENT_TEXT = "#FFFFFF"

WARN = "#A65423"
WARN_SOFT = "#FBF0E7"
DANGER = "#9B3535"
OK = "#2F6B45"
OK_SOFT = "#E8F2EC"

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

ICONS = Path(__file__).resolve().parent / "icons"


def _icon(name: str) -> str:
    return (ICONS / name).as_posix()


def stylesheet() -> str:
    check = _icon("check.svg")
    dot = _icon("dot.svg")
    down = _icon("chevron_down.svg")
    up = _icon("chevron_up.svg")
    right = _icon("chevron_right.svg")
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

QLabel#H1 {{ font-size: 20px; font-weight: 600; color: {TEXT}; }}
QLabel#H2 {{ font-size: 14.5px; font-weight: 600; color: {TEXT}; }}
QLabel#Sub {{ font-size: 12.5px; color: {TEXT_MUTED}; }}
QLabel#Hint {{ font-size: 12px; color: {TEXT_MUTED}; }}
QLabel#FieldLabel {{ font-size: 12.5px; color: {TEXT_MUTED}; }}
QLabel#Mono {{ font-family: {MONO_STACK}; font-size: 12px; color: {TEXT_MUTED}; }}
QLabel#Equation {{ color: {TEXT}; }}
QLabel#Warn {{
    background: {WARN_SOFT}; color: {WARN}; border: 1px solid #F0DCC9;
    border-radius: 7px; padding: 8px 10px; font-size: 12.5px;
}}
QLabel#Info {{
    background: {ACCENT_SOFT}; color: {TEXT}; border: 1px solid #D3E4E2;
    border-radius: 7px; padding: 8px 10px; font-size: 12.5px;
}}
QFrame#InfoBox {{ background: {ACCENT_SOFT}; border: 1px solid #D3E4E2; border-radius: 7px; }}
QFrame#WarnBox {{ background: {WARN_SOFT}; border: 1px solid #F0DCC9; border-radius: 7px; }}
QLabel#CalloutText {{ font-size: 12.5px; color: {TEXT}; background: transparent; border: none; }}
QFrame#WarnBox QLabel#CalloutText {{ color: {WARN}; }}
QLabel#Good {{ color: {OK}; font-size: 12.5px; }}
QLabel#Danger {{ color: {DANGER}; font-size: 12.5px; }}

/* ---- buttons: hover always keeps dark text on a light fill ---------------- */
QPushButton {{
    background: {SURFACE}; border: 1px solid {BORDER_STRONG};
    border-radius: 7px; padding: 6px 14px; color: {TEXT};
}}
QPushButton:hover {{ background: {ACCENT_SOFT}; border-color: {ACCENT}; color: {TEXT}; }}
QPushButton:pressed {{ background: {ACCENT_SELECTED}; color: {TEXT}; }}
QPushButton:disabled {{ color: {TEXT_FAINT}; background: {SURFACE_ALT}; border-color: {BORDER}; }}
QPushButton#Primary {{
    background: {ACCENT}; border: 1px solid {ACCENT}; color: {ACCENT_TEXT}; font-weight: 600;
    padding: 7px 20px;
}}
QPushButton#Primary:hover {{ background: {ACCENT_SELECTED}; border-color: {ACCENT}; color: {TEXT}; }}
QPushButton#Primary:pressed {{ background: {ACCENT_SOFT}; color: {TEXT}; }}
QPushButton#Primary:disabled {{ background: #B9CBCA; border-color: #B9CBCA; color: #EEF4F3; }}
QPushButton#Ghost {{ background: transparent; border: none; color: {ACCENT}; padding: 4px 8px; }}
QPushButton#Ghost:hover {{ background: {ACCENT_SOFT}; border-radius: 6px; color: {TEXT}; }}
QPushButton#Ghost:checked {{ background: {ACCENT_SELECTED}; border-radius: 6px; color: {TEXT}; }}
QPushButton#Ghost:disabled {{ color: {TEXT_FAINT}; background: transparent; }}
QToolButton {{ background: transparent; border: none; border-radius: 6px; padding: 3px; color: {TEXT}; }}
QToolButton:hover {{ background: {ACCENT_SOFT}; color: {TEXT}; }}
QToolButton:checked {{ background: {ACCENT_SELECTED}; }}

/* ---- inputs --------------------------------------------------------------- */
QLineEdit, QDoubleSpinBox, QSpinBox, QComboBox {{
    background: {SURFACE}; border: 1px solid {BORDER_STRONG}; border-radius: 7px;
    padding: 5px 9px; selection-background-color: {ACCENT_SELECTED}; selection-color: {TEXT};
    min-height: 20px;
}}
QLineEdit:hover, QDoubleSpinBox:hover, QSpinBox:hover, QComboBox:hover {{ border-color: {TEXT_FAINT}; }}
QLineEdit:focus, QDoubleSpinBox:focus, QSpinBox:focus, QComboBox:focus {{
    border-color: {ACCENT};
}}
QLineEdit:disabled, QDoubleSpinBox:disabled, QSpinBox:disabled, QComboBox:disabled {{
    color: {TEXT_FAINT}; background: {SURFACE_ALT}; border-color: {BORDER};
}}
QDoubleSpinBox[readOnly="true"] {{ background: {SURFACE_ALT}; }}
QComboBox::drop-down {{ border: none; width: 24px; }}
QComboBox::down-arrow {{ image: url("{down}"); width: 12px; height: 12px; margin-right: 6px; }}
QComboBox QAbstractItemView {{
    background: {SURFACE}; border: 1px solid {BORDER_STRONG};
    selection-background-color: {ACCENT_SELECTED}; selection-color: {TEXT}; outline: none;
    padding: 4px;
}}
QComboBox QAbstractItemView::item {{ padding: 5px 8px; min-height: 22px; }}
QComboBox QAbstractItemView::item:hover {{ background: {ACCENT_SOFT}; color: {TEXT}; }}
QComboBox QAbstractItemView::item:selected {{ background: {ACCENT_SELECTED}; color: {TEXT}; }}
QDoubleSpinBox::up-button, QSpinBox::up-button {{
    subcontrol-origin: border; subcontrol-position: top right;
    width: 18px; border: none; background: transparent; margin-top: 2px;
}}
QDoubleSpinBox::down-button, QSpinBox::down-button {{
    subcontrol-origin: border; subcontrol-position: bottom right;
    width: 18px; border: none; background: transparent; margin-bottom: 2px;
}}
QDoubleSpinBox::up-button:hover, QSpinBox::up-button:hover,
QDoubleSpinBox::down-button:hover, QSpinBox::down-button:hover {{ background: {ACCENT_SOFT}; }}
QDoubleSpinBox::up-arrow, QSpinBox::up-arrow {{ image: url("{up}"); width: 9px; height: 9px; }}
QDoubleSpinBox::down-arrow, QSpinBox::down-arrow {{ image: url("{down}"); width: 9px; height: 9px; }}
QDoubleSpinBox::up-arrow:disabled, QSpinBox::up-arrow:disabled,
QDoubleSpinBox::down-arrow:disabled, QSpinBox::down-arrow:disabled,
QDoubleSpinBox[readOnly="true"]::up-arrow, QDoubleSpinBox[readOnly="true"]::down-arrow {{ image: none; }}

/* ---- check boxes: a filled box with a visible tick ------------------------- */
QCheckBox, QRadioButton {{ spacing: 8px; color: {TEXT}; padding: 2px 0; }}
QCheckBox::indicator, QRadioButton::indicator,
QListView::indicator, QListWidget::indicator {{
    width: 16px; height: 16px; border: 1px solid {BORDER_STRONG};
    background: {SURFACE}; border-radius: 4px;
}}
QRadioButton::indicator {{ border-radius: 9px; }}
QCheckBox::indicator:hover, QRadioButton::indicator:hover,
QListView::indicator:hover, QListWidget::indicator:hover {{ border-color: {ACCENT}; }}
QCheckBox::indicator:checked, QListView::indicator:checked, QListWidget::indicator:checked {{
    background: {ACCENT}; border-color: {ACCENT}; image: url("{check}");
}}
QRadioButton::indicator:checked {{ background: {ACCENT}; border-color: {ACCENT}; image: url("{dot}"); }}
QCheckBox::indicator:disabled {{ background: {SURFACE_ALT}; border-color: {BORDER}; }}
QCheckBox::indicator:checked:disabled {{ background: #B9CBCA; border-color: #B9CBCA; }}
QCheckBox:disabled {{ color: {TEXT_FAINT}; }}

/* ---- lists and tables ------------------------------------------------------ */
QListWidget {{
    background: {SURFACE}; border: 1px solid {BORDER}; border-radius: 8px;
    outline: none; padding: 4px;
}}
QListWidget::item {{ padding: 7px 8px; border-radius: 6px; color: {TEXT}; border: 1px solid transparent; }}
QListWidget::item:hover {{ background: {ACCENT_SOFT}; color: {TEXT}; }}
QListWidget::item:selected, QListWidget::item:selected:!active {{
    background: {ACCENT_SELECTED}; color: {TEXT}; border: 1px solid #A9CBC7;
}}
QListWidget::item:selected:hover {{ background: {ACCENT_SELECTED}; color: {TEXT}; border: 1px solid {ACCENT}; }}
QTreeWidget {{
    background: {SURFACE}; border: 1px solid {BORDER}; border-radius: 8px;
    outline: none; padding: 4px;
}}
QTreeWidget::item {{ padding: 5px 6px; border-radius: 6px; color: {TEXT}; border: 1px solid transparent; }}
QTreeWidget::item:hover {{ background: {ACCENT_SOFT}; color: {TEXT}; }}
QTreeWidget::item:selected, QTreeWidget::item:selected:!active {{
    background: {ACCENT_SELECTED}; color: {TEXT}; border: 1px solid transparent;
}}
QTreeWidget::item:disabled {{ color: {TEXT_FAINT}; }}
QTreeView::branch {{ background: {SURFACE}; image: none; border-image: none; }}
QTreeView::branch:selected, QTreeView::branch:hover {{ background: {SURFACE}; image: none; border-image: none; }}
QTreeView::branch:has-children:closed {{ image: url({right}); }}
QTreeView::branch:has-children:open {{ image: url({down}); }}
QTableWidget {{
    background: {SURFACE}; border: 1px solid {BORDER}; border-radius: 8px; gridline-color: {BORDER};
}}
QHeaderView::section {{
    background: {SURFACE_ALT}; color: {TEXT_MUTED}; border: none;
    border-bottom: 1px solid {BORDER}; padding: 4px 6px; font-size: 12px;
}}

/* ---- rich text panes --------------------------------------------------------- */
QTextBrowser, QTextEdit, QPlainTextEdit {{
    background: {SURFACE}; border: 1px solid {BORDER}; border-radius: 10px;
    selection-background-color: {ACCENT_SELECTED}; selection-color: {TEXT};
}}

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
    font-size: 12px; color: {TEXT_MUTED}; padding: 4px 10px;
    border-radius: 13px; background: transparent;
}}
QLabel#StepDot:hover {{ background: {ACCENT_SOFT}; color: {TEXT}; }}
QLabel#StepDotActive {{
    font-size: 12px; font-weight: 600; color: {ACCENT_TEXT}; padding: 4px 10px;
    border-radius: 13px; background: {ACCENT};
}}
QLabel#StepDotDone {{
    font-size: 12px; color: {ACCENT}; padding: 4px 10px;
    border-radius: 13px; background: {ACCENT_SOFT};
}}
QLabel#StepDotDone:hover {{ background: {ACCENT_SELECTED}; color: {TEXT}; }}

/* ---- summary sidebar --------------------------------------------------------- */
QWidget#Summary {{ background: {SURFACE}; }}
QLabel#SummaryKey {{ font-size: 11.5px; color: {TEXT_FAINT}; }}
QLabel#SummaryVal {{ font-size: 12.5px; color: {TEXT}; }}
QLabel#SummaryHead {{
    font-size: 11px; font-weight: 600; color: {TEXT_MUTED};
    letter-spacing: 1px; padding-top: 4px;
}}

QProgressBar {{
    background: {BORDER}; border: none; border-radius: 4px; max-height: 7px; text-align: center;
    color: transparent;
}}
QProgressBar::chunk {{ background: {ACCENT}; border-radius: 4px; }}

QToolTip {{
    background: {SURFACE}; color: {TEXT}; border: 1px solid {BORDER_STRONG};
    padding: 5px 8px; border-radius: 6px;
}}
QMenuBar {{ background: {BG}; border-bottom: 1px solid {BORDER}; }}
QMenuBar::item {{ padding: 5px 11px; background: transparent; border-radius: 6px; color: {TEXT}; }}
QMenuBar::item:selected {{ background: {ACCENT_SOFT}; color: {TEXT}; }}
QMenu {{ background: {SURFACE}; border: 1px solid {BORDER_STRONG}; border-radius: 8px; padding: 5px; }}
QMenu::item {{ padding: 6px 22px; border-radius: 6px; color: {TEXT}; }}
QMenu::item:selected {{ background: {ACCENT_SOFT}; color: {TEXT}; }}
QSplitter::handle {{ background: {BORDER}; width: 1px; }}
QMessageBox QLabel {{ color: {TEXT}; }}
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
