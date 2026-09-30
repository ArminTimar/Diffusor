"""Background workers so the interface stays responsive.

Fitting one profile takes a second; a 1000-draw Monte Carlo takes minutes.
Both run on a QThread and report progress, and the Monte Carlo can be
cancelled.
"""
from __future__ import annotations

import time
import traceback
from typing import Dict, List, Optional, Sequence

import numpy as np
from PySide6.QtCore import QObject, QThread, Signal

from .. import updates
from ..fitting import (DiffusionModel, UncertaintyBudget, contributions, fit_time,
                       run_montecarlo)


class FitWorker(QObject):
    finished = Signal(object)
    failed = Signal(str)
    progress = Signal(str)

    def __init__(self, model: DiffusionModel, x, C, sigma, free_parameters,
                 t_min: float, t_max: float):
        super().__init__()
        self.args = (model, x, C, sigma, free_parameters, t_min, t_max)

    def run(self):
        try:
            model, x, C, sigma, free, t_min, t_max = self.args
            self.progress.emit("fitting...")
            res = fit_time(model, x, C, sigma, free, t_min, t_max)
            self.finished.emit(res)
        except Exception:
            self.failed.emit(traceback.format_exc())


class CompareWorker(QObject):
    """Fit the same profile with several diffusion coefficients."""
    finished = Signal(object)
    failed = Signal(str)
    progress = Signal(int, int, str)

    def __init__(self, models: Dict[str, DiffusionModel], x, C, sigma,
                 free_parameters, t_min: float, t_max: float):
        super().__init__()
        self.models = models
        self.rest = (x, C, sigma, free_parameters, t_min, t_max)
        self._abort = False

    def abort(self):
        self._abort = True

    def run(self):
        try:
            x, C, sigma, free, t_min, t_max = self.rest
            out = {}
            n = len(self.models)
            for i, (key, model) in enumerate(self.models.items()):
                if self._abort:
                    break
                self.progress.emit(i + 1, n, key)
                try:
                    out[key] = fit_time(model, x, C, sigma, free, t_min, t_max)
                except Exception as exc:
                    out[key] = exc
            self.finished.emit(out)
        except Exception:
            self.failed.emit(traceback.format_exc())


class MonteCarloWorker(QObject):
    """Monte Carlo with a live feed.

    Every successful draw is buffered and sent to the interface in batches, at
    most once a second, so the plot can grow while the run goes on
    without flooding the event loop.
    """
    finished = Signal(object)
    failed = Signal(str)
    progress = Signal(int, int)
    stage = Signal(str)
    draws = Signal(object)

    BATCH_SECONDS = 1.0

    def __init__(self, model: DiffusionModel, x, C, sigma, budget: UncertaintyBudget,
                 n_draws: int, seed: int, free_parameters, t_min: float, t_max: float,
                 do_contributions: bool = False, contribution_draws: int = 100,
                 workers: int = 1):
        super().__init__()
        self.args = (model, x, C, sigma, budget, n_draws, seed, free_parameters, t_min, t_max)
        self.do_contributions = do_contributions
        self.contribution_draws = contribution_draws
        self.workers = workers
        self._abort = False
        self._buffer: List[dict] = []
        self._last_emit = 0.0

    def abort(self):
        self._abort = True

    def _cb(self, i, n):
        self.progress.emit(i, n)
        return self._abort

    def _on_draw(self, info: dict):
        self._buffer.append(info)
        now = time.monotonic()
        if now - self._last_emit >= self.BATCH_SECONDS:
            self._flush()
            self._last_emit = now

    def _flush(self):
        if self._buffer:
            batch, self._buffer = self._buffer, []
            self.draws.emit(batch)

    def run(self):
        try:
            model, x, C, sigma, budget, n_draws, seed, free, t_min, t_max = self.args
            self.stage.emit("Monte Carlo")
            res = run_montecarlo(model, x, C, sigma, budget=budget, n_draws=n_draws,
                                 seed=seed, free_parameters=free, t_min=t_min, t_max=t_max,
                                 progress=self._cb, on_draw=self._on_draw,
                                 workers=self.workers)
            self._flush()
            if self.do_contributions and not self._abort:
                self.stage.emit("variance contributions")
                res.contributions = contributions(
                    model, x, C, sigma, budget=budget, n_draws=self.contribution_draws,
                    seed=seed + 1, free_parameters=free, progress=self._cb,
                    workers=self.workers)
            self.finished.emit(res)
        except Exception:
            self.failed.emit(traceback.format_exc())


class UpdateWorker(QObject):
    """Ask GitHub for the newest release without freezing the window."""
    finished = Signal(object)       # updates.ReleaseInfo
    failed = Signal(str)

    def __init__(self, current: str):
        super().__init__()
        self.current = current

    def run(self):
        try:
            self.finished.emit(updates.fetch_latest(self.current))
        except updates.UpdateError as exc:
            self.failed.emit(str(exc))
        except Exception:
            self.failed.emit("The update check failed unexpectedly.")


def start(worker: QObject) -> QThread:
    """Move a worker onto a fresh thread and start it."""
    thread = QThread()
    worker.moveToThread(thread)
    thread.started.connect(worker.run)
    worker.finished.connect(thread.quit)
    if hasattr(worker, "failed"):
        worker.failed.connect(thread.quit)
    thread.start()
    return thread
