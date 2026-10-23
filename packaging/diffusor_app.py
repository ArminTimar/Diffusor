"""Start script of the packaged application (the Diffusor program the installers ship).

``python -m diffusor`` is the start for a source checkout; this file plays that role
inside the bundle PyInstaller builds. ``--self-check [REPORT]`` builds the window off
screen, runs a short Monte Carlo fit on two worker processes, writes what happened to
REPORT and exits with 0 when both worked. The release build runs it on every bundle
before making the installer.
"""
import multiprocessing
import sys


def _self_check() -> int:
    import os
    os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
    import numpy as np
    from PySide6.QtWidgets import QApplication

    from diffusor.coefficients import Conditions, get
    from diffusor.datasets import EXAMPLES_DIR
    from diffusor.fitting import DiffusionModel, UncertaintyBudget, run_montecarlo
    from diffusor.gui.main_window import MainWindow
    from diffusor.solvers import InitialCondition, dirichlet

    app = QApplication(sys.argv[:1])
    win = MainWindow()
    win.close()
    del app
    if not (EXAMPLES_DIR / "validation" / "manifest.json").is_file():
        print(f"the examples are missing from {EXAMPLES_DIR}")
        return 1

    # the same set-up as the parallel test in tests/test_fitting.py
    model = DiffusionModel(coefficient=get("mt_Ti_vanorman_crispin2010"),
                           conditions=Conditions(T_K=950 + 273.15, P_Pa=1e8, log_fo2_bar=-11.0,
                                                 X={"xTi": 0.1}),
                           initial=InitialCondition("step", {"x0": 0.0, "C_left": 0.30,
                                                             "C_right": 0.18}),
                           bc_left=dirichlet(0.30), bc_right=dirichlet(0.18),
                           composition_dependent=False)
    x = np.linspace(-40, 40, 41)
    C = model.profile(8 * 86400.0, x) + np.random.default_rng(1).normal(0, 0.003, x.size)
    budget = UncertaintyBudget(sigma_T_K=20.0, buffer="NNO", delta_buffer=1.0,
                               sigma_delta_buffer=0.3)
    budget.sample_coefficient = False
    res = run_montecarlo(model, x, C, np.full(x.size, 0.003), budget=budget, n_draws=12,
                         seed=2, workers=2, parallel_threshold_s=0.0)
    stopped = [w for w in res.warnings if "worker processes" in w or "one core" in w]
    if res.workers != 2 or len(res.times) != 12 or stopped:
        print(f"the Monte Carlo run used {res.workers} worker(s) and kept "
              f"{len(res.times)} of 12 draws; {'; '.join(stopped)}")
        return 1
    print("self-check passed")
    return 0


if __name__ == "__main__":
    # Monte Carlo workers start this program again; freeze_support runs them as
    # workers instead of opening a second window.
    multiprocessing.freeze_support()
    if "--self-check" in sys.argv[1:]:
        # A windowed program has no console, so the outcome goes to the file named
        # after the flag, and an error must not open the bootloader's message box.
        import contextlib
        import io
        import traceback
        out = io.StringIO()
        try:
            with contextlib.redirect_stdout(out):
                code = _self_check()
        except BaseException:
            out.write(traceback.format_exc())
            code = 1
        args = sys.argv[1:]
        at = args.index("--self-check")
        if at + 1 < len(args):
            with open(args[at + 1], "w", encoding="utf8") as f:
                f.write(out.getvalue())
        else:
            sys.__stdout__ and sys.__stdout__.write(out.getvalue())
        raise SystemExit(code)
    from diffusor.gui.app import main
    raise SystemExit(main())
