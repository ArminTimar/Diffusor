from .gui.app import main

# The guard matters: parallel Monte Carlo workers re-import this module under
# another name, and without it each would open its own window.
if __name__ == "__main__":
    raise SystemExit(main())
