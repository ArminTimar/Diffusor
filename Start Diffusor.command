#!/bin/bash
# Double-click in Finder to start Diffusor on a Mac. On Linux, run it from a terminal.
# The first run sets up a private Python environment in the .venv folder here and
# installs what Diffusor needs. That takes a few minutes and an internet connection.
# Later runs start Diffusor straight away. Delete .venv to set it up again.

cd "$(dirname "$0")" || exit 1
VENV="$PWD/.venv"
READY="$VENV/diffusor-ready.txt"

venv_python() {
    if [ -x "$VENV/bin/python" ]; then
        echo "$VENV/bin/python"
    elif [ -x "$VENV/Scripts/python.exe" ]; then      # Git Bash on Windows
        echo "$VENV/Scripts/python.exe"
    fi
}

fail() {
    echo
    echo "$1"
    echo
    echo "Press Return to close this window."
    read -r _
    exit 1
}

if [ ! -f "$READY" ] || [ -z "$(venv_python)" ]; then
    echo "Setting up Diffusor for the first time."
    echo "This needs an internet connection and takes a few minutes. Later starts are quick."
    echo

    # Python 3.10 or newer. Versioned names first: the python3 that comes with
    # macOS is often 3.9, and python.org installs add python3.12 and the like.
    PY=""
    for cand in python3.14 python3.13 python3.12 python3.11 python3.10 python3 python; do
        if command -v "$cand" >/dev/null 2>&1 &&
           "$cand" -c "import sys; sys.exit(sys.version_info < (3, 10))" >/dev/null 2>&1; then
            PY="$cand"
            break
        fi
    done
    [ -n "$PY" ] || fail "Diffusor needs Python 3.10 or newer, and none was found.
Install it from https://www.python.org/downloads/ (the macOS installer),
then open this file again. The python3 that comes with macOS is too old."

    if [ -z "$(venv_python)" ]; then
        echo "Creating a private Python environment in the .venv folder ..."
        "$PY" -m venv "$VENV" || fail "Could not create the Python environment in .venv."
    fi
    VPY="$(venv_python)"
    echo "Installing Diffusor and the libraries it needs ..."
    "$VPY" -m pip install --disable-pip-version-check -e . ||
        fail "The setup did not finish. The messages above say why.
Check the internet connection and open this file again.
To start from scratch, delete the .venv folder first."
    "$VPY" -c "import diffusor.gui.main_window" ||
        fail "Diffusor was installed but does not start. The messages above say why."
    echo ready > "$READY"
    echo
    echo "Diffusor is set up."
fi

# start Diffusor on its own, so closing this Terminal window does not close it
nohup "$(venv_python)" -m diffusor >/dev/null 2>&1 &
disown 2>/dev/null
echo "Diffusor is starting. You can close this window."
exit 0
