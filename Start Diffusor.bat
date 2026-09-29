@echo off
rem Double-click to start Diffusor.
rem The first run sets up a private Python environment in the .venv folder here and
rem installs what Diffusor needs. That takes a few minutes and an internet connection.
rem Later runs start Diffusor straight away. Delete .venv to set it up again.
setlocal
title Diffusor
cd /d "%~dp0"
set "VENV=%~dp0.venv"
set "READY=%VENV%\diffusor-ready.txt"

if exist "%READY%" goto launch

echo Setting up Diffusor for the first time.
echo This needs an internet connection and takes a few minutes. Later starts are quick.
echo.

rem --- find Python 3.10 or newer: the py launcher first, then python on the PATH
set "PY="
py -3 -c "import sys; sys.exit(sys.version_info < (3, 10))" >nul 2>&1
if not errorlevel 1 set "PY=py -3"
if defined PY goto havepy
python -c "import sys; sys.exit(sys.version_info < (3, 10))" >nul 2>&1
if not errorlevel 1 set "PY=python"
if defined PY goto havepy
goto nopython

:havepy
if exist "%VENV%\Scripts\python.exe" goto install
echo Creating a private Python environment in the .venv folder ...
%PY% -m venv "%VENV%"
if errorlevel 1 goto failed

:install
echo Installing Diffusor and the libraries it needs ...
"%VENV%\Scripts\python.exe" -m pip install --disable-pip-version-check -e .
if errorlevel 1 goto failed
"%VENV%\Scripts\python.exe" -c "import diffusor.gui.main_window"
if errorlevel 1 goto failed
echo ready> "%READY%"
echo.
echo Diffusor is set up. Starting it now.

:launch
if not exist "%VENV%\Scripts\pythonw.exe" goto broken
rem a shortcut with Diffusor's icon beside this file, remade each time so it
rem still works after the folder is moved
powershell -NoProfile -Command "$s = (New-Object -ComObject WScript.Shell).CreateShortcut('%~dp0Diffusor.lnk'); $s.TargetPath = '%VENV%\Scripts\pythonw.exe'; $s.Arguments = '-m diffusor'; $s.WorkingDirectory = '%~dp0'; $s.IconLocation = '%~dp0diffusor\gui\icons\diffusor.ico,0'; $s.Description = 'Diffusor'; $s.Save()" >nul 2>&1
start "" "%VENV%\Scripts\pythonw.exe" -m diffusor
exit /b 0

:broken
del "%READY%" >nul 2>&1
echo The Diffusor environment in .venv is incomplete. Run this file again to repair it.
pause
exit /b 1

:nopython
echo Diffusor needs Python 3.10 or newer, and none was found.
echo Install it from https://www.python.org/downloads/ and tick "Add python.exe to PATH"
echo during the installation, then run this file again.
pause
exit /b 1

:failed
echo.
echo The setup did not finish. The messages above say why.
echo Check the internet connection and run this file again.
echo To start from scratch, delete the .venv folder first.
pause
exit /b 1
