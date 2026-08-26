@echo off
setlocal EnableExtensions
cd /d "%~dp0"

set "VENV_DIR=%~dp0.venv"
set "VENV_PYTHON=%VENV_DIR%\Scripts\python.exe"

echo [Voronoi Notch Simulation Launcher]
echo [STEP 1/2] Checking the virtual environment and required packages...

if not exist "%VENV_PYTHON%" (
    set "ERROR_MESSAGE=Project virtual environment not found. Run install.bat first."
    goto :error
)

"%VENV_PYTHON%" -c "import sys, struct; raise SystemExit(0 if (3, 10) <= sys.version_info[:2] < (3, 14) and struct.calcsize('P') * 8 == 64 else 1)"
if errorlevel 1 (
    set "ERROR_MESSAGE=The Python version in .venv is unsupported. Recreate it with install.bat."
    goto :error
)

"%VENV_PYTHON%" -m pip check
if errorlevel 1 (
    set "ERROR_MESSAGE=Package dependency validation failed. Run install.bat again."
    goto :error
)

"%VENV_PYTHON%" -c "import openpyxl, numpy, scipy, matplotlib, seaborn, pandas, tqdm, PyQt6.QtWidgets, PyQt6.Qsci, numba, h5py"
if errorlevel 1 (
    set "ERROR_MESSAGE=A required package is missing or cannot be imported. Run install.bat again."
    goto :error
)

echo [STEP 2/2] Starting the Voronoi Notch GUI...
"%VENV_PYTHON%" -m gui.main_window
set "APP_EXIT_CODE=%ERRORLEVEL%"

if not "%APP_EXIT_CODE%"=="0" (
    set "ERROR_MESSAGE=The application exited with code %APP_EXIT_CODE%."
    goto :app_error
)

exit /b 0

:app_error
echo.
echo [ERROR] %ERROR_MESSAGE%
if not defined VORONOI_NO_PAUSE pause
exit /b %APP_EXIT_CODE%

:error
echo.
echo [ERROR] %ERROR_MESSAGE%
if not defined VORONOI_NO_PAUSE pause
exit /b 1
