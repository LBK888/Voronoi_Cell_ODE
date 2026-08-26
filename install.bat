@echo off
setlocal EnableExtensions
cd /d "%~dp0"

set "VENV_DIR=%~dp0.venv"
set "VENV_PYTHON=%VENV_DIR%\Scripts\python.exe"
set "PYTHON_EXE="

echo [Voronoi Notch Simulation Installer]
echo [STEP 1/4] Looking for a supported Python installation...

call :find_python
if not defined PYTHON_EXE (
    set "ERROR_MESSAGE=No supported Python was found. Install 64-bit Python 3.10 through 3.13. Python 3.12 is recommended."
    goto :error
)

"%PYTHON_EXE%" -c "import sys, struct; raise SystemExit(0 if (3, 10) <= sys.version_info[:2] < (3, 14) and struct.calcsize('P') * 8 == 64 else 1)"
if errorlevel 1 (
    set "ERROR_MESSAGE=Python must be a 64-bit version from 3.10 through 3.13."
    goto :error
)

echo [INFO] Using this Python installation:
"%PYTHON_EXE%" -c "import sys; print(sys.executable); print(sys.version)"

echo [STEP 2/4] Creating or checking the project virtual environment...
if exist "%VENV_DIR%" if not exist "%VENV_PYTHON%" (
    set "ERROR_MESSAGE=.venv exists but is incomplete. Rename or remove .venv, then run install.bat again."
    goto :error
)

if not exist "%VENV_PYTHON%" (
    "%PYTHON_EXE%" -m venv "%VENV_DIR%"
    if errorlevel 1 (
        set "ERROR_MESSAGE=Could not create .venv. Verify that Python includes the venv module and that this folder is writable."
        goto :error
    )
) else (
    "%VENV_PYTHON%" -c "import sys, struct; raise SystemExit(0 if (3, 10) <= sys.version_info[:2] < (3, 14) and struct.calcsize('P') * 8 == 64 else 1)"
    if errorlevel 1 (
        set "ERROR_MESSAGE=The existing .venv uses an unsupported Python. Rename or remove .venv, then run install.bat again."
        goto :error
    )
)

echo [STEP 3/4] Updating pip and installing requirements.txt...
if not exist "%~dp0requirements.txt" (
    set "ERROR_MESSAGE=requirements.txt was not found."
    goto :error
)

"%VENV_PYTHON%" -m pip install --upgrade pip
if errorlevel 1 (
    set "ERROR_MESSAGE=Could not update pip. Check the network connection and Python installation."
    goto :error
)

"%VENV_PYTHON%" -m pip install --upgrade -r "%~dp0requirements.txt"
if errorlevel 1 (
    set "ERROR_MESSAGE=Could not install requirements.txt. Review the pip error above."
    goto :error
)

echo [STEP 4/4] Validating dependencies and native imports...
"%VENV_PYTHON%" -m pip check
if errorlevel 1 (
    set "ERROR_MESSAGE=Package dependency validation failed. Review the pip error above."
    goto :error
)

"%VENV_PYTHON%" -c "import openpyxl, numpy, scipy, matplotlib, seaborn, pandas, tqdm, PyQt6.QtWidgets, PyQt6.Qsci, numba, h5py"
if errorlevel 1 (
    set "ERROR_MESSAGE=Packages were installed, but at least one required package cannot be imported. Review the error above."
    goto :error
)

echo.
echo [DONE] Installation and validation succeeded. You can now run start.bat.
if not defined VORONOI_NO_PAUSE pause
exit /b 0

:find_python
where py >nul 2>nul
if errorlevel 1 goto :find_python_on_path

call :try_py_version 3.12
if defined PYTHON_EXE goto :eof
call :try_py_version 3.13
if defined PYTHON_EXE goto :eof
call :try_py_version 3.11
if defined PYTHON_EXE goto :eof
call :try_py_version 3.10
if defined PYTHON_EXE goto :eof

:find_python_on_path
where python >nul 2>nul
if errorlevel 1 goto :eof
for /f "delims=" %%P in ('python -c "import sys; print(sys.executable)" 2^>nul') do set "PYTHON_EXE=%%P"
if not defined PYTHON_EXE goto :eof
"%PYTHON_EXE%" -c "import sys, struct; raise SystemExit(0 if (3, 10) <= sys.version_info[:2] < (3, 14) and struct.calcsize('P') * 8 == 64 else 1)"
if errorlevel 1 set "PYTHON_EXE="
goto :eof

:try_py_version
for /f "delims=" %%P in ('py -%~1 -c "import sys; print(sys.executable)" 2^>nul') do set "PYTHON_EXE=%%P"
if not defined PYTHON_EXE goto :eof
"%PYTHON_EXE%" -c "import sys, struct; raise SystemExit(0 if sys.version_info[:2] == tuple(map(int, '%~1'.split('.'))) and struct.calcsize('P') * 8 == 64 else 1)"
if errorlevel 1 set "PYTHON_EXE="
goto :eof

:error
echo.
echo [ERROR] %ERROR_MESSAGE%
if not defined VORONOI_NO_PAUSE pause
exit /b 1
