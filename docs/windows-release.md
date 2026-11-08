# Voronoi Cell ODE — Windows 64-bit

## 直接執行 / Run without Python

下載 Release 的 `VoronoiCellODE.exe`，或解壓縮 Windows ZIP 後，雙擊 `VoronoiCellODE.exe`。不需要安裝 Python、venv 或 FFmpeg。第一次啟動需解壓縮內附的科學運算函式庫，可能需要稍候。

Download `VoronoiCellODE.exe` from the release, or extract the Windows ZIP and double-click the executable. Python and FFmpeg are bundled. The first startup can take a moment while libraries are extracted. This build targets 64-bit Windows 10/11 and is not digitally signed.

設定與日誌會保存在 `%LOCALAPPDATA%\VoronoiCellODE`，不會寫入 EXE 的暫存解壓縮位置，也不需要以系統管理員執行。更新 EXE 不會覆蓋既有設定。

Settings and logs persist in `%LOCALAPPDATA%\VoronoiCellODE`. Replacing the EXE keeps those settings. Templates are bundled; your selected export folder controls where simulation output is written.

## This release

- MP4 exports include labeled, per-frame colorbars.
- Save Frame defaults to editable PDF; SVG/PNG and right-click clipboard copying are available.
- Experiments follows Run, with two-axis parameter/initial-concentration sweeps, statistics, replicates and CSV/JSON export.
- Mouse-wheel scrolling no longer changes parameters, dropdown selections or slider values.
- The repository includes a paper-to-ODE agent skill and a validated teaching model.

The application executes user-authored Python model code. Use trusted models on your own computer, as described in the main project README.

## Rebuild from source

Use the runtime checkout on Windows with 64-bit Python 3.12:

```powershell
.\.venv\Scripts\python.exe -m pip install -r requirements.txt -r packaging\requirements-windows.txt -c packaging\constraints-windows.txt
.\.venv\Scripts\python.exe -m unittest discover -s tests -v
.\.venv\Scripts\python.exe packaging\build_windows.py --version 2026.11.08
```

Build outputs stay in `build/` and `dist/`; these folders are excluded from Git. The build runs a frozen-EXE smoke test covering the Qt window, NumPy/SciPy/Numba model execution, clipboard, colorbar MP4, PDF, Experiments, Excel, HDF5 and both editors. `SHA256SUMS.txt` identifies the release binaries. `third_party/` in the ZIP contains dependency licenses, build package versions and FFmpeg source/build links.

PyInstaller references: [spec files](https://pyinstaller.org/en/stable/spec-files.html), [runtime resource paths](https://pyinstaller.org/en/stable/runtime-information.html). FFmpeg binary distribution: [imageio-ffmpeg](https://github.com/imageio/imageio-ffmpeg), [Gyan builds and sources](https://www.gyan.dev/ffmpeg/builds/).
