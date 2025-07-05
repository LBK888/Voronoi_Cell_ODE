@echo off
chcp 65001 >nul
echo [Voronoi Notch Simulation 啟動器]
where python >nul 2>nul
if errorlevel 1 (
    echo [錯誤] 未偵測到 Python，請先安裝 Python 3.8 以上版本並加入 PATH。
    pause
    exit /b 1
)
echo [STEP1] 安裝/檢查必要套件...
python -m pip install -r requirements.txt
if errorlevel 1 (
    echo [錯誤] 依賴安裝失敗，請檢查網路或權限。
    pause
    exit /b 1
)
echo [STEP2] 啟動 Voronoi Notch GUI...
python -m gui.main_window
if errorlevel 1 (
    echo [錯誤] 程式執行失敗，請檢查錯誤訊息或 requirements.txt 是否完整。
    pause
    exit /b 1
)
pause 