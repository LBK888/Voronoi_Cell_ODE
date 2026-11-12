# Agent CLI

在主視窗選擇 **Agent → Enable Agent CLI**。CLI 使用相同使用者的設定，會在獨立程序執行模擬。將本說明及視窗上方的程式路徑交給 agent。模型可包含 Python 程式碼，請選用可信任的 JSON。

## 呼叫方式

EXE 是 Windows 視窗程式。從 PowerShell 呼叫時，使用 `Start-Process -Wait` 並指定 `--report`，以取得結束碼與 JSON。範例中的路徑請換成實際位置。使用原始碼時，將 `VoronoiCellODE.exe` 換成 `python windows_launcher.py`，可直接讀取標準輸出的 JSON。

```powershell
$exe = 'C:\Voronoi\VoronoiCellODE.exe'
$job = Start-Process -FilePath $exe -ArgumentList '--cli --report "C:\Voronoi\templates-report.json" list-templates' -Wait -PassThru -WindowStyle Hidden
Get-Content -LiteralPath 'C:\Voronoi\templates-report.json' -Raw
$job.ExitCode
```

Agent 以 subprocess 呼叫時，應傳入參數陣列，例如 `[exe, "--cli", "--report", report, "simulate", model, "--output", output]`。先等待程序結束，再讀取 report。使用新報告檔路徑，方便把回應對應到本次呼叫。成功回傳 `status: ok` 與結束碼 0，失敗回傳 `status: error`、`error` 與非零結束碼。把 `--report` 放在子指令前，報告檔放在輸出資料夾之外。

## 指令

| 子指令 | 用途 |
| --- | --- |
| `list-templates` | 列出模板名稱與絕對路徑。 |
| `inspect MODEL.json` | 讀取模型設定及檢查 Python 語法。 |
| `simulate MODEL.json --output NEW_DIR [--seed N]` | 執行單次模擬，輸出各時間點的 `cells.csv`、模型 JSON 與初始細胞位置。 |
| `experiment MODEL.json --output NEW_DIR [--settings SWEEP.json] [--seed N]` | 執行二維掃描，輸出熱圖 PDF、逐次結果 CSV、掃描設定與初始細胞位置。 |

每次模擬指定尚未存在的輸出資料夾。單次模擬預設 seed 0。掃描預設使用掃描設定的 seed，`--seed` 可覆寫兩者。掃描的報告有 `failed_runs`；每次失敗原因記錄在 CSV 的 `error` 欄。自行匯入的網格可用 **Save Current Model as Template** 保存 `cell_positions`，供 CLI 重用。

## 掃描設定範例

以下 JSON 可存成 `SWEEP.json` 搭配 `--settings`，或作為模型內的 `experiments` 物件。將 `a` 改成模型的參數名，`initial` 的 target 是從 0 開始的狀態變數索引。GUI 的 Experiments 設定也會隨模板儲存。

```json
{
  "x": {"kind": "parameter", "target": "a", "minimum": 0.1, "maximum": 1, "count": 5, "scale": "linear"},
  "y": {"kind": "initial", "target": 0, "minimum": 0, "maximum": 1, "count": 5, "scale": "linear"},
  "observable": 0,
  "metric": "final_mean",
  "repeats": 1,
  "seed": 0,
  "tail_percent": 25
}
```

可用指標為 `final_mean`、`final_std`、`final_cv`、`tail_mean`、`tail_amplitude`。開啟 JSON 後可用 **Generate Parameter Controls** 檢查與調整參數。GUI 與 CLI 使用同一模擬引擎。
