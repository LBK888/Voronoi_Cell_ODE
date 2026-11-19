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

X、Y 均可選 **Cell movement strength**。CLI 請設 `"kind": "movement"`，target 可填 `random_strength`、`away_strength`、`convergent_strength`、`repulsion_strength` 或 `vertex_strength`。掃描值直接套用到該次模擬，即使模板原先關閉該項移動也有效；0 表示關閉。

熱圖可選原有的 `final_mean`、`final_std`、`final_cv`、`tail_mean`、`tail_amplitude`，以及以下三項動態指標。每次重複實驗的 CSV 都包含 `final_mean` 與三項動態指標。`response`、`pixel_mean` 與 `pixel_std` 對應熱圖選定的指標。

| CSV 欄位 | 定義與單位 |
|---|---|
| `rise_slope` | 第一段明顯上升內，以連續五幀做局部直線回歸，取最大正斜率。單位為物種值／模型時間。 |
| `onset_time` | 將該段最陡切線回推至段首基線（前一明顯谷底或初值），限制於段首至切點之間，回報從模擬開始算起的時間。 |
| `frequency` | 相鄰明顯波峰之時間差的中位數取倒數，單位為週期數／模型時間。 |

三者都使用逐幀存活細胞的空間平均與該次模擬的 `dT`。異步細胞在空間平均時可能互相抵消。分裂或凋亡改變細胞組成時，也可能造成平均訊號的變化。無法辨識的指標在 CSV 留空，選為熱圖指標時顯示灰色。

斜率與起始時間採用最大斜率切線的幾何定義，將其應用到第一段明顯上升。峰谷偵測使用五幀移動平均，顯著度門檻取全段振幅的 10%、估計雜訊的五倍與浮點誤差下限的最大值。雜訊尺度以二階差分絕對值中位數除以 1.652 估算。[Kahm 等，grofit (2010)](https://www.jstatsoft.org/article/view/v033i07)

頻率採用平滑後的峰間距法。此實作用於掃描的判讀門檻為至少三個峰、週期至少五幀、所有峰距距離中位數在 25% 以內，且波形平移一個估計週期後的 Pearson 相關至少為 0.70。這些門檻可能排除具有長暫態、變頻或衰減的訊號。[Chemical modulation of circadian rhythms and assessment of cellular behavior via indirubin and derivatives](https://pmc.ncbi.nlm.nih.gov/articles/PMC8516081/)

### 隨附模板驗證

下表使用 seed 0、模板原有網格、時長與時間步長。列出第一個狀態變數的空間平均。這次七條曲線的頻率都未通過判讀門檻。

| 模板 | 觀測量 | 最終平均 | 上升斜率 | 起始時間 |
|---|---|---:|---:|---:|
| Noisy Oscillator | mRNA | 0.8052 | 0.4855 | 0 |
| Cell Cycle | Cyclin_CDK1 | 0.2830 | 0.1420 | 0.00035 |
| cardiac pulses | Voltage | −0.6080 | 0.2688 | 23.70 |
| Turing Pattern Zebrafish | Activator | 2.8700 | 11.2441 | 4.3014 |
| soluble Delta notochord | DeltaC | 10.9016 | 19.6994 | 0 |
| PSM w DAPT | Phase_Theta | 36.3768 | 2.0034 | 0.00464 |
| soluble Delta ionocyte | DeltaC | 11.2062 | 19.6449 | 0 |

心臟模板的刺激週期為 25。對左側 24 個直接受刺激細胞的平均電位獨立驗算，得到頻率 0.03987（預期 0.04000）。全組織平均出現大小交替的峰，約 12.5 的峰距只得到 0.53 的波形相關，因此頻率留空。其起始時間 23.70 表示通過顯著度門檻的第一段組織平均上升。左側細胞的起始估計為 0.061，約一個取樣間隔。另以已知斜率 2、延遲 4 的線性上升與週期 4 的正弦波測試，兩種時間步長（0.05、0.1）均重現預設值。
