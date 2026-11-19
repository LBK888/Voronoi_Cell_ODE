# Experiments：二維數值掃描與熱圖

## 使用流程

1. 在 Grid 產生格點，於 Dynamics 設定細胞動態，在 Model 輸入 ODE／Params。Run 的 `T` 決定每次試驗的模擬時間。
2. 選取 Run 後方的 **Experiments**。分頁會讀取目前模型；變更模型後也可按 **Refresh targets from current model** 重新整理。
3. 分別設定 X 與 Y。選擇 **Model parameter**、**Initial concentration (all cells)** 或 **Cell movement strength**，再設定目標、上下限、Samples 和 Linear／Logarithmic。
4. 選擇 Response 的 Species 與 Statistic，設定 Replicates / pixel、Random seed，以及需要時的 Tail window。
5. 按 **Run Experiment**。右側顯示進度與完成後的熱圖。可按 **Cancel Experiment**；取消後保留上一次完成的熱圖。
6. 按 **Save Heatmap...** 儲存 PDF（也可選 SVG／PNG）；按 **Export Data...** 儲存 CSV 與伴隨的 `_settings.json`。熱圖按右鍵可複製到剪貼簿。

## 兩軸的意義

Model parameter 可掃描 Params 中存在的有限數值 scalar；不包括布林、陣列、`n_var`。整數參數的掃描值也會是浮點數；如果參數代表離散個數，應在模型中明確處理其意義。`dT` 雖可掃描，但這是數值精度分析，不應直接當成生物參數。

Initial concentration 把指定 `Y[:, k]` 在 **t=0** 設為該軸的值，套用至所有細胞；其他物質維持原本初始設定。它不是在整段模擬中固定濃度，也不是改變觀察物質的顯示範圍。兩軸不可選同一個目標；同一物質的初始濃度與一個模型參數可以同時掃描。

Cell movement strength 可掃描隨機移動、遠離、匯聚、排斥或頂點模型的強度。掃描值會套用到每次試驗，即使模板原先關閉該種移動。強度 0 會關閉該種移動。

每軸 2–200 個點，最小值必須小於最大值。Logarithmic 要求兩端皆為正值。熱圖的欄對應 X，列對應 Y；較小的 Y 在下方。對數軸以對數座標呈現，色塊依實際取樣位置排列。每像素獨立重建格點並執行完整的 ODE、移動、分裂與凋亡設定，不會把前一像素的終點當成下一像素的初始條件。

總模擬次數 = X Samples × Y Samples × Replicates / pixel。建議先以 5×5 或 10×10、較短 T 試跑，再決定解析度。掃描期間暫停一般模擬／匯出入口與模型修改，避免共用亂數狀態互相影響。

## Response 定義

所有空間統計使用該影格**全部存活細胞（含邊界細胞）**，沒有依面積加權。下表的 S 指所選物質，`mean_cells` 為空間平均。

| Statistic | 單次試驗的數值 |
| --- | --- |
| Final spatial mean | 最後影格的 `mean_cells(S)` |
| Final spatial standard deviation | 最後影格細胞間標準差，`ddof=0` |
| Final spatial CV | 最後影格的 `std_cells(S) / abs(mean_cells(S))`；平均絕對值 ≤ 1e-12 時記為失敗 |
| Tail mean of spatial mean | 最後指定比例影格，各影格 `mean_cells(S)` 的時間平均 |
| Tail peak-to-peak of spatial mean | 同一末段時間窗內，空間平均時間序列的最大值減最小值；至少需要兩個影格 |
| First rise: maximum slope | 第一段明顯上升內的最大五幀回歸斜率 |
| First rise: tangent onset | 將該段最陡切線回推至上升前基線所得的時間 |
| Oscillation frequency | 規律波峰間距的中位數取倒數，單位為週期數／模型時間 |

Tail window 預設 25%，取最後 `ceil(影格數 × 比例)` 個影格。斜率、起始時間、頻率皆使用完整的空間平均時間序列。不同細胞反相振盪時，波形可能互相抵消。無法辨識的動態指標在 CSV 留空。演算法、門檻及七份模板的驗證結果見 [Agent CLI 說明](agent-cli.md)。

每個像素的顏色為所有 replicate response 的平均；CSV 的 `pixel_std` 是跨 replicate 的母體標準差（`ddof=0`），不是標準誤。任何 replicate 出錯，該像素整體顯示灰色 NaN，不會默默只平均成功試驗。CSV 的 `error` 欄保留原因。非有限值、無存活細胞、未定義 CV 等皆視為失敗。

## 重現與匯出

相同 replicate 在每個像素使用同一組 NumPy seed：`(seed + replicate_index) % 2**32`（index 從 0 起算），用於比較參數效應。每次 run 都建立獨立的初始狀態與格點；完成／取消後還原 NumPy 全域亂數狀態。這可重現使用 `np.random` 的模型；自行建立的 `default_rng()` 或外部隨機來源需要模型另外管理。

CSV 每列是一個 replicate，含 `x, y, replicate, seed, response, error, final_mean, rise_slope, onset_time, frequency, pixel_mean, pixel_std`。JSON 含兩軸、指標、T、原始 Params／ODE、實際起始細胞座標、移動／生命週期設定、覆寫值及模型警告。`config.json` 也保存 Experiments 控制項，重新開啟可恢復設定。掃描結果使用完成時的設定，不會因後續編輯控制項而改寫其座標或匯出內容。

## 數值解讀限制

Run 與 Experiments 使用相同的 explicit Euler 引擎。時間點為 `np.arange(0, T, dT)`，最後一個影格早於 T；Final 指最後儲存影格，並不宣稱已達穩態。生化導數在目前引擎中以新取樣時間評估；有時間相依刺激時請核對取樣定義。

濃度預設會截斷負值，過大的 dT 可能仍產生看似有限的錯誤結果。解讀熱圖前應比較 `dT` 與 `dT/2`、增加 T，並檢查模型警告、邊界細胞影響與 replicate 變異。此功能提供參數反應分析，不會自動判定 bifurcation、穩定性或重現論文的相圖。

## English quick reference

Open **Experiments** immediately after Run. Choose two distinct model parameters, species initial concentrations or cell movement strengths. Set the spacing and response statistic, then run in the background. Each pixel uses a fresh simulation with the current grid, dynamics, ODE and Run duration. Final statistics include all living boundary and interior cells. The replicate CSV includes the final spatial mean, rise slope, onset time and oscillation frequency. Replicates use common NumPy seeds across pixels. Any failed replicate makes the pixel NaN (gray). Export provides editable PDF/SVG, clipboard PNG, replicate CSV and reproducibility settings JSON.
