# Voronoi 細胞格模擬 + Simple ODE 圖形化介面

---

## ⚠️ 安全警告 (英文原文)

**此應用程式允許使用者在 GUI 中直接輸入並執行任意 Python 程式碼，用於 ODE 與顏色函數。**

- **請勿在任何公開、共用或不受信任的環境中執行此程式。**
- 任何在 ODE 或顏色函數文字框中輸入的程式碼，將與執行此程式的使用者擁有相同權限執行。
- 惡意或錯誤的程式碼可能：
  - 刪除或修改檔案
  - 存取敏感資料
  - 造成系統當機或凍結
  - 安裝或執行其他程式
- 目前版本 **未實作沙箱機制或程式碼限制**。

**在執行任何模擬之前，務必仔細檢閱您於 ODE 與顏色函數欄位中輸入的程式碼。**

> 本程式僅供在可受信任的本地環境中使用。若要公開部署，需自行新增程式碼沙箱與其他安全保護措施。

---

一個基於 PyQt6 的現代化互動 GUI，用於模擬與視覺化基於 Voronoi 鋪網的細胞 Notch 訊號模型。\\
此工具允許您設計、預覽並模擬各種細胞佈局與移動規則，並將結果匯出為影片與圖表。

![image](https://github.com/user-attachments/assets/5ef617fc-364e-4fae-ac55-20dfa61ec0b4)
---

## 功能介紹

- **單一視窗介面**：左側面板集中所有設定，右側面板即時預覽與結果。
- **彈性細胞網格生成**：蜂巢、隨機、規則網格、自訂或檔案匯入。
- **互動式預覽**：可用滑鼠拖曳、加入或刪除細胞；Voronoi 圖會即時更新。
- **多種細胞移動模式**：支援隨機移動、遠離中心、趨向延伸與鄰域排斥等組合。
- **ODE 與參數編輯**：可直接在 GUI 中修改 ODE 與模擬參數。
- **模擬與動畫**：執行完整模擬並預覽或匯出動畫與濃度曲線圖。
- **批次下載**：一次匯出所有結果（參數檔、影片、圖表）到指定資料夾。
- **狀態列與進度條**：隨時掌握當前狀態與進度。
- **健全的日誌系統**：所有操作與錯誤均會記錄於檔案並顯示在 GUI 中。

---

## 安裝方式

1. **複製此倉庫**

   ```bash
   git clone https://github.com/LBK888/Voronoi_Cell_ODE
   cd voronoi_Notch
   ```

2. **安裝相依套件**\
   建議使用虛擬環境：

   ```bash
   pip install -r requirements.txt
   ```

---

## 關於動畫匯出 (mp4)

**若要將動畫匯出為 .mp4，必須先於系統安裝 **[**FFMPEG**](https://ffmpeg.org/)**。**

- **Windows**：至官方網站下載 FFMPEG，並將 ffmpeg/bin 路徑加入系統 PATH。
- **macOS**：使用 Homebrew 安裝：`brew install ffmpeg`
- **Linux**：透過套件管理工具安裝，例如：`sudo apt install ffmpeg`

若未安裝或未設置於 PATH，matplotlib 將無法輸出 mp4，匯出時會出現錯誤。

---

## 如何啟動

1. **啟動 GUI**

   for Windows:  

   ```bash
   start.bat
   ```

   或：

   ```bash
   cd gui
   python -m gui.main_window
   ```

2. **首次啟動**

   - 主視窗將開啟，左側為所有設定，右側為預覽區。

---

## 使用手冊

### 1. **細胞網格設定**

- **Grid Mode（網格模式）**：可選蜂巢 (Honeycomb)、隨機 (Random)、規則網格 (Regular Grid)、自訂 (User-defined) 或匯入檔案 (Import from file)。
- **網格參數**：設定網格形狀、細胞間距與位置隨機度。

### 2. **細胞移動模式**

- 勾選一或多種移動規則：
  - 隨機移動 (Basic Random Movement)
  - 遠離中心 (Move Away From Center)
  - 趨向延伸 (Convergent Extension)
  - 鄰域排斥 (Repulsion: Neighbors, No Outer)
  - 細胞分裂 (Cell Division)（僅模擬中有效，預覽暫不支援互動）
  - 細胞凋亡 (Apoptosis)（同上）

### 3. **ODE 與參數**

- 直接於文字框中編輯 ODE 函數與模擬參數。

### 4. **模擬參數**

- 設定模擬時間 (`T`)、重複次數 (replicate) 與重複組數 (repeats)。

### 5. **預覽與編輯**

- 按下 **Preview** 生成並顯示目前的細胞佈局。
- 在預覽模式中可：
  - **拖曳** 細胞至新位置
  - **左鍵點擊** 新增細胞
  - **右鍵點擊** 刪除細胞
  - Voronoi 圖將即時更新

### 6. **套用移動**

- 按下 **Apply Movement**，將所選移動規則套用至目前細胞佈局。

### 7. **執行模擬**

- 按下 **Run Simulation**，以目前參數執行完整模擬。
- 模擬完成後會自動預覽動畫。

### 8. **下載結果**

- 按下 **Download Results** 匯出：
  - 參數檔 (`params.txt`)
  - 動畫影片 (`simulation.mp4`)
  - 濃度圖 (`concentration.pdf`)
- 選擇欲儲存的資料夾。

### 9. **狀態與進度**

- 狀態列（最下方）顯示當前操作、錯誤訊息與日誌。
- 進度條顯示模擬與匯出進度。

---

## 常見問題 (FAQ)

**Q: GUI 無法啟動或啟動即當機？**\
A: 請確認已安裝所有相依套件 (`pip install -r requirements.txt`)，且 Python 版本為 3.7 以上。

**Q: 如何匯入自訂細胞座標？**\
A: 在「Grid Mode」選擇「Import from file」，選擇含有兩欄 (x, y) 的 CSV 或 TXT 檔。

**Q: 是否可組合多種移動規則？**\
A: 可以！勾選多種移動模式，它們將按順序依次套用。

**Q: 細胞分裂與凋亡為何無法在預覽中互動？**\
A: 這些功能在模擬引擎中已支援，但尚未整合至互動預覽中，可於完整模擬中使用。

**Q: 日誌檔案儲存位置？**\
A: 所有日誌皆儲存於專案根目錄的 `app.log`。

**Q: 模擬速度過慢或 GUI 冻結？**\
A: 大型網格或長時間模擬會耗時，未來版本可能會新增背景執行支援。

**Q: 如何重置為預設 ODE 或參數？**\
A: 重新啟動 GUI 即可載入預設值。

---

## 授權條款

本專案採 CC BY-NC 4.0 授權。\
詳見 [LICENSE](https://creativecommons.org/licenses/by-nc/4.0/) 。

---

## 聯絡方式
- Email：liaobk@email.ntou.edu.tw
- Lab Web page： http://lbk.tw/
有任何問題、回報或功能建議，歡迎於 GitHub 開 issue 或聯絡維護者。

