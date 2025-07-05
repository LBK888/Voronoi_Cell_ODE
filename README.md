# Voronoi Cell Simulation + Simple ODE GUI  
  [中文說明](readme_tw.md)
---

## ⚠️ Security Warning (English)

**This application allows users to input and execute arbitrary Python code for ODE and color functions directly in the GUI.**

- **Do NOT use this program in any public, shared, or untrusted environment.**
- Any code entered in the ODE or color function text boxes will be executed with the same permissions as the user running the program.
- Malicious or incorrect code can:
  - Delete or modify files
  - Access sensitive data
  - Crash or freeze your system
  - Install or run other programs
- There is NO sandboxing or code restriction in the current version.

**Before running any simulation, always carefully review the code you enter in the ODE and color function fields.**

> This program is intended for local, private, and trusted use only. Public deployment is unsafe unless you implement additional code sandboxing and security measures.

---

A modern, interactive PyQt6-based GUI for simulating and visualizing cell-based Notch signaling models using Voronoi tessellation.  
This tool allows you to design, preview, and simulate various cell arrangements and movement rules, and export results as videos and figures.

---

## Features

- **Single-window GUI**: Left panel for all settings, right panel for interactive preview and results.
- **Flexible cell grid generation**: Honeycomb, random, regular grid, user-defined, or import from file.
- **Interactive preview**: Drag, add, or delete cells with your mouse; Voronoi diagram updates in real time.
- **Multiple cell movement modes**: Combine random movement, move away from center, convergent extension, and neighbor repulsion.
- **ODE and parameter editing**: Directly edit the ODE and simulation parameters in the GUI.
- **Simulation and animation**: Run full simulations and preview or export animations and concentration plots.
- **Batch download**: Export all results (parameters, video, figures) to a chosen folder.
- **Status bar and progress bar**: Always know the current state and progress.
- **Robust logger**: All actions and errors are logged to file and shown in the GUI.

---

## Installation

1. **Clone this repository**  
   ```bash
   git clone https://github.com/LBK888/Voronoi_Cell_ODE
   cd voronoi_Notch
   ```

2. **Install dependencies**  
   It is recommended to use a virtual environment.
   ```bash
   pip install -r requirements.txt
   ```

---

## Additional Notes on Animation Export (mp4)

**To export animation videos as .mp4 files, you must have [FFMPEG](https://ffmpeg.org/) installed on your system.**

- On Windows: Download FFMPEG from the [official website](https://ffmpeg.org/download.html) and add the ffmpeg/bin directory to your system PATH.
- On Mac: Install via Homebrew: `brew install ffmpeg`
- On Linux: Install via your package manager, e.g., `sudo apt install ffmpeg`

If FFMPEG is not installed or not in your PATH, matplotlib will not be able to save animations as mp4 files and you may encounter errors during export.

---

## How to Run

1. **Start the GUI**  
   
   For Windows OS:
   ```bash
   start.bat
   ```
   Or, 
   ```bash
   python -m gui.main_window
   ```

2. **First launch**  
   - The main window will open with all settings on the left and a preview area on the right.

---

## User Manual

### 1. **Cell Grid Settings**
- **Grid Mode**: Choose from Honeycomb, Random, Regular Grid, User-defined, or Import from file.
- **Grid Parameters**: Set grid shape, cell distance, and position randomness.

### 2. **Cell Movement Modes**
- Check one or more movement rules:
  - Basic Random Movement
  - Move Away From Center
  - Convergent Extension
  - Repulsion (Neighbors, No Outer)
  - Cell Division (checkbox, not yet interactive in preview)
  - Apoptosis (checkbox, not yet interactive in preview)

### 3. **ODE and Parameters**
- Edit the ODE function and simulation parameters directly in the provided text boxes.

### 4. **Simulation Parameters**
- Set simulation time (`T`), replicate, and repeats.

### 5. **Preview and Edit**
- Click **Preview** to generate and display the current cell arrangement.
- In preview mode, you can:
  - **Drag** cells to new positions.
  - **Left-click** to add a new cell.
  - **Right-click** to delete a cell.
  - The Voronoi diagram updates instantly.

### 6. **Apply Movement**
- Click **Apply Movement** to apply the selected movement rules to the current cell arrangement.

### 7. **Run Simulation**
- Click **Run Simulation** to execute the full simulation with the current settings.
- The animation will be previewed after simulation.

### 8. **Download Results**
- Click **Download Results** to export:
  - Simulation parameters (`params.txt`)
  - Animation video (`simulation.mp4`)
  - Concentration plot (`concentration.pdf`)
- Choose your target folder in the dialog.

### 9. **Status and Progress**
- The status bar (bottom) shows current actions, errors, and logger messages.
- The progress bar indicates simulation and export progress.

---

## FAQ

**Q: The GUI does not start or crashes on launch.**  
A: Make sure you have installed all dependencies (`pip install -r requirements.txt`). PyQt6 requires Python 3.7+.

**Q: How do I import my own cell coordinates?**  
A: Choose "Import from file" in Grid Mode, then select a CSV or TXT file with two columns (x, y).

**Q: Can I combine multiple movement rules?**  
A: Yes! You can check multiple movement modes and they will be applied in sequence.

**Q: Why are cell division and apoptosis not interactive in preview?**  
A: These features are available in the simulation engine, but not yet in the interactive preview. You can still use them in full simulations.

**Q: Where are logs saved?**  
A: All logs are saved to `app.log` in the project root.

**Q: The simulation is slow or the GUI freezes.**  
A: For large grids or long simulations, computation may take time. Future versions may support background threading.

**Q: How do I reset to default ODE or parameters?**  
A: Simply restart the GUI; default values will be loaded.

---

## License

MIT License.  
See [LICENSE](LICENSE) for details.

---

## Contact
- Email：liaobk@email.ntou.edu.tw
- Lab Web page： http://lbk.tw/

For questions, bug reports, or feature requests, please open an issue or contact the maintainer. 
