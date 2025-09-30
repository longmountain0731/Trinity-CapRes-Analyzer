# Trinity CapRes Analyzer  
*A GUI Tool for I–V / C–V Curve Visualization and CTLM Analysis*

---

## 📌 Overview
**Trinity CapRes Analyzer** is a graphical analysis tool designed for semiconductor electrical measurements (I–V, C–V, CTLM).  

It supports **Keysight B1500 CSV format**, automatically detects **Single / Double sweeps**, and provides:  
- I–V curve visualization and overlay  
- Differential R–V extraction  
- R₀ vs Spacing fitting (**Method-1 & Method-2**)  
- Specific contact resistivity (**ρc**) calculation  
- C–V curve plotting by frequency  
- High-resolution figure export (**PNG**)  

---

## 🚀 Usage

### 1. Launch
- **Windows**: run `TrinityCapResAnalyzer.exe`  
- **macOS**: run `Trinity CapRes Analyzer.app`  

### 2. Load Data
At startup, you will see a **Data Type Selection** window:  
- **Single CSV (multiple sweeps included)** → for B1500 data with multi-frequency C–V sweeps or double-sweep I–V.  
- **Multiple CSVs (one sweep per file)** → for batch loading separate files.  

The program will automatically parse and group sweeps.  

### 3. Interface
- **Global Width / Spacing**  
  - Enter or edit structure spacing (μm).  
  - Values will be applied to **R–Spacing fitting**.  

- **Sweeps Table**  
  - **Use**: enable/disable curve display  
  - **Follow**: sync with global spacing label  
  - **Label**: custom sweep name  
  - **Color / Line / Marker**: curve style  
  - **V Range / Y Range**: value ranges per sweep  

- **Right-Side Notebook Tabs**  
  - **I–V** → Current–Voltage curves  
  - **R–V** → Differential resistance curves  
  - **R–Spacing** → R₀ fitting & ρc extraction  
  - **R–Spacing Correlation** → Method-2 linearized correction  
  - **C–V** → Capacitance–Voltage curves (by frequency)  

### 4. Export
Click **Export** to save all selected sweeps with overlays.  

Generated files include:  
- `IV_overlay_selected.png`  
- `RV_overlay_selected.png`  
- `CV_overlay_selected.png`  
- `R_spacing_correlation.png`  
- `rho_summary.txt` (ρc fitting summary)  

---

## ⚙️ System Requirements
- Windows or macOS  
- **Python runtime** (for `.py` version) OR **standalone executable** (no Python required)  
- Recommended resolution: **1920×1080** or higher  

---

## 📖 Notes
- For **double sweep data**, curves are connected in acquisition order to preserve hysteresis.  
- All exported figures are **publication-ready (DPI ≥ 300)**.  
- If preview panels look distorted, adjust **Fig W / Fig H** or check **monitor scaling**.  

---

## ✨ Developer
- **Trinity CapRes Analyzer**  
- Version: **v1.0**  
- Author: **Chang-Shan Shen**  

---
