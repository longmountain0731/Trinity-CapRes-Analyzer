

from pathlib import Path
import tkinter as tk
from tkinter import ttk, filedialog, messagebox, colorchooser
from typing import List, Dict, Tuple, Optional
import numpy as np
import matplotlib
matplotlib.use("TkAgg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib import rcParams
from matplotlib.ticker import AutoMinorLocator

APP_TITLE = "Trinity CapRes Analyzer"

# ---- Matplotlib base style ----
rcParams['font.family'] = 'Times New Roman'
rcParams['axes.titleweight'] = 'bold'
rcParams['axes.labelweight'] = 'bold'

UM_TO_CM = 1e-4  # μm → cm

# ---------- Dark style + font scaling ----------
def apply_dark_style(root, base_font=("Segoe UI", 12)):
    root.configure(bg="#0f1216")
    style = ttk.Style(root)
    try:
        style.theme_use("clam")
    except Exception:
        pass
    # Global font for all ttk widgets
    style.configure(".", font=base_font)  # every widget
    # Colors
    fg = "#e9eef5"
    bg = "#0f1216"
    panel = "#161a20"
    accent = "#1f6feb"
    entry_bg = "#11161d"
    border = "#273043"
    style.configure("TFrame", background=bg)
    style.configure("TLabelframe", background=panel, foreground=fg)
    style.configure("TLabelframe.Label", background=panel, foreground=fg, font=(base_font[0], base_font[1], "bold"))
    style.configure("TLabel", background=bg, foreground=fg)
    style.configure("TCheckbutton", background=bg, foreground=fg)
    style.configure("TRadiobutton", background=bg, foreground=fg)
    style.configure("TEntry", fieldbackground=entry_bg, foreground=fg, background=entry_bg)
    style.configure("TCombobox", fieldbackground=entry_bg, foreground=fg, background=entry_bg)
    style.map("TCombobox", fieldbackground=[('readonly', entry_bg)])
    style.configure("TButton", background=panel, foreground=fg, bordercolor=border, focusthickness=1)
    style.map("TButton", background=[("active", accent)])
    style.configure("TNotebook", background=bg, tabmargins= [4, 4, 4, 0])
    style.configure("TNotebook.Tab", background=panel, foreground=fg, padding=[10,5], font=(base_font[0], base_font[1]-1, "bold"))
    style.map("TNotebook.Tab", background=[("selected", accent)])
    # Canvas default bg
    root.option_add("*Canvas.background", panel)
    root.option_add("*Text.background", panel)
    root.option_add("*Text.foreground", fg)
    root.option_add("*Listbox.background", panel)
    root.option_add("*Listbox.foreground", fg)

# ---------- Splash ----------
class Splash(tk.Toplevel):
    def __init__(self, master, *, secs=1.2):
        super().__init__(master)
        self.overrideredirect(True)
        self.configure(bg="#0f1216")
        w, h = 520, 260
        self.geometry(f"{w}x{h}+{self.winfo_screenwidth()//2 - w//2}+{self.winfo_screenheight()//2 - h//2}")
        # Branding
        title = tk.Label(self, text="TRINITY CAPRES ANALYZER", fg="#e9eef5", bg="#0f1216",
                         font=("Eurostile", 20, "bold"))
        subtitle = tk.Label(self, text="CTLM I–V · C–V  —  fast, precise, beautiful",
                            fg="#8aa5d8", bg="#0f1216", font=("Segoe UI", 11))
        title.pack(pady=(36,8))
        subtitle.pack()
        # Progress bar
        bar = ttk.Progressbar(self, mode="indeterminate", length=360)
        bar.pack(pady=28)
        bar.start(16)
        # Footer
        foot = tk.Label(self, text="Loading modules…", fg="#95a7c6", bg="#0f1216", font=("Segoe UI", 10))
        foot.pack(side="bottom", pady=10)
        # Auto close
        self.after(int(secs*1000), self._close)

    def _close(self):
        try:
            self.destroy()
        except Exception:
            pass

# ---------- 小工具 ----------
def rainbow_color(i: int, n: int):
    cmap = plt.get_cmap('rainbow')
    return matplotlib.colors.to_hex(cmap(i / max(n-1, 1)))

def safe_float(s, default=None):
    try:
        return float(s)
    except Exception:
        return default

def parse_numeric_from_label(lbl: str):
    import re
    m = re.search(r'[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?', lbl)
    return float(m.group(0)) if m else None

def _fmt_freq(f: float) -> str:
    if f >= 1e6:  return f"{f/1e6:g} MHz"
    if f >= 1e3:  return f"{f/1e3:g} kHz"
    return f"{f:g} Hz"

# ---------- CSV 解析（鍵在第2欄，值從第3欄起） ----------
def _find_header_params(lines):
    def first_float(tokens):
        for t in tokens:
            t = t.strip()
            try:
                return float(t)
            except Exception:
                pass
        return None

    locus = "single"; vstart=None; vstop=None; freqs=[]
    for ln in lines[:800]:
        parts = [p.strip() for p in ln.split(",")]
        if len(parts) < 2: continue
        key = parts[1].lower()
        if key == "measurement.primary.locus":
            for t in parts[2:]:
                if t:
                    locus = t.lower(); break
        elif key == "measurement.primary.start":
            vstart = first_float(parts[2:])
        elif key == "measurement.primary.stop":
            vstop = first_float(parts[2:])
        elif key == "measurement.secondary.frequency":
            freqs = []
            for t in parts[2:]:
                try: freqs.append(float(t))
                except Exception: pass
    return locus, vstart, vstop, (freqs if freqs else None)

def _find_dimension1_near(lines, dataname_line_index):
    for look in range(1, 60):
        j = dataname_line_index - look
        if j < 0: break
        parts = [p.strip() for p in lines[j].split(",")]
        if len(parts) >= 2 and parts[1].lower() == "dimension1":
            for t in parts[2:]:
                try: return int(float(t))
                except Exception: pass
            break
    return None

def _split_by_locus_or_wrap(V, locus, vstart, vstop, n_expected=None, npts_hint=None):
    if npts_hint and npts_hint > 2 and len(V) >= npts_hint:
        total = len(V); nseg = total // npts_hint
        if n_expected: nseg = min(nseg, n_expected)
        slices = [slice(k*npts_hint, (k+1)*npts_hint) for k in range(nseg)]
        if slices: return slices
    slices = []
    if vstart is not None and vstop is not None:
        span = abs(vstop - vstart); tol = max(1e-9, 0.01*span)
        if str(locus).startswith("double"):
            s = 0; seen_stop = False
            for k, v in enumerate(V):
                if not seen_stop and abs(v - vstop) <= tol: seen_stop = True
                elif seen_stop and abs(v - vstart) <= tol:
                    if k - s >= 3: slices.append(slice(s, k))
                    s = k; seen_stop = False
                    if n_expected and len(slices) >= n_expected: break
            if (not n_expected or len(slices) < n_expected) and len(V) - s >= 3:
                slices.append(slice(s, len(V)))
        else:
            s = 0
            for k in range(1, len(V)):
                if abs(V[k-1] - vstop) <= tol and abs(V[k] - vstart) <= tol:
                    if k - s >= 3: slices.append(slice(s, k))
                    s = k
                    if n_expected and len(slices) >= n_expected: break
            if (not n_expected or len(slices) < n_expected) and len(V) - s >= 3:
                slices.append(slice(s, len(V)))
    if not slices:
        if len(V) < 2: return [slice(0, len(V))]
        jumps = np.where(np.diff(V) < 0)[0]
        starts = np.r_[0, jumps + 1]; ends = np.r_[jumps + 1, len(V)]
        slices = [slice(s, e) for s, e in zip(starts, ends)]
    if n_expected and len(slices) > n_expected: slices = slices[:n_expected]
    return slices

def parse_b1500_csv_text(txt: str):
    lines = [ln.rstrip("\n") for ln in txt.splitlines() if ln.strip()]
    locus, vstart, vstop, freqs = _find_header_params(lines)
    curves = []; i = 0; sweep_idx = 0

    def find_idx(names_low, cands):
        for c in cands:
            if c in names_low: return names_low.index(c)
        return None

    while i < len(lines):
        line = lines[i].strip()
        if not line.lower().startswith("dataname"):
            i += 1; continue
        header = [h.strip() for h in line.split(",")]
        names = [h for h in header[1:] if h]
        names_low = [n.lower() for n in names]
        is_two_cols = (len(names) == 2)

        v_iv_idx = find_idx(names_low, ["vd","v_d","v drain","v"])
        i_iv_idx = find_idx(names_low, ["id","i_d","i drain","i"])
        v_cv_idx = find_idx(names_low, ["vbias","v","vd","vg","v gate"])
        c_cv_idx = find_idx(names_low, ["c","cap","capacitance"])

        i += 1
        block = []
        while i < len(lines) and lines[i].strip().lower().startswith("datavalue"):
            block.append([x.strip() for x in lines[i].split(",")][1:])
            i += 1
        if not block: continue
        arr = np.array(block, dtype=object)

        # I–V
        try:
            if not is_two_cols and v_iv_idx is not None and i_iv_idx is not None:
                V = np.asarray(arr[:, v_iv_idx], float); I = np.asarray(arr[:, i_iv_idx], float)
                sweep_idx += 1; curves.append(dict(label=f"Sweep_{sweep_idx}", V=V, I=I, type="iv")); continue
            elif is_two_cols and v_iv_idx is None and c_cv_idx is None:
                V = np.asarray(arr[:, 0], float); I = np.asarray(arr[:, 1], float)
                sweep_idx += 1; curves.append(dict(label=f"Sweep_{sweep_idx}", V=V, I=I, type="iv")); continue
        except Exception:
            pass

        # C–V
        cv_possible = False
        try:
            if not is_two_cols and v_cv_idx is not None and c_cv_idx is not None:
                V_all = np.asarray(arr[:, v_cv_idx], float); C_all = np.asarray(arr[:, c_cv_idx], float); cv_possible = True
            elif is_two_cols and c_cv_idx is not None:
                V_all = np.asarray(arr[:, 0], float); C_all = np.asarray(arr[:, 1], float); cv_possible = True
        except Exception:
            cv_possible = False
        if not cv_possible: continue

        n_expected = len(freqs) if freqs else None
        npts_hint = _find_dimension1_near(lines, i - len(block) - 1)
        slices = _split_by_locus_or_wrap(V_all, locus, vstart, vstop, n_expected, npts_hint)

        for k, sl in enumerate(slices):
            V = V_all[sl]; C = C_all[sl]
            if freqs is not None and k < len(freqs):
                f = freqs[k]; lbl = _fmt_freq(float(f))
            else:
                lbl = f"CV_{sweep_idx+1}"
            sweep_idx += 1
            curves.append(dict(label=lbl, V=V, C=C, type="cv"))
    return curves

def read_curves_from_file(path: Path):
    try:
        txt = path.read_text(encoding="utf-8")
    except Exception:
        txt = path.read_text(errors="ignore")
    cs = parse_b1500_csv_text(txt)
    if len(cs) == 1:
        cs[0]["label"] = path.stem
    return cs

def read_curves_single_csv(path: Path): return read_curves_from_file(path)
def read_curves_multi_files(paths: List[Path]):
    out = []
    for p in paths:
        cs = read_curves_from_file(p)
        for c in cs:
            if len(cs) == 1: c["label"] = p.stem
            out.append(c)
    return out



def make_scrollable(parent):
    """回傳 (canvas, inner_frame)。把 inner_frame 當成原本的 parent 用來 pack/grid。"""
    canvas = tk.Canvas(parent, highlightthickness=0)
    vsb = ttk.Scrollbar(parent, orient="vertical", command=canvas.yview)
    canvas.configure(yscrollcommand=vsb.set)
    vsb.pack(side="right", fill="y")
    canvas.pack(side="left", fill="both", expand=True)
    inner = ttk.Frame(canvas)
    canvas.create_window((0, 0), window=inner, anchor="nw")
    inner.bind("<Configure>", lambda e: canvas.configure(scrollregion=canvas.bbox("all")))

    # 滑鼠滾輪捲動
    def _on_mousewheel(event):
        canvas.yview_scroll(int(-1*(event.delta/120)), "units")
    canvas.bind_all("<MouseWheel>", _on_mousewheel)

    return canvas, inner


# --------------- Main App ---------------
class App(tk.Tk):
    def __init__(self, base_font=None):
        super().__init__()

        # --- DPI/解析度縮放設定 ---
        try:
            import ctypes
            ctypes.windll.shcore.SetProcessDpiAwareness(2)  # Windows: per-monitor DPI
        except Exception:
            pass

        scaling = self.winfo_fpixels('1i') / 72.0
        self.tk.call('tk', 'scaling', scaling)

        sw, sh = self.winfo_screenwidth(), self.winfo_screenheight()
        self.minsize(int(sw*0.75), int(sh*0.70))
        self.geometry(f"{int(sw*0.9)}x{int(sh*0.9)}+0+0")
        try:
            self.state('zoomed')  # Windows
        except Exception:
            pass

        # --- 全域字體（ttk + tk）---
        if base_font:
            style = ttk.Style(self)
            # ttk 常見元件
            style.configure(".", font=base_font)
            style.configure("TLabel", font=base_font)
            style.configure("TButton", font=base_font)
            style.configure("TCheckbutton", font=base_font)
            style.configure("TRadiobutton", font=base_font)
            style.configure("TEntry", font=base_font)
            style.configure("TCombobox", font=base_font)
            style.configure("TNotebook", font=base_font)
            style.configure("TNotebook.Tab", font=base_font)
            style.configure("Treeview", font=base_font)
            style.configure("Treeview.Heading", font=base_font)

            # tk 元件（Text、Canvas 內等）
            self.option_add("*Font", base_font)
            # Combobox 下拉清單字體
            self.option_add("*TCombobox*Listbox*Font", base_font)

        self.withdraw()  # 先隱藏，等 splash 結束
        apply_dark_style(self)
        self.title(APP_TITLE)

        self.curves: List[Dict] = []
        self.rows = []
        self.csv_path: Optional[Path] = None
        self.file_list: List[Path] = []
        self.outdir: Optional[Path] = None
        self.data_mode: Optional[str] = None

        Splash(self)
        self.after(1200, self._post_splash)

    def _post_splash(self):
        self.deiconify()
        self._build_ui()
        self.after(120, self.startup_wizard)

    # ----- 開場精靈 -----
    def startup_wizard(self):
        win = tk.Toplevel(self); win.title("Choose data type"); win.grab_set()
        ttk.Label(win, text="請選擇資料型態：", padding=10).pack()
        btns = ttk.Frame(win, padding=10); btns.pack(fill="x")
        ttk.Button(btns, text="單一 CSV（含多組 Sweep）",
                   command=lambda: self._choose_single_csv(win)).pack(fill="x", pady=4)
        ttk.Button(btns, text="多個 CSV（各一組 Sweep）",
                   command=lambda: self._choose_multi_files(win)).pack(fill="x", pady=4)

    def _choose_single_csv(self, win):
        p = filedialog.askopenfilename(title="選擇 CSV", filetypes=[("CSV","*.csv"),("All","*.*")])
        if not p: return
        win.destroy()
        self.data_mode = 'single_csv_multi'
        self.csv_path = Path(p)
        self.file_list = [self.csv_path]
        self.load_curves_from_selection()

    def _choose_multi_files(self, win):
        sub = tk.Toplevel(win); sub.title("選擇來源"); sub.grab_set()
        ttk.Label(sub, text="選擇方式：", padding=10).pack()
        fr = ttk.Frame(sub, padding=10); fr.pack(fill="x")
        def pick_files():
            files = filedialog.askopenfilenames(title="選擇多個 CSV", filetypes=[("CSV","*.csv"),("All","*.*")])
            if not files: return
            sub.destroy(); win.destroy()
            self.data_mode = 'multi_files_single'
            self.file_list = [Path(f) for f in files]
            self.csv_path = self.file_list[0]
            self.load_curves_from_selection()
        def pick_folder():
            d = filedialog.askdirectory(title="選擇資料夾")
            if not d: return
            files = sorted(Path(d).glob("*.csv"))
            if not files:
                messagebox.showwarning("提示", "該資料夾沒有 CSV 檔。"); return
            sub.destroy(); win.destroy()
            self.data_mode = 'multi_files_single'
            self.file_list = list(files)
            self.csv_path = self.file_list[0]
            self.load_curves_from_selection()
        ttk.Button(fr, text="多選檔案", command=pick_files).pack(fill="x", pady=4)
        ttk.Button(fr, text="選擇資料夾", command=pick_folder).pack(fill="x", pady=4)

    def load_curves_from_selection(self):
        try:
            if self.data_mode == 'single_csv_multi':
                curves = read_curves_single_csv(self.csv_path)
            else:
                curves = read_curves_multi_files(self.file_list)
        except Exception as e:
            messagebox.showerror("讀取失敗", str(e)); return
        self.curves = curves
        self._populate_rows()
        self.update_all_previews()
        self.log(f"Mode: {self.data_mode}, loaded {len(self.curves)} curves")

    # ----- UI -----
    def _build_ui(self):
        top = ttk.Frame(self, padding=8); top.pack(fill="x")
        self.src_var = tk.StringVar(value="未選擇")
        ttk.Label(top, text="Source:").pack(side="left")
        ttk.Label(top, textvariable=self.src_var, width=85).pack(side="left", padx=(0,6))
        ttk.Button(top, text="Choose again", command=self.startup_wizard).pack(side="left")
        ttk.Button(top, text="Output folder", command=self.on_choose_outdir).pack(side="left", padx=(6,0))
        ttk.Label(top, text="R2 (μm)").pack(side="left", padx=(12,4))
        self.r2_var = tk.StringVar(value="100")
        ttk.Entry(top, textvariable=self.r2_var, width=10).pack(side="left")

        main = ttk.Panedwindow(self, orient="horizontal"); main.pack(fill="both", expand=True, padx=8, pady=6)

        # left：全局 spacing + sweeps 表
        left = ttk.Frame(main); main.add(left, weight=1)

        gb = ttk.LabelFrame(left, text="Global Width/Spacing (1..9)", padding=6)
        gb.pack(fill="x", pady=(0,6))
        self.global_vars = [tk.StringVar(value=f"{i*10:.2f} um") for i in range(1,10)]
        row = ttk.Frame(gb); row.pack(fill="x")
        for i in range(9):
            ttk.Label(row, text=str(i+1)).grid(row=0, column=2*i, sticky="e")
            ttk.Entry(row, textvariable=self.global_vars[i], width=10).grid(row=0, column=2*i+1, padx=(2,8))
        ttk.Button(gb, text="Apply to table", command=self.update_all_previews).pack(anchor="e", pady=(6,0))

        tbl_frame = ttk.LabelFrame(left, text="Sweeps", padding=6)
        tbl_frame.pack(fill="both", expand=True)
        canvas = tk.Canvas(tbl_frame, highlightthickness=0)
        self.tbl = ttk.Frame(canvas)
        scroll = ttk.Scrollbar(tbl_frame, orient="vertical", command=canvas.yview)
        canvas.configure(yscrollcommand=scroll.set)
        scroll.pack(side="right", fill="y"); canvas.pack(side="left", fill="both", expand=True)
        self.tbl_id = canvas.create_window((0,0), window=self.tbl, anchor="nw")
        self.tbl.bind("<Configure>", lambda e: canvas.configure(scrollregion=canvas.bbox("all")))
        canvas.bind("<Configure>", lambda e: canvas.itemconfig(self.tbl_id, width=e.width))

        self.col_widths = [6, 4, 7, 18, 22, 12, 8, 8, 18, 18]
        hdr = ttk.Frame(self.tbl); hdr.grid(row=0, column=0, sticky="ew")
        heads = ["Use","#","Follow","Global#","Label","Color","Line","Marker","V Range","Y Range"]
        for j, w in enumerate(self.col_widths):
            ttk.Label(hdr, text=heads[j], width=w, anchor="w").grid(row=0, column=j, sticky="w")

        self.rows_container = ttk.Frame(self.tbl); self.rows_container.grid(row=1, column=0, sticky="nsew")
        for c, w in enumerate(self.col_widths):
            self.rows_container.grid_columnconfigure(c, minsize=int(w*8))

        # right notebook
        right_outer = ttk.Frame(main); main.add(right_outer, weight=1)
        _, right = make_scrollable(right_outer)  # right 變成可捲動的內層框
        nb = ttk.Notebook(right); nb.pack(fill="both", expand=True)


        # I–V
        self.tab_iv = ttk.Frame(nb); nb.add(self.tab_iv, text="I–V")
        self.preview_iv = self._build_preview_panel(self.tab_iv, "I–V Curves", "Voltage (V)", "Current (A)")

        # R–V
        self.tab_rv = ttk.Frame(nb); nb.add(self.tab_rv, text="R–V")
        self.preview_rv = self._build_preview_panel(self.tab_rv, "Differential Resistance R(V)", "Voltage (V)", "Resistance (Ω)")
        rv_cfg = ttk.LabelFrame(self.tab_rv, text="R0 fit window (|V| ≤ window)", padding=6)
        rv_cfg.pack(fill="x", padx=8, pady=(0,8))
        self.r0_window = tk.StringVar(value="0.5")
        ttk.Label(rv_cfg, text="window (V)").pack(side="left")
        ttk.Entry(rv_cfg, textvariable=self.r0_window, width=8).pack(side="left", padx=(4,10))
        ttk.Button(rv_cfg, text="Refresh previews", command=self.update_all_previews).pack(side="right")

        # R–Spacing
        self.tab_rs = ttk.Frame(nb); nb.add(self.tab_rs, text="R–Spacing")
        self.preview_rs = self._build_preview_panel(self.tab_rs, "R0 vs Spacing", "Spacing d (μm)", "Resistance R0 (Ω)")

        # --- Rt(R0) 清單（可滾動） ---
        rt_frame = ttk.LabelFrame(self.tab_rs, text="各點 Rt(R0) 清單", padding=6)
        rt_frame.pack(fill="both", padx=8, pady=(0,4))
        
        rt_wrap = ttk.Frame(rt_frame)
        rt_wrap.pack(fill="both", expand=True)
        
        rt_scroll = ttk.Scrollbar(rt_wrap, orient="vertical")
        self.rt_text = tk.Text(rt_wrap, height=8, wrap="none", yscrollcommand=rt_scroll.set)
        rt_scroll.config(command=self.rt_text.yview)
        
        self.rt_text.pack(side="left", fill="both", expand=True)
        rt_scroll.pack(side="right", fill="y")
        
        # 滑鼠滾輪綁定
        self.rt_text.bind("<MouseWheel>", lambda e: self.rt_text.yview_scroll(-int(e.delta/120), "units"))  # Windows/macOS
        self.rt_text.bind("<Button-4>", lambda e: self.rt_text.yview_scroll(-1, "units"))  # Linux up
        self.rt_text.bind("<Button-5>", lambda e: self.rt_text.yview_scroll( 1, "units"))  # Linux down

        self.result_panel = ttk.LabelFrame(self.tab_rs, text="Fits & ρc results", padding=6); self.result_panel.pack(fill="both", padx=8, pady=(4,4))
        wrap = ttk.Frame(self.result_panel); wrap.pack(fill="both", expand=True)
        self.result_text = tk.Text(wrap, height=8, wrap="none"); self.result_text.pack(side="left", fill="both", expand=True)
        scr = ttk.Scrollbar(wrap, orient="vertical", command=self.result_text.yview); scr.pack(side="right", fill="y")
        self.result_text.configure(yscrollcommand=scr.set)

        # Correlation
        self.tab_corr = ttk.Frame(nb); nb.add(self.tab_corr, text="R–Spacing Correlation")
        self.preview_corr = self._build_preview_panel(self.tab_corr, "Original vs Corrected Rt(d)", "Spacing d (μm)", "Resistance (Ω)")
        self.corr_text = tk.StringVar(value="")
        ttk.Label(self.tab_corr, textvariable=self.corr_text).pack(anchor="w", padx=8, pady=(0,8))

        # C–V
        self.tab_cv = ttk.Frame(nb); nb.add(self.tab_cv, text="C–V")
        self.preview_cv = self._build_preview_panel(self.tab_cv, "C–V (per frequency)", "Voltage (V)", "Capacitance (F)")

        # bottom
        btns = ttk.Frame(self, padding=8); btns.pack(fill="x")
        ttk.Button(btns, text="Select all", command=self.select_all).pack(side="left")
        ttk.Button(btns, text="Clear", command=self.clear_all).pack(side="left", padx=(6,0))
        ttk.Button(btns, text="Refresh previews", command=self.update_all_previews).pack(side="left", padx=(6,0))
        ttk.Button(btns, text="Export", command=self.export_all).pack(side="right")
        self.status = tk.Text(self, height=5); self.status.pack(fill="both", padx=8, pady=(0,8))

    def _build_preview_panel(self, parent, title, xl, yl):
        p = {}
        top = ttk.LabelFrame(parent, text="Figure settings", padding=6); top.pack(fill="x", padx=8, pady=(8,6))
        p['title'] = tk.StringVar(value=title); p['xlabel'] = tk.StringVar(value=xl); p['ylabel'] = tk.StringVar(value=yl)
        p['x_auto'] = tk.BooleanVar(value=True); p['y_auto'] = tk.BooleanVar(value=True)
        p['xmin'] = tk.StringVar(); p['xmax'] = tk.StringVar(); p['ymin'] = tk.StringVar(); p['ymax'] = tk.StringVar()
        p['xscale'] = tk.StringVar(value="linear"); p['yscale'] = tk.StringVar(value="linear")
        p['dpi'] = tk.StringVar(value="300"); p['figw'] = tk.StringVar(value="6"); p['figh'] = tk.StringVar(value="4")
        p['show_legend'] = tk.BooleanVar(value=True); p['legend_size'] = tk.StringVar(value="14"); p['legend_ncol'] = tk.StringVar(value="1")
        p['title_size'] = tk.StringVar(value="20"); p['label_size'] = tk.StringVar(value="18")
        r1 = ttk.Frame(top); r1.pack(fill="x")
        ttk.Label(r1, text="Title").grid(row=0,column=0,sticky="e"); ttk.Entry(r1,textvariable=p['title'],width=38).grid(row=0,column=1,sticky="w",padx=4)
        ttk.Label(r1, text="X label").grid(row=0,column=2,sticky="e"); ttk.Entry(r1,textvariable=p['xlabel'],width=18).grid(row=0,column=3,sticky="w",padx=4)
        ttk.Label(r1, text="Y label").grid(row=0,column=4,sticky="e"); ttk.Entry(r1,textvariable=p['ylabel'],width=18).grid(row=0,column=5,sticky="w",padx=4)
        r1b = ttk.Frame(top); r1b.pack(fill="x", pady=(4,0))
        ttk.Label(r1b, text="Title size").grid(row=0,column=0,sticky="e"); ttk.Entry(r1b,textvariable=p['title_size'],width=6).grid(row=0,column=1,sticky="w",padx=(2,12))
        ttk.Label(r1b, text="Label size").grid(row=0,column=2,sticky="e"); ttk.Entry(r1b,textvariable=p['label_size'],width=6).grid(row=0,column=3,sticky="w",padx=(2,12))
        r2 = ttk.Frame(top); r2.pack(fill="x", pady=(4,0))
        ttk.Checkbutton(r2, text="X auto", variable=p['x_auto']).grid(row=0,column=0,sticky="w")
        ttk.Label(r2, text="xmin").grid(row=0,column=1); ttk.Entry(r2,textvariable=p['xmin'],width=8).grid(row=0,column=2)
        ttk.Label(r2, text="xmax").grid(row=0,column=3); ttk.Entry(r2,textvariable=p['xmax'],width=8).grid(row=0,column=4)
        ttk.Label(r2, text="X scale").grid(row=0,column=5); ttk.Combobox(r2, values=["linear","log"], textvariable=p['xscale'], width=6, state="readonly").grid(row=0,column=6,padx=(2,8))
        ttk.Checkbutton(r2, text="Y auto", variable=p['y_auto']).grid(row=1,column=0,sticky="w")
        ttk.Label(r2, text="ymin").grid(row=1,column=1); ttk.Entry(r2,textvariable=p['ymin'],width=8).grid(row=1,column=2)
        ttk.Label(r2, text="ymax").grid(row=1,column=3); ttk.Entry(r2,textvariable=p['ymax'],width=8).grid(row=1,column=4)
        ttk.Label(r2, text="Y scale").grid(row=1,column=5); ttk.Combobox(r2, values=["linear","log"], textvariable=p['yscale'], width=6, state="readonly").grid(row=1,column=6,padx=(2,8))
        r3 = ttk.Frame(top); r3.pack(fill="x", pady=(4,0))
        ttk.Label(r3, text="DPI").grid(row=0,column=0); ttk.Entry(r3,textvariable=p['dpi'],width=6).grid(row=0,column=1,padx=(2,12))
        ttk.Label(r3, text="Fig W").grid(row=0,column=2); ttk.Entry(r3,textvariable=p['figw'],width=6).grid(row=0,column=3,padx=(2,12))
        ttk.Label(r3, text="Fig H").grid(row=0,column=4); ttk.Entry(r3,textvariable=p['figh'],width=6).grid(row=0,column=5,padx=(2,12))
        r4 = ttk.Frame(top); r4.pack(fill="x", pady=(4,0))
        ttk.Checkbutton(r4, text="Show legend", variable=p['show_legend']).grid(row=0,column=0,sticky="w")
        ttk.Label(r4, text="Legend size").grid(row=0,column=1,sticky="e"); ttk.Entry(r4,textvariable=p['legend_size'],width=6).grid(row=0,column=2,sticky="w",padx=(2,12))
        ttk.Label(r4, text="Legend ncol").grid(row=0,column=3,sticky="e"); ttk.Entry(r4,textvariable=p['legend_ncol'],width=6).grid(row=0,column=4,sticky="w",padx=(2,12))
        ttk.Label(parent, text="Preview", padding=(8,2)).pack(anchor="w")
        ttk.Label(parent, text="Preview", padding=(8,2)).pack(anchor="w")
        p['frame'] = ttk.Frame(parent, relief="groove", padding=4, height=400)
        p['frame'].pack(fill="both", expand=True, padx=8, pady=(0,8))
        p['frame'].pack_propagate(False)   # ★ 關閉自動縮放，保住高度

        return p

    def style_axes(self, ax, panel):
        ts = int(panel['title_size'].get() or 12)
        ls = int(panel['label_size'].get() or 12)
        ax.set_title(panel['title'].get(), fontsize=ts, fontweight='bold')
        ax.set_xlabel(panel['xlabel'].get(), fontsize=ls, fontweight='bold')
        ax.set_ylabel(panel['ylabel'].get(), fontsize=ls, fontweight='bold')
        ax.set_xscale(panel['xscale'].get()); ax.set_yscale(panel['yscale'].get())
        if not panel['x_auto'].get():
            xmin = safe_float(panel['xmin'].get()); xmax = safe_float(panel['xmax'].get())
            if xmin is not None and xmax is not None: ax.set_xlim([xmin, xmax])
        if not panel['y_auto'].get():
            ymin = safe_float(panel['ymin'].get()); ymax = safe_float(panel['ymax'].get())
            if ymin is not None and ymax is not None: ax.set_ylim([ymin, ymax])
        ax.tick_params(which='both', direction='in', length=8, width=2, top=False, right=False, bottom=True, left=True, labelsize=14)
        ax.tick_params(which='minor', length=4, width=1, direction='in', top=False, right=False, bottom=True, left=True)
        ax.xaxis.set_minor_locator(AutoMinorLocator()); ax.yaxis.set_minor_locator(AutoMinorLocator())
        ax.grid(which='major', linestyle='-', linewidth=0.5, color='darkgray')
        for s in ['top','right','bottom','left']:
            ax.spines[s].set_color('black'); ax.spines[s].set_linewidth(2)

    def _legend(self, ax, panel):
        if panel['show_legend'].get():
            fsz = int(panel['legend_size'].get() or 10)
            ncol = int(panel['legend_ncol'].get() or 2)
            leg = ax.legend(ncol=ncol, fancybox=True, loc='best', prop={'weight':'bold','size':fsz})
            leg.get_frame().set_edgecolor('black'); leg.get_frame().set_linewidth(2)

    # ----- sweeps 表 -----
    def _populate_rows(self):
        if self.data_mode == 'single_csv_multi':
            self.src_var.set(str(self.csv_path))
        else:
            self.src_var.set(f"{self.file_list[0]} ... ({len(self.file_list)} files)" if self.file_list else "未選擇")
        for ch in self.rows_container.winfo_children():
            ch.destroy()
        self.rows.clear()
        n = len(self.curves)
        for idx, c in enumerate(self.curves, start=1):
            var_use = tk.BooleanVar(value=True if idx <= 9 else False)
            ttk.Checkbutton(self.rows_container, variable=var_use, command=self.update_all_previews).grid(row=idx, column=0, sticky="w")
            ttk.Label(self.rows_container, text=str(idx), width=self.col_widths[1], anchor="w").grid(row=idx, column=1, sticky="w")
            var_follow = tk.BooleanVar(value=True)
            ttk.Checkbutton(self.rows_container, variable=var_follow, command=self.update_all_previews).grid(row=idx, column=2, sticky="w")
            var_g = tk.StringVar(value=str(min(idx,9)))
            cb = ttk.Combobox(self.rows_container, values=[str(i) for i in range(1,10)], textvariable=var_g, width=self.col_widths[3]-6, state="readonly")
            cb.grid(row=idx, column=3, sticky="w"); cb.bind("<<ComboboxSelected>>", lambda _e: self.update_all_previews())
            var_label = tk.StringVar(value=c["label"])
            ent = ttk.Entry(self.rows_container, textvariable=var_label, width=self.col_widths[4]); ent.grid(row=idx, column=4, sticky="w")
            ent.bind("<FocusOut>", lambda _e: self.update_all_previews())
            var_color = tk.StringVar(value=rainbow_color(idx-1, n))
            cf = ttk.Frame(self.rows_container); cf.grid(row=idx, column=5, sticky="w")
            ttk.Entry(cf, textvariable=var_color, width=10).pack(side="left")
            ttk.Button(cf, text="Pick", command=lambda v=var_color: self.pick_color(v)).pack(side="left", padx=(4,0))
            var_line = tk.BooleanVar(value=True); var_mark = tk.BooleanVar(value=False)
            ttk.Checkbutton(self.rows_container, variable=var_line, command=self.update_all_previews).grid(row=idx, column=6, sticky="w")
            ttk.Checkbutton(self.rows_container, variable=var_mark, command=self.update_all_previews).grid(row=idx, column=7, sticky="w")
            V = c["V"]
            if c.get("I") is not None:
                yr = f"[{np.min(c['I']):.3e}, {np.max(c['I']):.3e}]"
            else:
                yr = f"[{np.min(c['C']):.3e}, {np.max(c['C']):.3e}]"
            ttk.Label(self.rows_container, text=f"[{np.min(V):.3f}, {np.max(V):.3f}]", width=self.col_widths[8], anchor="w").grid(row=idx, column=8, sticky="w")
            ttk.Label(self.rows_container, text=yr, width=self.col_widths[9], anchor="w").grid(row=idx, column=9, sticky="w")
            self.rows.append((var_use, var_follow, var_g, var_label, var_color, var_line, var_mark, idx))

    def on_choose_outdir(self):
        d = filedialog.askdirectory(title="選擇輸出資料夾")
        if d:
            self.outdir = Path(d); self.log(f"Output folder: {self.outdir}")

    def pick_color(self, var_color):
        color = colorchooser.askcolor(title="Choose color", color=var_color.get())[1]
        if color:
            var_color.set(color); self.update_all_previews()

    def log(self, msg):
        self.status.insert("end", msg + "\n"); self.status.see("end")

    # ----- 選取 -----
    def current_selection(self):
        include = []; display = []
        globals_ = [v.get().strip() for v in self.global_vars]
        for i, row in enumerate(self.rows):
            use, follow, gidx, label, color, line_on, mark_on, idx = row
            if not use.get(): continue
            include.append(idx)
            if follow.get():
                try:
                    gname = globals_[int(gidx.get())-1] or f"W{gidx.get()}"
                except Exception:
                    gname = f"W{gidx.get()}"
                disp_label = gname
            else:
                disp_label = label.get().strip() or self.curves[i]["label"]
            clr = color.get().strip() or rainbow_color(i, len(self.rows))
            d = {
                "V": self.curves[i]["V"],
                "I": self.curves[i].get("I", None),
                "C": self.curves[i].get("C", None),
                "label": disp_label, "color": clr, "gidx": gidx.get(),
                "line": line_on.get(), "marker": mark_on.get()
            }
            display.append(d)
        return include, display

    # ----- 計算 -----
    def compute_rv(self, V, I):
        if I is None or len(V) < 3:
            return V, np.full_like(V, np.nan)
        dIdV = np.gradient(I, V, edge_order=2)
        with np.errstate(divide='ignore', invalid='ignore'):
            R = 1.0 / dIdV
        return V, R

    def compute_r0_at_zero(self, V, I, window=0.5):
        if I is None:
            return np.nan
        idx = np.where(np.abs(V) <= window)[0]
        if len(idx) < 3:
            order = np.argsort(np.abs(V)); idx = order[:max(3, min(7, len(V)))]
        v = V[idx]; i = I[idx]
        try:
            a, _b = np.polyfit(v, i, 1)
            if np.isclose(a, 0): return np.nan
            return 1.0 / a
        except Exception:
            return np.nan

    def build_rs_points(self, display):
        window = safe_float(self.r0_window.get(), 0.5)
        xs, ys, labs, clrs = [], [], [], []
        for d in display:
            if d.get("I") is None: continue
            R0 = self.compute_r0_at_zero(d["V"], d["I"], window=window)
            spacing_label = self.global_vars[int(d["gidx"])-1].get().strip() if d.get("gidx") else d["label"]
            spacing = parse_numeric_from_label(spacing_label) or parse_numeric_from_label(d["label"])
            if spacing is None or np.isnan(R0): continue
            xs.append(float(spacing)); ys.append(float(R0)); labs.append(d["label"]); clrs.append(d['color'])
        return xs, ys, labs, clrs

    # Model-1（內部長度用 cm）
    def rho_method1(self, xs: np.ndarray, ys: np.ndarray, R2_um: float):
        d_cm  = np.asarray(xs, float) * UM_TO_CM
        R2_cm = float(R2_um) * UM_TO_CM
        if np.any(d_cm <= 0) or np.any(d_cm >= R2_cm) or d_cm.size < 2:
            return None
        x1 = np.log(R2_cm / (R2_cm - d_cm))                 # 無因次
        x2 = 1.0/(R2_cm - d_cm) + 1.0/R2_cm                 # 1/cm
        X = np.column_stack([x1, x2])
        try:
            beta, *_ = np.linalg.lstsq(X, ys, rcond=None)
        except Exception:
            return None
        A, B = beta  # A: Ω, B: Ω·cm
        Rs = 2*np.pi*A                                   # Ω/□
        rhoc = (2*np.pi*B)**2 / Rs if Rs != 0 else np.nan  # Ω·cm²
        Lt_cm = np.sqrt(rhoc / Rs) if (np.isfinite(rhoc) and Rs > 0) else np.nan
        Lt_um = Lt_cm / UM_TO_CM if np.isfinite(Lt_cm) else np.nan
        return None if not np.isfinite(Lt_um) else (Rs, Lt_um, rhoc)

    # Model-2（correlation 線性化）
    def rho_method2(self, xs: np.ndarray, ys: np.ndarray, R2_um: float):
        d = np.asarray(xs, float); Rt = np.asarray(ys, float)
        if d.size < 2 or np.any(d <= 0) or np.any(d >= R2_um):
            return None
        C = (R2_um/d)*np.log(R2_um/(R2_um - d))
        ycorr = Rt / C
        X = np.column_stack([d, np.ones_like(d)])
        beta, *_ = np.linalg.lstsq(X, ycorr, rcond=None)
        m, c = float(beta[0]), float(beta[1])
        yfit = m*d + c
        ss_res = np.sum((ycorr - yfit)**2); ss_tot = np.sum((ycorr - ycorr.mean())**2)
        r2 = 1 - ss_res/ss_tot if ss_tot > 0 else np.nan
        Rs = m * 2*np.pi*R2_um
        Lt_um = c/(2*m) if m != 0 else np.nan
        if not np.isfinite(Lt_um):
            return None
        rhoc = Rs * (Lt_um*UM_TO_CM)**2
        return Rs, Lt_um, rhoc, (m, c, r2)

    # ----- 繪圖 -----
    def _set_blank(self, frame, text):
        for w in frame.winfo_children(): w.destroy()
        ttk.Label(frame, text=text).pack(anchor="center", expand=True)

    def _embed_figure_keep_ratio(self, panel, fig):
        """把 fig 以等比例縮放嵌入 panel['frame']，避免被擠壓變形。"""
        frm = panel['frame']
        for w in frm.winfo_children():
            w.destroy()
    
        fw = safe_float(panel['figw'].get(), 6.0)
        fh = safe_float(panel['figh'].get(), 4.0)
        aspect = max(fw / fh, 0.0001)
    
        container = ttk.Frame(frm)
        container.pack(fill="both", expand=True)
    
        canvas = FigureCanvasTkAgg(fig, master=container)
        widget = canvas.get_tk_widget()
        widget.place(x=0, y=0, relwidth=0, relheight=0)
    
        def on_resize(e):
            W, H = e.width, e.height
            tw = W
            th = int(tw / aspect)
            if th > H:
                th = H
                tw = int(th * aspect)
            x = (W - tw) // 2
            y = (H - th) // 2
            widget.place(x=x, y=y, width=tw, height=th)
            fig.set_size_inches(tw / fig.dpi, th / fig.dpi)
            canvas.draw_idle()
    
        container.bind("<Configure>", on_resize)
        return canvas




    def update_all_previews(self):
        if not self.curves:
            for p in [self.preview_iv, self.preview_rv, self.preview_rs, self.preview_corr, self.preview_cv]:
                self._set_blank(p['frame'], "請載入資料")
            return
        include, display = self.current_selection()
        self._draw_iv(self.preview_iv, display)
        self._draw_rv(self.preview_rv, display)
        self._draw_rs(self.preview_rs, display)
        self._draw_corr(self.preview_corr, display)
        self._draw_cv(self.preview_cv, display)

    def _draw_iv(self, panel, display):
        for w in panel['frame'].winfo_children(): w.destroy()
        items = [d for d in display if d.get("I") is not None]
        if not items:
            self._set_blank(panel['frame'], "No I–V curves"); return
        figw = safe_float(panel['figw'].get(), 6); figh = safe_float(panel['figh'].get(), 4)
        fig, ax = plt.subplots(figsize=(figw, figh), dpi=100)
        for d in items:
            lw = 1.5 if d["line"] else 0; ms = 4 if d["marker"] else 0
            ax.plot(d["V"], d["I"], label=d["label"], color=d["color"],
                    linewidth=lw, marker='o' if d["marker"] else None, markersize=ms,
                    markeredgewidth=1.0 if d["marker"] else 0, markeredgecolor='black' if d["marker"] else None,
                    markerfacecolor=d["color"] if d["marker"] else None)
        self.style_axes(ax, panel); self._legend(ax, panel)
        fig.tight_layout()
        self._embed_figure_keep_ratio(panel, fig)


    def _draw_rv(self, panel, display):
        for w in panel['frame'].winfo_children(): w.destroy()
        items = [d for d in display if d.get("I") is not None]
        if not items:
            self._set_blank(panel['frame'], "No I–V curves"); return
        figw = safe_float(panel['figw'].get(), 6); figh = safe_float(panel['figh'].get(), 4)
        fig, ax = plt.subplots(figsize=(figw, figh), dpi=100)
        for d in items:
            V, R = self.compute_rv(d["V"], d["I"])
            lw = 1.2 if d["line"] else 0; ms = 3 if d["marker"] else 0
            ax.plot(V, R, label=d["label"], color=d["color"],
                    linewidth=lw, marker='o' if d["marker"] else None, markersize=ms)
        self.style_axes(ax, panel); self._legend(ax, panel)
        fig.tight_layout()
        self._embed_figure_keep_ratio(panel, fig)

    def _draw_rs(self, panel, display):
        for w in panel['frame'].winfo_children(): w.destroy()
        items = [d for d in display if d.get("I") is not None]
        xs, ys, labs, clrs = self.build_rs_points(items)
        if not xs:
            self._set_blank(panel['frame'], "No R0 points"); self._update_rt_list([]); self.result_text.delete("1.0","end"); return
        figw = safe_float(panel['figw'].get(), 6); figh = safe_float(panel['figh'].get(), 4)
        fig, ax = plt.subplots(figsize=(figw, figh), dpi=100)
        for x, y, lab, clr in zip(xs, ys, labs, clrs):
            ax.scatter(x, y, label=lab, color=clr, edgecolors='black')
        # OLS: y = a x + b
        xarr = np.asarray(xs, float); yarr = np.asarray(ys, float)
        mask = np.isfinite(xarr) & np.isfinite(yarr); xarr = xarr[mask]; yarr = yarr[mask]
        fit_summary = "insufficient points"
        if xarr.size >= 2:
            X = np.column_stack([xarr, np.ones_like(xarr)])
            beta, *_ = np.linalg.lstsq(X, yarr, rcond=None)
            a, b = float(beta[0]), float(beta[1])
            xfit = np.linspace(xarr.min(), xarr.max(), 200); yfit = a*xfit + b
            ax.plot(xfit, yfit, color="black", linewidth=1.2, linestyle="--", label="fit")
            ss_res = np.sum((yarr - (a*xarr + b))**2); ss_tot = np.sum((yarr - yarr.mean())**2)
            r2 = 1 - ss_res/ss_tot if ss_tot > 0 else np.nan
            fit_summary = f"a={a:.6g} (Ω/μm), b={b:.6g} (Ω), R²={r2:.4f}"
            
        self.style_axes(ax, panel)
        self._legend(ax, panel)
        fig.tight_layout()
        self._embed_figure_keep_ratio(panel, fig)

    
        self._update_rt_list(list(zip(xs, ys, labs)))
        R2_um = safe_float(self.r2_var.get())
        lines = [f"[R0 vs Spacing] {fit_summary}"]
        if R2_um is None:
            lines.append("R2 invalid, ρc unavailable.")
        else:
            xarr = np.asarray(xs, float); yarr = np.asarray(ys, float)
            m1 = self.rho_method1(xarr, yarr, R2_um)
            if m1 is None: lines.append("Method-1: fail (check 0<d<R2 & points)")
            else:
                Rs1, Lt1_um, rhoc1 = m1
                lines.append(f"Method-1: Rs={Rs1:.6g} Ω/□, Lt={Lt1_um:.6g} μm, ρc={rhoc1:.6g} Ω·cm²")
            m2 = self.rho_method2(xarr, yarr, R2_um)
            if m2 is None: lines.append("Method-2: fail (check 0<d<R2 & points)")
            else:
                Rs2, Lt2_um, rhoc2, (m, c, r2c) = m2
                lines.append(f"Method-2: Rs={Rs2:.6g} Ω/□, Lt={Lt2_um:.6g} μm, ρc={rhoc2:.6g} Ω·cm²; y={m:.3g}x+{c:.3g}, R²={r2c:.4f}")
        self.result_text.delete("1.0","end"); self.result_text.insert("end", "\n".join(lines) + "\n")

    def _draw_corr(self, panel, display):
        for w in panel['frame'].winfo_children(): w.destroy()
        xs, ys, _, _ = self.build_rs_points([d for d in display if d.get("I") is not None])
        if len(xs) < 2:
            self._set_blank(panel['frame'], "Need at least two R–Spacing points"); self.corr_text.set(""); return
        R2_um = safe_float(self.r2_var.get())
        if R2_um is None:
            self._set_blank(panel['frame'], "Please input R2"); self.corr_text.set(""); return
        d = np.asarray(xs, float); Rt = np.asarray(ys, float)
        C = (R2_um/d)*np.log(R2_um/(R2_um - d)); Rt_corr = Rt / C
        X = np.column_stack([d, np.ones_like(d)])
        beta, *_ = np.linalg.lstsq(X, Rt_corr, rcond=None)
        m, c = float(beta[0]), float(beta[1])
        yfit = m*d + c
        ss_res = np.sum((Rt_corr - yfit)**2); ss_tot = np.sum((Rt_corr - Rt_corr.mean())**2)
        r2 = 1 - ss_res/ss_tot if ss_tot > 0 else np.nan
        figw = safe_float(panel['figw'].get(), 6); figh = safe_float(panel['figh'].get(), 4)
        fig, ax = plt.subplots(figsize=(figw, figh), dpi=100)
        ax.scatter(d, Rt, label="Original Rt(d)", marker='s')
        ax.scatter(d, Rt_corr, label="Corrected Rt/C(d)", marker='o')
        xfit = np.linspace(np.min(d), np.max(d), 200)
        ax.plot(xfit, m*xfit + c, linestyle='--', color='black', label="Linear fit on corrected")
        self.style_axes(ax, panel); self._legend(ax, panel)
        fig.tight_layout()
        self._embed_figure_keep_ratio(panel, fig)
        Rs = m * 2*np.pi*R2_um
        Lt_um = c/(2*m) if m != 0 else np.nan
        rhoc = Rs * (Lt_um*UM_TO_CM)**2 if np.isfinite(Lt_um) else np.nan
        self.corr_text.set(f"Linearized: m={m:.6g}, c={c:.6g}, R²={r2:.4f} | Rs={Rs:.6g} Ω/□, Lt={Lt_um:.6g} μm, ρc={rhoc:.6g} Ω·cm²")

    def _draw_cv(self, panel, display):
        for w in panel['frame'].winfo_children(): w.destroy()
        items = [d for d in display if d.get("C") is not None]
        if not items:
            self._set_blank(panel['frame'], "No C–V curves"); return
        figw = safe_float(panel['figw'].get(), 6); figh = safe_float(panel['figh'].get(), 4)
        fig, ax = plt.subplots(figsize=(figw, figh), dpi=100)
        for d in items:
            lw = 1.2 if d["line"] else 0; ms = 3 if d["marker"] else 0
            ax.plot(d["V"], d["C"], label=d["label"], color=d["color"],
                    linewidth=lw, marker='o' if d["marker"] else None, markersize=ms)
        self.style_axes(ax, panel); self._legend(ax, panel)
        fig.tight_layout()
        self._embed_figure_keep_ratio(panel, fig)

    # ----- Rt 清單 -----
    def _update_rt_list(self, items: List[Tuple[float, float, str]]):
        self.rt_text.delete("1.0", "end")
        if not items:
            self.rt_text.insert("end", "No data\n"); return
        self.rt_text.insert("end", f"{'Label/Spacing':<24}\tSpacing (μm)\tR0 (Ω)\n")
        self.rt_text.insert("end", "-"*60 + "\n")
        for spacing, R0, lab in items:
            self.rt_text.insert("end", f"{lab:<24}\t{spacing:.6g}\t{R0:.6g}\n")

    # ----- 批次與其他 -----
    def select_all(self):
        for (use, *_rest) in self.rows:
            use.set(True)
        self.update_all_previews()

    def clear_all(self):
        for (use, *_rest) in self.rows:
            use.set(False)
        self.update_all_previews()

    def _save_correlation_plot(self, d: np.ndarray, Rt: np.ndarray, R2_um: float, outfile: Path, panel_ref):
        if len(d) < 2 or np.any(d <= 0) or np.any(d >= R2_um):
            return None
        C = (R2_um/d)*np.log(R2_um/(R2_um - d)); Rt_corr = Rt / C
        X = np.column_stack([d, np.ones_like(d)]); beta, *_ = np.linalg.lstsq(X, Rt_corr, rcond=None)
        m, c = float(beta[0]), float(beta[1]); yfit = m*d + c
        ss_res = np.sum((Rt_corr - yfit)**2); ss_tot = np.sum((Rt_corr - Rt_corr.mean())**2)
        r2 = 1 - ss_res/ss_tot if ss_tot > 0 else np.nan
        Rs = m * 2*np.pi*R2_um; Lt_um = c/(2*m) if m != 0 else np.nan
        if not np.isfinite(Lt_um): return None
        rhoc = Rs * (Lt_um*UM_TO_CM)**2
        dpi = int(panel_ref['dpi'].get() or 300)
        figw = safe_float(panel_ref['figw'].get(), 6); figh = safe_float(panel_ref['figh'].get(), 4)
        fig, ax = plt.subplots(figsize=(figw, figh), dpi=dpi)
        ax.scatter(d, Rt, label="Original Rt(d)", marker='s')
        ax.scatter(d, Rt_corr, label="Corrected Rt/C(d)", marker='o')
        xfit = np.linspace(np.min(d), np.max(d), 200)
        ax.plot(xfit, m*xfit + c, linestyle='--', color='black', label="Linear fit on corrected")
        ax.set_title("Original vs Corrected Rt(d)", fontsize=int(panel_ref['title_size'].get() or 12), fontweight='bold')
        ax.set_xlabel("Spacing d (μm)", fontsize=int(panel_ref['label_size'].get() or 12), fontweight='bold')
        ax.set_ylabel("Resistance (Ω)", fontsize=int(panel_ref['label_size'].get() or 12), fontweight='bold')
        ax.tick_params(which='both', direction='in', length=8, width=2, top=False, right=False, bottom=True, left=True)
        ax.xaxis.set_minor_locator(AutoMinorLocator()); ax.yaxis.set_minor_locator(AutoMinorLocator())
        for s in ['top','right','bottom','left']:
            ax.spines[s].set_color('black'); ax.spines[s].set_linewidth(2)
        ax.grid(which='major', linestyle='-', linewidth=0.5, color='darkgray'); ax.legend()
        fig.tight_layout(); fig.savefig(outfile, dpi=dpi); plt.close(fig)
        return Rs, Lt_um, rhoc, (m, c, r2)

    def _save_panel_fig(self, panel, display, outfile: Path, kind: str):
        dpi = int(panel['dpi'].get() or 300)
        figw = safe_float(panel['figw'].get(), 6); figh = safe_float(panel['figh'].get(), 4)
        fig, ax = plt.subplots(figsize=(figw, figh), dpi=dpi)
        if kind == "iv":
            for d in display:
                if d.get("I") is None: continue
                lw = 1.5 if d["line"] else 0; ms = 4 if d["marker"] else 0
                ax.plot(d["V"], d["I"], label=d["label"], color=d["color"],
                        linewidth=lw, marker='o' if d["marker"] else None, markersize=ms)
        elif kind == "rv":
            for d in display:
                if d.get("I") is None: continue
                V, R = self.compute_rv(d["V"], d["I"])
                lw = 1.2 if d["line"] else 0; ms = 3 if d["marker"] else 0
                ax.plot(V, R, label=d["label"], color=d["color"],
                        linewidth=lw, marker='o' if d["marker"] else None, markersize=ms)
        elif kind == "cv":
            for d in display:
                if d.get("C") is None: continue
                lw = 1.2 if d["line"] else 0; ms = 3 if d["marker"] else 0
                ax.plot(d["V"], d["C"], label=d["label"], color=d["color"],
                        linewidth=lw, marker='o' if d["marker"] else None, markersize=ms)
        self.style_axes(ax, panel); self._legend(ax, panel)
        fig.tight_layout(); fig.savefig(outfile, dpi=dpi); plt.close(fig)

    
    def export_all(self):
        if not self.curves:
            messagebox.showwarning("提醒", "尚未載入資料。"); return
        include, display = self.current_selection()
        if not include:
            messagebox.showwarning("提醒", "請至少勾選一個 sweep。"); return
        outdir = self.outdir or (self.csv_path.parent if self.csv_path else Path.cwd()) / "trinity_capres_out"
        outdir.mkdir(parents=True, exist_ok=True)
        self._save_panel_fig(self.preview_iv, display, outdir / "IV_overlay_selected.png", "iv")
        self._save_panel_fig(self.preview_rv, display, outdir / "RV_overlay_selected.png", "rv")
        self._save_panel_fig(self.preview_cv, display, outdir / "CV_overlay_selected.png", "cv")
        xs, ys, _, _ = self.build_rs_points([d for d in display if d.get("I") is not None])
        R2_um = safe_float(self.r2_var.get())
        summary = []
        if xs and ys and R2_um is not None:
            xarr = np.asarray(xs, float); yarr = np.asarray(ys, float)
            m1 = self.rho_method1(xarr, yarr, R2_um)
            if m1:
                Rs, Lt_um, rhoc = m1
                summary.append(f"Method1: Rs={Rs:.9g} Ω/□, Lt={Lt_um:.9g} μm, rho_c={rhoc:.9g} Ω·cm²")
            else:
                summary.append("Method1: fail (check 0<d<R2 & points)")
            res2 = self._save_correlation_plot(np.asarray(xs, float), np.asarray(ys, float), R2_um, outdir / "R_spacing_correlation.png", self.preview_corr)
            if res2:
                Rs2, Lt2, rhoc2, (m, c, r2) = res2
                summary.append(f"Method2: Rs={Rs2:.9g} Ω/□, Lt={Lt2:.9g} μm, rho_c={rhoc2:.9g} Ω·cm²; y={m:.6g}x+{c:.6g}, R^2={r2:.6g}")
            else:
                summary.append("Method2: fail (check 0<d<R2 & points)")
        else:
            summary.append("Missing R0 or R2.")
        (outdir / "rho_summary.txt").write_text("\n".join(summary), encoding="utf-8")
        messagebox.showinfo("完成", f"輸出完成：\n{outdir}\n\n" + "\n".join(summary))

if __name__ == "__main__":
    # Tweak base font size here (global UI scaling)
    BASE_FONT = ("Segoe UI", 14)  # change 13 → 14/16/18 if you want bigger UI
    App(base_font=BASE_FONT).mainloop()

