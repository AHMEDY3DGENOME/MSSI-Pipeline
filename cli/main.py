import customtkinter as ctk
import threading
from tkinter import filedialog
from pathlib import Path
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from MSSI_Pipeline.loader     import MSSILoader
from MSSI_Pipeline.calculator import MSSICalculator
from MSSI_Pipeline.classifier import MSSIClassifier
from MSSI_Pipeline.visualizer import MSSIVisualizer
from MSSI_Pipeline.reporter   import MSSIReporter
from MSSI_Pipeline.risk_index import ResistanceRiskIndex, RISK_COLORS
from MSSI_Pipeline.mapper     import EgyptMapper, CACHE_FILE

ctk.set_appearance_mode("dark")
ctk.set_default_color_theme("blue")

STRESS_COLORS = {
    "very_high": "#E24B4A",
    "high":      "#EF9F27",
    "moderate":  "#378ADD",
    "low":       "#1D9E75",
}


class MSSIApp(ctk.CTk):

    def __init__(self):
        super().__init__()

        self.title("MSSI Pipeline v1.0.0")
        self.geometry("1100x800")
        self.minsize(1000, 700)

        self.filepath = None
        self.output_dir = None

        self.results = None
        self.risks = None
        self._mapper = None

        self.genes = []
        self.weight_vars = {}

        self._build_ui()

    def _build_ui(self):
        self.grid_columnconfigure(0, weight=1)
        self.grid_rowconfigure(1, weight=1)
        self._build_header()
        self._build_main()
        self._build_footer()

    def _build_header(self):
        header = ctk.CTkFrame(self, fg_color=("#5346B7", "#3C3489"), corner_radius=0)
        header.grid(row=0, column=0, sticky="ew")
        header.grid_columnconfigure(0, weight=1)
        ctk.CTkLabel(
            header, text="MSSI Pipeline",
            font=ctk.CTkFont(size=22, weight="bold"),
            text_color="white"
        ).grid(row=0, column=0, padx=20, pady=(12, 2), sticky="w")
        ctk.CTkLabel(
            header,
            text="Multi-HSP Stress Signature Index - Spodoptera frugiperda",
            font=ctk.CTkFont(size=12),
            text_color="#CECBF6"
        ).grid(row=1, column=0, padx=20, pady=(0, 12), sticky="w")

    def _build_main(self):
        main = ctk.CTkFrame(self, fg_color="transparent")
        main.grid(row=1, column=0, sticky="nsew", padx=16, pady=16)
        main.grid_columnconfigure(0, weight=1)
        main.grid_columnconfigure(1, weight=2)
        main.grid_rowconfigure(0, weight=1)
        self._build_left_panel(main)
        self._build_right_panel(main)

    def _build_left_panel(self, parent):
        left = ctk.CTkFrame(parent)
        left.grid(row=0, column=0, sticky="nsew", padx=(0, 8))
        left.grid_columnconfigure(0, weight=1)

        ctk.CTkLabel(left, text="Input File",
                     font=ctk.CTkFont(weight="bold")).grid(
            row=0, column=0, padx=16, pady=(16, 4), sticky="w")
        ctk.CTkButton(left, text="Browse Excel File",
                      command=self._browse_file).grid(
            row=1, column=0, padx=16, pady=4, sticky="ew")
        self.file_label = ctk.CTkLabel(
            left, text="No file selected",
            text_color="gray", wraplength=220)
        self.file_label.grid(row=2, column=0, padx=16, pady=4)

        ctk.CTkLabel(left, text="Output Directory",
                     font=ctk.CTkFont(weight="bold")).grid(
            row=3, column=0, padx=16, pady=(12, 4), sticky="w")
        ctk.CTkButton(left, text="Browse Output Folder",
                      command=self._browse_output).grid(
            row=4, column=0, padx=16, pady=4, sticky="ew")
        self.output_label = ctk.CTkLabel(
            left, text="outputs/", text_color="gray")
        self.output_label.grid(row=5, column=0, padx=16, pady=4)

        ctk.CTkLabel(left, text="Gene Weights",
                     font=ctk.CTkFont(weight="bold")).grid(
            row=6, column=0, padx=16, pady=(12, 4), sticky="w")
        self.weights_frame = ctk.CTkFrame(left, fg_color="transparent")
        self.weights_frame.grid(row=7, column=0, padx=16, pady=4, sticky="ew")
        self.weights_frame.grid_columnconfigure(0, weight=1)

        self.run_btn = ctk.CTkButton(
            left, text="Run Analysis",
            font=ctk.CTkFont(size=14, weight="bold"),
            fg_color="#1D9E75", hover_color="#0F6E56",
            height=40, command=self._run_analysis
        )
        self.run_btn.grid(row=8, column=0, padx=16, pady=16, sticky="ew")
        self.progress = ctk.CTkProgressBar(left)
        self.progress.grid(row=9, column=0, padx=16, pady=(0, 8), sticky="ew")
        self.progress.set(0)

    def _build_right_panel(self, parent):
        right = ctk.CTkFrame(parent)
        right.grid(row=0, column=1, sticky="nsew", padx=(8, 0))
        right.grid_columnconfigure(0, weight=1)
        right.grid_rowconfigure(1, weight=1)

        self.tab_view = ctk.CTkTabview(right)
        self.tab_view.grid(row=0, column=0, rowspan=2,
                           sticky="nsew", padx=16, pady=(8, 0))

        self.tab_view.add("MSSI Results")
        self.tab_view.add("Risk Index")
        self.tab_view.add("Egypt Map")

        self.results_frame = ctk.CTkScrollableFrame(
            self.tab_view.tab("MSSI Results"))
        self.results_frame.pack(fill="both", expand=True)
        self.results_frame.grid_columnconfigure(0, weight=1)

        self.risk_frame = ctk.CTkScrollableFrame(
            self.tab_view.tab("Risk Index"))
        self.risk_frame.pack(fill="both", expand=True)
        self.risk_frame.grid_columnconfigure(0, weight=1)

        self.map_frame = ctk.CTkFrame(
            self.tab_view.tab("Egypt Map"),
            fg_color="transparent")
        self.map_frame.pack(fill="both", expand=True)

        ctk.CTkLabel(right, text="Log",
                     font=ctk.CTkFont(weight="bold")).grid(
            row=2, column=0, padx=16, pady=(4, 2), sticky="w")
        self.log_box = ctk.CTkTextbox(right, height=120, state="disabled")
        self.log_box.grid(row=3, column=0, sticky="ew",
                          padx=16, pady=(0, 16))
        self._log("Ready. Please select a file to begin.")

    def _build_footer(self):
        footer = ctk.CTkFrame(self, fg_color="transparent", height=28)
        footer.grid(row=2, column=0, sticky="ew", padx=16)
        ctk.CTkLabel(
            footer,
            text="MSSI Pipeline v1.0.0  |  Ahmed Yassin  |  MIT License",
            text_color="gray", font=ctk.CTkFont(size=11)
        ).pack(side="right")

    def _log(self, msg: str):
        self.log_box.configure(state="normal")
        self.log_box.insert("end", msg + "\n")
        self.log_box.see("end")
        self.log_box.configure(state="disabled")

    def _browse_file(self):
        path = filedialog.askopenfilename(
            title="Select qPCR Input File",
            filetypes=[
                ("Supported Files", "*.csv *.xlsx *.xls"),
                ("CSV Files", "*.csv"),
                ("Excel Files", "*.xlsx *.xls"),
                ("All Files", "*.*"),
            ],
        )

        if not path:
            return

        # Store selected input file
        self.filepath = path

        # Create MSSI_Results beside the selected input file
        input_path = Path(path).expanduser().resolve()
        result_dir = input_path.parent / "MSSI_Results"

        result_dir.mkdir(
            parents=True,
            exist_ok=True,
        )

        # All pipeline outputs will now use this directory
        self.output_dir = str(result_dir)

        # Update GUI
        self.file_label.configure(
            text=input_path.name,
            text_color="white",
        )

        self.output_label.configure(
            text=str(result_dir),
            text_color="white",
        )

        # Log
        self._log(f"File loaded: {input_path}")
        self._log(f"Results will be saved to: {result_dir}")

        # Read genes from selected input file
        self._load_genes()
    def _browse_output(self):
        path = filedialog.askdirectory()
        if path:
            self.output_dir = path
            self.output_label.configure(text=path, text_color="white")
            self._log(f"Output directory: {path}")

    def _load_genes(self):
        try:
            loader     = MSSILoader(self.filepath)
            self.genes = loader.get_genes()
            for w in self.weights_frame.winfo_children():
                w.destroy()
            self.weight_vars = {}
            for i, gene in enumerate(self.genes):
                ctk.CTkLabel(self.weights_frame, text=gene).grid(
                    row=i*2, column=0, sticky="w", pady=(4, 0))
                var = ctk.DoubleVar(value=1.0)
                self.weight_vars[gene] = var
                ctk.CTkSlider(
                    self.weights_frame, from_=0.1, to=3.0,
                    variable=var, width=200
                ).grid(row=i*2+1, column=0, sticky="ew", pady=(0, 4))
            self._log(f"Genes found: {', '.join(self.genes)}")
        except Exception as e:
            self._log(f"Error: {e}")

    def _run_analysis(self):
        if not self.filepath:
            self._log("Please select a file first.")
            return
        threading.Thread(target=self._analysis_task, daemon=True).start()

    def _analysis_task(self):
        try:
            self.run_btn.configure(state="disabled")

            self.progress.set(0.1)
            self._log("[1/7] Loading data...")
            loader  = MSSILoader(self.filepath)
            data    = loader.load()
            weights = {g: self.weight_vars[g].get() for g in self.genes}

            self.progress.set(0.25)
            self._log("[2/7] Computing MSSI scores...")
            calculator = MSSICalculator(data, weights=weights)
            stats      = calculator.get_full_stats()
            normalized = calculator.normalize_scores()

            self.progress.set(0.40)
            self._log("[3/7] Classifying stress levels...")
            classifier   = MSSIClassifier(stats, normalized)
            self.results = classifier.classify_all()

            self.progress.set(0.55)
            self._log("[4/7] Computing Resistance Risk Index...")
            risk_engine  = ResistanceRiskIndex(self.results)
            self.risks   = risk_engine.compute_all()

            self.progress.set(0.65)
            self._log("[5/7] Generating plots...")
            visualizer = MSSIVisualizer(self.results, output_dir=self.output_dir)
            visualizer.plot_all()

            self.progress.set(0.78)
            first_time = not CACHE_FILE.exists()
            if first_time:
                self._log("[6/7] Downloading Egypt map (first time only)...")
            else:
                self._log("[6/7] Generating Egypt Risk Map...")
            self._mapper = EgyptMapper(self.risks, output_dir=self.output_dir)
            self._mapper.plot_risk_bar()
            self._log("[6/7] Map ready!")

            self.progress.set(0.90)
            self._log("[7/7] Exporting report...")
            reporter = MSSIReporter(self.results, output_dir=self.output_dir)
            reporter.export_all()

            self.progress.set(1.0)
            self._log(f"Done. Results saved to: {self.output_dir}/")
            self.after(0, self._show_results)
            self.after(0, self._show_risks)
            self.after(0, self._show_map)

        except Exception as e:
            self._log(f"Error: {e}")
            self.progress.set(0)
        finally:
            self.run_btn.configure(state="normal")

    def _show_map(self):
        for w in self.map_frame.winfo_children():
            w.destroy()
        if not self.risks:
            return

        loading = ctk.CTkLabel(
            self.map_frame,
            text="Rendering map...",
            font=ctk.CTkFont(size=14),
            text_color="gray"
        )
        loading.pack(expand=True)
        self.update()

        try:
            fig = self._mapper.plot_map(save=True)
            for w in self.map_frame.winfo_children():
                w.destroy()
            canvas = FigureCanvasTkAgg(fig, master=self.map_frame)
            canvas.draw()
            canvas.get_tk_widget().pack(fill="both", expand=True)
        except Exception as e:
            for w in self.map_frame.winfo_children():
                w.destroy()
            ctk.CTkLabel(
                self.map_frame,
                text=f"Map error: {e}",
                text_color="#E24B4A"
            ).pack(expand=True)

    def _show_results(self):
        for w in self.results_frame.winfo_children():
            w.destroy()
        if not self.results:
            return
        headers = ["Location", "MSSI Score", "Stress Level", "Confidence"]
        for col, h in enumerate(headers):
            ctk.CTkLabel(
                self.results_frame, text=h,
                font=ctk.CTkFont(weight="bold"),
                fg_color=("#3C3489", "#3C3489"),
                corner_radius=4, width=150
            ).grid(row=0, column=col, padx=4, pady=4, sticky="ew")
        for i, r in enumerate(self.results, start=1):
            color    = STRESS_COLORS[r.stress_level]
            row_data = [r.location, f"{r.mssi_score:.4f}",
                        r.stress_label, f"{int(r.confidence*100)}%"]
            for col, val in enumerate(row_data):
                tc = color if col == 2 else None
                ctk.CTkLabel(
                    self.results_frame, text=val,
                    text_color=tc, width=150
                ).grid(row=i, column=col, padx=4, pady=2, sticky="ew")
        high = max(self.results, key=lambda r: r.mssi_score)
        low  = min(self.results, key=lambda r: r.mssi_score)
        n    = len(self.results) + 2
        ctk.CTkLabel(
            self.results_frame,
            text=f"Highest stress: {high.location} ({high.mssi_score:.4f})",
            text_color=STRESS_COLORS["very_high"],
            font=ctk.CTkFont(weight="bold")
        ).grid(row=n, column=0, columnspan=4, pady=(12, 2), sticky="w")
        ctk.CTkLabel(
            self.results_frame,
            text=f"Lowest stress:  {low.location} ({low.mssi_score:.4f})",
            text_color=STRESS_COLORS["low"],
            font=ctk.CTkFont(weight="bold")
        ).grid(row=n+1, column=0, columnspan=4, pady=(2, 12), sticky="w")

    def _show_risks(self):
        for w in self.risk_frame.winfo_children():
            w.destroy()
        if not self.risks:
            return
        headers = ["Location", "Risk Index", "Risk Level", "Action"]
        for col, h in enumerate(headers):
            ctk.CTkLabel(
                self.risk_frame, text=h,
                font=ctk.CTkFont(weight="bold"),
                fg_color=("#3C3489", "#3C3489"),
                corner_radius=4,
                width=80 if col < 3 else 300
            ).grid(row=0, column=col, padx=4, pady=4, sticky="ew")
        for i, r in enumerate(self.risks, start=1):
            row_data = [
                r.location, f"{r.risk_index:.4f}",
                r.risk_level.replace("_", " ").title(),
                r.risk_action,
            ]
            for col, val in enumerate(row_data):
                tc = r.risk_color if col == 2 else None
                ctk.CTkLabel(
                    self.risk_frame, text=val,
                    text_color=tc,
                    width=80 if col < 3 else 300,
                    wraplength=280 if col == 3 else 0
                ).grid(row=i, column=col, padx=4, pady=2, sticky="ew")
        high = max(self.risks, key=lambda r: r.risk_index)
        n    = len(self.risks) + 2
        ctk.CTkLabel(
            self.risk_frame,
            text=f"Highest risk: {high.location} - {high.risk_label}",
            text_color=RISK_COLORS["critical"],
            font=ctk.CTkFont(weight="bold")
        ).grid(row=n, column=0, columnspan=4, pady=(12, 4), sticky="w")
        ctk.CTkLabel(
            self.risk_frame,
            text=f"Action: {high.risk_action}",
            text_color="#EF9F27",
            wraplength=500
        ).grid(row=n+1, column=0, columnspan=4, pady=(2, 12), sticky="w")


def main():
    app = MSSIApp()
    app.mainloop()


if __name__ == "__main__":
    main()