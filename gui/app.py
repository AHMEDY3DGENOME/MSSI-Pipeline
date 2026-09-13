import customtkinter as ctk
import matplotlib

matplotlib.use("Agg")

import threading
import json

from pathlib import Path
from tkinter import filedialog, messagebox
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg

from MSSI_Pipeline.loader import MSSILoader
from MSSI_Pipeline.calculator import MSSICalculator
from MSSI_Pipeline.classifier import MSSIClassifier
from MSSI_Pipeline.visualizer import MSSIVisualizer
from MSSI_Pipeline.reporter import MSSIReporter
from MSSI_Pipeline.risk_index import ResistanceRiskIndex
from MSSI_Pipeline.mapper import EgyptMapper


# ============================================================
# APPLICATION STYLE
# ============================================================

ctk.set_appearance_mode("dark")
ctk.set_default_color_theme("blue")


STRESS_COLORS = {
    "very_high": "#E24B4A",
    "high": "#EF9F27",
    "moderate": "#378ADD",
    "low": "#1D9E75",
}


# ============================================================
# MAIN APPLICATION
# ============================================================

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
        self._visualizer = None

        self.genes = []
        self.weight_vars = {}

        self.corrected_qpcr_mode = False

        self._build_ui()

    # ========================================================
    # MAIN UI
    # ========================================================

    def _build_ui(self):

        self.grid_columnconfigure(0, weight=1)
        self.grid_rowconfigure(1, weight=1)

        self._build_header()
        self._build_main()
        self._build_footer()

    # ========================================================
    # HEADER
    # ========================================================

    def _build_header(self):

        header = ctk.CTkFrame(
            self,
            fg_color=("#5346B7", "#3C3489"),
            corner_radius=0
        )

        header.grid(
            row=0,
            column=0,
            sticky="ew"
        )

        header.grid_columnconfigure(
            0,
            weight=1
        )

        ctk.CTkLabel(
            header,
            text="MSSI Pipeline",
            font=ctk.CTkFont(
                size=22,
                weight="bold"
            ),
            text_color="white"
        ).grid(
            row=0,
            column=0,
            padx=20,
            pady=(12, 2),
            sticky="w"
        )

        ctk.CTkLabel(
            header,
            text=(
                "Multi-HSP Stress Signature Index - "
                "Spodoptera frugiperda"
            ),
            font=ctk.CTkFont(size=12),
            text_color="#CECBF6"
        ).grid(
            row=1,
            column=0,
            padx=20,
            pady=(0, 12),
            sticky="w"
        )

    # ========================================================
    # MAIN BODY
    # ========================================================

    def _build_main(self):

        main = ctk.CTkFrame(
            self,
            fg_color="transparent"
        )

        main.grid(
            row=1,
            column=0,
            sticky="nsew",
            padx=16,
            pady=16
        )

        main.grid_columnconfigure(
            0,
            weight=1
        )

        main.grid_columnconfigure(
            1,
            weight=2
        )

        main.grid_rowconfigure(
            0,
            weight=1
        )

        self._build_left_panel(main)
        self._build_right_panel(main)

    # ========================================================
    # LEFT PANEL
    # ========================================================

    def _build_left_panel(self, parent):

        left = ctk.CTkFrame(parent)

        left.grid(
            row=0,
            column=0,
            sticky="nsew",
            padx=(0, 8)
        )

        left.grid_columnconfigure(
            0,
            weight=1
        )

        # Input file
        ctk.CTkLabel(
            left,
            text="Input File",
            font=ctk.CTkFont(
                weight="bold"
            )
        ).grid(
            row=0,
            column=0,
            padx=16,
            pady=(16, 4),
            sticky="w"
        )

        ctk.CTkButton(
            left,
            text="Browse qPCR / Excel File",
            command=self._browse_file
        ).grid(
            row=1,
            column=0,
            padx=16,
            pady=4,
            sticky="ew"
        )

        self.file_label = ctk.CTkLabel(
            left,
            text="No file selected",
            text_color="gray",
            wraplength=220
        )

        self.file_label.grid(
            row=2,
            column=0,
            padx=16,
            pady=4
        )

        # Output directory
        ctk.CTkLabel(
            left,
            text="Output Directory",
            font=ctk.CTkFont(
                weight="bold"
            )
        ).grid(
            row=3,
            column=0,
            padx=16,
            pady=(12, 4),
            sticky="w"
        )

        ctk.CTkButton(
            left,
            text="Browse Output Folder",
            command=self._browse_output
        ).grid(
            row=4,
            column=0,
            padx=16,
            pady=4,
            sticky="ew"
        )

        self.output_label = ctk.CTkLabel(
            left,
            text="Select an input file",
            text_color="gray"
        )

        self.output_label.grid(
            row=5,
            column=0,
            padx=16,
            pady=4
        )

        # Gene weights
        ctk.CTkLabel(
            left,
            text="Gene Weights",
            font=ctk.CTkFont(
                weight="bold"
            )
        ).grid(
            row=6,
            column=0,
            padx=16,
            pady=(12, 4),
            sticky="w"
        )

        self.weights_frame = ctk.CTkFrame(
            left,
            fg_color="transparent"
        )

        self.weights_frame.grid(
            row=7,
            column=0,
            padx=16,
            pady=4,
            sticky="ew"
        )

        self.weights_frame.grid_columnconfigure(
            0,
            weight=1
        )

        # Run button
        self.run_btn = ctk.CTkButton(
            left,
            text="Run Analysis",
            font=ctk.CTkFont(
                size=14,
                weight="bold"
            ),
            fg_color="#1D9E75",
            hover_color="#0F6E56",
            height=40,
            command=self._run_analysis
        )

        self.run_btn.grid(
            row=8,
            column=0,
            padx=16,
            pady=16,
            sticky="ew"
        )

        self.progress = ctk.CTkProgressBar(
            left
        )

        self.progress.grid(
            row=9,
            column=0,
            padx=16,
            pady=(0, 8),
            sticky="ew"
        )

        self.progress.set(0)

    # ========================================================
    # RIGHT PANEL
    # ========================================================

    def _build_right_panel(self, parent):

        right = ctk.CTkFrame(parent)

        right.grid(
            row=0,
            column=1,
            sticky="nsew",
            padx=(8, 0)
        )

        right.grid_columnconfigure(
            0,
            weight=1
        )

        right.grid_rowconfigure(
            1,
            weight=1
        )

        self.tab_view = ctk.CTkTabview(
            right
        )

        self.tab_view.grid(
            row=0,
            column=0,
            rowspan=2,
            sticky="nsew",
            padx=16,
            pady=(8, 0)
        )

        self.tab_view.add("MSSI Results")
        self.tab_view.add("Molecular Index")
        self.tab_view.add("Egypt Map")
        self.tab_view.add("Gene Plots")
        self.tab_view.add("Interpretation")
        self.tab_view.add("Report")

        # MSSI Results
        self.results_frame = ctk.CTkScrollableFrame(
            self.tab_view.tab(
                "MSSI Results"
            )
        )

        self.results_frame.pack(
            fill="both",
            expand=True
        )

        # Molecular Index
        self.risk_frame = ctk.CTkScrollableFrame(
            self.tab_view.tab(
                "Molecular Index"
            )
        )

        self.risk_frame.pack(
            fill="both",
            expand=True
        )

        # Egypt Map
        self.map_frame = ctk.CTkFrame(
            self.tab_view.tab(
                "Egypt Map"
            ),
            fg_color="transparent"
        )

        self.map_frame.pack(
            fill="both",
            expand=True
        )

        # Gene plots
        self.gene_frame = ctk.CTkFrame(
            self.tab_view.tab(
                "Gene Plots"
            ),
            fg_color="transparent"
        )

        self.gene_frame.pack(
            fill="both",
            expand=True
        )

        # Interpretation
        self.ai_advisory_frame = (
            ctk.CTkScrollableFrame(
                self.tab_view.tab(
                    "Interpretation"
                )
            )
        )

        self.ai_advisory_frame.pack(
            fill="both",
            expand=True
        )

        self._build_report_tab()

        # Log
        ctk.CTkLabel(
            right,
            text="Log",
            font=ctk.CTkFont(
                weight="bold"
            )
        ).grid(
            row=2,
            column=0,
            padx=16,
            pady=(4, 2),
            sticky="w"
        )

        self.log_box = ctk.CTkTextbox(
            right,
            height=120,
            state="disabled"
        )

        self.log_box.grid(
            row=3,
            column=0,
            sticky="ew",
            padx=16,
            pady=(0, 16)
        )

        self._log(
            "Ready. Please select a file to begin."
        )

    # ========================================================
    # REPORT TAB
    # ========================================================

    def _build_report_tab(self):

        tab = self.tab_view.tab(
            "Report"
        )

        container = ctk.CTkFrame(
            tab,
            fg_color="transparent"
        )

        container.pack(
            expand=True
        )

        ctk.CTkLabel(
            container,
            text="Export Final Analysis Reports",
            font=ctk.CTkFont(
                size=18,
                weight="bold"
            )
        ).pack(
            pady=20
        )

        self.pdf_btn = ctk.CTkButton(
            container,
            text="Save PDF Report to Desktop",
            fg_color="#D85A30",
            hover_color="#B44826",
            height=45,
            width=250,
            command=self._export_pdf_custom
        )

        self.pdf_btn.pack(
            pady=10
        )

        self.json_btn = ctk.CTkButton(
            container,
            text="Save JSON Data to Desktop",
            height=45,
            width=250,
            command=self._export_json_custom
        )

        self.json_btn.pack(
            pady=10
        )

        self.status_lbl = ctk.CTkLabel(
            container,
            text="Analysis must be completed first",
            text_color="gray"
        )

        self.status_lbl.pack(
            pady=10
        )

    # ========================================================
    # FOOTER
    # ========================================================

    def _build_footer(self):

        footer = ctk.CTkFrame(
            self,
            fg_color="transparent",
            height=28
        )

        footer.grid(
            row=2,
            column=0,
            sticky="ew",
            padx=16
        )

        ctk.CTkLabel(
            footer,
            text=(
                "MSSI Pipeline v1.0.0  |  "
                "Ahmed Yassin  |  MIT License"
            ),
            text_color="gray",
            font=ctk.CTkFont(
                size=11
            )
        ).pack(
            side="right"
        )

    # ========================================================
    # LOGGING
    # ========================================================

    def _log(self, message):

        self.log_box.configure(
            state="normal"
        )

        self.log_box.insert(
            "end",
            message + "\n"
        )

        self.log_box.see(
            "end"
        )

        self.log_box.configure(
            state="disabled"
        )

    # ========================================================
    # FILE SELECTION
    # ========================================================

    def _browse_file(self):

        path = filedialog.askopenfilename(
            title="Select MSSI input file",
            filetypes=[
                (
                    "Supported input files",
                    "*.csv *.xlsx *.xls"
                ),
                (
                    "Corrected qPCR CSV",
                    "*.csv"
                ),
                (
                    "Excel files",
                    "*.xlsx *.xls"
                ),
                (
                    "All files",
                    "*.*"
                )
            ]
        )

        if not path:
            return

        try:

            loader = MSSILoader(
                path
            )

            genes = loader.get_genes()

        except Exception as exc:

            messagebox.showerror(
                "Invalid input file",
                str(exc)
            )

            self._log(
                f"Input validation error: {exc}"
            )

            return

        # Store the selected input file.
        input_path = Path(path).expanduser().resolve()
        self.filepath = str(input_path)

        # Automatically create one results folder beside the input file.
        result_dir = input_path.parent / "MSSI_Results"
        result_dir.mkdir(
            parents=True,
            exist_ok=True
        )

        # Visualizer, mapper and reporter all receive this same directory
        # from _analysis_task().
        self.output_dir = str(result_dir)

        self.corrected_qpcr_mode = (
            input_path.suffix.lower()
            == ".csv"
        )

        self.file_label.configure(
            text=input_path.name,
            text_color="white"
        )

        self.output_label.configure(
            text=str(result_dir),
            text_color="white"
        )

        self._log(
            f"File loaded: {input_path}"
        )

        self._log(
            f"Results will be saved to: {result_dir}"
        )

        if self.corrected_qpcr_mode:

            self._log(
                "Corrected qPCR mode enabled."
            )

            self._log(
                "Biological replicate DeltaCt values "
                "will be used."
            )

            self._log(
                "Population fold changes are calculated "
                "from mean DeltaCt values."
            )

            self._log(
                "Technical replicates are not treated "
                "as independent observations."
            )

        self._load_genes(
            preloaded_genes=genes
        )

    # ========================================================
    # OUTPUT DIRECTORY
    # ========================================================

    def _browse_output(self):

        path = filedialog.askdirectory()

        if path:

            self.output_dir = path

            self.output_label.configure(
                text=path,
                text_color="white"
            )

            self._log(
                f"Output directory: {path}"
            )

    # ========================================================
    # GENE LOADING
    # ========================================================

    def _load_genes(
        self,
        preloaded_genes=None
    ):

        try:

            if preloaded_genes is None:

                loader = MSSILoader(
                    self.filepath
                )

                self.genes = (
                    loader.get_genes()
                )

            else:

                self.genes = list(
                    preloaded_genes
                )

            for widget in (
                self.weights_frame
                .winfo_children()
            ):
                widget.destroy()

            self.weight_vars = {}

            for i, gene in enumerate(
                self.genes
            ):

                ctk.CTkLabel(
                    self.weights_frame,
                    text=gene
                ).grid(
                    row=i * 3,
                    column=0,
                    sticky="w",
                    pady=(4, 0)
                )

                var = ctk.DoubleVar(
                    value=1.0
                )

                self.weight_vars[
                    gene
                ] = var

                slider = ctk.CTkSlider(
                    self.weights_frame,
                    from_=0.1,
                    to=3.0,
                    variable=var,
                    width=200
                )

                slider.grid(
                    row=i * 3 + 1,
                    column=0,
                    sticky="ew",
                    pady=(0, 2)
                )

                value_label = ctk.CTkLabel(
                    self.weights_frame,
                    text="Weight = 1.00",
                    text_color="gray"
                )

                value_label.grid(
                    row=i * 3 + 2,
                    column=0,
                    sticky="w",
                    pady=(0, 4)
                )

                slider.configure(
                    command=lambda value,
                    label=value_label:
                    label.configure(
                        text=(
                            f"Weight = "
                            f"{float(value):.2f}"
                        )
                    )
                )

            self._log(
                "Genes found: "
                + ", ".join(
                    self.genes
                )
            )

            self._log(
                "Default gene weights = 1.00 "
                "(equal weighting)."
            )

        except Exception as exc:

            self._log(
                f"Error loading genes: {exc}"
            )

            messagebox.showerror(
                "Gene loading error",
                str(exc)
            )

    # ========================================================
    # RUN ANALYSIS
    # ========================================================

    def _run_analysis(self):

        if not self.filepath:

            self._log(
                "Please select a file first."
            )

            return

        if not self.genes:

            self._log(
                "No genes were loaded."
            )

            return

        if not self.output_dir:

            input_path = Path(
                self.filepath
            ).expanduser().resolve()

            result_dir = (
                input_path.parent
                / "MSSI_Results"
            )

            result_dir.mkdir(
                parents=True,
                exist_ok=True
            )

            self.output_dir = str(
                result_dir
            )

            self.output_label.configure(
                text=self.output_dir,
                text_color="white"
            )

        threading.Thread(
            target=self._analysis_task,
            daemon=True
        ).start()

    # ========================================================
    # ANALYSIS TASK
    # ========================================================

    def _analysis_task(self):

        try:

            self.run_btn.configure(
                state="disabled"
            )

            self.progress.set(
                0.10
            )

            # ------------------------------------------------
            # 1. Load data
            # ------------------------------------------------

            self._log(
                "[1/7] Loading data..."
            )

            loader = MSSILoader(
                self.filepath
            )

            data = loader.load()

            weights = {
                gene:
                self.weight_vars[gene].get()
                for gene
                in self.genes
            }

            self._log(
                "Gene weights: "
                + ", ".join(
                    f"{gene}={weight:.2f}"
                    for gene, weight
                    in weights.items()
                )
            )

            # ------------------------------------------------
            # 2. MSSI
            # ------------------------------------------------

            self.progress.set(
                0.25
            )

            self._log(
                "[2/7] Computing MSSI scores..."
            )

            calculator = MSSICalculator(
                data,
                weights=weights
            )

            stats = (
                calculator.get_full_stats()
            )

            normalized = (
                calculator.normalize_scores()
            )

            # ------------------------------------------------
            # 3. Classification
            # ------------------------------------------------

            self.progress.set(
                0.40
            )

            self._log(
                "[3/7] Classifying "
                "molecular stress levels..."
            )

            classifier = MSSIClassifier(
                stats,
                normalized
            )

            self.results = (
                classifier.classify_all()
            )

            # ------------------------------------------------
            # 4. Exploratory molecular index
            # ------------------------------------------------

            self.progress.set(
                0.55
            )

            self._log(
                "[4/7] Computing exploratory "
                "molecular index..."
            )

            risk_engine = (
                ResistanceRiskIndex(
                    self.results
                )
            )

            self.risks = (
                risk_engine.compute_all()
            )

            if self.corrected_qpcr_mode:

                self._log(
                    "IMPORTANT: Molecular index is "
                    "exploratory only and must not be "
                    "interpreted as phenotypic "
                    "insecticide resistance."
                )

            # ------------------------------------------------
            # 5. Plots
            # ------------------------------------------------

            self.progress.set(
                0.65
            )

            self._log(
                "[5/7] Generating plots..."
            )

            self._visualizer = (
                MSSIVisualizer(
                    self.results,
                    output_dir=self.output_dir
                )
            )

            self._visualizer.plot_all()

            # ------------------------------------------------
            # 6. MSSI Egypt map
            # ------------------------------------------------

            self.progress.set(
                0.78
            )

            self._log(
                "[6/7] Generating MSSI Egypt map..."
            )

            # CRITICAL FIX:
            # Map uses MSSI results,
            # NOT resistance/risk results.
            self._mapper = EgyptMapper(
                self.results,
                output_dir=self.output_dir
            )

            self._mapper.plot_mssi_bar()

            # ------------------------------------------------
            # 7. Report
            # ------------------------------------------------

            self.progress.set(
                0.90
            )

            self._log(
                "[7/7] Exporting report..."
            )

            reporter = MSSIReporter(
                self.results,
                output_dir=self.output_dir
            )

            reporter.export_all()

            self.progress.set(
                1.0
            )

            self._log(
                "Done. Results saved to: "
                f"{self.output_dir}/"
            )

            self.after(
                0,
                self._show_results
            )

            self.after(
                0,
                self._show_risks
            )

            self.after(
                0,
                self._show_map
            )

            self.after(
                0,
                self._show_gene_plots
            )

            self.after(
                0,
                self._show_interpretation
            )

            self.after(
                0,
                lambda:
                self.status_lbl.configure(
                    text=(
                        "Ready for Desktop export"
                    ),
                    text_color="#1D9E75"
                )
            )

        except Exception as exc:

            self._log(
                f"Error: {exc}"
            )

            self.progress.set(
                0
            )

            error_message = str(exc)

            self.after(
                0,
                lambda message=error_message:
                messagebox.showerror(
                    "Analysis error",
                    message
                )
            )

        finally:

            self.run_btn.configure(
                state="normal"
            )

    # ========================================================
    # MSSI RESULTS
    # ========================================================

    def _show_results(self):

        for widget in (
            self.results_frame
            .winfo_children()
        ):
            widget.destroy()

        if not self.results:
            return

        headers = [
            "Location",
            "MSSI Score",
            "Stress Level",
            "Confidence"
        ]

        for col, header in enumerate(
            headers
        ):

            ctk.CTkLabel(
                self.results_frame,
                text=header,
                font=ctk.CTkFont(
                    weight="bold"
                ),
                fg_color="#3C3489",
                corner_radius=4,
                width=150
            ).grid(
                row=0,
                column=col,
                padx=4,
                pady=4,
                sticky="ew"
            )

        for row, result in enumerate(
            self.results,
            start=1
        ):

            stress_color = (
                STRESS_COLORS.get(
                    result.stress_level,
                    "white"
                )
            )

            row_data = [
                result.location,
                f"{result.mssi_score:.4f}",
                result.stress_label,
                f"{int(result.confidence * 100)}%"
            ]

            for col, value in enumerate(
                row_data
            ):

                text_color = (
                    stress_color
                    if col == 2
                    else None
                )

                ctk.CTkLabel(
                    self.results_frame,
                    text=value,
                    text_color=text_color,
                    width=150
                ).grid(
                    row=row,
                    column=col,
                    padx=4,
                    pady=2,
                    sticky="ew"
                )

        highest = max(
            self.results,
            key=lambda r:
            r.mssi_score
        )

        lowest = min(
            self.results,
            key=lambda r:
            r.mssi_score
        )

        row = len(
            self.results
        ) + 2

        ctk.CTkLabel(
            self.results_frame,
            text=(
                f"Highest MSSI: "
                f"{highest.location} "
                f"({highest.mssi_score:.4f})"
            ),
            text_color="#E24B4A",
            font=ctk.CTkFont(
                weight="bold"
            )
        ).grid(
            row=row,
            column=0,
            columnspan=4,
            pady=(12, 2),
            sticky="w"
        )

        ctk.CTkLabel(
            self.results_frame,
            text=(
                f"Lowest MSSI: "
                f"{lowest.location} "
                f"({lowest.mssi_score:.4f})"
            ),
            text_color="#1D9E75",
            font=ctk.CTkFont(
                weight="bold"
            )
        ).grid(
            row=row + 1,
            column=0,
            columnspan=4,
            pady=(2, 12),
            sticky="w"
        )

    # ========================================================
    # EXPLORATORY MOLECULAR INDEX
    # ========================================================

    def _show_risks(self):

        for widget in (
            self.risk_frame
            .winfo_children()
        ):
            widget.destroy()

        if not self.risks:
            return

        warning = (
            "Exploratory molecular index only. "
            "No insecticide susceptibility bioassays "
            "were performed. These values must not be "
            "interpreted as probabilities or confirmation "
            "of phenotypic insecticide resistance."
        )

        ctk.CTkLabel(
            self.risk_frame,
            text=warning,
            text_color="#EF9F27",
            wraplength=650,
            justify="left"
        ).grid(
            row=0,
            column=0,
            columnspan=3,
            padx=8,
            pady=(8, 12),
            sticky="w"
        )

        headers = [
            "Location",
            "Index",
            "Level"
        ]

        for col, header in enumerate(
            headers
        ):

            ctk.CTkLabel(
                self.risk_frame,
                text=header,
                font=ctk.CTkFont(
                    weight="bold"
                ),
                fg_color="#3C3489",
                corner_radius=4,
                width=150
            ).grid(
                row=1,
                column=col,
                padx=4,
                pady=4,
                sticky="ew"
            )

        for row, result in enumerate(
            self.risks,
            start=2
        ):

            values = [
                result.location,
                f"{result.risk_index:.4f}",
                result.risk_level
                .replace("_", " ")
                .title()
            ]

            for col, value in enumerate(
                values
            ):

                ctk.CTkLabel(
                    self.risk_frame,
                    text=value,
                    width=150
                ).grid(
                    row=row,
                    column=col,
                    padx=4,
                    pady=2,
                    sticky="ew"
                )

    # ========================================================
    # EGYPT MAP
    # ========================================================

    def _show_map(self):

        for widget in (
            self.map_frame
            .winfo_children()
        ):
            widget.destroy()

        # IMPORTANT:
        # map depends on MSSI results.
        if not self.results:
            return

        if not self._mapper:
            return

        loading = ctk.CTkLabel(
            self.map_frame,
            text="Rendering MSSI map...",
            font=ctk.CTkFont(
                size=14
            ),
            text_color="gray"
        )

        loading.pack(
            expand=True
        )

        self.update()

        try:

            fig = self._mapper.plot_map(
                save=True
            )

            for widget in (
                self.map_frame
                .winfo_children()
            ):
                widget.destroy()

            canvas = FigureCanvasTkAgg(
                fig,
                master=self.map_frame
            )

            canvas.draw()

            canvas.get_tk_widget().pack(
                fill="both",
                expand=True
            )

        except Exception as exc:

            for widget in (
                self.map_frame
                .winfo_children()
            ):
                widget.destroy()

            ctk.CTkLabel(
                self.map_frame,
                text=(
                    f"Map error: {exc}"
                ),
                text_color="#E24B4A"
            ).pack(
                expand=True
            )

    # ========================================================
    # GENE PLOTS
    # ========================================================

    def _show_gene_plots(self):

        for widget in (
            self.gene_frame
            .winfo_children()
        ):
            widget.destroy()

        if (
            not self.results
            or not self._visualizer
        ):
            return

        try:

            fig = (
                self._visualizer
                .plot_publication_heatmap(
                    save=True
                )
            )

            canvas = FigureCanvasTkAgg(
                fig,
                master=self.gene_frame
            )

            canvas.draw()

            canvas.get_tk_widget().pack(
                fill="both",
                expand=True
            )

        except Exception as exc:

            ctk.CTkLabel(
                self.gene_frame,
                text=(
                    f"Plot error: {exc}"
                ),
                text_color="#E24B4A"
            ).pack(
                expand=True
            )

    # ========================================================
    # INTERPRETATION
    # ========================================================

    def _show_interpretation(self):

        for widget in (
            self.ai_advisory_frame
            .winfo_children()
        ):
            widget.destroy()

        if not self.results:
            return

        ctk.CTkLabel(
            self.ai_advisory_frame,
            text=(
                "Exploratory Molecular "
                "Stress Interpretation"
            ),
            font=ctk.CTkFont(
                size=16,
                weight="bold"
            ),
            text_color="#378ADD"
        ).pack(
            pady=10
        )

        warning = (
            "MSSI summarizes the relative expression "
            "patterns of the analyzed sHSP genes. "
            "It represents a molecular stress-response "
            "signature within this dataset and must not "
            "be interpreted as evidence or probability "
            "of phenotypic insecticide resistance."
        )

        ctk.CTkLabel(
            self.ai_advisory_frame,
            text=warning,
            wraplength=650,
            justify="left",
            text_color="#EF9F27"
        ).pack(
            padx=20,
            pady=15,
            anchor="w"
        )

        for result in self.results:

            card = ctk.CTkFrame(
                self.ai_advisory_frame
            )

            card.pack(
                fill="x",
                padx=15,
                pady=8
            )

            ctk.CTkLabel(
                card,
                text=result.location,
                font=ctk.CTkFont(
                    weight="bold"
                )
            ).pack(
                anchor="w",
                padx=10,
                pady=(10, 2)
            )

            ctk.CTkLabel(
                card,
                text=(
                    f"MSSI = "
                    f"{result.mssi_score:.4f}"
                )
            ).pack(
                anchor="w",
                padx=10
            )

            ctk.CTkLabel(
                card,
                text=(
                    f"Stress tier: "
                    f"{result.stress_label}"
                )
            ).pack(
                anchor="w",
                padx=10,
                pady=(2, 10)
            )

    # ========================================================
    # PDF EXPORT
    # ========================================================

    def _export_pdf_custom(self):

        if not self.results:

            messagebox.showwarning(
                "Warning",
                "Run analysis first."
            )

            return

        path = filedialog.asksaveasfilename(
            defaultextension=".pdf",
            initialfile="MSSI_Report.pdf",
            title="Save PDF Report"
        )

        if not path:
            return

        reporter = MSSIReporter(
            self.results,
            output_dir=Path(path).parent
        )

        reporter.export_pdf(
            filename=Path(path).name
        )

        messagebox.showinfo(
            "Success",
            f"PDF report saved at:\n{path}"
        )

    # ========================================================
    # JSON EXPORT
    # ========================================================

    def _export_json_custom(self):

        if not self.results:

            messagebox.showwarning(
                "Warning",
                "Run analysis first."
            )

            return

        path = filedialog.asksaveasfilename(
            defaultextension=".json",
            initialfile="MSSI_Data.json",
            title="Save JSON Results"
        )

        if not path:
            return

        export_data = {
            "input_file": str(
                self.filepath
            ),
            "input_mode": (
                "corrected_qpcr_delta_ct"
                if self.corrected_qpcr_mode
                else "legacy_excel"
            ),
            "gene_weights": {
                gene:
                self.weight_vars[
                    gene
                ].get()
                for gene
                in self.genes
            },
            "mssi_results": [
                result.__dict__
                for result
                in self.results
            ],
            "note": (
                "MSSI represents an integrated "
                "molecular stress-response signature "
                "and not phenotypic insecticide "
                "resistance."
            )
        }

        if self.risks:

            export_data[
                "exploratory_molecular_index"
            ] = [
                result.__dict__
                for result
                in self.risks
            ]

        with open(
            path,
            "w",
            encoding="utf-8"
        ) as file:

            json.dump(
                export_data,
                file,
                indent=4,
                ensure_ascii=False
            )

        messagebox.showinfo(
            "Success",
            f"JSON data saved at:\n{path}"
        )


# ============================================================
# ENTRY POINT
# ============================================================

def main():

    app = MSSIApp()
    app.mainloop()


if __name__ == "__main__":
    main()