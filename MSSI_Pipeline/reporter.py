import json
from datetime import datetime
from pathlib import Path
from typing import Dict, List
from fpdf import FPDF
from .classifier import StressResult

STRESS_COLORS_RGB = {
    "very_high": (226, 75,  74),
    "high":      (239, 159, 39),
    "moderate":  (55,  138, 221),
    "low":       (29,  158, 117),
}

class MSSIReporter:

    def __init__(self, results: List[StressResult], output_dir: str = "outputs"):
        self.results    = results
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(exist_ok=True)
        self.timestamp  = datetime.now().strftime("%Y-%m-%d %H:%M")

    def _add_header(self, pdf: FPDF):
        pdf.set_fill_color(83, 74, 183)
        pdf.rect(0, 0, 210, 30, style="F")
        pdf.set_text_color(255, 255, 255)
        pdf.set_font("Helvetica", "B", 16)
        pdf.set_xy(10, 8)
        pdf.cell(0, 8, "MSSI Pipeline - Stress Analysis Report", ln=True)
        pdf.set_font("Helvetica", size=9)
        pdf.set_x(10)
        pdf.cell(0, 6, f"Generated: {self.timestamp}  |  Species: Spodoptera frugiperda")
        pdf.set_text_color(0, 0, 0)
        pdf.ln(12)

    def _add_summary_table(self, pdf: FPDF):
        pdf.set_font("Helvetica", "B", 12)
        pdf.cell(0, 8, "Summary Results", ln=True)
        pdf.ln(2)

        headers = ["Location", "MSSI Score", "Stress Level", "Confidence"]
        widths  = [50, 35, 75, 30]

        pdf.set_fill_color(240, 240, 248)
        pdf.set_font("Helvetica", "B", 10)
        for h, w in zip(headers, widths):
            pdf.cell(w, 8, h, border=1, fill=True)
        pdf.ln()

        pdf.set_font("Helvetica", size=10)
        for r in self.results:
            rgb = STRESS_COLORS_RGB[r.stress_level]
            pdf.set_fill_color(*rgb)
            pdf.set_text_color(255, 255, 255)
            row = [
                r.location,
                f"{r.mssi_score:.4f}",
                r.stress_label,
                f"{int(r.confidence * 100)}%"
            ]
            for val, w in zip(row, widths):
                pdf.cell(w, 8, val, border=1, fill=True)
            pdf.ln()
            pdf.set_text_color(0, 0, 0)
        pdf.ln(6)

    def _add_gene_table(self, pdf: FPDF):
        pdf.set_font("Helvetica", "B", 12)
        pdf.cell(0, 8, "Gene Expression Details", ln=True)
        pdf.ln(2)

        genes   = list(self.results[0].gene_profile.keys())
        headers = ["Location"] + genes
        widths  = [50] + [46] * len(genes)

        pdf.set_fill_color(240, 240, 248)
        pdf.set_font("Helvetica", "B", 10)
        for h, w in zip(headers, widths):
            pdf.cell(w, 8, h, border=1, fill=True)
        pdf.ln()

        pdf.set_font("Helvetica", size=10)
        for r in self.results:
            pdf.cell(50, 8, r.location, border=1)
            for gene in genes:
                fc = r.gene_profile[gene]["fold_change"]
                pdf.cell(46, 8, f"{fc:.2f}x", border=1)
            pdf.ln()
        pdf.ln(6)

    def _add_images(self, pdf: FPDF):
        images = [
            ("mssi_bar.png",      "MSSI Score by Location"),
            ("heatmap.png",        "Gene Expression Heatmap"),
            ("publication_heatmap.png", "Publication Quality Heatmap"),
            ("gene_profiles.png",  "Gene Profiles by Location"),
            ("confidence.png",     "Classification Confidence"),
        ]
        for filename, title in images:
            path = self.output_dir / filename
            if path.exists():
                pdf.add_page()
                pdf.set_font("Helvetica", "B", 12)
                pdf.cell(0, 10, title, ln=True)
                pdf.image(str(path), x=10, w=185)

    def export_pdf(self, filename: str = "mssi_report.pdf") -> Path:
        pdf = FPDF()
        pdf.add_page()
        pdf.set_auto_page_break(auto=True, margin=15)
        self._add_header(pdf)
        self._add_summary_table(pdf)
        self._add_gene_table(pdf)
        self._add_images(pdf)
        out = self.output_dir / filename
        pdf.output(str(out))
        return out

    def export_json(self, filename: str = "mssi_results.json") -> Path:
        data = [
            {
                "location":          r.location,
                "mssi_score":        r.mssi_score,
                "normalized_score":  r.normalized_score,
                "stress_level":      r.stress_level,
                "stress_label":      r.stress_label,
                "confidence":        r.confidence,
                "gene_profile":      r.gene_profile,
            }
            for r in self.results
        ]
        out = self.output_dir / filename
        with open(out, "w") as f:
            json.dump(data, f, indent=2)
        return out

    def export_all(self):
        self.export_pdf()
        self.export_json()