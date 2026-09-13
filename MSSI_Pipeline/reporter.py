import json
from datetime import datetime
from pathlib import Path
from typing import List, Optional

from fpdf import FPDF

from .classifier import StressResult


STRESS_COLORS_RGB = {
    "very_high": (226, 75, 74),
    "high": (239, 159, 39),
    "moderate": (55, 138, 221),
    "low": (29, 158, 117),
}

STRESS_LABELS = {
    "very_high": "Critical",
    "high": "High",
    "moderate": "Moderate",
    "low": "Low",
}


class MSSIReporter:

    def __init__(
        self,
        results: List[StressResult],
        input_file: Optional[str] = None,
        output_dir: Optional[str] = None,
    ):
        if not results:
            raise ValueError("No StressResult objects were supplied to MSSIReporter.")

        self.results = results
        self.timestamp = datetime.now().strftime("%Y-%m-%d %H:%M")

        # -------------------------------------------------------------
        # Output directory logic
        #
        # Priority:
        # 1) Explicit output_dir
        # 2) MSSI_Results folder beside selected input file
        # 3) Local fallback folder only if neither is supplied
        # -------------------------------------------------------------
        if output_dir:
            self.output_dir = Path(output_dir).expanduser().resolve()

        elif input_file:
            input_path = Path(input_file).expanduser().resolve()
            self.output_dir = input_path.parent / "MSSI_Results"

        else:
            self.output_dir = Path("MSSI_Results").resolve()

        self.output_dir.mkdir(parents=True, exist_ok=True)

    # =================================================================
    # PDF HELPERS
    # =================================================================

    def _safe_stress_label(self, result: StressResult) -> str:
        """
        Return publication-facing MSSI tier label.
        Keeps compatibility with classifier values such as 'very_high'.
        """

        if hasattr(result, "stress_label") and result.stress_label:
            label = str(result.stress_label).strip()

            if label.lower() == "very high":
                return "Critical"

            return label

        return STRESS_LABELS.get(
            result.stress_level,
            str(result.stress_level).replace("_", " ").title(),
        )

    def _add_header(self, pdf: FPDF):
        pdf.set_fill_color(83, 74, 183)
        pdf.rect(0, 0, 210, 30, style="F")

        pdf.set_text_color(255, 255, 255)
        pdf.set_font("Helvetica", "B", 16)

        pdf.set_xy(10, 7)
        pdf.cell(
            0,
            8,
            "MSSI Pipeline - Molecular Expression Analysis Report",
            ln=True,
        )

        pdf.set_font("Helvetica", size=9)
        pdf.set_x(10)

        pdf.cell(
            0,
            6,
            f"Generated: {self.timestamp}  |  "
            "Species: Spodoptera frugiperda",
        )

        pdf.set_text_color(0, 0, 0)
        pdf.ln(13)

    def _add_interpretation_note(self, pdf: FPDF):
        """
        Add a short interpretation disclaimer consistent with the manuscript.
        """

        pdf.set_font("Helvetica", "B", 10)
        pdf.cell(
            0,
            6,
            "Interpretation note",
            ln=True,
        )

        pdf.set_font("Helvetica", size=9)

        note = (
            "MSSI scores and tiers are exploratory computational summaries "
            "of the selected sHSP expression profiles. They should not be "
            "interpreted as validated measures or probabilities of "
            "insecticide resistance, phenotypic stress tolerance, or "
            "environmental adaptation."
        )

        pdf.multi_cell(
            0,
            5,
            note,
        )

        pdf.ln(4)

    def _add_summary_table(self, pdf: FPDF):
        pdf.set_font("Helvetica", "B", 12)

        pdf.cell(
            0,
            8,
            "Summary Results",
            ln=True,
        )

        pdf.ln(2)

        headers = [
            "Location",
            "MSSI Score",
            "MSSI Tier",
            "Computational Confidence",
        ]

        widths = [
            45,
            32,
            48,
            65,
        ]

        pdf.set_fill_color(240, 240, 248)
        pdf.set_font("Helvetica", "B", 9)

        for header, width in zip(headers, widths):
            pdf.cell(
                width,
                8,
                header,
                border=1,
                fill=True,
                align="C",
            )

        pdf.ln()

        pdf.set_font("Helvetica", size=9)

        for result in self.results:
            rgb = STRESS_COLORS_RGB.get(
                result.stress_level,
                (128, 128, 128),
            )

            pdf.set_fill_color(*rgb)
            pdf.set_text_color(255, 255, 255)

            stress_label = self._safe_stress_label(result)

            confidence = getattr(
                result,
                "confidence",
                0.0,
            )

            row = [
                str(result.location),
                f"{float(result.mssi_score):.2f}",
                stress_label,
                f"{float(confidence) * 100:.0f}%",
            ]

            for value, width in zip(row, widths):
                pdf.cell(
                    width,
                    8,
                    value,
                    border=1,
                    fill=True,
                    align="C",
                )

            pdf.ln()

            pdf.set_text_color(0, 0, 0)

        pdf.ln(6)

    def _add_gene_table(self, pdf: FPDF):
        pdf.set_font("Helvetica", "B", 12)

        pdf.cell(
            0,
            8,
            "Gene Expression Details",
            ln=True,
        )

        pdf.ln(2)

        genes = list(
            self.results[0].gene_profile.keys()
        )

        headers = [
            "Location"
        ] + genes

        number_of_genes = len(genes)

        if number_of_genes == 0:
            return

        # Keep table within A4 printable width.
        location_width = 42
        remaining_width = 190 - location_width
        gene_width = remaining_width / number_of_genes

        widths = [
            location_width
        ] + [
            gene_width
        ] * number_of_genes

        pdf.set_fill_color(
            240,
            240,
            248,
        )

        pdf.set_font(
            "Helvetica",
            "B",
            8,
        )

        for header, width in zip(
            headers,
            widths,
        ):
            pdf.cell(
                width,
                8,
                str(header),
                border=1,
                fill=True,
                align="C",
            )

        pdf.ln()

        pdf.set_font(
            "Helvetica",
            size=9,
        )

        for result in self.results:
            pdf.cell(
                location_width,
                8,
                str(result.location),
                border=1,
                align="C",
            )

            for gene in genes:
                profile = result.gene_profile[gene]

                fold_change = float(
                    profile.get(
                        "fold_change",
                        0.0,
                    )
                )

                pdf.cell(
                    gene_width,
                    8,
                    f"{fold_change:.2f}",
                    border=1,
                    align="C",
                )

            pdf.ln()

        pdf.ln(4)

        pdf.set_font(
            "Helvetica",
            size=8,
        )

        pdf.multi_cell(
            0,
            5,
            "Values represent population-level relative expression "
            "(2^-DeltaDeltaCt) relative to the control calibrator.",
        )

        pdf.ln(4)

    def _add_images(self, pdf: FPDF):
        """
        Add publication figures.

        New publication names are checked first.
        Legacy names are retained as fallback for GUI compatibility.
        """

        image_groups = [
            (
                [
                    "Figure1_relative_mRNA_expression.png",
                    "gene_profiles.png",
                ],
                "Figure 1. Relative mRNA Expression",
            ),
            (
                [
                    "Figure2_expression_heatmap.png",
                    "publication_heatmap.png",
                    "heatmap.png",
                ],
                "Figure 2. Relative Expression Heatmap",
            ),
            (
                [
                    "Figure3_MSSI_scores.png",
                    "mssi_bar.png",
                ],
                "Figure 3. Multi-HSP Stress Signature Index",
            ),
            (
                [
                    "egypt_mssi_map.png",
                ],
                "Figure 4. Geographic Distribution of MSSI Signatures",
            ),
        ]

        for candidate_names, title in image_groups:
            selected_path = None

            for filename in candidate_names:
                candidate = self.output_dir / filename

                if candidate.exists():
                    selected_path = candidate
                    break

            if selected_path is None:
                continue

            pdf.add_page()

            pdf.set_font(
                "Helvetica",
                "B",
                12,
            )

            pdf.cell(
                0,
                10,
                title,
                ln=True,
            )

            # Leave enough margins for publication plots.
            pdf.image(
                str(selected_path),
                x=10,
                w=185,
            )

    # =================================================================
    # EXPORT PDF
    # =================================================================

    def export_pdf(
        self,
        filename: str = "mssi_report.pdf",
    ) -> Path:

        pdf = FPDF(
            orientation="P",
            unit="mm",
            format="A4",
        )

        pdf.set_auto_page_break(
            auto=True,
            margin=15,
        )

        pdf.add_page()

        self._add_header(pdf)
        self._add_interpretation_note(pdf)
        self._add_summary_table(pdf)
        self._add_gene_table(pdf)
        self._add_images(pdf)

        output_path = (
            self.output_dir /
            filename
        )

        pdf.output(
            str(output_path)
        )

        return output_path

    # =================================================================
    # EXPORT JSON
    # =================================================================

    def export_json(
        self,
        filename: str = "mssi_results.json",
    ) -> Path:

        data = []

        for result in self.results:
            result_data = {
                "location": result.location,
                "mssi_score": float(
                    result.mssi_score
                ),
                "stress_level_internal": result.stress_level,
                "mssi_tier": self._safe_stress_label(
                    result
                ),
                "computational_confidence": float(
                    getattr(
                        result,
                        "confidence",
                        0.0,
                    )
                ),
                "gene_profile": result.gene_profile,
            }

            # Retain normalized score only if available.
            if hasattr(
                result,
                "normalized_score",
            ):
                result_data[
                    "normalized_score"
                ] = float(
                    result.normalized_score
                )

            data.append(
                result_data
            )

        output_path = (
            self.output_dir /
            filename
        )

        with open(
            output_path,
            "w",
            encoding="utf-8",
        ) as file:
            json.dump(
                data,
                file,
                indent=2,
                ensure_ascii=False,
            )

        return output_path

    # =================================================================
    # EXPORT ALL
    # =================================================================

    def export_all(self):
        """
        Export final PDF and JSON analysis outputs.

        Returns:
            Dictionary containing generated paths.
        """

        pdf_path = self.export_pdf()
        json_path = self.export_json()

        return {
            "pdf": pdf_path,
            "json": json_path,
            "output_dir": self.output_dir,
        }