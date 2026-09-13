import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import seaborn as sns
import numpy as np
import pandas as pd

from typing import List, Optional, Tuple
from pathlib import Path

from .classifier import StressResult


# ---------------------------------------------------------------------
# Internal classifier colors
# Keep "very_high" internally for compatibility with the classifier,
# but display it as "Critical" in publication figures.
# ---------------------------------------------------------------------
STRESS_COLORS = {
    "very_high": "#E24B4A",
    "high": "#EF9F27",
    "moderate": "#378ADD",
    "low": "#1D9E75",
}

STRESS_LABELS = {
    "very_high": "Critical",
    "high": "High",
    "moderate": "Moderate",
    "low": "Low",
}

# Order used in publication figures
LOCATION_ORDER = [
    "Beni Suef",
    "Sohag",
    "Giza",
    "Qalyubia",
]

GENE_ORDER = [
    "Hsp 19.74",
    "Hsp 20.7",
    "Hsp 19.07",
]

GENE_DISPLAY_NAMES = {
    "Hsp 19.74": "sHsp19.74",
    "Hsp 20.7": "sHsp20.7",
    "Hsp 19.07": "sHsp19.07",
    "sHsp19.74": "sHsp19.74",
    "sHsp20.7": "sHsp20.7",
    "sHsp19.07": "sHsp19.07",
}


class MSSIVisualizer:
    """
    Publication-quality visualization module for the MSSI Pipeline.

    Expected gene_profile structure for each StressResult:

        result.gene_profile[gene]["fold_change"]

    Figure 1 uses a point-based display rather than bars or error bars.
    The two independent biological-replicate relative-expression values
    are shown directly, and a short horizontal segment marks the
    population-level fold change calculated from mean DeltaCt.

    Where raw relative-expression replicate values are not directly
    available in gene_profile, the two values are reconstructed exactly
    from the population fold change and DeltaCt SD for n = 2 biological
    replicates. This is possible because, for two observations, the
    sample SD uniquely determines their symmetric distance from the mean
    on the DeltaCt scale.

    The y-axis is logarithmic (base 2), which is appropriate for
    multiplicative 2^(-DeltaDeltaCt) values and prevents extreme
    biological variability from compressing the rest of the data.
    """

    def __init__(
        self,
        results: List[StressResult],
        output_dir: str = "outputs",
    ):
        if not results:
            raise ValueError("No StressResult objects were supplied to MSSIVisualizer.")

        self.results = results
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)

        # Publication-style plotting
        sns.set_theme(
            style="whitegrid",
            context="paper",
            font_scale=1.15,
        )

        plt.rcParams.update(
            {
                "figure.facecolor": "white",
                "axes.facecolor": "white",
                "savefig.facecolor": "white",
                "axes.edgecolor": "black",
                "axes.labelcolor": "black",
                "xtick.color": "black",
                "ytick.color": "black",
                "text.color": "black",
                "font.family": "DejaVu Sans",
                "axes.titleweight": "bold",
                "axes.titlesize": 12,
                "axes.labelsize": 11,
                "xtick.labelsize": 9,
                "ytick.labelsize": 9,
                "legend.fontsize": 9,
            }
        )

    # =================================================================
    # GENERAL HELPERS
    # =================================================================

    def _ordered_results(self) -> List[StressResult]:
        """
        Return results in manuscript order whenever possible.
        """

        result_map = {
            str(r.location).strip(): r
            for r in self.results
        }

        ordered = []

        for location in LOCATION_ORDER:
            if location in result_map:
                ordered.append(result_map[location])

        # Append any additional locations not included in LOCATION_ORDER.
        existing = {r.location for r in ordered}

        for result in self.results:
            if result.location not in existing:
                ordered.append(result)

        return ordered

    def _get_colors(self) -> List[str]:
        return [
            STRESS_COLORS.get(r.stress_level, "#808080")
            for r in self._ordered_results()
        ]

    def _get_gene_names(self) -> List[str]:
        """
        Preserve the actual gene keys stored by the pipeline while ordering
        the three target genes consistently.
        """

        available = list(self.results[0].gene_profile.keys())

        normalized_lookup = {}

        for gene in available:
            cleaned = (
                str(gene)
                .replace("sHsp", "Hsp ")
                .replace("Hsp  ", "Hsp ")
                .strip()
            )
            normalized_lookup[cleaned] = gene

        ordered = []

        for target in GENE_ORDER:
            if target in normalized_lookup:
                ordered.append(normalized_lookup[target])

        # If naming differs unexpectedly, keep remaining genes.
        for gene in available:
            if gene not in ordered:
                ordered.append(gene)

        return ordered

    def _display_gene_name(self, gene: str) -> str:
        cleaned = (
            str(gene)
            .replace("sHsp", "Hsp ")
            .replace("Hsp  ", "Hsp ")
            .strip()
        )

        return GENE_DISPLAY_NAMES.get(
            gene,
            GENE_DISPLAY_NAMES.get(cleaned, str(gene)),
        )

    def _get_matrix(self) -> pd.DataFrame:
        genes = self._get_gene_names()
        ordered_results = self._ordered_results()

        locations = [r.location for r in ordered_results]

        matrix = np.array(
            [
                [
                    float(r.gene_profile[gene]["fold_change"])
                    for gene in genes
                ]
                for r in ordered_results
            ],
            dtype=float,
        )

        display_genes = [
            self._display_gene_name(gene)
            for gene in genes
        ]

        return pd.DataFrame(
            matrix,
            index=locations,
            columns=display_genes,
        )

    def _extract_fc_errors(
        self,
        result: StressResult,
        gene: str,
    ) -> Tuple[float, float, str]:
        """
        Extract asymmetric fold-change error lengths for Figure 1.

        Preferred corrected-qPCR fields:
            fc_error_lower_sd
            fc_error_upper_sd

        These values are derived from biological-replicate DeltaCt
        variability and transformed to the fold-change scale.

        Legacy fallback:
            symmetric SD > SEM > SE

        Returns
        -------
        lower_error, upper_error, label
        """

        profile = result.gene_profile[gene]

        lower_sd = profile.get("fc_error_lower_sd")
        upper_sd = profile.get("fc_error_upper_sd")

        if lower_sd is not None and upper_sd is not None:
            lower_value = float(lower_sd)
            upper_value = float(upper_sd)

            if (
                np.isfinite(lower_value)
                and np.isfinite(upper_value)
                and lower_value >= 0
                and upper_value >= 0
            ):
                return (
                    lower_value,
                    upper_value,
                    "SD transformed from DeltaCt"
                )

        lower_se = profile.get("fc_error_lower_se")
        upper_se = profile.get("fc_error_upper_se")

        if lower_se is not None and upper_se is not None:
            lower_value = float(lower_se)
            upper_value = float(upper_se)

            if (
                np.isfinite(lower_value)
                and np.isfinite(upper_value)
                and lower_value >= 0
                and upper_value >= 0
            ):
                return (
                    lower_value,
                    upper_value,
                    "SEM transformed from DeltaCt"
                )

        # ---------------------------------------------------------
        # Legacy fallback for historical files that do not contain
        # corrected DeltaCt-derived asymmetric bounds.
        # ---------------------------------------------------------

        for key, label in (
            ("sd", "SD"),
            ("sem", "SEM"),
            ("se", "SEM"),
        ):
            value = profile.get(key)

            if value is None:
                continue

            value = float(value)

            if np.isfinite(value) and value >= 0:
                fc = float(profile["fold_change"])

                # Prevent the lower visual bound from crossing zero.
                lower_error = min(value, fc)
                upper_error = value

                return (
                    lower_error,
                    upper_error,
                    f"{label} (legacy fallback)"
                )

        return (
            np.nan,
            np.nan,
            "Not available"
        )

    def _extract_biological_points(
        self,
        result: StressResult,
        gene: str,
    ) -> np.ndarray:
        """
        Return biological-replicate relative-expression values for Figure 1.

        Preferred source
        ----------------
        If the classifier/gene profile contains an explicit ``replicates``
        field, those values are used directly.

        Corrected qPCR fallback
        -----------------------
        The current corrected pipeline stores:

            fold_change
            delta_ct_sd
            n_biological

        but older StressResult objects may not carry the original replicate
        array. For exactly n = 2 biological replicates, the two DeltaCt
        observations can be reconstructed from their mean and sample SD:

            x1, x2 = mean ± SD / sqrt(2)

        Because population fold change is:

            FC = 2^[-(mean DeltaCt_population - mean DeltaCt_control)]

        the corresponding individual relative-expression values are:

            FC * 2^(-SD/sqrt(2))
            FC * 2^(+SD/sqrt(2))

        This reconstruction is exact for n = 2, apart from any rounding
        already present in the stored SD.
        """

        profile = result.gene_profile[gene]

        # ---------------------------------------------------------
        # Preferred: explicit biological-replicate expression values
        # ---------------------------------------------------------

        explicit = profile.get("replicates")

        if explicit is not None:

            values = np.asarray(
                explicit,
                dtype=float,
            )

            values = values[
                np.isfinite(values)
            ]

            if len(values) > 0:
                return values

        # ---------------------------------------------------------
        # Exact n = 2 reconstruction from DeltaCt SD
        # ---------------------------------------------------------

        n_biological = profile.get(
            "n_biological"
        )

        delta_ct_sd = profile.get(
            "delta_ct_sd"
        )

        fold_change = profile.get(
            "fold_change"
        )

        if (
            n_biological == 2
            and delta_ct_sd is not None
            and fold_change is not None
        ):

            sd = float(
                delta_ct_sd
            )

            fc = float(
                fold_change
            )

            if (
                np.isfinite(sd)
                and np.isfinite(fc)
                and sd >= 0
                and fc > 0
            ):

                delta = (
                    sd
                    / np.sqrt(2.0)
                )

                lower_point = (
                    fc
                    * (2.0 ** (-delta))
                )

                upper_point = (
                    fc
                    * (2.0 ** delta)
                )

                return np.asarray(
                    [
                        lower_point,
                        upper_point,
                    ],
                    dtype=float,
                )

        return np.asarray(
            [],
            dtype=float,
        )

    def _save_figure(
            self,
            fig: plt.Figure,
            basename: str,
    ) -> None:
        """
        Save publication figures in high-quality formats.

        PNG:
            High-resolution 600 dpi image suitable for direct
            insertion into Microsoft Word.

        TIFF:
            High-resolution 600 dpi image suitable for journal
            submission when required.

        PDF:
            Vector-format copy retained for publication/archive use.
        """

        # ---------------------------------------------------------
        # Output paths
        # ---------------------------------------------------------

        png_path = (
                self.output_dir
                / f"{basename}.png"
        )

        tiff_path = (
                self.output_dir
                / f"{basename}.tiff"
        )

        pdf_path = (
                self.output_dir
                / f"{basename}.pdf"
        )

        # ---------------------------------------------------------
        # PNG
        # Main image for Microsoft Word
        # ---------------------------------------------------------

        fig.savefig(
            png_path,
            format="png",
            dpi=600,
            bbox_inches="tight",
            pad_inches=0.15,
            facecolor="white",
            edgecolor="none",
        )

        # ---------------------------------------------------------
        # TIFF
        # High-resolution journal copy
        # ---------------------------------------------------------

        fig.savefig(
            tiff_path,
            format="tiff",
            dpi=600,
            bbox_inches="tight",
            pad_inches=0.15,
            facecolor="white",
            edgecolor="none",
        )

        # ---------------------------------------------------------
        # PDF
        # Vector archive/publication copy
        # ---------------------------------------------------------

        fig.savefig(
            pdf_path,
            format="pdf",
            bbox_inches="tight",
            pad_inches=0.15,
            facecolor="white",
            edgecolor="none",
        )

    # =================================================================
    # FIGURE 1
    # RELATIVE mRNA EXPRESSION
    # =================================================================

    def plot_gene_profiles(
            self,
            save: bool = True,
    ) -> plt.Figure:
        """
        Figure 1.

        Publication-ready relative mRNA expression figure.

        Final Word-friendly layout
        --------------------------
        - Three panels are arranged vertically (3 rows × 1 column).
        - Two biological replicates are shown as individual open circles.
        - A short horizontal segment marks the population-level fold change
          calculated from mean DeltaCt.
        - No bars are used because the population fold change is not the
          arithmetic mean of the transformed replicate values.
        - No error bars are used because n = 2 provides an unstable estimate
          of dispersion after exponential transformation.
        - The y-axis is logarithmic (base 2), matching the multiplicative
          nature of 2^(-DeltaDeltaCt) fold-change values.
        - The dashed horizontal line at 1.0 represents the control calibrator.

        Formal statistical inference remains based on DeltaCt values.
        """

        # =========================================================
        # DATA
        # =========================================================

        genes = self._get_gene_names()
        ordered_results = self._ordered_results()

        if len(genes) != 3:
            raise ValueError(
                "Figure 1 expects exactly three sHSP genes, "
                f"but {len(genes)} genes were found: {genes}"
            )

        locations = [
            result.location
            for result in ordered_results
        ]

        x_positions = np.arange(
            len(locations),
            dtype=float,
        )

        # =========================================================
        # WORD-FRIENDLY FIGURE
        #
        # Vertical layout prevents the three panels from becoming
        # too small when the figure is inserted into a portrait page.
        # =========================================================

        fig, axes = plt.subplots(
            3,
            1,
            figsize=(8.5, 12.5),
            constrained_layout=False,
        )

        fig.patch.set_facecolor(
            "white"
        )

        panel_labels = [
            "A",
            "B",
            "C",
        ]

        gene_colors = [
            "#D55E00",
            "#0072B2",
            "#009E73",
        ]

        # Small x offsets keep the two biological replicates visible.
        point_offsets = np.asarray(
            [-0.10, 0.10],
            dtype=float,
        )

        # =========================================================
        # PANELS
        # =========================================================

        for panel_idx, (
                ax,
                gene,
                marker_color,
        ) in enumerate(
            zip(
                axes,
                genes,
                gene_colors,
            )
        ):

            all_positive_values = [1.0]

            # -----------------------------------------------------
            # Biological replicate points + population FC
            # -----------------------------------------------------

            for idx_result, result in enumerate(
                    ordered_results
            ):

                profile = result.gene_profile[
                    gene
                ]

                fold_change = float(
                    profile["fold_change"]
                )

                if (
                        not np.isfinite(fold_change)
                        or fold_change <= 0
                ):
                    raise ValueError(
                        f"Non-positive or invalid fold change for "
                        f"{gene} / {result.location}: {fold_change}"
                    )

                all_positive_values.append(
                    fold_change
                )

                # -------------------------------------------------
                # Individual biological replicates
                # -------------------------------------------------

                replicate_values = (
                    self._extract_biological_points(
                        result,
                        gene,
                    )
                )

                replicate_values = np.asarray(
                    replicate_values,
                    dtype=float,
                )

                replicate_values = replicate_values[
                    np.isfinite(
                        replicate_values
                    )
                    & (
                            replicate_values > 0
                    )
                    ]

                if len(replicate_values) > 0:

                    all_positive_values.extend(
                        replicate_values.tolist()
                    )

                    if len(replicate_values) == 1:

                        x_rep = np.asarray(
                            [
                                x_positions[
                                    idx_result
                                ]
                            ],
                            dtype=float,
                        )

                    elif len(replicate_values) == 2:

                        x_rep = (
                                x_positions[
                                    idx_result
                                ]
                                + point_offsets
                        )

                    else:

                        x_rep = np.linspace(
                            x_positions[
                                idx_result
                            ]
                            - 0.14,
                            x_positions[
                                idx_result
                            ]
                            + 0.14,
                            len(
                                replicate_values
                            ),
                        )

                    ax.scatter(
                        x_rep,
                        replicate_values,
                        s=75,
                        facecolors="white",
                        edgecolors=marker_color,
                        linewidths=1.8,
                        zorder=5,
                    )

                # -------------------------------------------------
                # Population-level FC from mean DeltaCt
                # -------------------------------------------------

                ax.hlines(
                    y=fold_change,
                    xmin=(
                            x_positions[
                                idx_result
                            ]
                            - 0.23
                    ),
                    xmax=(
                            x_positions[
                                idx_result
                            ]
                            + 0.23
                    ),
                    color="black",
                    linewidth=2.6,
                    zorder=6,
                )

                # -------------------------------------------------
                # Numerical FC label
                # -------------------------------------------------

                ax.annotate(
                    f"{fold_change:.2f}",
                    xy=(
                        x_positions[
                            idx_result
                        ],
                        fold_change,
                    ),
                    xytext=(
                        0,
                        8,
                    ),
                    textcoords="offset points",
                    ha="center",
                    va="bottom",
                    fontsize=11,
                    fontweight="bold",
                    zorder=7,
                )

            # =====================================================
            # CONTROL CALIBRATOR
            # =====================================================

            ax.axhline(
                1.0,
                color="gray",
                linestyle="--",
                linewidth=1.5,
                zorder=1,
            )

            # =====================================================
            # LOG2 SCALE
            # =====================================================

            ax.set_yscale(
                "log",
                base=2,
            )

            positive_values = np.asarray(
                all_positive_values,
                dtype=float,
            )

            data_min = float(
                np.min(
                    positive_values
                )
            )

            data_max = float(
                np.max(
                    positive_values
                )
            )

            lower_exp = int(
                np.floor(
                    np.log2(
                        data_min
                    )
                )
            )

            upper_exp = int(
                np.ceil(
                    np.log2(
                        data_max
                    )
                )
            )

            y_min = (
                            2.0 ** lower_exp
                    ) / 1.20

            y_max = (
                            2.0 ** upper_exp
                    ) * 1.20

            ax.set_ylim(
                y_min,
                y_max,
            )

            # -----------------------------------------------------
            # Y-axis ticks
            # -----------------------------------------------------

            tick_exponents = np.arange(
                lower_exp,
                upper_exp + 1,
                dtype=int,
            )

            tick_values = (
                    2.0
                    ** tick_exponents
            )

            ax.set_yticks(
                tick_values
            )

            def _format_fc_tick(value):

                if value >= 1:
                    return f"{value:g}"

                if value >= 0.1:
                    return f"{value:.2g}"

                return f"{value:.3g}"

            ax.set_yticklabels(
                [
                    _format_fc_tick(
                        float(value)
                    )
                    for value in tick_values
                ],
                fontsize=10,
            )

            # =====================================================
            # PANEL FORMATTING
            # =====================================================

            ax.set_title(
                self._display_gene_name(
                    gene
                ),
                fontsize=15,
                fontweight="bold",
                pad=10,
            )

            ax.set_xlabel(
                "Population",
                fontsize=12,
                fontweight="bold",
                labelpad=8,
            )

            ax.set_ylabel(
                "Relative mRNA expression\n"
                r"($2^{-\Delta\Delta Ct}$, log$_2$ scale)",
                fontsize=11,
                fontweight="bold",
                labelpad=10,
            )

            ax.set_xticks(
                x_positions
            )

            ax.set_xticklabels(
                locations,
                rotation=0,
                ha="center",
                fontsize=11,
            )

            ax.tick_params(
                axis="y",
                labelsize=10,
                width=1.2,
                length=4,
            )

            ax.tick_params(
                axis="x",
                width=1.2,
                length=4,
            )

            # -----------------------------------------------------
            # Grid
            # -----------------------------------------------------

            ax.grid(
                axis="y",
                which="major",
                linestyle="--",
                linewidth=0.9,
                alpha=0.35,
                zorder=0,
            )

            ax.grid(
                axis="y",
                which="minor",
                visible=False,
            )

            ax.grid(
                axis="x",
                visible=False,
            )

            # -----------------------------------------------------
            # Panel borders
            # -----------------------------------------------------

            for spine in ax.spines.values():
                spine.set_linewidth(
                    1.2
                )

            # -----------------------------------------------------
            # Panel label
            # -----------------------------------------------------

            ax.text(
                -0.07,
                1.04,
                panel_labels[
                    panel_idx
                ],
                transform=ax.transAxes,
                fontsize=16,
                fontweight="bold",
                va="top",
                ha="left",
            )

            # No significance stars are shown because none of the
            # population-versus-control comparisons remained
            # significant after Holm correction.

        # =========================================================
        # MAIN TITLE
        # =========================================================

        fig.suptitle(
            "Relative expression of sHsp genes in "
            r"$Spodoptera\ frugiperda$ populations",
            fontsize=17,
            fontweight="bold",
            y=0.975,
        )

        # =========================================================
        # FIGURE NOTE
        # =========================================================

        fig.text(
            0.5,
            0.025,
            "Open circles represent individual biological replicates "
            "(n = 2 per population); short horizontal lines indicate "
            "population-level relative expression calculated from mean ΔCt. "
            "The y-axis is shown on a log₂ scale and the dashed line at 1.0 "
            "represents the control calibrator. "
            "No population-versus-control comparison was significant after "
            "Holm correction (adjusted P > 0.05).",
            ha="center",
            va="bottom",
            fontsize=9.5,
            wrap=True,
        )

        # =========================================================
        # LAYOUT
        # =========================================================

        fig.subplots_adjust(
            left=0.15,
            right=0.97,
            bottom=0.10,
            top=0.92,
            hspace=0.62,
        )

        # =========================================================
        # SAVE
        # =========================================================

        if save:
            self._save_figure(
                fig,
                "Figure1_relative_mRNA_expression",
            )

            # Legacy filename retained for GUI compatibility.
            fig.savefig(
                self.output_dir
                / "gene_profiles.png",
                format="png",
                dpi=600,
                bbox_inches="tight",
                pad_inches=0.15,
                facecolor="white",
                edgecolor="none",
            )

        return fig

    # =================================================================
    # FIGURE 2
    # PUBLICATION HEATMAP
    # =================================================================

    def plot_heatmap(
        self,
        save: bool = True,
    ) -> plt.Figure:
        """
        Standard heatmap.
        """

        df = self._get_matrix()

        fig, ax = plt.subplots(
            figsize=(8.5, 5.5)
        )

        sns.heatmap(
            df,
            annot=True,
            fmt=".2f",
            cmap="RdYlGn",
            center=1.0,
            linewidths=0.8,
            linecolor="white",
            ax=ax,
            cbar_kws={
                "label":
                    r"Relative expression ($2^{-\Delta\Delta Ct}$)"
            },
        )

        ax.set_title(
            "sHsp Gene Expression Heatmap",
            fontsize=14,
            fontweight="bold",
        )

        ax.set_xlabel("sHsp Gene")
        ax.set_ylabel("Population")

        ax.tick_params(
            axis="x",
            rotation=0,
        )

        ax.tick_params(
            axis="y",
            rotation=0,
        )

        plt.tight_layout()

        if save:
            self._save_figure(
                fig,
                "heatmap",
            )

        return fig

    def plot_publication_heatmap(
            self,
            save: bool = True,
    ) -> plt.Figure:
        """
        Figure 2.

        Publication-ready heatmap of population-level relative-expression
        values.

        Cell annotations show the original 2^(-DeltaDeltaCt) fold-change
        values, while cell colors represent log2-transformed fold change.
        This provides a symmetric visual scale around the no-change value:

            Fold change = 1.0  ->  log2FC = 0
            Fold change > 1.0  ->  positive log2FC
            Fold change < 1.0  ->  negative log2FC

        The transformation affects visualization only and does not alter
        the reported relative-expression values.
        """

        # ---------------------------------------------------------
        # Original population-level relative-expression matrix
        # ---------------------------------------------------------

        df = self._get_matrix()

        # ---------------------------------------------------------
        # Validate fold-change values
        #
        # log2 transformation requires strictly positive values.
        # qPCR 2^(-DeltaDeltaCt) values should always be > 0.
        # ---------------------------------------------------------

        if (df <= 0).any().any():
            raise ValueError(
                "Figure 2 cannot calculate log2 fold change because "
                "one or more relative-expression values are <= 0."
            )

        # ---------------------------------------------------------
        # Transform only for the heatmap color scale
        #
        # Keep df unchanged so the annotations continue to show
        # the original 2^(-DeltaDeltaCt) values.
        # ---------------------------------------------------------

        log2_df = np.log2(
            df.astype(float)
        )

        # ---------------------------------------------------------
        # Symmetric color limits around log2FC = 0
        #
        # This ensures equivalent up- and down-regulation have
        # comparable color intensity.
        # ---------------------------------------------------------

        max_abs_log2 = float(
            np.nanmax(
                np.abs(
                    log2_df.to_numpy(
                        dtype=float
                    )
                )
            )
        )

        if not np.isfinite(max_abs_log2):
            raise ValueError(
                "Figure 2 contains invalid values after log2 transformation."
            )

        # Avoid a zero-width color scale if every value equals 1.
        if max_abs_log2 == 0:
            max_abs_log2 = 1.0

        # ---------------------------------------------------------
        # Figure
        # ---------------------------------------------------------

        fig, ax = plt.subplots(
            figsize=(9.5, 6.2)
        )

        fig.patch.set_facecolor(
            "white"
        )

        # ---------------------------------------------------------
        # Diverging palette
        #
        # Blue  = down-regulation
        # White = no change
        # Red   = up-regulation
        # ---------------------------------------------------------

        cmap = sns.diverging_palette(
            240,
            10,
            s=80,
            l=45,
            n=11,
            as_cmap=True,
        )

        # ---------------------------------------------------------
        # Heatmap
        #
        # Data used for colors:
        #     log2 fold change
        #
        # Annotation labels:
        #     original relative expression / fold change
        # ---------------------------------------------------------

        sns.heatmap(
            log2_df,
            annot=df,
            fmt=".2f",
            cmap=cmap,
            center=0.0,
            vmin=-max_abs_log2,
            vmax=max_abs_log2,
            linewidths=1.2,
            linecolor="white",
            ax=ax,
            annot_kws={
                "size": 11,
                "weight": "bold",
            },
            cbar_kws={
                "label": r"log$_2$(Fold Change)",
                "shrink": 0.82,
                "aspect": 20,
            },
        )

        # ---------------------------------------------------------
        # Title
        # ---------------------------------------------------------

        ax.set_title(
            "Relative Expression Profiles of sHsp Genes\n"
            r"in $Spodoptera\ frugiperda$ Populations",
            fontsize=14,
            fontweight="bold",
            pad=16,
        )

        # ---------------------------------------------------------
        # Axis labels
        # ---------------------------------------------------------

        ax.set_xlabel(
            "sHsp Gene",
            fontsize=11,
            fontweight="bold",
            labelpad=10,
        )

        ax.set_ylabel(
            "Population",
            fontsize=11,
            fontweight="bold",
            labelpad=10,
        )

        # ---------------------------------------------------------
        # Tick formatting
        # ---------------------------------------------------------

        ax.tick_params(
            axis="x",
            labelsize=10,
            rotation=0,
        )

        ax.tick_params(
            axis="y",
            labelsize=10,
            rotation=0,
        )

        # ---------------------------------------------------------
        # Figure note
        #
        # Important:
        # Numbers remain original 2^(-DeltaDeltaCt) values.
        # Only the color scale is log2 transformed.
        # ---------------------------------------------------------

        ax.text(
            1.0,
            -0.13,
            "Cell labels show population-level relative expression "
            r"($2^{-\Delta\Delta Ct}$); colors represent log$_2$-transformed "
            "fold change relative to the control calibrator.",
            transform=ax.transAxes,
            ha="right",
            va="top",
            fontsize=9,
            color="#333333",
            style="italic",
        )

        # ---------------------------------------------------------
        # Layout
        # ---------------------------------------------------------

        plt.tight_layout()

        # ---------------------------------------------------------
        # Save publication outputs
        # ---------------------------------------------------------

        if save:
            self._save_figure(
                fig,
                "Figure2_expression_heatmap",
            )

            # Keep previous output filename for GUI compatibility.
            fig.savefig(
                self.output_dir
                / "publication_heatmap.png",
                dpi=600,
                bbox_inches="tight",
                facecolor="white",
            )

        return fig

    # =================================================================
    # FIGURE 3
    # MSSI SCORE
    # =================================================================

    def plot_mssi_bar(
        self,
        save: bool = True,
    ) -> plt.Figure:
        """
        Figure 3.

        MSSI score by population using the corrected MSSI values.
        """

        ordered_results = self._ordered_results()

        locations = [
            r.location
            for r in ordered_results
        ]

        scores = [
            float(r.mssi_score)
            for r in ordered_results
        ]

        colors = [
            STRESS_COLORS.get(
                r.stress_level,
                "#808080",
            )
            for r in ordered_results
        ]

        fig, ax = plt.subplots(
            figsize=(9, 5.5)
        )

        fig.patch.set_facecolor("white")

        bars = ax.bar(
            locations,
            scores,
            color=colors,
            edgecolor="black",
            linewidth=0.8,
        )

        for bar, score in zip(
            bars,
            scores,
        ):
            ax.text(
                bar.get_x() +
                bar.get_width() / 2,
                bar.get_height() + 0.05,
                f"{score:.2f}",
                ha="center",
                va="bottom",
                fontsize=10,
                fontweight="bold",
            )

        # MSSI classification thresholds
        ax.axhline(
            y=0.8,
            color="gray",
            linestyle=":",
            linewidth=0.9,
        )

        ax.axhline(
            y=1.5,
            color="gray",
            linestyle=":",
            linewidth=0.9,
        )

        ax.axhline(
            y=2.5,
            color="gray",
            linestyle=":",
            linewidth=0.9,
        )

        ax.set_title(
            "Multi-HSP Stress Signature Index (MSSI) by Population",
            fontsize=14,
            fontweight="bold",
            pad=12,
        )

        ax.set_xlabel(
            "Population"
        )

        ax.set_ylabel(
            "MSSI Score"
        )

        legend_patches = [
            mpatches.Patch(
                color=STRESS_COLORS["low"],
                label="Low",
            ),
            mpatches.Patch(
                color=STRESS_COLORS["moderate"],
                label="Moderate",
            ),
            mpatches.Patch(
                color=STRESS_COLORS["high"],
                label="High",
            ),
            mpatches.Patch(
                color=STRESS_COLORS["very_high"],
                label="Critical",
            ),
        ]

        ax.legend(
            handles=legend_patches,
            title="MSSI Tier",
            loc="upper left",
            frameon=True,
        )

        ax.grid(
            axis="y",
            linestyle="--",
            alpha=0.30,
        )

        ax.grid(
            axis="x",
            visible=False,
        )

        plt.tight_layout()

        if save:
            self._save_figure(
                fig,
                "Figure3_MSSI_scores",
            )

            # Keep old filename for GUI compatibility.
            fig.savefig(
                self.output_dir /
                "mssi_bar.png",
                dpi=600,
                bbox_inches="tight",
                facecolor="white",
            )

        return fig

    # =================================================================
    # OPTIONAL CONFIDENCE GRAPH
    # Retained for GUI compatibility.
    # This is not intended as a manuscript resistance-probability figure.
    # =================================================================

    def plot_confidence(
        self,
        save: bool = True,
    ) -> plt.Figure:

        ordered_results = self._ordered_results()

        locations = [
            r.location
            for r in ordered_results
        ]

        confidences = [
            float(r.confidence) * 100
            for r in ordered_results
        ]

        colors = [
            STRESS_COLORS.get(
                r.stress_level,
                "#808080",
            )
            for r in ordered_results
        ]

        fig, ax = plt.subplots(
            figsize=(8, 4.5)
        )

        bars = ax.barh(
            locations,
            confidences,
            color=colors,
            edgecolor="black",
            linewidth=0.7,
        )

        for bar, conf in zip(
            bars,
            confidences,
        ):
            ax.text(
                bar.get_width() + 1,
                bar.get_y() +
                bar.get_height() / 2,
                f"{conf:.0f}%",
                va="center",
                fontsize=9,
            )

        ax.set_xlim(
            0,
            115,
        )

        ax.set_title(
            "Computational Classification Confidence",
            fontsize=13,
            fontweight="bold",
        )

        ax.set_xlabel(
            "Computational confidence (%)"
        )

        ax.set_ylabel(
            "Population"
        )

        plt.tight_layout()

        if save:
            self._save_figure(
                fig,
                "classification_confidence",
            )

        return fig

    # =================================================================
    # GENERATE ALL VISUALIZER OUTPUTS
    # =================================================================

    def plot_all(self):
        """
        Generate all visual outputs.

        Manuscript figures:
            Figure 1 = Relative mRNA expression
            Figure 2 = Expression heatmap
            Figure 3 = MSSI scores

        Figure 4 (geographic MSSI map) remains generated by mapper.py.
        """

        self.plot_gene_profiles()
        self.plot_publication_heatmap()
        self.plot_mssi_bar()

        # Auxiliary plots retained for pipeline/GUI compatibility.
        self.plot_heatmap()
        self.plot_confidence()