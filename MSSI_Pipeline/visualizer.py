import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import seaborn as sns
import numpy as np
import pandas as pd
from typing import Dict, List
from pathlib import Path
from .classifier import StressResult

STRESS_COLORS = {
    "very_high": "#E24B4A",
    "high": "#EF9F27",
    "moderate": "#378ADD",
    "low": "#1D9E75",
}


class MSSIVisualizer:

    def __init__(self, results: List[StressResult], output_dir: str = "outputs"):
        self.results = results
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(exist_ok=True)
        sns.set_theme(style="whitegrid", font_scale=1.1)

    def _get_colors(self) -> List[str]:
        return [STRESS_COLORS[r.stress_level] for r in self.results]

    def _get_matrix(self) -> pd.DataFrame:
        genes = list(self.results[0].gene_profile.keys())
        locations = [r.location for r in self.results]
        matrix = np.array([
            [r.gene_profile[gene]["fold_change"] for gene in genes]
            for r in self.results
        ])
        return pd.DataFrame(matrix, index=locations, columns=genes)

    def plot_mssi_bar(self, save: bool = True) -> plt.Figure:
        fig, ax = plt.subplots(figsize=(9, 5))
        locations = [r.location for r in self.results]
        scores = [r.mssi_score for r in self.results]
        colors = self._get_colors()

        bars = ax.bar(locations, scores, color=colors,
                      edgecolor="white", linewidth=0.8)
        for bar, score in zip(bars, scores):
            ax.text(
                bar.get_x() + bar.get_width() / 2,
                bar.get_height() + 0.05,
                f"{score:.2f}",
                ha="center", va="bottom",
                fontsize=10, fontweight="bold"
            )
        ax.axhline(y=1.0, color="gray", linestyle="--",
                   linewidth=1, label="Baseline (1.0)")
        ax.set_title("MSSI Score by Location", fontsize=14, fontweight="bold")
        ax.set_xlabel("Location")
        ax.set_ylabel("MSSI Score")
        legend_patches = [
            mpatches.Patch(color=color, label=level.replace("_", " ").title())
            for level, color in STRESS_COLORS.items()
        ]
        ax.legend(handles=legend_patches, loc="upper left", fontsize=9)
        plt.tight_layout()
        if save:
            fig.savefig(self.output_dir / "mssi_bar.png", dpi=300)
        return fig

    def plot_heatmap(self, save: bool = True) -> plt.Figure:
        df = self._get_matrix()
        fig, ax = plt.subplots(figsize=(8, 5))
        sns.heatmap(
            df, annot=True, fmt=".2f", cmap="RdYlGn",
            center=1.0, linewidths=0.5, ax=ax,
            cbar_kws={"label": "Fold Change"}
        )
        ax.set_title("sHsp Gene Expression Heatmap",
                     fontsize=14, fontweight="bold")
        ax.set_xlabel("Gene")
        ax.set_ylabel("Location")
        plt.tight_layout()
        if save:
            fig.savefig(self.output_dir / "heatmap.png", dpi=300)
        return fig

    def plot_publication_heatmap(self, save: bool = True) -> plt.Figure:
        df = self._get_matrix()

        fig, ax = plt.subplots(figsize=(10, 7))
        fig.patch.set_facecolor("white")

        # High-contrast diverging palette for academic publishing
        cmap = sns.diverging_palette(240, 10, s=80, l=40, n=9, as_cmap=True)

        hm = sns.heatmap(
            df,
            annot=True,
            fmt=".2f",
            cmap=cmap,
            center=1.0,
            robust=True,
            linewidths=1.5,
            linecolor="white",
            ax=ax,
            annot_kws={"size": 12, "weight": "bold"},
            cbar_kws={
                "label": "Fold Change (2^-ΔΔCt)",
                "shrink": 0.8,
                "aspect": 20
            }
        )

        ax.set_title(
            "sHsp Gene Expression Profile across Egyptian Governorates\n"
            "Spodoptera frugiperda (Fall Armyworm)",
            fontsize=14,
            fontweight="bold",
            pad=20,
            color="black"
        )

        ax.set_xlabel("sHsp Gene", fontsize=12, fontweight="bold", labelpad=12)
        ax.set_ylabel("Governorate", fontsize=12, fontweight="bold", labelpad=12)

        ax.tick_params(axis="x", labelsize=11, rotation=0)
        ax.tick_params(axis="y", labelsize=11, rotation=0)

        # Footnote with statistical context
        ax.text(
            1.0, -0.15,
            "Values represent mean fold change relative to control. "
            "Red = Up-regulated (>1.0), Blue = Down-regulated (<1.0).",
            transform=ax.transAxes,
            ha="right",
            va="top",
            fontsize=9,
            color="#333333",
            style="italic"
        )

        plt.tight_layout()

        if save:
            output_path = self.output_dir / "publication_heatmap.png"
            fig.savefig(
                output_path,
                dpi=300,
                bbox_inches="tight",
                facecolor="white"
            )

        return fig

    def plot_gene_profiles(self, save: bool = True) -> plt.Figure:
        genes = list(self.results[0].gene_profile.keys())
        locations = [r.location for r in self.results]
        x = np.arange(len(locations))
        width = 0.25

        fig, ax = plt.subplots(figsize=(10, 5))
        gene_colors = ["#534AB7", "#1D9E75", "#D85A30"]

        for i, (gene, color) in enumerate(zip(genes, gene_colors)):
            values = [r.gene_profile[gene]["fold_change"] for r in self.results]
            errors = [r.gene_profile[gene]["se"] for r in self.results]
            ax.bar(x + i * width, values, width,
                   label=gene, color=color,
                   yerr=errors, capsize=4,
                   error_kw={"elinewidth": 1.2})

        ax.axhline(y=1.0, color="gray", linestyle="--", linewidth=1)
        ax.set_title("sHsp Gene Expression by Location",
                     fontsize=14, fontweight="bold")
        ax.set_xlabel("Location")
        ax.set_ylabel("Fold Change")
        ax.set_xticks(x + width)
        ax.set_xticklabels(locations)
        ax.legend()
        plt.tight_layout()
        if save:
            fig.savefig(self.output_dir / "gene_profiles.png", dpi=300)
        return fig

    def plot_confidence(self, save: bool = True) -> plt.Figure:
        locations = [r.location for r in self.results]
        confidences = [r.confidence * 100 for r in self.results]
        colors = self._get_colors()

        fig, ax = plt.subplots(figsize=(8, 4))
        bars = ax.barh(locations, confidences, color=colors, edgecolor="white")
        for bar, conf in zip(bars, confidences):
            ax.text(
                bar.get_width() + 1,
                bar.get_y() + bar.get_height() / 2,
                f"{conf:.0f}%", va="center", fontsize=10
            )
        ax.set_xlim(0, 115)
        ax.set_title("Classification Confidence by Location",
                     fontsize=14, fontweight="bold")
        ax.set_xlabel("Confidence (%)")
        plt.tight_layout()
        if save:
            fig.savefig(self.output_dir / "confidence.png", dpi=300)
        return fig

    def plot_all(self):
        self.plot_mssi_bar()
        self.plot_heatmap()
        self.plot_gene_profiles()
        self.plot_confidence()
        self.plot_publication_heatmap()
