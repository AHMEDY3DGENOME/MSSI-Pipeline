import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import urllib.request
from pathlib import Path
from typing import List
from .risk_index import RiskResult, RISK_COLORS

try:
    import geopandas as gpd
    GEOPANDAS_AVAILABLE = True
except ImportError:
    GEOPANDAS_AVAILABLE = False

SHAPEFILE_URL = (
    "https://geodata.ucdavis.edu/gadm/gadm4.1/json/gadm41_EGY_1.json.zip"
)
CACHE_DIR  = Path.home() / ".mssi_pipeline"
CACHE_FILE = CACHE_DIR / "egypt_gov.zip"

NAME_MAP = {
    "AlQahirah":      "Cairo",
    "AlJizah":        "Giza",
    "AlQalyubiyah":   "Qalyubia",
    "AlIskandariyah": "Alexandria",
    "AdDaqahliyah":   "Dakahlia",
    "AshSharqiyah":   "Sharqia",
    "AlGharbiyah":    "Gharbia",
    "AlMinufiyah":    "Menoufia",
    "AlBuhayrah":     "Beheira",
    "KafrashShaykh":  "Kafr El Sheikh",
    "Dumyat":         "Damietta",
    "BurSa`id":       "Port Said",
    "AlIsma`iliyah":  "Ismailia",
    "AsSuways":       "Suez",
    "BaniSuwayf":     "Beni Suef",
    "AlFayyum":       "Faiyum",
    "AlMinya":        "Minya",
    "Asyut":          "Asyut",
    "Suhaj":          "Sohag",
    "Qina":           "Qena",
    "AlUqsur":        "Luxor",
    "Aswan":          "Aswan",
    "AlBahralAhmar":  "Red Sea",
    "ShamalSina`":    "North Sinai",
    "JanubSina`":     "South Sinai",
    "Matrouh":        "Matrouh",
    "AlWadialJadid":  "New Valley",
}

GOVERNORATE_COORDS = {
    "Cairo":          (31.2357, 30.0444),
    "Giza":           (31.2089, 29.9870),
    "Qalyubia":       (31.2000, 30.3300),
    "Alexandria":     (29.9187, 31.2001),
    "Dakahlia":       (31.4167, 31.0000),
    "Sharqia":        (31.5500, 30.7333),
    "Gharbia":        (31.0333, 30.8667),
    "Menoufia":       (30.9833, 30.5833),
    "Beheira":        (30.3480, 30.8480),
    "Kafr El Sheikh": (30.9408, 31.1108),
    "Damietta":       (31.8133, 31.4167),
    "Port Said":      (32.2667, 31.2500),
    "Ismailia":       (32.2667, 30.6000),
    "Suez":           (32.5333, 29.9667),
    "Beni Suef":      (31.0833, 29.0667),
    "Faiyum":         (30.8418, 29.3084),
    "Minya":          (30.7500, 28.1000),
    "Asyut":          (31.1667, 27.1667),
    "Sohag":          (31.6948, 26.5569),
    "Qena":           (32.7167, 26.1667),
    "Luxor":          (32.6396, 25.6872),
    "Aswan":          (32.9000, 24.0889),
    "Red Sea":        (33.8333, 26.0000),
    "North Sinai":    (33.7667, 30.2833),
    "South Sinai":    (33.6333, 28.5000),
    "Matrouh":        (27.2333, 31.3500),
    "New Valley":     (28.5000, 25.4500),
}

BG_COLOR       = "#FFFFFF"
DEFAULT_COLOR  = "#EEEEEE"
BOUNDARY_COLOR = "#333333"


def _get_egypt_geodata():
    if not GEOPANDAS_AVAILABLE:
        raise ImportError("geopandas not installed.")
    if CACHE_FILE.exists():
        try:
            gdf = gpd.read_file(CACHE_FILE)
            if len(gdf) > 0:
                gdf["location"] = gdf["NAME_1"].map(NAME_MAP)
                return gdf
        except Exception:
            CACHE_FILE.unlink(missing_ok=True)
    CACHE_DIR.mkdir(parents=True, exist_ok=True)
    urllib.request.urlretrieve(SHAPEFILE_URL, CACHE_FILE)
    gdf = gpd.read_file(CACHE_FILE)
    gdf["location"] = gdf["NAME_1"].map(NAME_MAP)
    return gdf


class EgyptMapper:

    def __init__(self, risk_results: List[RiskResult], output_dir: str = "outputs"):
        self.risk_results = risk_results
        self.output_dir   = Path(output_dir)
        self.output_dir.mkdir(exist_ok=True)
        self.risk_map     = {r.location: r for r in risk_results}

    def plot_map(self, save: bool = True) -> plt.Figure:
        if not GEOPANDAS_AVAILABLE:
            return self._plot_map_dots(save)
        try:
            gdf = _get_egypt_geodata()
            return self._plot_map_geo(gdf, save)
        except Exception:
            return self._plot_map_dots(save)

    def _plot_map_geo(self, gdf, save: bool) -> plt.Figure:
        fig, ax = plt.subplots(figsize=(10, 12))
        ax.set_facecolor(BG_COLOR)
        fig.patch.set_facecolor(BG_COLOR)

        gdf["color"] = gdf["location"].apply(
            lambda loc: self.risk_map[loc].risk_color
            if loc in self.risk_map else DEFAULT_COLOR
        )
        gdf.plot(ax=ax, color=gdf["color"],
                 edgecolor=BOUNDARY_COLOR, linewidth=0.6)

        for _, row in gdf.iterrows():
            loc = row["location"]
            if loc in self.risk_map:
                r        = self.risk_map[loc]
                centroid = row.geometry.centroid
                ax.annotate(
                    f"{loc}\n{r.risk_index:.2f}",
                    xy=(centroid.x, centroid.y),
                    ha="center", va="center",
                    fontsize=8, fontweight="bold",
                    color="white",
                    bbox=dict(boxstyle="round,pad=0.2",
                              fc=r.risk_color, alpha=0.85, ec="none")
                )

        self._add_legend(ax)
        self._style_ax(ax)
        plt.tight_layout()
        if save:
            fig.savefig(self.output_dir / "egypt_risk_map.png",
                        dpi=300, bbox_inches="tight",
                        facecolor=BG_COLOR)
        return fig

    def _plot_map_dots(self, save: bool) -> plt.Figure:
        fig, ax = plt.subplots(figsize=(10, 12))
        ax.set_facecolor(BG_COLOR)
        fig.patch.set_facecolor(BG_COLOR)
        ax.set_xlim(24, 38)
        ax.set_ylim(21, 32)

        for gov, (lon, lat) in GOVERNORATE_COORDS.items():
            if gov in self.risk_map:
                r    = self.risk_map[gov]
                size = 300 + r.risk_index * 500
                ax.scatter(lon, lat, c=r.risk_color, s=size,
                           zorder=5, alpha=0.9,
                           edgecolors="#333333", linewidths=1.5)
                ax.annotate(
                    f"{gov}\n{r.risk_index:.2f}",
                    (lon, lat),
                    textcoords="offset points",
                    xytext=(8, 4),
                    fontsize=9, fontweight="bold",
                    color=r.risk_color
                )
            else:
                ax.scatter(lon, lat, c="#CCCCCC",
                           s=80, zorder=3, alpha=0.6)
                ax.annotate(gov, (lon, lat),
                            textcoords="offset points",
                            xytext=(5, 3),
                            fontsize=7, color="#666666")

        self._add_legend(ax)
        self._style_ax(ax)
        ax.text(0.98, 0.02, "Larger circle = Higher risk",
                transform=ax.transAxes, color="#666666",
                fontsize=8, ha="right")
        plt.tight_layout()
        if save:
            fig.savefig(self.output_dir / "egypt_risk_map.png",
                        dpi=300, bbox_inches="tight",
                        facecolor=BG_COLOR)
        return fig

    def _add_legend(self, ax):
        legend_patches = [
            mpatches.Patch(color=RISK_COLORS[level],
                           label=level.replace("_", " ").title())
            for level in RISK_COLORS
        ]
        legend = ax.legend(
            handles=legend_patches, loc="lower left",
            facecolor="white", edgecolor="#333333",
            labelcolor="black", fontsize=10,
            title="Risk Level", title_fontsize=11,
        )
        legend.get_title().set_color("black")

    def _style_ax(self, ax):
        ax.set_title(
            "FAW Resistance Risk Map - Egypt\nMSSI Pipeline v1.0.0",
            color="black", fontsize=14,
            fontweight="bold", pad=15
        )
        ax.tick_params(colors="black")
        for spine in ax.spines.values():
            spine.set_color("#333333")
        ax.set_xlabel("Longitude", color="black", fontsize=11)
        ax.set_ylabel("Latitude",  color="black", fontsize=11)

    def plot_risk_bar(self, save: bool = True) -> plt.Figure:
        fig, ax = plt.subplots(figsize=(9, 5))
        locations = [r.location   for r in self.risk_results]
        indices   = [r.risk_index for r in self.risk_results]
        colors    = [r.risk_color for r in self.risk_results]

        bars = ax.bar(locations, indices, color=colors,
                      edgecolor="#333333", linewidth=0.8)
        for bar, idx in zip(bars, indices):
            ax.text(
                bar.get_x() + bar.get_width() / 2,
                bar.get_height() + 0.01,
                f"{idx:.3f}",
                ha="center", va="bottom",
                fontsize=10, fontweight="bold"
            )
        ax.axhline(y=0.6, color="#FC8D59", linestyle="--",
                   linewidth=1, label="High risk threshold")
        ax.axhline(y=0.8, color="#D73027", linestyle="--",
                   linewidth=1, label="Critical risk threshold")
        ax.set_ylim(0, 1.1)
        ax.set_title("Resistance Risk Index by Location",
                     fontsize=14, fontweight="bold")
        ax.set_xlabel("Location", fontsize=11)
        ax.set_ylabel("Risk Index (0-1)", fontsize=11)
        ax.legend(fontsize=9)
        plt.tight_layout()
        if save:
            fig.savefig(self.output_dir / "risk_bar.png",
                        dpi=300)
        return fig

    def plot_all(self):
        self.plot_map()
        self.plot_risk_bar()