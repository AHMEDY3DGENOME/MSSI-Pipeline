import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import urllib.request

from pathlib import Path
from typing import List, Any


# ============================================================
# OPTIONAL GEOPANDAS SUPPORT
# ============================================================

try:
    import geopandas as gpd
    GEOPANDAS_AVAILABLE = True
except ImportError:
    GEOPANDAS_AVAILABLE = False


# ============================================================
# GEOGRAPHIC DATA
# ============================================================

SHAPEFILE_URL = (
    "https://geodata.ucdavis.edu/gadm/"
    "gadm4.1/json/gadm41_EGY_1.json.zip"
)

CACHE_DIR = Path.home() / ".mssi_pipeline"
CACHE_FILE = CACHE_DIR / "egypt_gov.zip"


# ============================================================
# GOVERNORATE NAME NORMALIZATION
# ============================================================

NAME_MAP = {
    "AlQahirah": "Cairo",
    "AlJizah": "Giza",
    "AlQalyubiyah": "Qalyubia",
    "AlIskandariyah": "Alexandria",
    "AdDaqahliyah": "Dakahlia",
    "AshSharqiyah": "Sharqia",
    "AlGharbiyah": "Gharbia",
    "AlMinufiyah": "Menoufia",
    "AlBuhayrah": "Beheira",
    "KafrashShaykh": "Kafr El Sheikh",
    "Dumyat": "Damietta",
    "BurSa`id": "Port Said",
    "AlIsma`iliyah": "Ismailia",
    "AsSuways": "Suez",
    "BaniSuwayf": "Beni Suef",
    "AlFayyum": "Faiyum",
    "AlMinya": "Minya",
    "Asyut": "Asyut",
    "Suhaj": "Sohag",
    "Qina": "Qena",
    "AlUqsur": "Luxor",
    "Aswan": "Aswan",
    "AlBahralAhmar": "Red Sea",
    "ShamalSina`": "North Sinai",
    "JanubSina`": "South Sinai",
    "Matrouh": "Matrouh",
    "AlWadialJadid": "New Valley",
}


# ============================================================
# FALLBACK GOVERNORATE COORDINATES
# ============================================================

GOVERNORATE_COORDS = {
    "Cairo": (31.2357, 30.0444),
    "Giza": (31.2089, 29.9870),
    "Qalyubia": (31.2000, 30.3300),
    "Alexandria": (29.9187, 31.2001),
    "Dakahlia": (31.4167, 31.0000),
    "Sharqia": (31.5500, 30.7333),
    "Gharbia": (31.0333, 30.8667),
    "Menoufia": (30.9833, 30.5833),
    "Beheira": (30.3480, 30.8480),
    "Kafr El Sheikh": (30.9408, 31.1108),
    "Damietta": (31.8133, 31.4167),
    "Port Said": (32.2667, 31.2500),
    "Ismailia": (32.2667, 30.6000),
    "Suez": (32.5333, 29.9667),
    "Beni Suef": (31.0833, 29.0667),
    "Faiyum": (30.8418, 29.3084),
    "Minya": (30.7500, 28.1000),
    "Asyut": (31.1667, 27.1667),
    "Sohag": (31.6948, 26.5569),
    "Qena": (32.7167, 26.1667),
    "Luxor": (32.6396, 25.6872),
    "Aswan": (32.9000, 24.0889),
    "Red Sea": (33.8333, 26.0000),
    "North Sinai": (33.7667, 30.2833),
    "South Sinai": (33.6333, 28.5000),
    "Matrouh": (27.2333, 31.3500),
    "New Valley": (28.5000, 25.4500),
}


# ============================================================
# DISPLAY COLORS
# ============================================================

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

BG_COLOR = "#FFFFFF"
DEFAULT_COLOR = "#EEEEEE"
BOUNDARY_COLOR = "#333333"


# ============================================================
# LABEL OFFSETS
# ============================================================

# These offsets move only the annotation boxes.
# They do NOT change the geographic position of governorates.

LABEL_OFFSETS = {
    "Giza": (-45, 20),
    "Beni Suef": (48, -12),
    "Qalyubia": (8, 18),
    "Sohag": (35, 0),
}


# ============================================================
# LOAD EGYPT GOVERNORATE DATA
# ============================================================

def _get_egypt_geodata():

    if not GEOPANDAS_AVAILABLE:
        raise ImportError(
            "geopandas is not installed."
        )

    # --------------------------------------------------------
    # Use cached geographic data if available
    # --------------------------------------------------------

    if CACHE_FILE.exists():

        try:

            gdf = gpd.read_file(
                CACHE_FILE
            )

            if len(gdf) > 0:

                gdf["location"] = (
                    gdf["NAME_1"]
                    .map(NAME_MAP)
                )

                return gdf

        except Exception:

            CACHE_FILE.unlink(
                missing_ok=True
            )

    # --------------------------------------------------------
    # Download geographic data if needed
    # --------------------------------------------------------

    CACHE_DIR.mkdir(
        parents=True,
        exist_ok=True
    )

    urllib.request.urlretrieve(
        SHAPEFILE_URL,
        CACHE_FILE
    )

    gdf = gpd.read_file(
        CACHE_FILE
    )

    gdf["location"] = (
        gdf["NAME_1"]
        .map(NAME_MAP)
    )

    return gdf


# ============================================================
# MSSI MAP CLASS
# ============================================================

class EgyptMapper:
    """
    Visualize MSSI molecular stress-response results across
    sampled Egyptian governorates.

    This mapper uses MSSI results directly.

    It does NOT map:
        - insecticide resistance,
        - resistance probability,
        - phenotypic susceptibility,
        - RRI values.

    Expected result objects:
        location
        mssi_score
        stress_level
        stress_label
    """

    def __init__(
        self,
        mssi_results: List[Any],
        output_dir: str = "outputs"
    ):

        if not mssi_results:
            raise ValueError(
                "No MSSI results were provided to EgyptMapper."
            )

        self.mssi_results = mssi_results

        self.output_dir = Path(
            output_dir
        )

        self.output_dir.mkdir(
            parents=True,
            exist_ok=True
        )

        required_attributes = [
            "location",
            "mssi_score",
            "stress_level",
        ]

        for result in self.mssi_results:

            for attribute in required_attributes:

                if not hasattr(
                    result,
                    attribute
                ):

                    raise ValueError(
                        "Invalid MSSI result object: "
                        f"missing '{attribute}'."
                    )

        self.mssi_map = {
            result.location: result
            for result in self.mssi_results
        }

    # ========================================================
    # PUBLIC MAP METHOD
    # ========================================================

    def plot_map(
        self,
        save: bool = True
    ) -> plt.Figure:
        """
        Plot geographic distribution of MSSI scores.

        Polygon governorate boundaries are used when
        geopandas is available. Otherwise, a coordinate-based
        fallback map is generated.
        """

        if not GEOPANDAS_AVAILABLE:

            return self._plot_map_dots(
                save
            )

        try:

            gdf = _get_egypt_geodata()

            return self._plot_map_geo(
                gdf,
                save
            )

        except Exception as exc:

            print(
                "Geographic map unavailable; "
                f"using coordinate fallback: {exc}"
            )

            return self._plot_map_dots(
                save
            )

    # ========================================================
    # POLYGON MAP
    # ========================================================

    def _plot_map_geo(
        self,
        gdf,
        save: bool
    ) -> plt.Figure:

        fig, ax = plt.subplots(
            figsize=(10, 12)
        )

        ax.set_facecolor(
            BG_COLOR
        )

        fig.patch.set_facecolor(
            BG_COLOR
        )

        # ----------------------------------------------------
        # Assign stress-tier colors
        # ----------------------------------------------------

        def get_color(location):

            if location not in self.mssi_map:
                return DEFAULT_COLOR

            result = self.mssi_map[
                location
            ]

            return STRESS_COLORS.get(
                result.stress_level,
                DEFAULT_COLOR
            )

        gdf = gdf.copy()

        gdf["color"] = (
            gdf["location"]
            .apply(get_color)
        )

        gdf.plot(
            ax=ax,
            color=gdf["color"],
            edgecolor=BOUNDARY_COLOR,
            linewidth=0.6
        )

        # ----------------------------------------------------
        # Label sampled governorates
        # ----------------------------------------------------

        for _, row in gdf.iterrows():

            location = row[
                "location"
            ]

            if (
                location
                not in self.mssi_map
            ):
                continue

            result = self.mssi_map[
                location
            ]

            centroid = (
                row.geometry.centroid
            )

            color = STRESS_COLORS.get(
                result.stress_level,
                "#666666"
            )

            offset = LABEL_OFFSETS.get(
                location,
                (0, 0)
            )

            annotation_kwargs = {
                "text": (
                    f"{location}\n"
                    f"MSSI = {result.mssi_score:.2f}"
                ),
                "xy": (
                    centroid.x,
                    centroid.y
                ),
                "xytext": offset,
                "textcoords": "offset points",
                "ha": "center",
                "va": "center",
                "fontsize": 8,
                "fontweight": "bold",
                "color": "white",
                "bbox": dict(
                    boxstyle="round,pad=0.25",
                    fc=color,
                    alpha=0.90,
                    ec="none"
                )
            }

            # Add connector line if label is moved.
            if offset != (0, 0):

                annotation_kwargs[
                    "arrowprops"
                ] = dict(
                    arrowstyle="-",
                    color=color,
                    linewidth=0.9,
                    alpha=0.85
                )

            ax.annotate(
                **annotation_kwargs
            )

        self._add_legend(ax)
        self._style_ax(ax)

        plt.tight_layout()

        if save:

            fig.savefig(
                self.output_dir
                / "egypt_mssi_map.png",
                dpi=300,
                bbox_inches="tight",
                facecolor=BG_COLOR
            )

        return fig

    # ========================================================
    # FALLBACK COORDINATE MAP
    # ========================================================

    def _plot_map_dots(
        self,
        save: bool
    ) -> plt.Figure:

        fig, ax = plt.subplots(
            figsize=(10, 12)
        )

        ax.set_facecolor(
            BG_COLOR
        )

        fig.patch.set_facecolor(
            BG_COLOR
        )

        ax.set_xlim(
            24,
            38
        )

        ax.set_ylim(
            21,
            32
        )

        scores = [
            float(
                result.mssi_score
            )
            for result
            in self.mssi_results
        ]

        max_score = max(
            scores
        )

        if max_score <= 0:
            max_score = 1.0

        for governorate, (
            lon,
            lat
        ) in GOVERNORATE_COORDS.items():

            if (
                governorate
                in self.mssi_map
            ):

                result = self.mssi_map[
                    governorate
                ]

                score = float(
                    result.mssi_score
                )

                color = STRESS_COLORS.get(
                    result.stress_level,
                    "#666666"
                )

                normalized_score = (
                    score / max_score
                )

                size = (
                    300
                    + normalized_score
                    * 500
                )

                ax.scatter(
                    lon,
                    lat,
                    c=color,
                    s=size,
                    zorder=5,
                    alpha=0.9,
                    edgecolors="#333333",
                    linewidths=1.5
                )

                offset = LABEL_OFFSETS.get(
                    governorate,
                    (8, 4)
                )

                ax.annotate(
                    (
                        f"{governorate}\n"
                        f"MSSI = {score:.2f}"
                    ),
                    xy=(
                        lon,
                        lat
                    ),
                    xytext=offset,
                    textcoords="offset points",
                    fontsize=9,
                    fontweight="bold",
                    color=color,
                    arrowprops=dict(
                        arrowstyle="-",
                        color=color,
                        linewidth=0.8,
                        alpha=0.8
                    )
                )

            else:

                ax.scatter(
                    lon,
                    lat,
                    c="#CCCCCC",
                    s=80,
                    zorder=3,
                    alpha=0.6
                )

                ax.annotate(
                    governorate,
                    (
                        lon,
                        lat
                    ),
                    textcoords="offset points",
                    xytext=(5, 3),
                    fontsize=7,
                    color="#666666"
                )

        self._add_legend(ax)
        self._style_ax(ax)

        ax.text(
            0.98,
            0.02,
            "Larger circle = higher MSSI score",
            transform=ax.transAxes,
            color="#666666",
            fontsize=8,
            ha="right"
        )

        plt.tight_layout()

        if save:

            fig.savefig(
                self.output_dir
                / "egypt_mssi_map.png",
                dpi=300,
                bbox_inches="tight",
                facecolor=BG_COLOR
            )

        return fig

    # ========================================================
    # LEGEND
    # ========================================================

    def _add_legend(
        self,
        ax
    ):

        order = [
            "low",
            "moderate",
            "high",
            "very_high",
        ]

        legend_patches = [

            mpatches.Patch(
                color=STRESS_COLORS[
                    level
                ],
                label=STRESS_LABELS[
                    level
                ]
            )

            for level
            in order
        ]

        legend = ax.legend(
            handles=legend_patches,
            loc="lower left",
            facecolor="white",
            edgecolor="#333333",
            labelcolor="black",
            fontsize=10,
            title="MSSI Stress Tier",
            title_fontsize=11
        )

        legend.get_title().set_color(
            "black"
        )

    # ========================================================
    # MAP STYLE
    # ========================================================

    def _style_ax(
        self,
        ax
    ):

        ax.set_title(
            (
                "Geospatial Distribution of "
                "Molecular Stress-Response Signatures\n"
                "Multi-HSP Stress Signature Index (MSSI)"
            ),
            color="black",
            fontsize=14,
            fontweight="bold",
            pad=15
        )

        ax.tick_params(
            colors="black"
        )

        for spine in (
            ax.spines.values()
        ):

            spine.set_color(
                "#333333"
            )

        ax.set_xlabel(
            "Longitude",
            color="black",
            fontsize=11
        )

        ax.set_ylabel(
            "Latitude",
            color="black",
            fontsize=11
        )

    # ========================================================
    # MSSI BAR PLOT
    # ========================================================

    def plot_mssi_bar(
        self,
        save: bool = True
    ) -> plt.Figure:

        fig, ax = plt.subplots(
            figsize=(9, 5)
        )

        locations = [
            result.location
            for result
            in self.mssi_results
        ]

        scores = [
            float(
                result.mssi_score
            )
            for result
            in self.mssi_results
        ]

        colors = [
            STRESS_COLORS.get(
                result.stress_level,
                "#666666"
            )
            for result
            in self.mssi_results
        ]

        bars = ax.bar(
            locations,
            scores,
            color=colors,
            edgecolor="#333333",
            linewidth=0.8
        )

        for bar, score in zip(
            bars,
            scores
        ):

            ax.text(
                (
                    bar.get_x()
                    + bar.get_width()
                    / 2
                ),
                (
                    bar.get_height()
                    + 0.03
                ),
                f"{score:.3f}",
                ha="center",
                va="bottom",
                fontsize=10,
                fontweight="bold"
            )

        # ----------------------------------------------------
        # MSSI classification thresholds
        # ----------------------------------------------------

        ax.axhline(
            y=0.8,
            linestyle="--",
            linewidth=1,
            label=(
                "Moderate threshold "
                "(MSSI = 0.8)"
            )
        )

        ax.axhline(
            y=1.5,
            linestyle="--",
            linewidth=1,
            label=(
                "High threshold "
                "(MSSI = 1.5)"
            )
        )

        ax.axhline(
            y=2.5,
            linestyle="--",
            linewidth=1,
            label=(
                "Critical threshold "
                "(MSSI = 2.5)"
            )
        )

        ymax = max(
            max(scores) * 1.20,
            3.0
        )

        ax.set_ylim(
            0,
            ymax
        )

        ax.set_title(
            (
                "Multi-HSP Stress Signature "
                "Index by Population"
            ),
            fontsize=14,
            fontweight="bold"
        )

        ax.set_xlabel(
            "Population",
            fontsize=11
        )

        ax.set_ylabel(
            "MSSI Score",
            fontsize=11
        )

        ax.legend(
            fontsize=8
        )

        plt.tight_layout()

        if save:

            fig.savefig(
                self.output_dir
                / "mssi_bar.png",
                dpi=300,
                bbox_inches="tight"
            )

        return fig

    # ========================================================
    # BACKWARD-COMPATIBLE METHOD
    # ========================================================

    def plot_risk_bar(
        self,
        save: bool = True
    ) -> plt.Figure:
        """
        Backward-compatible alias.

        Historical GUI versions may still call
        plot_risk_bar(). It now returns the MSSI bar plot.
        """

        return self.plot_mssi_bar(
            save=save
        )

    # ========================================================
    # GENERATE ALL OUTPUTS
    # ========================================================

    def plot_all(self):

        self.plot_map(
            save=True
        )

        self.plot_mssi_bar(
            save=True
        )