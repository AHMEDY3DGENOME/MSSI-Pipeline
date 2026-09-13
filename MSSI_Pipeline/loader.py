import pandas as pd
import numpy as np

from pathlib import Path
from typing import Dict, List


# ============================================================
# CONSTANTS
# ============================================================

LOCATION_MAP = {
    "B": "Beni Suef",
    "S": "Sohag",
    "G": "Giza",
    "Q": "Qalyubia",
}


EXPECTED_GENES = [
    "Hsp 19.74",
    "Hsp 20.7",
    "Hsp 19.07",
]


EXPECTED_LOCATIONS = [
    "Beni Suef",
    "Sohag",
    "Giza",
    "Qalyubia",
]


ALL_GROUPS = [
    "Control",
    "Beni Suef",
    "Sohag",
    "Giza",
    "Qalyubia",
]


# ============================================================
# LOADER
# ============================================================

class MSSILoader:
    """
    Load MSSI input data.

    Supported input formats
    -----------------------
    1. Corrected qPCR CSV:
       01_qPCR_DeltaCt_individual_samples.csv

       Required columns:
           Gene
           Population
           DeltaCt

       For corrected qPCR CSV input, the loader calculates:

           - mean Control DeltaCt for each gene
           - DeltaDeltaCt for each biological replicate
           - individual relative expression = 2^(-DeltaDeltaCt)
           - population fold change from mean DeltaCt

       The loader also retains the original biological-replicate
       DeltaCt values for each population. These values are required
       for downstream statistical summaries and for calculating
       fold-change uncertainty on the correct qPCR scale.

    2. Legacy Excel input:
       .xlsx / .xls

       Legacy support is retained for backward compatibility.

       Corrected manuscript analyses should use the qPCR CSV
       derived from the raw Ct data.
    """

    # ========================================================
    # INITIALIZATION
    # ========================================================

    def __init__(self, filepath: str):

        self.filepath = (
            Path(filepath)
            .expanduser()
            .resolve()
        )

        self._validate_file()

    # ========================================================
    # FILE VALIDATION
    # ========================================================

    def _validate_file(self):

        if not self.filepath.exists():

            raise FileNotFoundError(
                f"File not found: {self.filepath}"
            )

        suffix = self.filepath.suffix.lower()

        if suffix not in [
            ".xlsx",
            ".xls",
            ".csv",
        ]:

            raise ValueError(
                "Unsupported input format. "
                "Expected .csv, .xlsx, or .xls"
            )

    # ========================================================
    # CSV VALIDATION
    # ========================================================

    @staticmethod
    def _validate_qpcr_csv(df: pd.DataFrame):

        required_columns = {
            "Gene",
            "Population",
            "DeltaCt",
        }

        missing = (
            required_columns
            - set(df.columns)
        )

        if missing:

            raise ValueError(
                "Invalid qPCR CSV. "
                "Missing required columns: "
                + ", ".join(
                    sorted(missing)
                )
            )

        if df.empty:

            raise ValueError(
                "The qPCR CSV contains no data."
            )

    # ========================================================
    # NORMALIZE GENE NAMES
    # ========================================================

    @staticmethod
    def _normalize_gene_name(
        gene: str,
    ) -> str:
        """
        Normalize common gene-name variations.
        """

        text = str(gene).strip()

        aliases = {
            "hsp 19.74": "Hsp 19.74",
            "hsp19.74": "Hsp 19.74",
            "shsp 19.74": "Hsp 19.74",
            "shsp19.74": "Hsp 19.74",

            "hsp 20.7": "Hsp 20.7",
            "hsp20.7": "Hsp 20.7",
            "shsp 20.7": "Hsp 20.7",
            "shsp20.7": "Hsp 20.7",

            "hsp 19.07": "Hsp 19.07",
            "hsp19.07": "Hsp 19.07",
            "shsp 19.07": "Hsp 19.07",
            "shsp19.07": "Hsp 19.07",
        }

        return aliases.get(
            text.lower(),
            text,
        )

    # ========================================================
    # NORMALIZE POPULATION NAMES
    # ========================================================

    @staticmethod
    def _normalize_population(
        population: str,
    ) -> str:

        text = (
            str(population)
            .strip()
        )

        aliases = {
            "control": "Control",

            "b": "Beni Suef",
            "beni suef": "Beni Suef",
            "bani suef": "Beni Suef",

            "s": "Sohag",
            "sohag": "Sohag",

            "g": "Giza",
            "giza": "Giza",

            "q": "Qalyubia",
            "qalyubia": "Qalyubia",
            "qalubia": "Qalyubia",
            "qaliubiya": "Qalyubia",
        }

        return aliases.get(
            text.lower(),
            text,
        )

    # ========================================================
    # LOAD CORRECTED qPCR CSV
    # ========================================================

    def _load_qpcr_csv(
        self,
    ) -> Dict:

        df = pd.read_csv(
            self.filepath
        )

        self._validate_qpcr_csv(
            df
        )

        df = df.copy()

        # ----------------------------------------------------
        # Normalize gene names
        # ----------------------------------------------------

        df["Gene"] = (
            df["Gene"]
            .astype(str)
            .map(
                self._normalize_gene_name
            )
        )

        # ----------------------------------------------------
        # Normalize population names
        # ----------------------------------------------------

        df["Population"] = (
            df["Population"]
            .astype(str)
            .map(
                self._normalize_population
            )
        )

        # ----------------------------------------------------
        # Validate DeltaCt values
        # ----------------------------------------------------

        df["DeltaCt"] = pd.to_numeric(
            df["DeltaCt"],
            errors="coerce",
        )

        if df["DeltaCt"].isna().any():

            bad_rows = df[
                df["DeltaCt"].isna()
            ]

            raise ValueError(
                "Non-numeric or missing DeltaCt values found "
                f"in {len(bad_rows)} row(s)."
            )

        # ----------------------------------------------------
        # Validate genes
        # ----------------------------------------------------

        genes_found = list(
            dict.fromkeys(
                df["Gene"].tolist()
            )
        )

        missing_genes = [
            gene
            for gene in EXPECTED_GENES
            if gene not in genes_found
        ]

        if missing_genes:

            raise ValueError(
                "Missing expected gene(s): "
                + ", ".join(
                    missing_genes
                )
            )

        # ----------------------------------------------------
        # Validate biological replicate counts
        #
        # Current corrected manuscript dataset contains
        # exactly n = 2 biological replicates per group.
        # ----------------------------------------------------

        for gene in EXPECTED_GENES:

            gene_df = df[
                df["Gene"] == gene
            ]

            for population in ALL_GROUPS:

                n = len(
                    gene_df[
                        gene_df[
                            "Population"
                        ]
                        == population
                    ]
                )

                if n != 2:

                    raise ValueError(
                        f"{gene} / "
                        f"{population}: "
                        "expected 2 biological replicates, "
                        f"found {n}."
                    )

        # ----------------------------------------------------
        # Build pipeline data
        # ----------------------------------------------------

        data = {}

        for gene in EXPECTED_GENES:

            gene_df = df[
                df["Gene"] == gene
            ].copy()

            # ------------------------------------------------
            # Control DeltaCt biological replicates
            # ------------------------------------------------

            control_delta_ct = (
                gene_df.loc[
                    gene_df[
                        "Population"
                    ]
                    == "Control",
                    "DeltaCt",
                ]
                .to_numpy(
                    dtype=float
                )
            )

            control_mean_delta_ct = float(
                np.mean(
                    control_delta_ct
                )
            )

            # ------------------------------------------------
            # DeltaDeltaCt for every biological replicate
            # ------------------------------------------------

            gene_df[
                "DeltaDeltaCt"
            ] = (
                gene_df["DeltaCt"]
                - control_mean_delta_ct
            )

            # ------------------------------------------------
            # Individual relative expression
            #
            # Biological-replicate expression values:
            #
            # 2^(-DeltaDeltaCt)
            # ------------------------------------------------

            gene_df[
                "RelativeExpression"
            ] = (
                2.0
                ** (
                    -gene_df[
                        "DeltaDeltaCt"
                    ]
                )
            )

            # ------------------------------------------------
            # Store relative-expression biological replicates
            # ------------------------------------------------

            replicates = {}

            for population in EXPECTED_LOCATIONS:

                values = (
                    gene_df.loc[
                        gene_df[
                            "Population"
                        ]
                        == population,
                        "RelativeExpression",
                    ]
                    .to_numpy(
                        dtype=float
                    )
                )

                replicates[
                    population
                ] = values.tolist()

            # Control relative-expression replicates
            control_expression = (
                gene_df.loc[
                    gene_df[
                        "Population"
                    ]
                    == "Control",
                    "RelativeExpression",
                ]
                .to_numpy(
                    dtype=float
                )
            )

            replicates[
                "Control"
            ] = (
                control_expression
                .tolist()
            )

            # ------------------------------------------------
            # Store original DeltaCt biological replicates
            #
            # IMPORTANT:
            # Statistical inference is performed on DeltaCt,
            # not on transformed fold-change values.
            #
            # These values are also required to calculate
            # asymmetric fold-change error bars correctly.
            # ------------------------------------------------

            delta_ct_replicates = {}

            for population in ALL_GROUPS:

                values = (
                    gene_df.loc[
                        gene_df[
                            "Population"
                        ]
                        == population,
                        "DeltaCt",
                    ]
                    .to_numpy(
                        dtype=float
                    )
                )

                delta_ct_replicates[
                    population
                ] = values.tolist()

            # ------------------------------------------------
            # Population fold change
            #
            # IMPORTANT:
            #
            # Population fold change is calculated from
            # the MEAN DeltaCt:
            #
            # FC =
            # 2^[-(
            #     mean DeltaCt population
            #     -
            #     mean DeltaCt control
            # )]
            #
            # Do NOT calculate population FC as the
            # arithmetic mean of individual FC values.
            # ------------------------------------------------

            fold_change = {}

            mean_delta_ct = {}

            for population in EXPECTED_LOCATIONS:

                population_delta_ct = (
                    gene_df.loc[
                        gene_df[
                            "Population"
                        ]
                        == population,
                        "DeltaCt",
                    ]
                    .to_numpy(
                        dtype=float
                    )
                )

                population_mean_delta_ct = float(
                    np.mean(
                        population_delta_ct
                    )
                )

                mean_delta_ct[
                    population
                ] = (
                    population_mean_delta_ct
                )

                delta_delta_ct = (
                    population_mean_delta_ct
                    - control_mean_delta_ct
                )

                fold_change[
                    population
                ] = float(
                    2.0
                    ** (
                        -delta_delta_ct
                    )
                )

            # ------------------------------------------------
            # Store complete gene-level structure
            # ------------------------------------------------

            data[gene] = {

                # Biological replicate values on the
                # transformed expression scale.
                "replicates":
                    replicates,

                # Biological replicate values on the
                # original DeltaCt scale.
                "delta_ct_replicates":
                    delta_ct_replicates,

                # Population-level FC calculated from
                # mean DeltaCt.
                "fold_change":
                    fold_change,

                # Control biological DeltaCt values.
                "control_delta_ct_replicates":
                    control_delta_ct.tolist(),

                # Mean Control DeltaCt calibrator.
                "control_mean_delta_ct":
                    control_mean_delta_ct,

                # Population mean DeltaCt values.
                "mean_delta_ct":
                    mean_delta_ct,

                # Number of independent biological
                # replicates per group.
                "n_biological":
                    2,
            }

        return data

    # ========================================================
    # LEGACY EXCEL SUPPORT
    # ========================================================

    def _extract_replicates_excel(
        self,
        df: pd.DataFrame,
    ) -> Dict:

        results = {}

        header_rows = []

        for idx, row in df.iterrows():

            values = (
                row
                .dropna()
                .values
            )

            if "Control" in str(
                values
            ):

                header_rows.append(
                    idx
                )

        if not header_rows:

            raise ValueError(
                "No valid data structure "
                "found in Excel sheet."
            )

        first_header = (
            header_rows[0]
        )

        cols = [
            "Control",
            "B",
            "S",
            "G",
            "Q",
        ]

        data_rows = []

        # Legacy Excel format historically
        # contained three rows.
        #
        # This block is retained only for
        # backward compatibility.

        for idx in range(
            first_header + 1,
            first_header + 4,
        ):

            if idx in df.index:

                row = (
                    df.loc[idx]
                    .dropna()
                )

                if len(row) >= 5:

                    data_rows.append(
                        row.values[
                            :5
                        ].tolist()
                    )

        if not data_rows:

            raise ValueError(
                "No replicate data found "
                "in Excel sheet."
            )

        for i, col in enumerate(
            cols
        ):

            label = LOCATION_MAP.get(
                col,
                col,
            )

            values = []

            for row in data_rows:

                try:

                    values.append(
                        float(
                            row[i]
                        )
                    )

                except (
                    ValueError,
                    TypeError,
                    IndexError,
                ):

                    continue

            results[
                label
            ] = values

        return results

    # ========================================================
    # LEGACY FOLD CHANGE
    # ========================================================

    def _extract_fold_change_excel(
        self,
        df: pd.DataFrame,
    ) -> Dict:

        fold_row = None

        for _, row in df.iterrows():

            if "fold" in str(
                row.values
            ).lower():

                fold_row = row

                break

        if fold_row is None:

            raise ValueError(
                "No fold-change row "
                "found in Excel sheet."
            )

        numeric = []

        for value in (
            fold_row
            .dropna()
        ):

            if isinstance(
                value,
                (
                    int,
                    float,
                    np.integer,
                    np.floating,
                ),
            ):

                numeric.append(
                    float(value)
                )

        cols = [
            "B",
            "S",
            "G",
            "Q",
        ]

        if len(numeric) < 4:

            raise ValueError(
                "Fold-change row does not contain "
                "four population values."
            )

        return {
            LOCATION_MAP[col]:
                float(
                    numeric[i]
                )

            for i, col
            in enumerate(
                cols
            )
        }

    # ========================================================
    # LOAD LEGACY EXCEL
    # ========================================================

    def _load_excel(
        self,
    ) -> Dict:

        sheets = pd.read_excel(
            self.filepath,
            sheet_name=None,
            header=None,
        )

        data = {}

        for (
            gene_name,
            df,
        ) in sheets.items():

            normalized_gene = (
                self._normalize_gene_name(
                    gene_name
                )
            )

            data[
                normalized_gene
            ] = {

                "replicates":
                    self._extract_replicates_excel(
                        df
                    ),

                "fold_change":
                    self._extract_fold_change_excel(
                        df
                    ),

                # Legacy Excel input does not contain
                # the corrected biological DeltaCt
                # structure used by the current qPCR CSV.
                "delta_ct_replicates":
                    None,

                "control_delta_ct_replicates":
                    None,

                "control_mean_delta_ct":
                    None,

                "mean_delta_ct":
                    None,

                "n_biological":
                    None,
            }

        return data

    # ========================================================
    # PUBLIC LOAD METHOD
    # ========================================================

    def load(
        self,
    ) -> Dict:

        suffix = (
            self.filepath
            .suffix
            .lower()
        )

        if suffix == ".csv":

            return (
                self._load_qpcr_csv()
            )

        if suffix in [
            ".xlsx",
            ".xls",
        ]:

            return (
                self._load_excel()
            )

        raise ValueError(
            "Unsupported file format: "
            f"{suffix}"
        )

    # ========================================================
    # GET GENES
    # ========================================================

    def get_genes(
        self,
    ) -> List[str]:

        suffix = (
            self.filepath
            .suffix
            .lower()
        )

        if suffix == ".csv":

            df = pd.read_csv(
                self.filepath
            )

            self._validate_qpcr_csv(
                df
            )

            genes = [
                self._normalize_gene_name(
                    gene
                )
                for gene
                in df[
                    "Gene"
                ].unique()
            ]

            return genes

        sheets = pd.read_excel(
            self.filepath,
            sheet_name=None,
            header=None,
        )

        return [
            self._normalize_gene_name(
                name
            )
            for name
            in sheets.keys()
        ]

    # ========================================================
    # GET LOCATIONS
    # ========================================================

    def get_locations(
        self,
    ) -> List[str]:

        return (
            EXPECTED_LOCATIONS
            .copy()
        )