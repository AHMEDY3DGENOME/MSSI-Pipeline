import numpy as np
from typing import Dict, List, Optional, Tuple


class MSSICalculator:

    def __init__(
        self,
        data: Dict,
        weights: Optional[Dict] = None
    ):
        if not data:
            raise ValueError("No input data supplied to MSSICalculator.")

        self.data = data
        self.genes = list(data.keys())

        if not self.genes:
            raise ValueError("No genes found in input data.")

        first_gene = self.genes[0]

        if "fold_change" not in data[first_gene]:
            raise ValueError(
                "Input data is missing the 'fold_change' structure."
            )

        self.locations = list(
            data[first_gene]["fold_change"].keys()
        )

        self.weights = (
            weights
            if weights is not None
            else {
                gene: 1.0
                for gene in self.genes
            }
        )

        for gene in self.genes:
            if gene not in self.weights:
                raise ValueError(
                    f"Missing weight for gene: {gene}"
                )

    # ============================================================
    # BASIC DATA ACCESS
    # ============================================================

    def _get_fold_changes(
        self,
        location: str
    ) -> List[float]:

        return [
            float(
                self.data[gene]["fold_change"][location]
            )
            for gene in self.genes
        ]

    def _get_replicates(
        self,
        gene: str,
        location: str
    ) -> List[float]:
        """
        Return individual relative-expression biological replicates.

        For corrected qPCR CSV input, these are:
            2^(-DeltaDeltaCt)

        This method is retained for backward compatibility.
        """

        return [
            float(value)
            for value in
            self.data[gene]["replicates"][location]
        ]

    def _get_delta_ct_replicates(
        self,
        gene: str,
        location: str
    ) -> Optional[List[float]]:
        """
        Return biological-replicate DeltaCt values.

        Available for corrected qPCR CSV input.
        Legacy Excel input may not contain them.
        """

        delta_ct_data = self.data[gene].get(
            "delta_ct_replicates"
        )

        if not delta_ct_data:
            return None

        if location not in delta_ct_data:
            return None

        values = delta_ct_data[location]

        if values is None:
            return None

        return [
            float(value)
            for value in values
        ]

    # ============================================================
    # MSSI
    # ============================================================

    def compute_mssi(
        self,
        location: str
    ) -> float:

        fold_changes = self._get_fold_changes(
            location
        )

        weighted_sum = sum(
            self.weights[gene]
            * fold_changes[i]
            for i, gene in enumerate(
                self.genes
            )
        )

        return round(
            weighted_sum
            / len(self.genes),
            4
        )

    def compute_all(
        self
    ) -> Dict:

        return {
            location:
                self.compute_mssi(
                    location
                )
            for location
            in self.locations
        }

    # ============================================================
    # LEGACY RELATIVE-EXPRESSION SUMMARY
    # ============================================================

    def compute_mean(
        self,
        gene: str,
        location: str
    ) -> float:
        """
        Arithmetic mean of individual relative-expression values.

        Retained for compatibility with existing pipeline code.

        NOTE:
        This value is not used as the primary population fold change
        in the corrected qPCR workflow.
        """

        replicates = self._get_replicates(
            gene,
            location
        )

        return round(
            float(
                np.mean(replicates)
            ),
            4
        )

    def compute_sd(
        self,
        gene: str,
        location: str
    ) -> float:
        """
        SD of individual relative-expression values.

        Retained for backward compatibility only.

        Figure 1 should use the asymmetric fold-change bounds
        calculated from DeltaCt values instead.
        """

        replicates = self._get_replicates(
            gene,
            location
        )

        if len(replicates) < 2:
            return 0.0

        return round(
            float(
                np.std(
                    replicates,
                    ddof=1
                )
            ),
            4
        )

    def compute_se(
        self,
        gene: str,
        location: str
    ) -> float:
        """
        SEM of individual relative-expression values.

        Retained for compatibility only.
        """

        replicates = self._get_replicates(
            gene,
            location
        )

        if len(replicates) < 2:
            return 0.0

        sd = float(
            np.std(
                replicates,
                ddof=1
            )
        )

        return round(
            sd
            / np.sqrt(
                len(replicates)
            ),
            4
        )

    # ============================================================
    # DeltaCt SUMMARY
    # ============================================================

    def compute_delta_ct_mean(
        self,
        gene: str,
        location: str
    ) -> Optional[float]:

        replicates = (
            self._get_delta_ct_replicates(
                gene,
                location
            )
        )

        if not replicates:
            return None

        return float(
            np.mean(
                replicates
            )
        )

    def compute_delta_ct_sd(
        self,
        gene: str,
        location: str
    ) -> Optional[float]:

        replicates = (
            self._get_delta_ct_replicates(
                gene,
                location
            )
        )

        if not replicates:
            return None

        if len(replicates) < 2:
            return 0.0

        return float(
            np.std(
                replicates,
                ddof=1
            )
        )

    def compute_delta_ct_se(
        self,
        gene: str,
        location: str
    ) -> Optional[float]:

        replicates = (
            self._get_delta_ct_replicates(
                gene,
                location
            )
        )

        if not replicates:
            return None

        if len(replicates) < 2:
            return 0.0

        sd = float(
            np.std(
                replicates,
                ddof=1
            )
        )

        return (
            sd
            / np.sqrt(
                len(replicates)
            )
        )

    # ============================================================
    # CORRECT FOLD-CHANGE ERROR BOUNDS
    # ============================================================

    def compute_fold_change_bounds(
        self,
        gene: str,
        location: str,
        error_type: str = "sd"
    ) -> Tuple[float, float]:
        """
        Calculate asymmetric fold-change bounds from DeltaCt values.

        Population fold change is:

            FC =
            2^[-(
                mean DeltaCt_population
                -
                mean DeltaCt_control
            )]

        Variability is calculated on the DeltaCt scale and then
        transformed to the fold-change scale.

        For SD:

            lower =
            2^[-(
                (mean DeltaCt_population + SD)
                -
                mean DeltaCt_control
            )]

            upper =
            2^[-(
                (mean DeltaCt_population - SD)
                -
                mean DeltaCt_control
            )]

        For SEM, SD is replaced by SEM.

        Returns
        -------
        tuple
            (lower_bound, upper_bound)

        For legacy Excel input without DeltaCt metadata,
        symmetric relative-expression SD/SE fallback bounds are used,
        with the lower bound constrained to zero.
        """

        fold_change = float(
            self.data[gene][
                "fold_change"
            ][location]
        )

        delta_ct_replicates = (
            self._get_delta_ct_replicates(
                gene,
                location
            )
        )

        control_mean_delta_ct = (
            self.data[gene].get(
                "control_mean_delta_ct"
            )
        )

        # --------------------------------------------------------
        # Corrected qPCR pathway
        # --------------------------------------------------------

        if (
            delta_ct_replicates
            and
            control_mean_delta_ct
            is not None
        ):

            mean_delta_ct = float(
                np.mean(
                    delta_ct_replicates
                )
            )

            if len(
                delta_ct_replicates
            ) < 2:
                spread = 0.0

            else:
                sd_delta_ct = float(
                    np.std(
                        delta_ct_replicates,
                        ddof=1
                    )
                )

                if error_type.lower() in (
                    "se",
                    "sem"
                ):
                    spread = (
                        sd_delta_ct
                        / np.sqrt(
                            len(
                                delta_ct_replicates
                            )
                        )
                    )

                else:
                    spread = (
                        sd_delta_ct
                    )

            lower_delta_ct = (
                mean_delta_ct
                + spread
            )

            upper_delta_ct = (
                mean_delta_ct
                - spread
            )

            lower_bound = float(
                2.0
                ** (
                    -(
                        lower_delta_ct
                        - float(
                            control_mean_delta_ct
                        )
                    )
                )
            )

            upper_bound = float(
                2.0
                ** (
                    -(
                        upper_delta_ct
                        - float(
                            control_mean_delta_ct
                        )
                    )
                )
            )

            return (
                lower_bound,
                upper_bound
            )

        # --------------------------------------------------------
        # Legacy fallback
        # --------------------------------------------------------

        if error_type.lower() in (
            "se",
            "sem"
        ):
            error = self.compute_se(
                gene,
                location
            )
        else:
            error = self.compute_sd(
                gene,
                location
            )

        lower_bound = max(
            0.0,
            fold_change - error
        )

        upper_bound = (
            fold_change + error
        )

        return (
            float(
                lower_bound
            ),
            float(
                upper_bound
            )
        )

    def compute_fold_change_errors(
        self,
        gene: str,
        location: str,
        error_type: str = "sd"
    ) -> Tuple[float, float]:
        """
        Return asymmetric error lengths required by matplotlib.

        Returns
        -------
        tuple
            (lower_error, upper_error)
        """

        fold_change = float(
            self.data[gene][
                "fold_change"
            ][location]
        )

        lower_bound, upper_bound = (
            self.compute_fold_change_bounds(
                gene,
                location,
                error_type=error_type
            )
        )

        lower_error = max(
            0.0,
            fold_change
            - lower_bound
        )

        upper_error = max(
            0.0,
            upper_bound
            - fold_change
        )

        return (
            float(
                lower_error
            ),
            float(
                upper_error
            )
        )

    # ============================================================
    # FULL STATS
    # ============================================================

    def get_full_stats(
        self
    ) -> Dict:

        stats = {}

        for location in self.locations:

            stats[location] = {
                "mssi_score":
                    self.compute_mssi(
                        location
                    ),
                "genes": {}
            }

            for gene in self.genes:

                fold_change = float(
                    self.data[gene][
                        "fold_change"
                    ][location]
                )

                delta_ct_mean = (
                    self.compute_delta_ct_mean(
                        gene,
                        location
                    )
                )

                delta_ct_sd = (
                    self.compute_delta_ct_sd(
                        gene,
                        location
                    )
                )

                delta_ct_se = (
                    self.compute_delta_ct_se(
                        gene,
                        location
                    )
                )

                fc_lower_sd, fc_upper_sd = (
                    self.compute_fold_change_bounds(
                        gene,
                        location,
                        error_type="sd"
                    )
                )

                (
                    fc_error_lower_sd,
                    fc_error_upper_sd
                ) = (
                    self.compute_fold_change_errors(
                        gene,
                        location,
                        error_type="sd"
                    )
                )

                fc_lower_se, fc_upper_se = (
                    self.compute_fold_change_bounds(
                        gene,
                        location,
                        error_type="se"
                    )
                )

                (
                    fc_error_lower_se,
                    fc_error_upper_se
                ) = (
                    self.compute_fold_change_errors(
                        gene,
                        location,
                        error_type="se"
                    )
                )

                stats[
                    location
                ][
                    "genes"
                ][
                    gene
                ] = {

                    # ------------------------------------------------
                    # Existing compatibility fields
                    # ------------------------------------------------

                    "mean":
                        self.compute_mean(
                            gene,
                            location
                        ),

                    "sd":
                        self.compute_sd(
                            gene,
                            location
                        ),

                    "se":
                        self.compute_se(
                            gene,
                            location
                        ),

                    "fold_change":
                        fold_change,

                    # ------------------------------------------------
                    # Corrected DeltaCt statistics
                    # ------------------------------------------------

                    "delta_ct_mean":
                        (
                            round(
                                delta_ct_mean,
                                4
                            )
                            if delta_ct_mean
                            is not None
                            else None
                        ),

                    "delta_ct_sd":
                        (
                            round(
                                delta_ct_sd,
                                4
                            )
                            if delta_ct_sd
                            is not None
                            else None
                        ),

                    "delta_ct_se":
                        (
                            round(
                                delta_ct_se,
                                4
                            )
                            if delta_ct_se
                            is not None
                            else None
                        ),

                    # ------------------------------------------------
                    # Asymmetric SD fold-change bounds
                    # ------------------------------------------------

                    "fc_lower_sd":
                        round(
                            fc_lower_sd,
                            6
                        ),

                    "fc_upper_sd":
                        round(
                            fc_upper_sd,
                            6
                        ),

                    "fc_error_lower_sd":
                        round(
                            fc_error_lower_sd,
                            6
                        ),

                    "fc_error_upper_sd":
                        round(
                            fc_error_upper_sd,
                            6
                        ),

                    # ------------------------------------------------
                    # Asymmetric SEM fold-change bounds
                    # ------------------------------------------------

                    "fc_lower_se":
                        round(
                            fc_lower_se,
                            6
                        ),

                    "fc_upper_se":
                        round(
                            fc_upper_se,
                            6
                        ),

                    "fc_error_lower_se":
                        round(
                            fc_error_lower_se,
                            6
                        ),

                    "fc_error_upper_se":
                        round(
                            fc_error_upper_se,
                            6
                        ),

                    # ------------------------------------------------
                    # Biological replicate count
                    # ------------------------------------------------

                    "n_biological":
                        self.data[gene].get(
                            "n_biological"
                        ),
                }

        return stats

    # ============================================================
    # RANKING
    # ============================================================

    def rank_locations(
        self
    ) -> List[tuple]:

        scores = self.compute_all()

        return sorted(
            scores.items(),
            key=lambda x: x[1],
            reverse=True
        )

    # ============================================================
    # NORMALIZATION
    # ============================================================

    def normalize_scores(
        self
    ) -> Dict:

        scores = (
            self.compute_all()
        )

        max_score = max(
            scores.values()
        )

        min_score = min(
            scores.values()
        )

        if max_score == min_score:

            return {
                location: 0.0
                for location
                in scores
            }

        return {
            location:
                round(
                    (
                        score
                        - min_score
                    )
                    /
                    (
                        max_score
                        - min_score
                    ),
                    4
                )

            for location, score
            in scores.items()
        }