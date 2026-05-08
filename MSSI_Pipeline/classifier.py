from typing import Dict, List, Tuple
from dataclasses import dataclass


@dataclass
class StressResult:
    location: str
    mssi_score: float
    normalized_score: float
    stress_level: str
    stress_label: str
    confidence: float
    gene_profile: Dict


STRESS_THRESHOLDS = {
    "very_high": (2.5, float("inf")),
    "high":      (1.5, 2.5),
    "moderate":  (0.8, 1.5),
    "low":       (0.0, 0.8),
}

STRESS_LABELS = {
    "very_high": "Critical Environmental Stress",
    "high":      "High Environmental Stress",
    "moderate":  "Moderate Environmental Stress",
    "low":       "Low Environmental Stress",
}


class MSSIClassifier:

    def __init__(self, stats: Dict, normalized_scores: Dict):
        self.stats = stats
        self.normalized_scores = normalized_scores

    def _classify_stress_level(self, mssi_score: float) -> Tuple[str, str]:
        for level, (low, high) in STRESS_THRESHOLDS.items():
            if low <= mssi_score < high:
                return level, STRESS_LABELS[level]
        return "low", STRESS_LABELS["low"]

    def _compute_confidence(self, gene_profile: Dict) -> float:
        fold_values = [
            v["fold_change"]
            for v in gene_profile.values()
        ]
        all_up   = all(f > 1.0 for f in fold_values)
        all_down = all(f < 1.0 for f in fold_values)
        if all_up or all_down:
            return 1.0
        agreeing = sum(1 for f in fold_values if f > 1.0)
        return round(max(agreeing, len(fold_values) - agreeing) / len(fold_values), 2)

    def _build_gene_profile(self, location: str) -> Dict:
        return self.stats[location]["genes"]

    def classify_location(self, location: str) -> StressResult:
        mssi_score       = self.stats[location]["mssi_score"]
        normalized_score = self.normalized_scores[location]
        gene_profile     = self._build_gene_profile(location)
        stress_level, stress_label = self._classify_stress_level(mssi_score)
        confidence       = self._compute_confidence(gene_profile)

        return StressResult(
            location        = location,
            mssi_score      = mssi_score,
            normalized_score= normalized_score,
            stress_level    = stress_level,
            stress_label    = stress_label,
            confidence      = confidence,
            gene_profile    = gene_profile,
        )

    def classify_all(self) -> List[StressResult]:
        return [
            self.classify_location(location)
            for location in self.stats.keys()
        ]

    def get_highest_stress(self) -> StressResult:
        results = self.classify_all()
        return max(results, key=lambda r: r.mssi_score)

    def get_lowest_stress(self) -> StressResult:
        results = self.classify_all()
        return min(results, key=lambda r: r.mssi_score)

    def summary(self) -> List[Dict]:
        return [
            {
                "location":         r.location,
                "mssi_score":       r.mssi_score,
                "stress_level":     r.stress_level,
                "stress_label":     r.stress_label,
                "confidence":       f"{int(r.confidence * 100)}%",
            }
            for r in self.classify_all()
        ]