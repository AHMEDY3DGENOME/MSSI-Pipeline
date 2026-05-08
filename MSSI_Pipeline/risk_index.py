from dataclasses import dataclass
from typing import Dict, List
from .classifier import StressResult

RISK_THRESHOLDS = {
    "critical": (0.80, 1.00),
    "high":     (0.60, 0.80),
    "moderate": (0.40, 0.60),
    "low":      (0.00, 0.40),
}

RISK_LABELS = {
    "critical": "Critical Resistance Risk",
    "high":     "High Resistance Risk",
    "moderate": "Moderate Resistance Risk",
    "low":      "Low Resistance Risk",
}

RISK_ACTIONS = {
    "critical": "Immediate intervention required. Rotate insecticide classes urgently.",
    "high":     "Monitor resistance closely. Consider insecticide rotation.",
    "moderate": "Apply preventive measures. Continue standard monitoring.",
    "low":      "Standard pest management protocols sufficient.",
}

RISK_STRATEGIES = {
    "critical": (
        "STRICT MoA ROTATION: Suspend use of current insecticide class. "
        "Switch to IRAC Group 28 (Diamides) or biologicals (Bacillus thuringiensis). "
        "Implement non-chemical cultural controls immediately."
    ),
    "high": (
        "RESISTANCE MANAGEMENT: Implement a 30-day 'window' strategy. "
        "Rotate between Group 5 (Spinosyns) and Group 22 (Indoxacarb). "
        "Increase pheromone trap monitoring frequency."
    ),
    "moderate": (
        "PREVENTIVE CARE: Utilize Integrated Pest Management (IPM). "
        "Incorporate neem-based products and support natural predators. "
        "Scout fields weekly to detect early population spikes."
    ),
    "low": (
        "ROUTINE MONITORING: Continue standard localized treatment. "
        "Follow label rates strictly to prevent selection pressure. "
        "Maintain baseline susceptibility records."
    ),
}

RISK_COLORS = {
    "critical": "#D73027",
    "high":     "#FC8D59",
    "moderate": "#91BFDB",
    "low":      "#4DAC26",
}


@dataclass
class RiskResult:
    location:         str
    mssi_score:       float
    normalized_score: float
    risk_index:       float
    risk_level:       str
    risk_label:       str
    risk_action:      str
    risk_strategy:    str
    risk_color:       str
    gene_profile:     Dict


class ResistanceRiskIndex:

    def __init__(self, stress_results: List[StressResult]):
        self.stress_results = stress_results
        self.max_mssi       = max(r.mssi_score for r in stress_results)
        self.min_mssi       = min(r.mssi_score for r in stress_results)

    def _compute_risk_index(self, result: StressResult) -> float:
        norm       = result.normalized_score
        confidence = result.confidence
        up_genes   = sum(
            1 for g in result.gene_profile.values()
            if g["fold_change"] > 1.0
        )
        gene_ratio = up_genes / len(result.gene_profile)
        risk = (
            0.50 * norm      +
            0.30 * gene_ratio +
            0.20 * confidence
        )
        return round(min(risk, 1.0), 4)

    def _classify_risk(self, risk_index: float) -> str:
        for level, (low, high) in RISK_THRESHOLDS.items():
            if low <= risk_index <= high:
                return level
        return "low"

    def compute_risk(self, result: StressResult) -> RiskResult:
        risk_index = self._compute_risk_index(result)
        risk_level = self._classify_risk(risk_index)
        return RiskResult(
            location         = result.location,
            mssi_score       = result.mssi_score,
            normalized_score = result.normalized_score,
            risk_index       = risk_index,
            risk_level       = risk_level,
            risk_label       = RISK_LABELS[risk_level],
            risk_action      = RISK_ACTIONS[risk_level],
            risk_strategy    = RISK_STRATEGIES[risk_level],
            risk_color       = RISK_COLORS[risk_level],
            gene_profile     = result.gene_profile,
        )

    def compute_all(self) -> List[RiskResult]:
        return [self.compute_risk(r) for r in self.stress_results]

    def get_highest_risk(self) -> RiskResult:
        return max(self.compute_all(), key=lambda r: r.risk_index)

    def get_lowest_risk(self) -> RiskResult:
        return min(self.compute_all(), key=lambda r: r.risk_index)

    def summary(self) -> List[Dict]:
        return [
            {
                "location":    r.location,
                "mssi_score":  r.mssi_score,
                "risk_index":  r.risk_index,
                "risk_level":  r.risk_level,
                "risk_label":  r.risk_label,
                "risk_action": r.risk_action,
                "strategy":    r.risk_strategy
            }
            for r in self.compute_all()
        ]