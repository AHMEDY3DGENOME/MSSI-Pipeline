import numpy as np
from typing import Dict, List, Optional


class MSSICalculator:

    def __init__(self, data: Dict, weights: Optional[Dict] = None):
        self.data = data
        self.genes = list(data.keys())
        self.locations = list(data[self.genes[0]]["fold_change"].keys())
        self.weights = weights or {gene: 1.0 for gene in self.genes}

    def _get_fold_changes(self, location: str) -> List[float]:
        return [
            self.data[gene]["fold_change"][location]
            for gene in self.genes
        ]

    def _get_replicates(self, gene: str, location: str) -> List[float]:
        return self.data[gene]["replicates"][location]

    def compute_mssi(self, location: str) -> float:
        fold_changes = self._get_fold_changes(location)
        weighted_sum = sum(
            self.weights[gene] * fold_changes[i]
            for i, gene in enumerate(self.genes)
        )
        return round(weighted_sum / len(self.genes), 4)

    def compute_all(self) -> Dict:
        return {
            location: self.compute_mssi(location)
            for location in self.locations
        }

    def compute_se(self, gene: str, location: str) -> float:
        replicates = self._get_replicates(gene, location)
        return round(np.std(replicates, ddof=1) / np.sqrt(len(replicates)), 4)

    def compute_sd(self, gene: str, location: str) -> float:
        replicates = self._get_replicates(gene, location)
        return round(float(np.std(replicates, ddof=1)), 4)

    def compute_mean(self, gene: str, location: str) -> float:
        replicates = self._get_replicates(gene, location)
        return round(float(np.mean(replicates)), 4)

    def get_full_stats(self) -> Dict:
        stats = {}
        for location in self.locations:
            stats[location] = {
                "mssi_score": self.compute_mssi(location),
                "genes": {}
            }
            for gene in self.genes:
                stats[location]["genes"][gene] = {
                    "mean": self.compute_mean(gene, location),
                    "sd": self.compute_sd(gene, location),
                    "se": self.compute_se(gene, location),
                    "fold_change": self.data[gene]["fold_change"][location],
                }
        return stats

    def rank_locations(self) -> List[tuple]:
        scores = self.compute_all()
        return sorted(
            scores.items(),
            key=lambda x: x[1],
            reverse=True
        )

    def normalize_scores(self) -> Dict:
        scores = self.compute_all()
        max_score = max(scores.values())
        min_score = min(scores.values())
        return {
            loc: round((score - min_score) / (max_score - min_score), 4)
            for loc, score in scores.items()
        }