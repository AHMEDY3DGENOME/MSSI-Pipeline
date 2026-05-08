import pandas as pd
from pathlib import Path
from typing import Dict, List


LOCATION_MAP = {
    "B": "Beni Suef",
    "S": "Sohag",
    "G": "Giza",
    "Q": "Qalyubia",
}

EXPECTED_GENES = ["Hsp 19.74", "Hsp 20.7", "Hsp 19.07"]


class MSSILoader:

    def __init__(self, filepath: str):
        self.filepath = Path(filepath)
        self._validate_file()

    def _validate_file(self):
        if not self.filepath.exists():
            raise FileNotFoundError(
                f"File not found: {self.filepath}"
            )
        if self.filepath.suffix not in [".xlsx", ".xls"]:
            raise ValueError(
                "File must be .xlsx or .xls format"
            )

    def _extract_replicates(self, df: pd.DataFrame) -> Dict:
        results = {}
        header_rows = []
        for idx, row in df.iterrows():
            values = row.dropna().values
            if "Control" in str(values):
                header_rows.append(idx)

        if not header_rows:
            raise ValueError("No valid data structure found in sheet")

        first_header = header_rows[0]
        cols = ["Control", "B", "S", "G", "Q"]
        data_rows = []

        for idx in range(first_header + 1, first_header + 4):
            if idx in df.index:
                row = df.iloc[idx].dropna()
                if len(row) >= 5:
                    data_rows.append(row.values[:5].tolist())

        for i, col in enumerate(cols):
            label = LOCATION_MAP.get(col, col)
            results[label] = [float(row[i]) for row in data_rows]

        return results

    def _extract_fold_change(self, df: pd.DataFrame) -> Dict:
        fold_row = None
        for idx, row in df.iterrows():
            if "fold" in str(row.values).lower():
                fold_row = row
                break

        if fold_row is None:
            raise ValueError("No fold change row found")

        numeric = fold_row.dropna()
        numeric = [x for x in numeric if isinstance(x, (int, float))]
        cols = ["B", "S", "G", "Q"]

        return {
            LOCATION_MAP[cols[i]]: float(numeric[i])
            for i in range(len(cols))
            if i < len(numeric)
        }

    def load(self) -> Dict:
        sheets = pd.read_excel(self.filepath, sheet_name=None, header=None)
        data = {}

        for gene_name, df in sheets.items():
            data[gene_name] = {
                "replicates": self._extract_replicates(df),
                "fold_change": self._extract_fold_change(df),
            }

        return data

    def get_genes(self) -> List[str]:
        sheets = pd.read_excel(self.filepath, sheet_name=None, header=None)
        return list(sheets.keys())

    def get_locations(self) -> List[str]:
        return list(LOCATION_MAP.values())