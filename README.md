# 🧬 MSSI Pipeline

**Multi-HSP Stress Signature Index — Computational Pipeline for Stress Classification in *Spodoptera frugiperda***

[![Python](https://img.shields.io/badge/Python-3.9%2B-blue)](https://www.python.org/)
[![License](https://img.shields.io/badge/License-MIT-green)](LICENSE)
[![PyPI](https://img.shields.io/badge/PyPI-mssi--pipeline-orange)](https://pypi.org/project/mssi-pipeline/)

---

## Overview

MSSI Pipeline is an open-source Python tool that transforms raw qRT-PCR fold change data from multiple small heat shock protein (sHsp) genes into a single, interpretable **Multi-HSP Stress Signature Index (MSSI)**. It enables researchers to classify environmental stress levels across geographic locations and generate publication-ready reports automatically.

Developed using expression data from *Spodoptera frugiperda* (Fall Armyworm) collected from four Egyptian governorates, but designed to work with **any insect species** and **any geographic locations** worldwide.

---

## Features

- Accepts multi-sheet Excel files with qRT-PCR fold change data
- Computes weighted MSSI scores per location
- Classifies stress levels: Low / Moderate / High / Critical
- Generates confidence scores based on gene agreement
- Produces heatmaps, bar charts, and gene profile plots
- Exports full PDF reports and JSON results
- Available as both a Command Line Interface (CLI) and a Streamlit GUI
- Fully open-source and pip-installable

---

## Installation

```bash
pip install mssi-pipeline
```

Or install from source:

```bash
git clone https://github.com/YourUsername/MSSI-Pipeline.git
cd MSSI-Pipeline
pip install -e .
```

---

## Quick Start

### Command Line

```bash
mssi --input data.xlsx --output results --summary
```

### Python API

```python
from MSSI_Pipeline.loader     import MSSILoader
from MSSI_Pipeline.calculator import MSSICalculator
from MSSI_Pipeline.classifier import MSSIClassifier
from MSSI_Pipeline.visualizer import MSSIVisualizer
from MSSI_Pipeline.reporter   import MSSIReporter

loader     = MSSILoader("data.xlsx")
data       = loader.load()

calculator = MSSICalculator(data)
stats      = calculator.get_full_stats()
normalized = calculator.normalize_scores()

classifier = MSSIClassifier(stats, normalized)
results    = classifier.classify_all()

visualizer = MSSIVisualizer(results, output_dir="outputs")
visualizer.plot_all()

reporter   = MSSIReporter(results, output_dir="outputs")
reporter.export_all()
```

### GUI

```bash
mssi-gui
```

---

## Input Format

The input must be an `.xlsx` file with one sheet per sHsp gene. Each sheet must contain:

| Column | Description |
|--------|-------------|
| Control | Control group values (3 replicates) |
| B | Location 1 values |
| S | Location 2 values |
| G | Location 3 values |
| Q | Location 4 values |

A row labeled `fold` containing the pre-calculated fold change values must be present in each sheet.

See `data/example.xlsx` for a template.

---

## CLI Options

```
usage: mssi --input FILE [options]

required:
  --input,  -i     Path to Excel file (.xlsx)

optional:
  --output, -o     Output directory (default: outputs)
  --weights        Gene weights space-separated (e.g. 1.0 1.5 2.0)
  --summary        Print results to terminal
  --no-plots       Skip plot generation
  --no-pdf         Skip PDF report
  --json-only      Export JSON results only
```

---

## MSSI Formula

```
MSSI = Σ (w_i × FoldChange_i) / n
```

Where:
- `w_i` = weight of gene i (default = 1.0)
- `FoldChange_i` = fold change of gene i at the location
- `n` = total number of genes

---

## Stress Classification

| MSSI Score | Stress Level | Label |
|------------|-------------|-------|
| ≥ 2.5 | very_high | Critical Environmental Stress |
| 1.5 – 2.5 | high | High Environmental Stress |
| 0.8 – 1.5 | moderate | Moderate Environmental Stress |
| < 0.8 | low | Low Environmental Stress |

---

## Output Files

| File | Description |
|------|-------------|
| `mssi_report.pdf` | Full PDF report with tables and plots |
| `mssi_results.json` | Machine-readable results |
| `mssi_bar.png` | MSSI scores bar chart |
| `heatmap.png` | Gene expression heatmap |
| `gene_profiles.png` | Per-gene grouped bar chart |
| `confidence.png` | Classification confidence chart |

---

## Case Study: Fall Armyworm in Egypt

Expression of three sHsp genes (sHsp19.74, sHsp20.7, sHsp19.07) was measured in *S. frugiperda* collected from four Egyptian governorates using qRT-PCR with Actin as housekeeping gene.

| Location | MSSI Score | Stress Level |
|----------|-----------|-------------|
| Qalyubia | 3.07 | Critical Environmental Stress |
| Giza | 1.87 | High Environmental Stress |
| Beni Suef | 0.62 | Low Environmental Stress |
| Sohag | 0.36 | Low Environmental Stress |

---

## Project Structure

```
MSSI-Pipeline/
├── MSSI_Pipeline/
│   ├── __init__.py
│   ├── loader.py
│   ├── calculator.py
│   ├── classifier.py
│   ├── visualizer.py
│   └── reporter.py
├── gui/
│   └── app.py
├── cli/
│   └── main.py
├── tests/
│   └── test_calculator.py
├── data/
│   └── example.xlsx
├── pyproject.toml
├── requirements.txt
├── README.md
└── LICENSE
```

---

## Dependencies

```
pandas >= 2.0
numpy >= 1.24
matplotlib >= 3.7
seaborn >= 0.12
fpdf2 >= 2.7
openpyxl >= 3.1
streamlit >= 1.30
```

---

## Contributing

Contributions are welcome. Please open an issue or submit a pull request on GitHub.

---

## Citation

If you use MSSI Pipeline in your research, please cite:

```
Your Name et al. (2026). MSSI Pipeline: A Computational Tool for
Multi-HSP Stress Signature Index in Insects.
GitHub: https://github.com/YourUsername/MSSI-Pipeline
```

---

## License

This project is licensed under the MIT License. See [LICENSE](LICENSE) for details.

---

## Contact

For questions or collaboration: your@email.com