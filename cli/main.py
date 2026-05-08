import argparse
import sys
from pathlib import Path
from MSSI_Pipeline.loader     import MSSILoader
from MSSI_Pipeline.calculator import MSSICalculator
from MSSI_Pipeline.classifier import MSSIClassifier
from MSSI_Pipeline.visualizer import MSSIVisualizer
from MSSI_Pipeline.reporter   import MSSIReporter


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="mssi",
        description="MSSI Pipeline — Multi-HSP Stress Signature Index",
    )
    parser.add_argument("--input",    "-i", required=True,  help="Path to Excel file (.xlsx)")
    parser.add_argument("--output",   "-o", default="outputs", help="Output directory")
    parser.add_argument("--no-plots",       action="store_true", help="Skip plot generation")
    parser.add_argument("--no-pdf",         action="store_true", help="Skip PDF report")
    parser.add_argument("--json-only",      action="store_true", help="Export JSON only")
    parser.add_argument("--weights", nargs="+", type=float,      help="Gene weights (space-separated)")
    parser.add_argument("--summary",        action="store_true", help="Print summary to terminal")
    return parser


def print_banner():
    print("""
╔══════════════════════════════════════════════╗
║         MSSI Pipeline  v1.0.0                ║
║  Multi-HSP Stress Signature Index            ║
║  Target: Spodoptera frugiperda (FAW)         ║
╚══════════════════════════════════════════════╝
""")


def print_summary(results):
    print(f"\n{'Location':<20} {'MSSI':>8} {'Stress Level':<30} {'Confidence':>10}")
    print("-" * 72)
    for r in results:
        print(
            f"{r.location:<20} {r.mssi_score:>8.4f} "
            f"{r.stress_label:<30} {int(r.confidence*100):>9}%"
        )
    print()


def run(args):
    print_banner()

    print("[1/5] Loading data...")
    loader = MSSILoader(args.input)
    data   = loader.load()
    genes  = loader.get_genes()

    weights = None
    if args.weights:
        if len(args.weights) != len(genes):
            print(f"Error: expected {len(genes)} weights, got {len(args.weights)}")
            sys.exit(1)
        weights = dict(zip(genes, args.weights))

    print("[2/5] Computing MSSI scores...")
    calculator = MSSICalculator(data, weights=weights)
    stats      = calculator.get_full_stats()
    normalized = calculator.normalize_scores()
    ranked     = calculator.rank_locations()

    print("[3/5] Classifying stress levels...")
    classifier = MSSIClassifier(stats, normalized)
    results    = classifier.classify_all()

    if args.summary:
        print_summary(results)

    if not args.no_plots:
        print("[4/5] Generating plots...")
        visualizer = MSSIVisualizer(results, output_dir=args.output)
        visualizer.plot_all()
    else:
        print("[4/5] Skipping plots.")

    print("[5/5] Exporting report...")
    reporter = MSSIReporter(results, output_dir=args.output)
    if args.json_only:
        reporter.export_json()
    elif args.no_pdf:
        reporter.export_json()
    else:
        reporter.export_all()

    print(f"Done. Results saved to: {args.output}/")
    print(f"Highest stress: {classifier.get_highest_stress().location}")
    print(f"Lowest stress:  {classifier.get_lowest_stress().location}")


def main():
    parser = build_parser()
    args   = parser.parse_args()
    try:
        run(args)
    except FileNotFoundError as e:
        print(f"Error: {e}")
        sys.exit(1)
    except ValueError as e:
        print(f"Error: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()