"""
Command-line entry point for MSSI Pipeline.

Usage:
    mssi run
    mssi version
    mssi --help
"""

import argparse
import sys


VERSION = "1.0.0"


def run_gui():
    """
    Launch the MSSI graphical user interface.
    """
    try:
        from gui.app import main as gui_main
    except ImportError as exc:
        print("Error: Could not import the MSSI GUI.")
        print("Make sure the package is installed correctly using:")
        print("  pip install -e .")
        print(f"\nDetails: {exc}")
        sys.exit(1)

    gui_main()


def main():
    parser = argparse.ArgumentParser(
        prog="mssi",
        description="MSSI Pipeline - Multi-HSP Stress Signature Index"
    )

    parser.add_argument(
        "command",
        nargs="?",
        default="run",
        choices=["run", "version"],
        help="Command to execute. Use 'run' to launch the GUI."
    )

    parser.add_argument(
        "--version",
        action="version",
        version=f"MSSI Pipeline v{VERSION}"
    )

    args = parser.parse_args()

    if args.command == "run":
        run_gui()

    elif args.command == "version":
        print(f"MSSI Pipeline v{VERSION}")


if __name__ == "__main__":
    main()