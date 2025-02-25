# Standard libraries
import argparse


def parse_arguments():
    """
    Parsing arguments given to script
    """

    parser = argparse.ArgumentParser()

    parser.add_argument(
            "input",
            type=str,
            nargs="*",
            help="Input folder to compute homo-lumo gaps in")

    parser.add_argument(
            "-t", "--type",
            choices=["calculations", "sequence"],
            default="calculations",
            help="Input folder can be a calculation folder or a sequence of calculation folders")

    parser.add_argument(
            "-l", "--labels",
            type=str,
            nargs="*",
            help="Labels of folders with calculations or sequences of folders with calculations")

    parser.add_argument(
            "-w", "--what",
            choices=["avg", "all"],
            default="all",
            help="Plot average gaps as function of N H2O or all values against frames")

    parser.add_argument(
            "-ho", "--how",
            choices=["together", "separate"],
            default="together",
            help="Plot data in together in one graph or in separate sub graphs (only for -w: all)")

    parser.add_argument(
            "-s", "--separate-states",
            action="store_true",
            help="Separate data for state where fragment configuration has changed")

    parser.add_argument(
            "-ma", "--moving-average",
            type=int,
            default=None,
            help="Npoints to include on both sides for moving averages (only for -w: all)")

    parser.add_argument(
            "-e", "--include-energies",
            action="store_true",
            help="Include total energies in the plot(s)")

    parser.add_argument(
            "-lp", "--legend-position",
            type=str,
            default="best",
            help="The position of the legend")

    parser.add_argument(
            "-lpc", "--legend-position-coords",
            type=float,
            nargs=2,
            default=None,
            help="Coordinates of legend position for better control (combine with -lp)")

    parser.add_argument(
            "-lc", "--legend-columns",
            type=int,
            default=1,
            help="Number of columns the legend should have")

    parser.add_argument(
            "-li", "--legend-index",
            type=int,
            default=0,
            help="In which plot to include the legend (in case of separate plots)")

    parser.add_argument(
            "-sp", "--plot-name",
            type=str,
            help="Save plot as plot name")

    return parser.parse_args()
