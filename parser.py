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
            "-l", "--labels",
            type=str,
            nargs="*",
            help="Labels of folders with calculations or sequences of folders with calculations")

    parser.add_argument(
            "-p", "--select_points",
            type=int,
            help="Take only the p points with the lowest total energy for each sequence")

    parser.add_argument(
            "-avg", "--average",
            action="store_true",
            help="Plot average data as function of N H2O instead of all data against frames")

    parser.add_argument(
            "-m", "--multi-plot",
            action="store_true",
            help="Enforce plottting data in multiple subgraphs instead of the default type")

    parser.add_argument(
            "-v", "--vertical",
            action="store_true",
            help="Transpose columns and rows of subplots, making them extend vertically")

    parser.add_argument(
            "-s", "--separate-states",
            action="store_true",
            help="Separate data for states where the fragment configuration has changed")

    parser.add_argument(
            "-ma", "--moving-average",
            type=int,
            default=None,
            help="Npoints to include on both sides for moving averages (not of -avg or -p)")

    parser.add_argument(
            "-e", "--include-energies",
            action="store_true",
            help="Include total energies in the plot(s)")

    parser.add_argument(
            "-a", "--ascending-energies",
            action="store_true",
            help="Gaps and energies are reordered in ascending energy order (to find correlation)")

    parser.add_argument(
            "-ln", "--legend-none",
            action="store_true",
            help="Do not include legend in plot")

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
