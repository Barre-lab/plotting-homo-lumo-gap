#!/usr/bin/env amspython

# Local libraries
from parser import parse_arguments
from support import find_input_type
from plotting import plot_all_together


def main(args):
    """
    Main part of the code
    """

    # Defining colors to plot with
    colors = ["blue", "orange", "limegreen", "red"]

    # Defining symbols to indicate different states
    symbols = ["o", "v", "s", "^"]

    # Accepted output files
    accepted_files = ["single-point.out", "out"]

    input_type = find_input_type(args.input)

    if input_type == "sequence" and not args.multi_plot:
        plot_all_together(args, input_type, accepted_files, colors, symbols)
    elif input_type == "collection" and args.multi_plot:
        plot_all_together(args, input_type, accepted_files, colors, symbols)
    else:
        print("Not plotting anything")


if __name__ == '__main__':
    arguments = parse_arguments()
    main(arguments)
