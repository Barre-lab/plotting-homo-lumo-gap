#!/usr/bin/env amspython

# Local libraries
from parser import parse_arguments
from support import find_input_type
from plotting import plot_all_together
from settings import Settings


def main(args):
    """
    Main part of the code
    """

    # Accepted output files
    accepted_files = ["single-point.out", "out"]

    # Input type: either "sequence" or "collection"
    input_type = find_input_type(args.input)

    # Plot settings (a separate class)
    plot_settings = Settings()

    if input_type == "sequence" and not args.multi_plot:
        plot_all_together(args, input_type, accepted_files, plot_settings)
    elif input_type == "collection" and args.multi_plot:
        plot_all_together(args, input_type, accepted_files, plot_settings)
    else:
        print("Not plotting anything")


if __name__ == '__main__':
    arguments = parse_arguments()
    main(arguments)
