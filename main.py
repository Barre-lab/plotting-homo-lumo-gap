#!/usr/bin/env amspython

# Local libraries
from parser import parse_arguments
from support import find_input_type
from plotting import plot_all_together
from plotting import plot_all_separate
from plotting import plot_averages
from plotting import plot_distribution
from plotting import plot_multiple_distributions
from settings import Settings


def main(args):
    """
    Main part of the code
    """

    # Accepted output files
    accepted_files = ["geometry_optimization.out", "geometry-optimization.out",
                      "single-point.out", "out"]

    # Input type: either "sequence" or "collection"
    input_type = find_input_type(args.input)

    # Plot settings (a separate class)
    plot_settings = Settings()

    # Making two arguments mutually exclusive (should be done in parser)
    if args.moving_average and args.separate_states:
        raise ValueError("Combination of moving average and separating states is not supported")

    if args.moving_average and args.include_energies:
        raise ValueError("Combination of moving average and including energies is not supported")

    # Plotting
    if args.distribution and input_type == "sequence":
        plot_distribution(args, input_type, accepted_files, plot_settings)
    elif args.distribution and input_type == "collection":
        plot_multiple_distributions(args, input_type, accepted_files, plot_settings)
    elif input_type == "sequence" and not args.multi_plot:
        plot_all_together(args, input_type, accepted_files, plot_settings)
    elif input_type == "collection" and args.average:
        plot_averages(args, input_type, accepted_files, plot_settings)
    elif input_type == "collection" and args.multi_plot:
        plot_all_together(args, input_type, accepted_files, plot_settings)
    elif input_type == "collection":
        plot_all_separate(args, input_type, accepted_files, plot_settings)
    elif input_type == "sequence" and args.multi_plot:
        plot_all_separate(args, input_type, accepted_files, plot_settings)
    else:
        print("Not plotting anything")


if __name__ == '__main__':
    arguments = parse_arguments()
    main(arguments)
