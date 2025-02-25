#!/usr/bin/env amspython

# Third-party libraries
import matplotlib.pyplot as plt
import numpy as np

# Local libraries
from parser import parse_arguments
from support import find_input_type


def main(args):
    """
    Main part of the code
    """

    # Defining colors to plot with
    colors = ["blue", "orange", "limegreen", "red"]

    input_type = find_input_type(args.input)


if __name__ == '__main__':
    arguments = parse_arguments()
    main(arguments)
