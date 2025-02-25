# standard libraries
import os
import time
import argparse

# Third-part libraries
import matplotlib.pyplot as plt
import matplotlib.pyplot

# Local libraries
from sequence import calculation_sequence
from support import get_moving_average


def plot_all_together(args: argparse.Namespace, input_type: str, colors: list) -> None:
    """
    Plots homo-lumo gaps and energies as a function of the frame number for all input
    sequences together in one graph
    """

    # In case input is a collection of calculation sequences
    if input_type == "collection":
        if len(args.input) > 1:
            raise ValueError("Cannot plot multiple collections together in one graph")

        sequences = [args.input[0] + direc for direc in os.listdir(args.input[0]) if
                     direc.startswith("calculations")]
    elif input_type == "sequence":
        sequences = args.input
    else:
        raise ValueError(f"Unknown input type: {input_type}")
    
    for index, sequence in enumerate(sequences):



