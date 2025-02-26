# standard libraries
import os
import time
import argparse
from typing import List

# Third-part libraries
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np

# Local libraries
from sequence import calculation_sequence
from support import get_moving_average


def plot_all_together(args: argparse.Namespace, input_type: str, accepted_files: List[str],
                      colors: list, symbols: list) -> None:
    """
    Plots homo-lumo gaps and energies as a function of the frame number for all input
    sequences together in one graph
    """

    # Collecting calculation sequences from input
    if input_type == "collection":
        if len(args.input) > 1:
            raise ValueError("Cannot plot multiple collections together in one graph")

        sequences = [args.input[0] + direc for direc in os.listdir(args.input[0]) if
                     direc.startswith("calculations")]
    elif input_type == "sequence":
        sequences = args.input
    else:
        raise ValueError(f"Unknown input type: {input_type}")

    # Initializing legends for different states and sequences
    handles_states = []
    handles_sequences = []

    # Initializing figure
    plt.figure(figsize=(4.5, 3.4))

    # Looping over all the sequence with calculations to plot them
    for index, calc_sequence in enumerate(sequences):

        sequence = calculation_sequence(calc_sequence, accepted_files)

        # Defining the label for sequence
        if args.labels:
            label = args.labels[index]
        elif args.moving_average:
            label = f"{sequence.nwater} H2O (MA)"
        else:
            label = f"{sequence.nwater} H2O"

        # Creating sequence handle for legend
        sequence_handle = Line2D([0], [0], marker="s", label=label,
                                 color=colors[index], ms=5, ls="")
        handles_sequences.append(sequence_handle)

        # Moving average
        if args.moving_average and not args.separate_states:
            averaged_gaps = get_moving_average(sequence.gaps, args.moving_average)

            plt.plot(sequence.frames, sequence.gaps, color=colors[index], alpha=0.1, lw=0.8)
            plt.plot(sequence.frames, averaged_gaps, color=colors[index], lw=1.2,
                     label=label+" (MA)")

        # Separating states
        elif args.separate_states and not args.moving_average:

            # Plotting the connecting lines
            plt.plot(sequence.frames, sequence.gaps, color=colors[index], lw=0.4, zorder=1)

            # Looping over each state to plot
            fillcolors = np.array(["white", colors[index]])

            for subindex, state in enumerate(sequence.states):
                fillcolor = fillcolors[subindex % 2]

                plt.plot(sequence.frames[state], sequence.gaps[state], color=colors[index],
                         marker=symbols[subindex], ms=5, ls="", markerfacecolor=fillcolor)

                if len(sequences) == 1:
                    state_handle = Line2D([0], [0], marker=symbols[subindex], color=colors[index],
                                          label=f"State {subindex + 1}", markerfacecolor=fillcolor,
                                          ls="", ms=5)
                    handles_states.append(state_handle)

        # No moving average and not separating states
        elif not args.separate_states and not args.moving_average:
            plt.plot(sequence.frames, sequence.gaps, color=colors[index], lw=0.8, label=label)

        else:
            raise ValueError("Combination of moving average and separating states is not supported")

    # Naming plot axis
    plt.xlabel("Frame number")
    plt.ylabel("Homo-lumo gap [eV]")

    # Including legend
    if handles_states:
        plt.legend(handles=handles_states, fontsize=9, loc=args.legend_position, frameon=False)
    elif len(sequences) > 1:
        plt.legend(handles=handles_sequences, fontsize=9, loc=args.legend_position, frameon=False)

    plt.tight_layout()

    # Saving plot
    if args.plot_name:
        plt.savefig(args.plot_name)

    plt.show()
