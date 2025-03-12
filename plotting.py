# standard libraries
import os
import re
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
                      settings: classmethod) -> None:
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
    plt.figure(figsize=settings.figsize)

    # Looping over all the sequence with calculations to plot them
    for index, calc_sequence in enumerate(sequences):

        sequence = calculation_sequence(calc_sequence, accepted_files)

        color = settings.colors[index]

        # Defining the label for sequence
        if args.labels:
            label = args.labels[index]
        elif args.moving_average:
            label = f"{sequence.nwater} H2O (MA)"
        else:
            label = f"{sequence.nwater} H2O"

        # Creating sequence handle for legend
        sequence_handle = Line2D([0], [0], marker="s", label=label,
                                 color=settings.colors[index], ms=5, ls="")
        handles_sequences.append(sequence_handle)

        # Moving average
        if args.moving_average:
            averaged_gaps = get_moving_average(sequence.gaps, args.moving_average)

            plt.plot(sequence.frames, sequence.gaps, color=color, alpha=0.1, lw=0.8)
            plt.plot(sequence.frames, averaged_gaps, color=color, lw=1.2, label=label+" (MA)")

        # Separating states
        elif args.separate_states:

            # Plotting the connecting lines
            plt.plot(sequence.frames, sequence.gaps, color=color, lw=0.4, zorder=1)

            # Looping over each state to plot
            fillcolors = np.array(["white", color])

            for subindex, state in enumerate(sequence.states):
                fillcolor = fillcolors[subindex % 2]
                symbol = settings.symbols[subindex]

                plt.plot(sequence.frames[state], sequence.gaps[state], color=color, marker=symbol,
                         ms=5, ls="", markerfacecolor=fillcolor)

                if len(sequences) == 1:
                    state_handle = Line2D([0], [0], marker=symbol, label=f"State {subindex + 1}",
                                          color=color, ls="", ms=5, markerfacecolor=fillcolor)
                    handles_states.append(state_handle)

        # No moving average and not separating states
        else:
            plt.plot(sequence.frames, sequence.gaps, color=color, lw=0.8, label=label)

    # Defining font size for axis and ticks
    fs = settings.axes_size
    ts = settings.tick_size

    # Naming plot axis
    plt.xlabel("Frame number", fontsize=fs)
    plt.ylabel("Homo-lumo gap [eV]", fontsize=fs)

    # Specifying axis ticks
    plt.tick_params(axis="both", labelsize=ts)

    # Including legend
    if handles_states:
        plt.legend(handles=handles_states, fontsize=fs, loc=args.legend_position, frameon=False,
                   ncol=args.legend_columns, bbox_to_anchor=args.legend_position_coords)
    elif len(sequences) > 1:
        plt.legend(handles=handles_sequences, fontsize=fs, loc=args.legend_position, frameon=False,
                   ncol=args.legend_columns, bbox_to_anchor=args.legend_position_coords)

    plt.tight_layout()

    # Saving plot
    if args.plot_name:
        plt.savefig(args.plot_name)

    plt.show()


def plot_all_separate(args: argparse.Namespace, input_type: str, accepted_files: List[str],
                      settings: classmethod) -> None:
    """
    Plots homo-lumo gaps and energies as a function of the frame number for all input
    sequences together in sub-graphs based on the number of water molecules
    """

    def _plot_sequence_data_in_subplot(sequence, axs, all_water, color) -> None:
        """
        Finds the correct sub-graph and plots sequence there (based on number of water molecules)
        """

        # Initializing sequence class (which gets all data)
        sequence = calculation_sequence(sequence, accepted_files)

        # Finding index of axs to plot in
        ax_index = all_water.index(sequence.nwater)
        if axs.ndim > 1:
            ax_index = (ax_index // 2, ax_index % 2)

        # Plotting data in correct sub-graph
        if args.moving_average:
            averaged_gaps = get_moving_average(sequence.gaps, args.moving_average)

            axs[ax_index].plot(sequence.frames, sequence.gaps, color=color, lw=0.8, alpha=0.1)
            axs[ax_index].plot(sequence.frames, averaged_gaps, color=color, lw=1.0)

        elif args.separate_states:
            # Plotting the connecting lines
            axs[ax_index].plot(sequence.frames, sequence.gaps, color=color, lw=0.2, zorder=1)

            # Looping over each state
            fillcolors = np.array(["white", color])
            for subindex, state in enumerate(sequence.states):
                fillcolor = fillcolors[subindex % 2]
                symbol = settings.symbols[subindex]

                axs[ax_index].plot(sequence.frames[state], sequence.gaps[state], color=color,
                                   marker=symbol, ms=3, ls="", markerfacecolor=fillcolor)

        else:
            axs[ax_index].plot(sequence.frames, sequence.gaps, color=color, lw=0.8)

    # In case a collection of calculation sequences is provided
    if input_type == "collection":

        # Getting a list with all the number of water molecules in the collections
        all_water = []
        for collection in args.input:
            water_collection = [re.search(r"\_(.*?)\_", folder).group(1) for folder in
                                os.listdir(collection) if folder.startswith("calculations")]
            for nwater in water_collection:
                if int(nwater) not in all_water:
                    all_water.append(int(nwater))

        # Sorting list with all waters in ascending order
        all_water.sort()

        # Initializing figure with sub-figures (one sub-figure for each water amount)
        nrows = -(-len(all_water) // 2)
        fig, axs = plt.subplots(nrows=nrows, ncols=2, figsize=(6, 2*nrows))

        # Initializing handles for legend
        handles_collections = []

        # Looping over each collection
        for index, collection in enumerate(args.input):

            color = settings.colors[index]

            # Getting list of all calculation sequences
            sequences = [os.path.join(collection, folder) for folder
                         in os.listdir(collection) if folder.startswith("calculations")]

            # Plotting the data of each sequence in the correct subplot
            for sequence in sequences:
                _plot_sequence_data_in_subplot(sequence, axs, all_water, color)

            # Creating legend for collection
            if args.labels:
                label = args.labels[index]
            else:
                basename = os.path.basename(os.path.normpath(collection))
                label = basename.replace("-", " ").replace("_", " ")

            if args.moving_average:
                label += " (MA)"

            handles_collections.append(Line2D([0], [0], marker="s", label=label, ms=5, ls="",
                                              color=settings.colors[index]))

    elif input_type == "sequence":

        # Initializing figure with sub-figures
        nrows = -(-len(args.input) // 2)
        fig, axs = plt.subplots(nrows=nrows, ncols=2, figsize=(6, 2*nrows))

        # Finding all number of waters in the input sequences
        all_water = []
        for sequence in args.input:
            all_water.append(int(re.search(r"\_(.*?)\_", sequence).group(1)))
        all_water.sort()

        for sequence in sequences:
            _plot_sequence_data_in_subplot(sequence, axs, all_water, color=settings.colors[0])

    else:
        raise ValueError(f"Unknown input type: {input_type}")

    # Defining sub-plot names, axis names and ticks and getting y-ranges
    fs = settings.axes_size
    ts = settings.tick_size
    min_yrange, max_yrange = np.inf, 0
    for ax, nwater in zip(axs.flat, all_water):

        ax.set_ylabel("Homo-lumo gap [eV]", fontsize=fs)
        ax.set_xlabel("Frame number", fontsize=fs)
        ax.label_outer()

        ax.xaxis.set_tick_params(labelbottom=True, labelsize=ts)
        ax.yaxis.set_tick_params(labelleft=True, labelsize=ts)

        if ax.lines:
            ax.set_title(f"{nwater} H2O", fontsize=10)

            min_yrange = min(ax.get_ylim()[0], min_yrange)
            max_yrange = max(ax.get_ylim()[1], max_yrange)
        else:
            ax.set_visible(False)

    # Setting equal y-range for all subplots
    for ax in axs.flat:
        ax.set_ylim(min_yrange, max_yrange)

    # Defining spacing between subplots
    plt.tight_layout()
    fig.subplots_adjust(hspace=0.4)
    fig.subplots_adjust(wspace=0.2)

    # Adding legend
    if axs.ndim > 1:
        ax_index = (args.legend_index // 2, args.legend_index % 2)
    else:
        ax_index = args.legend_index

    if input_type == "collection" and len(args.input) > 1:
        axs[ax_index].legend(handles=handles_collections, fontsize=fs, loc=args.legend_position,
                             frameon=False, ncol=args.legend_columns,
                             bbox_to_anchor=args.legend_position_coords)

    if args.plot_name:
        plt.savefig(args.plot_name)

    plt.show()


def plot_averages(args: argparse.Namespace, input_type: str, accepted_file: List[str],
                  settings: classmethod) -> None:
    """
    Plots average homo-lumo gaps and energies as a function of the number water molecules for each
    sequence in all the input collections
    """

    # Making sure input type is collection
    if input_type != "collection":
        raise ValueError("Cannot plot averages for input type: {input_type}")

    # Initializing plot
    plt.figure(figsize=settings.figsize)

    # Initializing legends for different states and sequences
    handles_collections = []

    # Looping over the collections to plot the average of each sequence
    for index, collection in enumerate(args.input):

        color = settings.color[index]

        
        sequence = calculation_sequence(collection, accepted_file)
