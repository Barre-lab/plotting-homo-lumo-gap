# standard libraries
import os
import re
import argparse
from typing import List

# Third-part libraries
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib import ticker
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
                     direc.startswith(("calculations", "configurations"))]
    elif input_type == "sequence":
        sequences = args.input
        if len(sequences) > 1 and args.include_energies:
            raise ValueError("Cannot plot multiple sequences with energies in one graph")
    else:
        raise ValueError(f"Unknown input type: {input_type}")

    # Initializing legend handles
    handles = []

    # Initializing figure
    fig, ax1 = plt.subplots(figsize=settings.figsize)

    # Looping over all the sequence with calculations to plot them
    for index, calc_sequence in enumerate(sequences):

        sequence = calculation_sequence(calc_sequence, accepted_files, args.ascending_energies,
                                        args.separate_states, args.select_points)

        color = settings.colors[index]

        # Defining the label for sequence
        if args.labels:
            label = args.labels[index]
        elif args.moving_average:
            label = f"{sequence.nwater} H2O (MA)"
        else:
            label = f"{sequence.nwater} H2O"

        # Creating sequence handle for legend
        if len(sequences) > 1:
            handle = Line2D([0], [0], marker="s", label=label,
                            color=settings.colors[index], ms=5, ls="")
            handles.append(handle)

        # Moving average
        if args.moving_average:
            averaged_gaps = get_moving_average(sequence.gaps, args.moving_average)

            plt.plot(sequence.frames, sequence.gaps, color=color, alpha=0.1, lw=0.8)
            plt.plot(sequence.frames, averaged_gaps, color=color, lw=1.2, label=label+" (MA)")

        # Separating states
        elif args.separate_states:

            # Plotting the connecting lines
            ax1.plot(sequence.frames, sequence.gaps, color=color, lw=0.4, zorder=1)
            if args.include_energies:
                ax2 = ax1.twinx()
                color_energies = settings.colors[1]
                ax2.plot(sequence.frames, sequence.energies, color=color_energies, lw=0.4, zorder=1)

                fillcolors_energies = np.array(["white", color_energies])

            # Looping over each state to plot
            fillcolors_gaps = np.array(["white", color])

            for subindex, state in enumerate(sequence.states):
                fillcolor = fillcolors_gaps[subindex % 2]
                symbol = settings.symbols[subindex]

                ax1.plot(sequence.frames[state], sequence.gaps[state], color=color, marker=symbol,
                         ms=5, ls="", markerfacecolor=fillcolor)

                # Including energies
                if args.include_energies:
                    fillcolor_energies = fillcolors_energies[subindex % 2]
                    ax2.plot(sequence.frames[state], sequence.energies[state], marker=symbol,
                             color=color_energies, ms=5, ls="", markerfacecolor=fillcolor_energies)

                # Legend for states
                elif len(sequences) == 1:
                    handle = Line2D([0], [0], marker=symbol, label=f"State {subindex + 1}",
                                    color=color, ls="", ms=5, markerfacecolor=fillcolor)
                    handles.append(handle)

        # Including the energies
        elif args.include_energies and not args.separate_states:

            plt.plot(sequence.frames, sequence.gaps, settings.colors[0], lw=0.8, label=label)
            ax2 = ax1.twinx()
            ax2.plot(sequence.frames, sequence.energies, settings.colors[1], lw=0.8, label=label)

        # No moving average and not separating states
        else:
            plt.plot(sequence.frames, sequence.gaps, color=color, lw=0.8, label=label)

    # Legend for Homo-Lumo gap and energies
    if args.include_energies:
        handles.append(Line2D([0], [0], marker="s", label="Gap", c=settings.colors[0], ls="", ms=5))
        handles.append(Line2D([0], [0], marker="s", label="Energy", c=settings.colors[1], ls="", ms=5))

    # Defining font size for axis and ticks
    fs = settings.axes_size
    ts = settings.tick_size

    # Naming plot axis
    ax1.set_xlabel("Frame number", fontsize=fs)
    ax1.set_ylabel("Homo-lumo gap [eV]", fontsize=fs)
    if args.include_energies:
        ax2.set_ylabel("Total energy [eV]", fontsize=fs)
        offset = ax2.get_yticks()[1]
        ax2.ticklabel_format(axis="y", style="plain", useOffset=offset, useMathText=True)
        ax2.yaxis.get_offset_text().set(fontsize=ts)

    # Specifying axis tickss
    ax1.tick_params(axis="both", labelsize=ts)
    if args.include_energies:
        ax2.tick_params(axis="both", labelsize=ts)

    # Removing ticks and label from x-axis in case of ascending gaps
    if args.ascending_energies:
        ax1.tick_params(bottom=False, labelbottom=False)
        ax1.set_xlabel("")

    # Including legend
    plt.legend(handles=handles, fontsize=fs, loc=args.legend_position, frameon=False,
               ncol=args.legend_columns, bbox_to_anchor=args.legend_position_coords)

    fig.tight_layout()

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

    def _plot_sequence_data_in_subplot(sequence, axs, axs2, all_water, color) -> None:
        """
        Finds the correct sub-graph and plots sequence there (based on number of water molecules)
        """

        # Initializing sequence class (which gets all data)
        sequence = calculation_sequence(sequence, accepted_files, args.ascending_energies,
                                        args.separate_states, args.select_points)

        # Finding index of axs to plot in
        ax_index = all_water.index(sequence.nwater)
        ncols = -(-len(all_water) // 2)
        if axs.ndim > 1:
            if args.vertical:
                ax_index = (ax_index // ncols, ax_index % ncols)
            else:
                ax_index = (ax_index // 2, ax_index % 2)

        # Plotting data in correct sub-graph
        if args.moving_average:
            averaged_gaps = get_moving_average(sequence.gaps, args.moving_average)

            axs[ax_index].plot(sequence.frames, sequence.gaps, color=color, lw=0.8, alpha=0.1)
            axs[ax_index].plot(sequence.frames, averaged_gaps, color=color, lw=1.0)

        elif args.separate_states:
            # Plotting the connecting lines
            axs[ax_index].plot(sequence.frames, sequence.gaps, color=color, lw=0.2, zorder=1)
            if args.include_energies:
                axs2[ax_index] = axs[ax_index].twinx()
                color_ener = settings.colors[1]
                axs2[ax_index].plot(sequence.frames, sequence.energies, color=color_ener, lw=0.2, zorder=1)

                fillcolors_ener = np.array(["white", color_ener])

            # Looping over each state
            fillcolors = np.array(["white", color])
            for subindex, state in enumerate(sequence.states):
                fillcolor = fillcolors[subindex % 2]
                symbol = settings.symbols[subindex]

                axs[ax_index].plot(sequence.frames[state], sequence.gaps[state], color=color,
                                   marker=symbol, ms=3, ls="", markerfacecolor=fillcolor)

                if args.include_energies:
                    fillcolor_ener = fillcolors_ener[subindex % 2]
                    axs2[ax_index].plot(sequence.frames[state], sequence.energies[state], color=color_ener,
                                        marker=symbol, ms=3, ls="", markerfacecolor=fillcolor_ener)

        elif args.include_energies and not args.separate_states:
            axs[ax_index].plot(sequence.frames, sequence.gaps, settings.colors[0], lw=0.8)
            axs2[ax_index] = axs[ax_index].twinx()
            axs2[ax_index].plot(sequence.frames, sequence.energies, settings.colors[1], lw=0.8)

        else:
            axs[ax_index].plot(sequence.frames, sequence.gaps, color=color, lw=0.8)

    # In case a collection of calculation sequences is provided
    if input_type == "collection":

        if args.include_energies and len(args.input) > 1:
            raise ValueError("Cannot plot energies for more than 1 sequence per subplot")

        # Getting a list with all the number of water molecules in the collections
        all_water = []
        for collection in args.input:
            try:
                water_collection = [re.search(r"\_(.*?)\_", folder).group(1) for folder in
                                    os.listdir(collection) if
                                    folder.startswith(("calculations", "configurations"))]
            except AttributeError:
                water_collection = [re.search(r"\-(.*?)\-", folder).group(1) for folder in
                                    os.listdir(collection) if
                                    folder.startswith(("calculations", "configurations"))]
            for nwater in water_collection:
                if int(nwater) not in all_water:
                    all_water.append(int(nwater))

        # Sorting list with all waters in ascending order
        all_water.sort()

        # Initializing figure with sub-figures (one sub-figure for each water amount)
        if args.vertical:
            ncols = -(-len(all_water) // 2)
            fig, axs = plt.subplots(nrows=2, ncols=ncols, figsize=(3*ncols, 4))
        else:
            nrows = -(-len(all_water) // 2)
            fig, axs = plt.subplots(nrows=nrows, ncols=2, figsize=(6, 2*nrows))

        axs2 = np.empty(axs.shape, dtype=object)

        # Initializing handles for legend
        handles = []

        # Looping over each collection
        for index, collection in enumerate(args.input):

            color = settings.colors[index]

            # Getting list of all calculation sequences
            sequences = [os.path.join(collection, folder) for folder in os.listdir(collection)
                         if folder.startswith(("calculations", "configurations"))]

            # Plotting the data of each sequence in the correct subplot
            for sequence in sequences:
                _plot_sequence_data_in_subplot(sequence, axs, axs2, all_water, color)

            # Creating legend for collection
            if args.labels:
                label = args.labels[index]
            else:
                basename = os.path.basename(os.path.normpath(collection))
                label = basename.replace("-", " ").replace("_", " ")

            if args.moving_average:
                label += " (MA)"

            # Appending to legend
            if not args.include_energies and len(args.input) > 1:
                handles.append(Line2D([0], [0], marker="s", label=label, ms=5, ls="",
                                      color=settings.colors[index]))

    elif input_type == "sequence":

        # Initializing figure with sub-figures
        if args.vertical:
            ncols = -(-len(all_water) // 2)
            fig, axs = plt.subplots(nrows=2, ncols=ncols, figsize=(3*ncols, 4))
        else:
            nrows = -(-len(all_water) // 2)
            fig, axs = plt.subplots(nrows=nrows, ncols=2, figsize=(6, 2*nrows))

        axs2 = np.empty(axs.shape, dtype=object)

        # Finding all number of waters in the input sequences
        all_water = []
        for sequence in args.input:
            all_water.append(int(re.search(r"\_(.*?)\_", sequence).group(1)))
        all_water.sort()

        sequences = args.input

        for sequence in sequences:
            _plot_sequence_data_in_subplot(sequence, axs, axs2, all_water, color=settings.colors[0])

    else:
        raise ValueError(f"Unknown input type: {input_type}")

    # Adding legend handles for gaps and energies
    if args.include_energies:
        handles.append(Line2D([0], [0], marker="s", label="Gap", c=settings.colors[0], ls="", ms=5))
        handles.append(Line2D([0], [0], marker="s", label="Energy", c=settings.colors[1], ls="", ms=5))

    # Defining sub-plot names, axis names and ticks and getting y-ranges
    fs = settings.axes_size
    ts = settings.tick_size
    min_yrange, max_yrange = np.inf, 0
    for ax, ax2, nwater in zip(axs.flat, axs2.flat, all_water):

        ax.set_ylabel("Homo-lumo gap [eV]", fontsize=fs)
        ax.set_xlabel("Frame number", fontsize=fs)
        ax.label_outer()

        ax.xaxis.set_tick_params(labelbottom=True, labelsize=ts)

        if ax.lines:
            if args.include_energies:
                ax.set_title(f"{nwater} H2O", fontsize=10, loc="left")
            else:
                ax.set_title(f"{nwater} H2O", fontsize=10)

            min_yrange = min(ax.get_ylim()[0], min_yrange)
            max_yrange = max(ax.get_ylim()[1], max_yrange)
        else:
            ax.set_visible(False)
            ax2.set_visible(False)

        if args.include_energies:
            ax2.set_ylabel("Total energy [eV]", fontsize=fs)
            ax2.label_outer()
            ax2.yaxis.set_tick_params(labelright=True, labelsize=ts)

            offset = ax2.get_yticks()[1]
            ax2.ticklabel_format(axis="y", style="plain", useOffset=offset, useMathText=True)
            ax2.yaxis.get_offset_text().set(fontsize=fs-2)
        else:
            ax.yaxis.set_tick_params(labelleft=True, labelsize=ts)

        if args.ascending_energies:
            ax.tick_params(bottom=False, labelbottom=False)
            ax.set_xlabel("")

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

    if handles and not args.legend_none:
        axs[ax_index].legend(handles=handles, fontsize=fs, loc=args.legend_position,
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

    if args.include_energies and len(args.input) > 1:
        raise ValueError("Cannot include average energies for more than 1 collection")

    # Initializing plot
    fig, ax1 = plt.subplots(figsize=settings.figsize)

    # Initializing legends for different states and sequences
    handles = []

    # Legend for Homo-Lumo gap and energies
    if args.include_energies:
        handles.append(Line2D([0], [0], marker="s", label="Gap", c=settings.colors[0], ls="", ms=5))
        handles.append(Line2D([0], [0], marker="s", label="Energy", c=settings.colors[1], ls="", ms=5))

    # Looping over the collections to plot the average of each sequence
    for index, collection in enumerate(args.input):

        color = settings.colors[index]

        # Finding each sequence in collection
        sequences = [os.path.join(collection, direc) for direc in os.listdir(collection) if
                     direc.startswith(("calculations", "configurations"))]

        # Initializing data arrays
        nsequences = len(sequences)
        waters = np.zeros(nsequences, dtype=int)
        avg_gaps = np.zeros(nsequences)
        std_gaps = np.zeros(nsequences)
        avg_energies = np.zeros(nsequences)
        std_energies = np.zeros(nsequences)

        # Looping over each sequence in collection
        for sub_index, calc_sequence in enumerate(sequences):
            sequence = calculation_sequence(calc_sequence, accepted_file,
                                            select_points=args.select_points)
            waters[sub_index] = sequence.nwater
            avg_gaps[sub_index] = np.average(sequence.gaps)
            std_gaps[sub_index] = np.std(sequence.gaps)
            avg_energies[sub_index] = np.average(sequence.energies)
            std_energies[sub_index] = np.std(sequence.energies)

        # Ordering waters and avg gaps in ascending water order
        order = np.argsort(waters, kind="heapsort")
        waters = waters[order]
        avg_gaps = avg_gaps[order]
        std_gaps = std_gaps[order]
        avg_energies = avg_energies[order]
        std_energies = std_energies[order]

        # Plotting the average gaps with connecting lines
        ax1.plot(waters, avg_gaps, color=color, marker="o", ls="dashed", lw=1, ms=3)

        # Plotting error bars in the average gaps
        ax1.errorbar(waters, avg_gaps, yerr=std_gaps, color=color, alpha=0.75, lw=0.4,
                     ls="", capsize=2, markeredgewidth=0.4, zorder=1)

        # Plotting average energies and standard deviations when requested
        if args.include_energies:
            ax2 = ax1.twinx()

            ax2.plot(waters, avg_energies, color=settings.colors[1], marker="o",
                     ls="dashed", lw=1, ms=3)

            ax2.errorbar(waters, avg_energies, yerr=std_energies, color=settings.colors[1],
                         alpha=0.75, lw=0.4, ls="", capsize=2, markeredgewidth=0.4, zorder=1)

        # Appending collection label to legend
        if not args.include_energies:
            if args.labels:
                label = args.labels[index]
            else:
                basename = os.path.basename(os.path.normpath(collection))
                label = basename.replace("-", " ").replace("_", " ")

            handles.append(Line2D([0], [0], marker="s", label=label, c=color, ls="", ms=5))

    # Naming plot axis
    fs = settings.axes_size
    ts = settings.tick_size

    ax1.set_xlabel("Number of H2O molecules", fontsize=fs)
    ax1.set_ylabel("Average Homo-Lumo gap [eV]", fontsize=fs)
    if args.include_energies:
        ax2.set_ylabel("Average total energy [eV]", fontsize=fs)
        offset = ax2.get_yticks()[1]
        ax2.ticklabel_format(axis="y", style="plain", useOffset=offset, useMathText=True)
        ax2.yaxis.get_offset_text().set(fontsize=ts)

    # Specifying axis ticks
    ax1.tick_params(axis="both", which="major", labelsize=ts)
    if args.include_energies:
        ax2.tick_params(axis="both", labelsize=ts)
    axes = plt.gca()
    axes.xaxis.set_major_locator(ticker.MaxNLocator(integer=True))

    # Including legend
    if (not args.legend_none and len(args.input)) > 1 or args.include_energies:
        plt.legend(handles=handles, fontsize=fs, loc=args.legend_position, frameon=False,
                   ncol=args.legend_columns, bbox_to_anchor=args.legend_position_coords)

    plt.tight_layout()

    # Saving plot
    if args.plot_name:
        plt.savefig(args.plot_name)

    plt.show()
