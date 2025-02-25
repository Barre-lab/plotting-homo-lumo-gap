#!/usr/bin/env amspython

# Standard libraries
import os
import re
import time
import argparse
from typing import Tuple

# Third-party libraries
import scm.plams
import numpy as np
import matplotlib.pyplot as plt
import matplotlib


def parse_arguments():
    """
    Parsing arguments given to script
    """

    parser = argparse.ArgumentParser()

    parser.add_argument(
            "inputfolders",
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


def get_data(inputfolder: str, separate_states=True) -> Tuple[np.ndarray]:
    """
    Gets frame number, homo-lumo gap and total energy for each calculation in the input folder

    If separate_states is True, the indices of calculations where the configuration of
    the fragment has changed w.r.t the fragment in the first frame are found and returned
    """

    def _find_paths(directory, filenames, str_fragment=".xyz"):
        path_to_output = []
        path_to_fragment = []

        for dir_path, _, file_names in os.walk(directory):
            for file in file_names:
                if file in filenames:
                    path_to_output.append(os.path.join(dir_path, file))
                if file.endswith(str_fragment) and not file.startswith("output"):
                    path_to_fragment.append(dir_path)

        if not path_to_output or not path_to_fragment:
            raise ValueError(f"No path to outputfile or fragment in {directory}")
        if len(path_to_output) > 1 or len(path_to_fragment) > 1:
            raise ValueError(f"More than 1 file or frame found: {outputfiles}, {frames}")

        return path_to_output[0].replace(directory, ""), path_to_fragment[0].replace(directory, "")

    def get_connections(fragment_file):
        """
        Finds connection table of fragment file: List with neighbors for each fragment atom
        """

        molecule = scm.plams.Molecule(fragment_file)
        molecule.guess_bonds()

        components = molecule.separate()
        fragment = scm.plams.Molecule()
        for component in components:
            if component.get_formula() != "H2O":
                fragment += component

        return fragment.get_connection_table()

    # List of possible output files
    outputfiles = ["single-point.out", "out"]

    # Making sure each folder ends with /
    if not inputfolder.endswith("/"):
        inputfolder += "/"

    # Getting list of directories with calculations
    directories = [direc for direc in os.listdir(inputfolder) if os.path.isdir(inputfolder + direc)]
    directories.sort()

    # Creating arrays for the standard variables
    gaps = np.zeros(len(directories))
    energies = np.zeros(len(directories))
    frames = np.zeros(len(directories), dtype=int)

    # Creating list for indices in which the fragment configuration has changed
    if separate_states:
        newstates = []

    file_path, fragment_path = None, None

    # Looping over each folder within the given input folder
    for index, direc in enumerate(directories):

        # Finding paths to output files and fragment files, and initial configuration
        if index == 0:
            file_path, fragment_path = _find_paths(inputfolder + direc, outputfiles)

            # Finding configuration of the first fragment
            for file in os.listdir(inputfolder + direc + fragment_path):
                if file.endswith(".xyz"):
                    full_fragment_path = inputfolder + direc + fragment_path + "/"
                    initial_connections = get_connections(full_fragment_path + file)
                    break

            if not initial_connections:
                raise ValueError("Unable to find initial fragment connection table")

        # Finding HOMO-LUMO gap and total energy (in eV)
        with open(inputfolder + direc + file_path, "r") as file:
            gap, energy = False, False
            for line in file:
                if line.startswith("HOMO-LUMO"):
                    gaps[index] = line.split()[3]
                    gap = True
                if line.startswith("Total Energy (eV)"):
                    energies[index] = line.split()[3]
                    energy = True
                if gap and energy:
                    break

        # Finding frame number
        for file in os.listdir(inputfolder + direc + fragment_path):
            if file.endswith(".xyz"):
                frames[index] = re.sub(r'[^\d]+', '', file)

                # Comparing fragment configuration with initial configuration
                if separate_states:
                    full_fragment_path = inputfolder + direc + fragment_path + "/" + file
                    if get_connections(full_fragment_path) != initial_connections:
                        newstates.append(index)
                break

    # Ordering the arrays in ascending frame order
    order = np.argsort(frames, kind="heapsort")
    if order[0] != 0:
        raise ValueError("Initial fragment does not correspond to first frame")

    frames = frames[order]
    gaps = gaps[order]
    energies = energies[order]

    # Updating the list of new state indices using the new (sorted) frame order
    if separate_states:
        newstates = np.argsort(order, kind="heapsort")[newstates]
    else:
        newstates = None

    return frames, gaps, energies, newstates


def get_average_data(inputsequence: str, separate_states=False) -> Tuple[np.ndarray]:
    """
    Computes averages and standard deviations of homo-lumo gaps and energies for all calculation
    folders in calculation sequence
    """

    # Getting list with all calculation folders in input sequence
    calculation_folders = [direc for direc in os.listdir(inputsequence)
                           if direc.startswith("calculations")]

    # Creating arrays for the variables
    nwaters = np.zeros(len(calculation_folders), dtype=int)
    avg_gaps_init = np.zeros(len(calculation_folders))
    avg_gaps_new = np.zeros(len(calculation_folders))
    stdv_gaps_init = np.zeros(len(calculation_folders))
    stdv_gaps_new = np.zeros(len(calculation_folders))

    # Looping over each calculation folder within the sequence to obtain averages
    for index, folder in enumerate(calculation_folders):

        # Getting number of water molecules
        nwaters[index] = re.search(r"\_(.*?)\_", folder).group(1)

        # Getting homo-lumo gaps and energies of all calculations
        _, gaps, _, newstates = get_data(inputsequence + folder)

        # Separating the states
        if separate_states:

            # The initial state
            avg_gaps_init[index] = np.average(np.delete(gaps, newstates))
            stdv_gaps_init[index] = np.std(np.delete(gaps, newstates))

            # The new state
            if len(newstates) < 1:
                avg_gaps_new[index], stdv_gaps_new[index] = np.nan, np.nan
            else:
                avg_gaps_new[index] = np.average(gaps[newstates])
                stdv_gaps_new[index] = np.std(gaps[newstates])

        else:
            avg_gaps_init[index] = np.average(gaps)
            stdv_gaps_init[index] = np.std(gaps)

    # Ordering the arrays in ascending order for number of water molecules
    order = np.argsort(nwaters, kind="heapsort")

    nwaters = nwaters[order]
    avg_gaps_init = avg_gaps_init[order]
    stdv_gaps_init = stdv_gaps_init[order]

    if separate_states:
        avg_gaps_new = avg_gaps_new[order]
        stdv_gaps_new = stdv_gaps_new[order]
    else:
        avg_gaps_new, stdv_gaps_new = None, None

    return nwaters, avg_gaps_init, stdv_gaps_init, avg_gaps_new, stdv_gaps_new


def get_moving_average(values: np.ndarray, window: int) -> np.ndarray:
    """
    Calculates and return moving average of homo-lumo gaps, energies or other quantity
    as a function of the frame numbers
    """

    # Initializing array with nan's
    moving_average = np.zeros(len(values))
    moving_average.fill(np.nan)

    # Fill array with values of moving average
    for index in range(len(values) - window + 1):
        min_index = max(0, index - window)
        max_index = min(len(values), index + window + 1)
        moving_average[index] = np.mean(values[min_index:max_index])

    return moving_average


def plot_averages(args: argparse.Namespace, colors: list) -> None:
    """
    Plots average homo-lumo gaps and energies of all calculation folders in given calculation
    sequences against the corresponding number of water molecules
    """

    # Defining input sequences
    inputsequences = args.inputfolders

    # Initializing plot
    plt.figure(figsize=(4.5, 3.4))

    # Defining the labels
    if args.labels:
        labels = [" ".join(label.split("_")) for label in args.labels]
    else:
        labels = [" ".join(sequence.split("-")).replace("/", "") for sequence in inputsequences]

    # Creating handles for legend
    handles = []

    # Adding labels for initial and new states to legend
    if args.separate_states:
        label_initial = matplotlib.lines.Line2D([0], [0], marker="", label="Initial state",
                                                ms=3, c="black", ls="dashed")
        label_new = matplotlib.lines.Line2D([0], [0], marker="", label="New state",
                                            ms=3, c="black", ls="dotted")

        handles.append(label_initial)
        handles.append(label_new)

    # Looping over each given sequence of calculation folders to plot
    for sequence, label, color in zip(inputsequences, labels, colors):

        # Making sure the directory name ends with /
        if not sequence.endswith("/"):
            sequence += "/"

        # Getting the average data
        nwaters, avg_gaps_init, stdv_gaps_init, avg_gaps_new, stdv_gaps_new = get_average_data(
                sequence, args.separate_states)

        # Appending label for sequence to legend
        if len(args.inputfolders) > 1:
            label_sequence = matplotlib.lines.Line2D([0], [0], marker="s", label=label, ms=5,
                                                     color=color, ls="")
            handles.append(label_sequence)

        if args.separate_states:

            # Plotting the average gaps with connecting lines
            plt.plot(nwaters, avg_gaps_init, color=color, marker="o", ls="dashed", lw=1, ms=3)

            plt.plot(nwaters[~np.isnan(avg_gaps_new)], avg_gaps_new[~np.isnan(avg_gaps_new)],
                     color=color, marker="v", ls="dotted", lw=1, ms=3)

            # Plotting error bars in average gaps
            plt.errorbar(nwaters, avg_gaps_init, yerr=stdv_gaps_init, color=color,
                         alpha=0.75, lw=0.4, ls="", capsize=2, markeredgewidth=0.4)

            plt.errorbar(nwaters, avg_gaps_new, yerr=stdv_gaps_new, color=color,
                         alpha=0.75, lw=0.4, ls="", capsize=2, markeredgewidth=0.4)
        else:

            # Plotting the average gaps with connecting lines
            plt.plot(nwaters, avg_gaps_init, color=color, marker="o", ls="dashed", lw=1, ms=3)

            # Plotting error bars in average gaps
            plt.errorbar(nwaters, avg_gaps_init, yerr=stdv_gaps_init, color=color,
                         alpha=0.75, lw=0.4, ls="", capsize=2, markeredgewidth=0.4)

    # Naming plot axis
    plt.xlabel("Number of H2O molecules", fontsize=9)
    plt.ylabel("Average Homo-Lumo gap [eV]", fontsize=9)

    # Specifying axis ticks
    plt.tick_params(axis="both", which="major", labelsize=8)
    axes = plt.gca()
    axes.xaxis.set_major_locator(matplotlib.ticker.MaxNLocator(integer=True))

    # Including legend
    plt.legend(handles=handles, fontsize=9, loc=args.legend_position, frameon=False,
               ncol=args.legend_columns, bbox_to_anchor=args.legend_position_coords)

    plt.tight_layout()

    # Saving plot
    if args.plot_name:
        plt.savefig(args.plot_name)

    plt.show()


def plot_all_together(args: argparse.Namespace, colors: list) -> None:
    """
    Plots homo-lumo gaps and energies as a function of the frame numbers for all input calculation
    folders together in one graph
    """

    # In case input is a sequence, all calculation inside are taken
    if args.type == "sequence":
        if len(args.inputfolders) > 1:
            raise ValueError("More than 1 sequence of calculation folders given, cannot plot all")

        calculation_folders = [args.inputfolders[0] + direc for direc in
                               os.listdir(args.inputfolders[0]) if direc.startswith("calculations")]
    elif args.type == "calculations":
        calculation_folders = args.inputfolders
    else:
        raise ValueError(f"Unknown input type: {args.type}")

    # Plotting all data in a single graph
    plt.figure(figsize=(4.5, 3.4))

    # Defining the labels
    if args.labels:
        labels = [" ".join(label.split("_")) for label in args.labels]
    else:
        if args.type == "sequence":
            labels = [folder.replace(args.inputfolders[0], "").replace("calculations_", "").
                      replace("_", " ") for folder in calculation_folders]
        else:
            labels = ["".join(folder.rsplit("/", 1)).replace("calculations_", " ").replace("_", " ")
                      for folder in calculation_folders]

    for folder, label, color in zip(calculation_folders, labels, colors):

        # Getting data from all calculations in folder
        frames, gaps, energies, newstates = get_data(folder, args.separate_states)

        if args.moving_average:
            averaged_gaps = get_moving_average(gaps, args.moving_average)

            # Plotting both moving averages and original data
            plt.plot(frames, gaps, color=color, alpha=0.1, lw=0.8)
            plt.plot(frames, averaged_gaps, color=color, lw=0.8, label=label + " (MA)")

        elif args.separate_states:
            # Plotting the connecting lines
            plt.plot(frames, gaps, color=color, lw=0.2, zorder=1)

            # Plotting data of the initial state
            plt.plot(np.delete(frames, newstates), np.delete(gaps, newstates), color=color, ls="",
                     fillstyle="none", marker="o", ms=2, label=label + " initial state", alpha=0.5)

            # Plotting data of the new state
            plt.plot(frames[newstates], gaps[newstates], color=color, ls="", marker="o",
                     ms=4, fillstyle="full", label=label + " new state", zorder=2)

        else:
            plt.plot(frames, gaps, color=color, lw=0.8, label=label)

    # Naming plot axis
    plt.xlabel("Frame number")
    plt.ylabel("Homo-Lumo gap [eV]")

    # Including legend
    if args.separate_states and len(args.inputfolders) == 1:
        handles_states = []

        point1 = matplotlib.lines.Line2D([0], [0], marker="o", label="Initial state",
                                         ms=1, c=colors[0], fillstyle="none", ls="")
        point2 = matplotlib.lines.Line2D([0], [0], marker="o", label="New state",
                                         c=colors[0], ls="", ms=3)
        handles_states.append(point1)
        handles_states.append(point2)

        plt.legend(handles=handles_states, fontsize=9, loc=args.legend_position, frameon=False)

    elif len(args.inputfolders) > 1:
        plt.legend(loc=args.legend_position, frameon=False)

    plt.tight_layout()

    # Saving plot
    if args.plot_name:
        plt.savefig(args.plot_name)

    plt.show()


def plot_all_separate(args: argparse.Namespace, colors) -> None:
    """
    Plots homo-lumo gaps and energies as a function of the frame numbers for all input calculation
    folders in separate graphs
    """

    def _plot_folder_data_in_subplots(folder, axs, list_all_waters, color, sequence="") -> None:

        # Finding number of water molecules and the corresponding plot index
        nwater = int(re.search(r"\_(.*?)\_", folder).group(1))
        ax_index = list_all_waters.index(nwater)

        # Getting data
        frames, gaps, energies, newstates = get_data(sequence + folder)

        # Correcting subplot indices in case axis are 2D
        if axs.ndim > 1:
            ax_index = (ax_index // 2, ax_index % 2)

        # Plotting gaps
        if args.moving_average:
            averaged_gaps = get_moving_average(gaps, args.moving_average)

            axs[ax_index].plot(frames, gaps, color=color, lw=0.8, alpha=0.1)
            axs[ax_index].plot(frames, averaged_gaps, color=color, lw=1.0)

        elif args.separate_states:
            # Plotting the connecting lines
            axs[ax_index].plot(frames, gaps, color=color, lw=0.4, zorder=1, alpha=0.6)

            # Plotting data of the initial state
            axs[ax_index].plot(np.delete(frames, newstates), np.delete(gaps, newstates),
                               color=color, ls="", marker="", fillstyle="none", ms=1, alpha=0.5)

            # Plotting data of the new state
            axs[ax_index].plot(frames[newstates], gaps[newstates], color=color, ls="",
                               marker="o", ms=3, fillstyle="full", zorder=2)

        else:
            axs[ax_index].plot(frames, gaps, color=color, lw=0.8)

    # In case a collection of sequences of calculation folders is provided
    if args.type == "sequence":

        list_all_waters = []

        # Looping over all sequences
        for sequence in args.inputfolders:

            # Finding the list of all number of water molecules in the sequence
            list_waters = [re.search(r"\_(.*?)\_", folder).group(1) for folder in
                           os.listdir(sequence) if folder.startswith("calculations")]

            # Appending all amounts of water if not already in the list
            for nwater in list_waters:
                if int(nwater) not in list_all_waters:
                    list_all_waters.append(int(nwater))

        # Sorting list with all waters in ascending order
        list_all_waters.sort()

        # Defining the labels
        if args.labels:
            labels = args.labels
        else:
            labels = [sequence.replace("/", "").replace("-", " ") for sequence in args.inputfolders]
        if args.moving_average:
            labels = [label + " (MA)" for label in labels]

        # Initializing figure with sub-figures
        nrows = -(-len(list_all_waters) // 2)
        fig, axs = plt.subplots(nrows=nrows, ncols=2, figsize=(6, 2*nrows))

        # Creating handles for legend
        handles = []

        # Looping over each sequence to plot the data of all the folders inside
        for sequence, label, color in zip(args.inputfolders, labels, colors):

            # Getting list of all calculation folders
            folders = [folder for folder in os.listdir(sequence)
                       if folder.startswith("calculations")]

            # Plotting data of each folder in correct subplot
            for folder in folders:
                _plot_folder_data_in_subplots(folder, axs, list_all_waters, color, sequence)

            # Constructing legend
            point = matplotlib.lines.Line2D([0], [0], marker="s", label=label, ms=5, color=color, ls="")
            handles.append(point)

        # Adding legend to plot
        if axs.ndim > 1:
            ax_index = (args.legend_index // 2, args.legend_index % 2)
        else:
            ax_index = args.legend_index

        # Adding legend for the separated states
        if args.separate_states and len(args.inputfolders) == 1:
            handles_states = []

            point1 = matplotlib.lines.Line2D([0], [0], marker="o", label="Initial state",
                                             ms=1, c=colors[0], fillstyle="none", ls="")
            point2 = matplotlib.lines.Line2D([0], [0], marker="o", label="New state",
                                             c=colors[0], ls="", ms=3)
            handles_states.append(point1)
            handles_states.append(point2)

            axs[ax_index].legend(handles=handles_states, fontsize=9, loc=args.legend_position,
                                 frameon=False)

        # Adding legend for the different sequences
        elif len(args.inputfolders) > 1:
            axs[ax_index].legend(handles=handles, fontsize=9, loc=args.legend_position,
                                 frameon=False)

    # In case a collection of calculation folders is provided
    elif args.type == "calculations":

        # Finding all calculation folders from input
        folders = [folder for folder in args.inputfolders if folder.startswith("calculations")]

        # Initializing figure with sub-figures
        nrows = -(-len(folders) // 2)
        fig, axs = plt.subplots(nrows=nrows, ncols=2, figsize=(6, 2*nrows))

        list_all_waters = []
        for folder in folders:
            if folder.count("/") > 1:
                raise ValueError("Plotting folders separatly only works within 1 sequence")

            list_all_waters.append(int(re.search(r"\_(.*?)\_", folder).group(1)))

        # Sorting list with all waters in ascending order
        list_all_waters.sort()

        # Looping over the folders to plot
        for folder in folders:
            _plot_folder_data_in_subplots(folder, axs, list_all_waters, colors[0], sequence="")

    else:
        raise ValueError(f" Unknown input type {args.type}")

    # Defining sub-plot names, axis names and ticks and getting y-ranges
    min_yrange, max_yrange = np.inf, 0
    for ax, nwater in zip(axs.flat, list_all_waters):

        ax.set_ylabel("Homo-Lumo gap [ev]", fontsize=9)
        ax.set_xlabel("Frame number", fontsize=9)
        ax.label_outer()

        ax.xaxis.set_tick_params(labelbottom=True, labelsize=8)
        ax.yaxis.set_tick_params(labelleft=True, labelsize=8)

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

    if args.plot_name:
        plt.savefig(args.plot_name)

    plt.show()


def main(args):
    """
    Main part of the script
    """

    colors = ["blue", "orange", "limegreen", "red"]

    if args.type == "sequence" and args.what == "avg":
        plot_averages(args, colors)

    elif args.what == "all" and args.how == "together":
        plot_all_together(args, colors)

    elif args.what == "all" and args.how == "separate":
        plot_all_separate(args, colors)

    else:
        print("Combination of given arguments is (currently) not supported")


if __name__ == '__main__':
    arguments = parse_arguments()
    main(arguments)
