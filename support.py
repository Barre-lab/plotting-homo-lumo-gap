# Standard libraries
import os
import re
import time

# Third-party libraries
import scm.plams
import numpy as np


def find_paths(directory, filenames, str_fragment=".xyz"):
    """
    Finds the paths to the out file of a single-point calculation (any file in 'filenames')
    and the fragment file within in the input directory.

    Additionally, the function will try to find an engine specific rkf file. If it can not
    find a fragment file, it will instead return this rkf file from which the molecule
    geometry can be extracted.
    """

    path_to_output = []
    path_to_fragment = []
    path_to_rkf = []

    for dir_path, _, file_names in os.walk(directory):
        for file in file_names:
            if file in filenames:
                path_to_output.append(os.path.join(dir_path, file))
            if file.endswith(str_fragment) and not file.startswith("output"):
                path_to_fragment.append(dir_path)
            if file in ["dftb.rkf", "adf.rkf"]:
                path_to_rkf.append(os.path.join(dir_path, file))

    if not path_to_output:
        raise ValueError(f"No path to outputfile found in: {directory}")
    if not path_to_fragment:
        if path_to_rkf:
            path_to_fragment.append(path_to_rkf[0])
        else:
            raise ValueError(f"No path to fragment file or rkf file in: {directory}")

    if len(path_to_output) > 1:
        for output_file in path_to_output.copy():
            if "optimization" in output_file:
                path_to_output.remove(output_file)

        if len(path_to_output) > 1:
            raise ValueError(f"More than 1 path to output file found: {path_to_output}")

    if len(path_to_fragment) > 1:
        raise ValueError(f"More than 1 path to fragment file found: {path_to_fragment}")

    return path_to_output[0].replace(directory, ""), path_to_fragment[0].replace(directory, "")


def get_connections(fragment_path):
    """
    Finds connection table of scm.plams.Molecule (excluding water) in fragment path.
    The fragment path can either point to an rkf or an xyz file.

    Connection table: list with list of neighbors for each fragment atom.
    """

    fragment_file = None

    # Finding the rkf or xyz file
    if os.path.isfile(fragment_path) and fragment_path.endswith(".rkf"):
        job = scm.plams.AMSJob.load_external(fragment_path)
        molecule = job.results.get_main_molecule()
    else:
        for file in os.listdir(fragment_path):
            if file.endswith(".xyz"):
                fragment_file = os.path.join(fragment_path, file)
                molecule = scm.plams.Molecule(fragment_file)

    if not molecule:
        raise FileNotFoundError(f"Could not find fragment file in: {fragment_path}")

    # Removing any water molecules from the fragment
    molecule.guess_bonds()
    components = molecule.separate()
    fragment = scm.plams.Molecule()

    for component in components:
        if component.get_formula() != "H2O":
            fragment += component

    return fragment.get_connection_table()


def get_data(output_file):
    """
    Gets the homo-lumo gap and total energy from the given output file.
    The output file must be a relative path to (and including) the file.
    """

    gap, energy = None, None

    with open(output_file) as file:
        found_gap, found_energy = False, False
        for line in file:
            if line.startswith("HOMO-LUMO"):
                gap = line.split()[3]
                found_gap = True
            if line.startswith("Energy (hartree)"):
                energy = float(line.split()[2]) * 27.211407953
                found_energy = True
            if found_gap and found_energy:
                break

    if not gap or not energy:
        raise NameError(f" Could not find gap or energy in {output_file}")

    return gap, energy


def get_frame_number(fragment_path):
    """
    Finds the frame number from the fragment file that is found with
    the given fragment path.
    """

    frame = None

    if os.path.isfile(fragment_path) and fragment_path.endswith(".rkf"):
        raise NameError(f"Cannot get frame number from .rkf: {fragment_path}")

    for file in os.listdir(fragment_path):
        if file.endswith(".xyz"):
            frame = re.sub(r'[^\d]+', '', file)

    if not frame:
        raise NameError(f"Could not find frame from fragment: {fragment_path}")

    return frame


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


def find_input_type(args_input: list):
    """
    Determines if the input consists of sequences of calculations or collections
    of these sequences.
    """

    check_folders = np.zeros(len(args_input))
    for index, folder in enumerate(args_input):
        basename = os.path.basename(os.path.normpath(folder))

        if basename.startswith(("calculations", "configurations")):
            check_folders[index] = True
        else:
            check_folders[index] = False

    if len(set(check_folders)) != 1:
        raise ValueError("Input is a mix of calculation sequences and collection of sequences")

    if all(check_folders):
        input_type = "sequence"
    else:
        input_type = "collection"

    return input_type
