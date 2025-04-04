# Standard libraries
import os
import re
import time
from typing import List

# Third-party libraries
import numpy as np

# Local libraries
from support import find_paths
from support import get_connections
from support import get_data
from support import get_frame_number


class calculation_sequence:
    """
    This class represents a sequence of single-point calculations on frames of a single
    MD trajectory

    This class is created for plotting purposes. This means it's purpose is to store data
    in a clear manner.

    Instances of this class have the following attributes:

    * nwater -- the number of water molecules present in each single-point calculation
    * frames -- list of MD trajectory frame numbers (ascending order)
    * gaps -- list of homo-lumo gaps (corresponding to frames)
    * energies -- list of total energies (corresponding to frames)
    * connections -- list of connection tables
    * states -- list of indices (of the ordered frames) for each connection table
    * avg_gaps -- list of average homo-lumo gaps (for each state)
    * avg_energies -- list of average total energies (for each state)
    * std_gaps -- list of standard deviations in homo-lumo gaps (for each state)
    * std_energies -- list of standard deviations in total energies (for each state)
    """

    def __init__(self, foldername: str, accepted_outputfiles: List[str], ascending_energies=False,
                 separate_states=False, select_points=None):

        if not foldername.endswith("/"):
            foldername += "/"

        # Finding the number of water molecules
        basename = os.path.basename(os.path.normpath(foldername))
        try:
            self.nwater = int(re.search(r"\_(.*?)\_", basename).group(1))
        except AttributeError:
            self.nwater = int(re.search(r"\-(.*?)\-", basename).group(1))

        # Getting list of directories with calculations (ascending integer order)
        calculations = np.array([direc for direc in os.listdir(foldername) if
                                 os.path.isdir(foldername + direc) and direc.isnumeric()])

        int_calculations = np.array([int(direc) for direc in calculations])
        order_calculations = np.argsort(int_calculations, kind="heapsort")
        calculations = calculations[order_calculations]

        # Initializing frames, gaps and energies
        ncalculations = len(calculations)
        self.frames = np.zeros(ncalculations, dtype=int)
        self.gaps = np.zeros(ncalculations)
        self.energies = np.zeros(ncalculations)

        # Initializing connections and states
        self.connections = []
        self.states = []

        # Looping over every calculation the sequence
        output_file, fragment_path = None, None
        for index, calculation in enumerate(calculations):

            # Finding the single-point output file and path to fragment
            if index == 0:
                output_file, fragment_path = find_paths(foldername + calculation, accepted_outputfiles)

            # Getting the homo-lumo gap and energy
            gap, energy = get_data(foldername + calculation + output_file)
            self.gaps[index] = gap
            self.energies[index] = energy

            # Getting the frame number
            try:
                self.frames[index] = get_frame_number(foldername + calculation + fragment_path)
            except NameError:
                self.frames[index] = int(calculation)

            if separate_states:
                # Finding the connection table of the fragment
                connection_table = get_connections(foldername + calculation + fragment_path)
                if connection_table not in self.connections:
                    self.connections.append(connection_table)
                    self.states.append([])

                # classifying the fragment state based on connection table
                state_index = self.connections.index(connection_table)
                self.states[state_index].append(index)

        # Ordering the (frames), gaps and energies in ascending frame or ascending energy order
        if ascending_energies or select_points:
            order = np.argsort(self.energies, kind="heapsort")
            self.frames = np.linspace(0, 1, num=len(self.energies))
        else:
            order = np.argsort(self.frames, kind="heapsort")
            self.frames = self.frames[order]

        self.gaps = self.gaps[order]
        self.energies = self.energies[order]

        # Keeping only the number of selected points with the lowest energies
        if select_points:
            self.frames = self.frames[:select_points]
            self.gaps = self.gaps[:select_points]
            self.energies = self.energies[:select_points]

        if separate_states:
            # Updating the indices corresponding to each state using the sorted order
            for index, state in enumerate(self.states):
                ordered_state = np.argsort(order, kind="heapsort")[state]
                self.states[index] = ordered_state

            # Calculating the averages and standard deviations for the gaps and energies in each state
            nstates = len(self.states)
            self.avg_gaps = np.zeros(nstates)
            self.avg_energies = np.zeros(nstates)
            self.std_gaps = np.zeros(nstates)
            self.std_energies = np.zeros(nstates)

            for index, state in enumerate(self.states):
                self.avg_gaps[index] = np.average(self.gaps[state])
                self.avg_energies[index] = np.average(self.energies[state])
                self.std_gaps[index] = np.std(self.gaps[state])
                self.std_energies[index] = np.std(self.energies[state])
