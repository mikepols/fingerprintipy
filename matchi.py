#!/usr/bin/env python3

import argparse

import numpy as np
from scipy import sparse

import ase
from ase.io import read
import ase.neighborlist as nl

def main():
    """
    Main function for analyzing molecular fingerprint
    """
    args = parse_command_line_arguments()
    input_file, output_file = args.input, args.output
    cutoff_tolerance = args.tolerance
    N_mps = args.N_mps

    # Get structural and molecular fingerprint
    s = read(input_file) # read structure
    active_idxs = np.array([atom.index for atom in s if atom.symbol in args.active]) # indices of active elements
    s_active = s[active_idxs] # active selection of structure
    n_molecules, molecular_idxs, molecular_fps = get_molecular_information(s_active, N_mps, cutoff_tolerance) # get molecular fingerprint
    symbols = s.get_chemical_symbols() # get elements

    # Selected indices matching fingerprint
    rows, cols = np.where(molecular_fps == args.fp) # match fingerprints
    matched_idxs = molecular_idxs[rows, cols] # selected indices
    output_matched_indices(matched_idxs, output_file) # write matched indices to file
 
def parse_command_line_arguments():
    """
    Get command-line arguments
    """
    parser = argparse.ArgumentParser(description="Analyze the chirality of metal halide perovskite structures or trajectories.")

    # Main arguments
    parser.add_argument('--fp', help='Atomic fingerprint to match in the structure', required=True, type=str)
    parser.add_argument('-i', '--input', help='Structure to analyze the atomic fingerprint of', required=True, type=str)
    parser.add_argument('-o', '--output', help='File to output the molecular fingerprint to', required=False, type=str, default='matched_indices.fp')

    # Optional arguments
    parser.add_argument('--N_mps', help='Number of message passing steps through the molecule', required=False, type=int, default=10)
    parser.add_argument('--tolerance', help='Structure to analyze the atomic fingerprint of', required=False, type=float, default=0.1)
    parser.add_argument('--active', help='Active elements for organic compounds', required=False, type=str, nargs='+', default=['H', 'C', 'N', 'O'])
    parser.add_argument('-v', '--verbose', help='A boolean switch for verbose output', action='store_true', default=False)

    return parser.parse_args()

def output_matched_indices(matched_idxs, output_file):
    """
    Output zero-based indices that match the atomic fingerprint to a file
    """
    np.savetxt(output_file, matched_idxs, fmt='%.0f') # output column vector of indices to file

def get_molecular_information(s, N_mps, cutoff_tolerance):
    """
    Construct arrays with molecular indices and fingerprints from a structure
    """
    cutoffs = np.array(nl.natural_cutoffs(s)) + cutoff_tolerance # generate standard cutoffs
    nlist = nl.build_neighbor_list(s, cutoffs=cutoffs, self_interaction=False) # create neighborlist

    cmat = nlist.get_connectivity_matrix(sparse=False) # connectivity matrix
    dmat = nl.get_distance_matrix(cmat, limit=N_mps).toarray() # distance matrix
    n_molecules, idxs_molecules = sparse.csgraph.connected_components(cmat) # determine molecules from graph

    # Iterate over all molecules
    molecular_idxs, molecular_fps = [], []
    for idx in range(n_molecules):
        mol_idxs = np.array([i for i in range(len(idxs_molecules)) if idxs_molecules[i] == idx], dtype=int) # indices of species in molecule
        molecular_idxs.append(mol_idxs) # store molecular indices
        molecular_fps.append(get_molecular_fingerprint(s, dmat, mol_idxs)) # store molecular fingerprints
    molecular_idxs, molecular_fps = np.array(molecular_idxs), np.array(molecular_fps)

    return n_molecules, molecular_idxs, molecular_fps

def get_molecular_fingerprint(s, dmat, atom_idxs):
    """
    Determine a string-based molecular fingerprint
    """
    elements = np.array(s.get_chemical_symbols())[atom_idxs] # elements in structure
    unique_elements = np.array(list(set(elements))) # unique elements in structure

    molecular_dmat = dmat[atom_idxs, :][:, atom_idxs]

    atomic_fingerprint = []
    for idx in range(len(atom_idxs)):
        distances = molecular_dmat[idx, :]
        fp = {}
        for element in unique_elements:
            active_idxs = np.where(elements == element)[0]
            active_distances = distances[active_idxs]
            fp[element] = np.sum(active_distances)
        atomic_fingerprint.append(fingerprint2string(fp))

    return atomic_fingerprint

def fingerprint2string(fp):
    """
    Convert the fingerprint dictionary to a string
    """
    ptable = ase.data.chemical_symbols
    
    s = ''
    for atom in ptable:
        if atom in fp:
            s += '{:}{:d}'.format(atom, fp[atom])

    return s

if __name__ == "__main__":
    main()
