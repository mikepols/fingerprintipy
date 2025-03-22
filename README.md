# FingerprintiPy

A Python tool to fingerprint atoms in molecules, with the purpose of identifying specific atoms within a larger structure.

## Requirements

A working installation of [`ASE`](https://wiki.fysik.dtu.dk/ase/index.html) is required. The tool has been tested and developed with version 3.22.1.

## Installation

The script can be automatically installed through
```
./install.sh
```
which adds both scripts to the `$PATH`.

## Usage

The utility consists of two scripts: 
1. `fingerprinti.py`: used to obtain atom fingerprints for a molecule using specific parameters
2. `matchi.py`: used to obtain the indices of the atoms with the same atomic fingerprints in larger structures

## Example

Some example files are provided in the (`example`) to demonstrate the usage of the included scripts. Navigate into this folder and execute the following command:
```
fingerprinti.py -i molecule.xyz -o molecule.fp -v
```
The fingerprint of all atoms in the molecule (`-i molecule.xyz`) has been output to both the screen (`-v`) and a file (`-o molecule.fp`). To find all atoms at the tail end of organic compounds, we now match the `H36C6N4` atomic fingerprint to all species, the atoms that match this fingerprint are considered equivalent. To do this we run the following command, where we select the species that can be part of the molecules (`--active C H N`):
```
matchi.py --fp H36C6N4 -i bulk.xyz -o matched_indices.fp --active C H N
```
The relevant indices are ouput to a file (`-o matched_indices.fp`).

## References

The utility was developed for use in the following publication: M. Pols *et al.*, *J. Phys. Chem. Lett.*, 15, 8057-8064 (2024), DOI:  [`10.1021/acs.jpclett.4c01629`](https://doi.org/10.1021/acs.jpclett.4c01629).

## License

The code is available as open source under the terms of the [MIT License](LICENSE).
