MaterialsPES
============

A set of Python scripts for generating a potential energy surface (PES) for the shearing of a material using the PWSCF program within Quantum ESPRESSO. It does this by defining two slabs in the simulation cell and calculates the energy of a grid of points representing one slab moving relative to the other. The outputs from PWSCF are read and converted into a format gnuplot can read to make a graphic representing the PES.

pes_builder.py - This script takes an input file containing a set of atomic coordinates and lattice vectors, definition of the potential energy surface region and resolution, and some information required by Quantum ESPRESSO. It calculates the movement of the atoms and lattice vectors for each point on the PES and writes a PWSCF input file for each point. Also, writes an XYZ file for easy visualization of each point on the PES. An example of the input file, calculating the shear of an alumina cell, can be found in aluminaPES_example.

results_grabber_pes.py - This script reads the Quantum ESPRESSO output files generated from running the PWSCF calculations set up by the previous script. Calculates the energies of each point relative to the lowest energy point and generates a formatted table of energies and an input file for gnuplot to make a graphic of the PES.

Compatible with Python 2.6/2.7
