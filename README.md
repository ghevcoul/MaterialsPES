MaterialsPES
============

A set of Python scripts for generating a potential energy surface (PES) for the shearing of a material using the PWSCF program within Quantum ESPRESSO. It does this by defining two slabs in the simulation cell and calculates the energy of a grid of points representing one slab moving relative to the other. The outputs from PWSCF are read and converted into a format gnuplot can read to make a graphic representing the PES.
