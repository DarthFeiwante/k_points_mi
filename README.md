
The program is intended for calculations of k-points density in the reciprocal space, 
per reciprocal angstrom, from given grid or k-points grid from given density.
Vectors of the reciprocal space are calculated from ones of the real space.

Requirements: python3

Program uses following moduli:

sympy.matrices, math

There is class Reciprocal_Lattice, which include following methods:
1) k_points(rprim_list, density_kpt_string)
  Its input parameters:
  a) rprim_list - list of vectors (each vector is list of dimensionless coordinates) of the elemental cell 
  (e.g. [[0.0, 0.5, 0.5], [0.5, 0.0, 0.5], [0.5, 0.5, 0.0]]).

  b) density_kpt_string - string of the required k-points densities for each vector of reciprocal space per \AA$^{-1}$ 
  (e.g. '8 8 8')
  
  Results of the method are: 
    1) k-points grid parameters, corresponding to the given k-points density, written in attribute 'self.ngkpt'
    2) coordinates of the reciprocal space vectors, written in attribute 'self.recip_vecs'
    3) volume of real space elemental cell, written in attribute 'self.volume'

2) k_density(rprim_list, ngkpt_string)}
  Its input parameters:
  a) rprim_list - list of vectors (each vector is list of dimensionless coordinates) of the elemental cell 
  (e.g. [[0.0, 0.5, 0.5], [0.5, 0.0, 0.5], [0.5, 0.5, 0.0]]).

  b) ngkpt_string - string of k-points grid, e.g. '8 8 8', that means 8x8x8 grid in reciprocal space
        
  Results of the method are: 
    1) list of reciprocal vectors lengths in reciprocal angstrom units, written in attribute 'self.recip_len_list',
    2) list of k-points densities for each vector per reciprocal angstrom, written in attribute 'self.list_den'
    \end{itemize}

Input file k_points_mi.in (see example) contains following parameters:
1) rprim - dimensionless real space lattice vectors of elemental cell
2) acell - lattice parameters in angstroms (multipliers for dimensionless vectors)
3) regime - type of program work: 'density' - calculations of k-points density; 'ngkpt' - calculations
of the k-points grid

Program can be used as executable or importable 
You should input additional data corresponding to prompts if program is used as executable 
