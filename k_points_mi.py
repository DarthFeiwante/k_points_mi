#!/usr/bin/env python3
'''
\section{k\_points\_mi.py}

\\hyperlink{content}{Table of content}

The program is intended for calculations of k-points density in the reciprocal space, 
per reciprocal angstrom, from given grid or k-points grid from given density.
Vectors of the reciprocal space are calculated from ones of the real space.

Requirements: python3

Program uses following moduli:

sympy.matrices, math

There is class Reciprocal\_Lattice, which include following methods:

    \\begin{itemize}

    \item \\textcolor{green}{k\_points(rprim\_list, density\_kpt\_string)}

        Its input parameters:

        \\begin{itemize}

        \item rprim\_list - list of vectors (each vector is list of dimensionless coordinates) of the elemental cell (e.g. [[0.0, 0.5, 0.5], [0.5, 0.0, 0.5], [0.5, 0.5, 0.0]]).

        \item density\_kpt\_string - string of the required k-points densities for each vector of reciprocal space per \AA$^{-1}$ (e.g. '8 8 8')

        \end{itemize}

        Results of the method are: 

        k-points grid parameters, corresponding to the given k-points density, written in attribute 'self.ngkpt'

        coordinates of the reciprocal space vectors, written in attribute 'self.recip\_vecs'

        volume of real space elemental cell, written in attribute 'self.volume'

    \item \\textcolor{green}{k\_density(rprim\_list, ngkpt\_string)}

        Its input parameters:

        \\begin{itemize}

        \item rprim\_list - list of vectors (each vector is list of dimensionless coordinates) of the elemental cell (e.g. [[0.0, 0.5, 0.5], [0.5, 0.0, 0.5], [0.5, 0.5, 0.0]]).

        \item ngkpt\_string - string of k-points grid, e.g. '8 8 8', that means 8x8x8 grid in reciprocal space
        
        \end{itemize}

        Results of the method are: 

        list of reciprocal vectors lengths in reciprocal angstrom units, written in attribute 'self.recip\_len\_list',

        list of k-points densities for each vector per reciprocal angstrom, written in attribute 'self.list\_den'
    \end{itemize}

Input file k\_points\_mi.in (see example) contains following parameters:

\\begin{itemize}

\item rprim - dimensionless real space lattice vectors of elemental cell

\item acell - lattice parameters in angstroms (multipliers for dimensionless vectors)

\item regime - type of program work: 'density' - calculations of k-points density; 'ngkpt' - calculations
of the k-points grid

\end{itemize}

Program can be used as executable or importable 
You should input additional data corresponding to prompts if program is used as executable 

\\newpage
'''

from sympy.matrices import *
from math import pi
class Reciprocal_Lattice:
    def k_points(self, rprim_list, density_kpt_string):
        density_kpt_string_list = density_kpt_string.split()
        density_kpt_list = [int(i) for i in density_kpt_string_list]
        dkl = density_kpt_list
        m = Matrix((rprim_list[0], rprim_list[1], rprim_list[2]))
        v = m.det()
        f = (2*pi/v)
        b1 = [f*(m[4]*m[8]-m[7]*m[5]), f*(m[6]*m[5]-m[3]*m[8]), f*(m[3]*m[7]-m[4]*m[6])]
        b2 = [f*(m[7]*m[2]-m[8]*m[1]), f*(m[8]*m[0]-m[6]*m[2]), f*(m[6]*m[1]-m[7]*m[0])]
        b3 = [f*(m[1]*m[5]-m[2]*m[4]), f*(m[3]*m[2]-m[5]*m[0]), f*(m[0]*m[4]-m[1]*m[3])]
        len_b1 = (b1[0] ** 2 + b1[1] ** 2 + b1[2] ** 2) ** 0.5
        len_b2 = (b2[0] ** 2 + b2[1] ** 2 + b2[2] ** 2) ** 0.5
        len_b3 = (b3[0] ** 2 + b3[1] ** 2 + b3[2] ** 2) ** 0.5
        ngkpt_i = [int(dkl[0] * len_b1 + 1), int(dkl[1] * len_b2 + 1), int(dkl[2] * len_b3 + 1)]
        ngkpt_l = [str(i) for i in ngkpt_i]
        ngkpt = ''
        for i in ngkpt_l:
            ngkpt += i + ' '
        self.ngkpt = ngkpt
        self.density_list = dkl
        self.volume = v
        self.recip_vecs = [b1, b2, b3]
    def k_density(self, rprim_list, ngkpt_string):
        m = Matrix((rprim_list[0], rprim_list[1], rprim_list[2]))
        v = m.det()
        f = (2*pi/v)
        b1 = [f*(m[4]*m[8]-m[7]*m[5]), f*(m[6]*m[5]-m[3]*m[8]), f*(m[3]*m[7]-m[4]*m[6])]
        b2 = [f*(m[7]*m[2]-m[8]*m[1]), f*(m[8]*m[0]-m[6]*m[2]), f*(m[6]*m[1]-m[7]*m[0])]
        b3 = [f*(m[1]*m[5]-m[2]*m[4]), f*(m[3]*m[2]-m[5]*m[0]), f*(m[0]*m[4]-m[1]*m[3])]
        len_b1 = (b1[0] ** 2 + b1[1] ** 2 + b1[2] ** 2) ** 0.5
        len_b2 = (b2[0] ** 2 + b2[1] ** 2 + b2[2] ** 2) ** 0.5
        len_b3 = (b3[0] ** 2 + b3[1] ** 2 + b3[2] ** 2) ** 0.5
        recip_len_list = [len_b1, len_b2, len_b3]
        ngkpt_list = ngkpt_string.split()
        ngkpt_list_i = [int(i) for i in ngkpt_list]
        den_1 = int(round(ngkpt_list_i[0]/len_b1))
        den_2 = int(round(ngkpt_list_i[1]/len_b2))
        den_3 = int(round(ngkpt_list_i[2]/len_b3))
        self.ngkpt_list_i = ngkpt_list_i
        self.recip_len_list = recip_len_list
        self.list_den = [den_1, den_2, den_3]

if __name__ == '__main__':
    from k_points_mi import Reciprocal_Lattice as RL
    f = open('k_points_mi.in')
    rprim = []
    for i in f:
        if len(i) == 1 or i[0] == '#': continue
        elif 'rprim' == i.split()[0]: 
            rprim.append([float(i.split()[1]), float(i.split()[2]), float(i.split()[3])])
            ii = f.readline()
            rprim.append([float(ii.split()[0]), float(ii.split()[1]), float(ii.split()[2])])
            iii = f.readline()
            rprim.append([float(iii.split()[0]), float(iii.split()[1]), float(iii.split()[2])])
#        elif 'acell' == i.split()[0] and 'angstrom'== i.split()[4]: acell_list = [float(j) for j in i.split()[1:4]]
        elif 'acell' == i.split()[0]: acell = [float(j) for j in i.split()[1:]]
        elif 'regime' == i.split()[0]: regime = i.split()[1]

    rprimd_list = [[i*acell[0] for i in rprim[0]], [i*acell[1] for i in rprim[1]], [i*acell[2] for i in rprim[2]]]
    print(rprimd_list)
    a = RL()
    if regime == 'density': 
        ngkpt_string = input('Input value of the ngkpt parameter ')
        a.k_density(rprimd_list, ngkpt_string)
        print('Density of k-points along each vector equals to'+str(a.list_den[0])+' '+str(a.list_den[1])+' '+str(a.list_den[2]))
        print(a.recip_len_list)
    if regime == 'ngkpt': 
        dkl_raw = input('input value of the k-points density per reciprocal angstrom along each vector ')
        dkl = str(dkl_raw)+' '+str(dkl_raw)+' '+str(dkl_raw)
        a.k_points(rprimd_list, dkl)
        print('ngkpt '+a.ngkpt)
















	
