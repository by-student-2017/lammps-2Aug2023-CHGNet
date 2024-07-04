Customized LAMMPS(2Aug2023) for Neural Network Potential, by AdvanceSoft Corp. <http://www.advancesoft.jp>.

Comments -- By STUDENT ---

If you want to calculate an interface with an applied electric field, we recommend that you use QEq to calculate with the original state (advancesoftcorp/lammps).

Since the calculation of "d3" is very heavy, there is also a method to use momb(=-D2) implemented in Lammps. The parameters are described in the original paper. I prepared a folder named momb as a reference example. When used in conjunction with "QEq", "CHGNet" and "M3GNet" will enable a fairly wide range of calculations.

momb (=DFT-D2): make yes-EXTRA-PAIR

we found that the accuracy of the DFT-D2 method of Grimme is comparable to that of the DFT+vdWsurf method. [Ref; https://doi.org/10.1021/jp412098n]

QEq: make yes-QEQ

GCMC: make yes-MC

Currently, the following improvements do not work well with calculations that include "d3". 

A "va0_q_magmom" is a prototype that can output MAGMOM (magnetic moment) as q (=charge) on Lammps. I don't have the capacity to convert MAGMOM values ​​to electric charge (q), so I'm leaving it as it is. Please be careful.

A "va0" version gives a charge from Magmom. Currently, the provisional formula is given, but after calculating the same structure in the same structure for the calculation system, the charge is calculated so that the same value is calculated from the MAGMOM value. " pair_chgnet.cpp "should be added. 

The "va0" version has been turned off by warning by setting "Boundary". If you want to set "Boundary" to "f", please prepare a vacuum layer or solvent with 10 Angstrom or higher for the axis to "f". In practical use, it is rare to calculate with 10 Angstrom or less in Lammps, so there is no problem.

"Magnetic interactions in molecules and an analysis of molecular electronic charge distribution from magnetic parameters"(https://pubs.acs.org/doi/pdf/10.1021/cr60292a003)

"Geometric, electronic, and magnetic structure of Co2⁢FeSi: Curie temperature and magnetic moment measurements and calculations"(https://doi.org/10.1103/PhysRevB.72.184434)

I've managed to assign the magnetic moment value from CHGNet to the q part where the charge value is stored in Lammps, so I've shown this below. From here, rewrite line 123 of chgnet_driver.py so that the charge is calculated from the magnetic moment value, or rewrite the "this->magmoms[iatom]" part on lines 235 or 729 of pair_chgnet.cpp. Simply write into the code an equation that uses "this->magmoms[iatom]" and atomic information to calculate the charge. (https://github.com/by-student-2017/lammps-2Aug2023-CHGNet)


Or you could use "chgnet_driver.py" and put the charges into magmoms, if there is code in python or other languages ​​that predicts the charges from the structure, magnetic moments, etc.

----- Note -----

Failed：Born matrix: make yes-EXTRA-COMPUTE, make yes-python

----- Problems solved after operation check -----

One thing that is bothering me is that the force and virial are set to "+=", so I think I'll have to check that later. It would be fine if they were processed with "+=" after they reach 0 each time. I'll make time to check that later too. I am very worried because PLUMUD uses "=". I compared the force at the 1000th cycle of a simple system and the force at the initial structure before and after changing the code, and the values ​​were almost the same.

In the "va0" version, "+=" is changed to "=" in "void PairCHGNet::performGNN()".
If you use "va0", you do so at your own risk. Either delete "va0" from the file name or replace it with a file that does not have "va0".
