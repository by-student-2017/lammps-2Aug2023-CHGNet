Customized LAMMPS(2Aug2023) for Neural Network Potential, by AdvanceSoft Corp. <http://www.advancesoft.jp>.

This is a prototype that can output MAGMOM (magnetic moment) as q (=charge) on Lammps. I don't have the capacity to convert MAGMOM values ​​to electric charge (q), so I'm leaving it as it is. Please be careful.

"va0" Version gives a charge from Magmom. Currently, the provisional formula is given, but after calculating the same structure in the same structure for the calculation system, the charge is calculated so that the same value is calculated from the MAGMOM value. " Pair_chgnet.cpp "should be added.
