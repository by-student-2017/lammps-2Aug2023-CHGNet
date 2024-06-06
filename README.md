Customized LAMMPS(2Aug2023) for Neural Network Potential, by AdvanceSoft Corp. <http://www.advancesoft.jp>.

Comment --By STUDENT ---

This is a prototype that can output MAGMOM (magnetic moment) as q (=charge) on Lammps. I don't have the capacity to convert MAGMOM values ​​to electric charge (q), so I'm leaving it as it is. Please be careful.

A "va0" version gives a charge from Magmom. Currently, the provisional formula is given, but after calculating the same structure in the same structure for the calculation system, the charge is calculated so that the same value is calculated from the MAGMOM value. " pair_chgnet.cpp "should be added. 

The "va0" version has been turned off by warning by setting "Boundary". If you want to set "Boundary" to "f", please prepare a vacuum layer or solvent with 10 Angstrom or higher for the axis to "f". In practical use, it is rare to calculate with 10 Angstrom or less in Lammps, so there is no problem.

"Magnetic interactions in molecules and an analysis of molecular electronic charge distribution from magnetic parameters"(https://pubs.acs.org/doi/pdf/10.1021/cr60292a003)

"Geometric, electronic, and magnetic structure of Co2⁢FeSi: Curie temperature and magnetic moment measurements and calculations"(https://doi.org/10.1103/PhysRevB.72.184434)

I've managed to assign the magnetic moment value from CHGNet to the q part where the charge value is stored in Lammps, so I've shown this below. From here, rewrite line 123 of chgnet_driver.py so that the charge is calculated from the magnetic moment value, or rewrite the "this->magmoms[iatom]" part on lines 235 or 729 of pair_chgnet.cpp. Simply write into the code an equation that uses "this->magmoms[iatom]" and atomic information to calculate the charge. (https://github.com/by-student-2017/lammps-2Aug2023-CHGNet)

Or you could use "chgnet_driver.py" and put the charges into magmoms, if there is code in python or other languages ​​that predicts the charges from the structure, magnetic moments, etc.

https://gitlab.com/jmargraf/qpac

[Schrödinger Python API 2018](https://content.schrodinger.com/Docs/r2018-2/python_api/api/schrodinger.structure.html)

https://docs.eyesopen.com/toolkits/python/quacpactk/index.html

https://pypi.org/project/pyeqeq/

https://www.biotite-python.org/apidoc/biotite.structure.partial_charges.html

https://www.rdkit.org/docs/GettingStartedInPython.html

https://github.com/sb-ncbr/AtomicChargeCalculator2

https://github.com/MergunFrimen/molstar-partial-charges

One thing that is bothering me is that the force and virial are set to "+=", so I think I'll have to check that later. It would be fine if they were processed with "+=" after they reach 0 each time. I'll make time to check that later too. I am very worried because PLUMUD uses "=".
