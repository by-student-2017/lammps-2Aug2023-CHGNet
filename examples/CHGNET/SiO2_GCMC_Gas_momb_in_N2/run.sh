#!/bin/bash

#lammps_adress=/mnt/d/lammps
lammps_adress=$HOME/lammps

mkdir cfg

${lammps_adress}/src/lmp_serial -in in.lmp
