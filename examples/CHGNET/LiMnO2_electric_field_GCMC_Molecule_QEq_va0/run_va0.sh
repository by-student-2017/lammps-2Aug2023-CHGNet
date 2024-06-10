#!/bin/bash

export OMP_NUM_THREADS=1

mkdir cfg

/mnt/d/lammps/src/lmp_serial -in in_va0.lmp

