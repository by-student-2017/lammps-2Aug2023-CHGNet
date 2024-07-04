# NOTE: This script can be modified for different pair styles 
# See in.elastic for more info.

variable elem string "Co Ni Cr Fe Mn"

#-------------------- Force field --------------------------------------------------------
pair_style      m3gnet /mnt/d/lammps/potentials/M3GNET
#pair_style      m3gnet/d3 /mnt/d/lammps/potentials/M3GNET

#pair_coeff      * *  MP-2021.2.8-EFS ${elem} # M3GNet <https://github.com/materialsvirtuallab/m3gnet> will be called
pair_coeff      * *  M3GNet-MP-2021.2.8-PES ${elem} # MatGL <https://github.com/materialsvirtuallab/matgl> will be called

# Setup neighbor style
neighbor 1.0 nsq
neigh_modify once no every 1 delay 0 check yes

# Setup minimization style
min_style	     cg
min_modify	     dmax ${dmax} line quadratic

# Setup output
thermo		1
thermo_style custom step temp pe press pxx pyy pzz pxy pxz pyz lx ly lz vol
thermo_modify norm no
