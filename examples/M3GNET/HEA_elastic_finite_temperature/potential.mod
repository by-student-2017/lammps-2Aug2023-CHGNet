# NOTE: This script can be modified for different pair styles 
# See in.elastic for more info.

# we must undefine any fix ave/* fix before using reset_timestep
if "$(is_defined(fix,avp))" then "unfix avp"
reset_timestep 0

variable elem string "Co Ni Cr Fe Mn"

#-------------------- Force field --------------------------------------------------------
pair_style      m3gnet /mnt/d/lammps/potentials/M3GNET
#pair_style      m3gnet/d3 /mnt/d/lammps/potentials/M3GNET

#pair_coeff      * *  MP-2021.2.8-EFS ${elem} # M3GNet <https://github.com/materialsvirtuallab/m3gnet> will be called
pair_coeff      * *  M3GNet-MP-2021.2.8-PES ${elem} # MatGL <https://github.com/materialsvirtuallab/matgl> will be called

# Setup neighbor style
neighbor 1.0 nsq
neigh_modify once no every 1 delay 0 check yes

# Setup output

fix avp all ave/time  ${nevery} ${nrepeat} ${nfreq} c_thermo_press mode vector
thermo		${nthermo}
thermo_style custom step temp pe press f_avp[1] f_avp[2] f_avp[3] f_avp[4] f_avp[5] f_avp[6]
thermo_modify norm no

# Setup MD

timestep ${timestep}
fix 4 all nve
if "${thermostat} == 1" then &
   "fix 5 all langevin ${temp} ${temp} ${tdamp} ${seed}"


