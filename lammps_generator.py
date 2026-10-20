import os
import sys
sys.dont_write_bytecode = True



def generate_bcp_input_file(temperature, timestep, epsab=0.91, filename="bcp.inp", run_dir="."):
    """generating LAMMPS input file"""
    filename = os.path.join(run_dir, filename)
    
    content = f"""variable   T equal {temperature}             
variable   t equal {timestep}

###Initialization

units			lj
atom_style		full
#processors      4 4 4
boundary        p p p
read_data	    bulk.lammps
neighbor		1.0 bin
neigh_modify            every 1 delay 0 check yes
neigh_modify            one 4000
thermo_style	        custom step temp vol press ke pe ebond eimp evdwl ecoul elong
thermo			10000
timestep		0.005

###Force Field 1: Intramolecular interactions turned on and pushoff

bond_style		fene
bond_coeff      1 30     1.5   1.0    1.0
bond_coeff      2 30     1.5   {epsab}    1.0

special_bonds    	fene
group polymer type 1 2
group substrate type 3 4
neigh_modify exclude group substrate substrate

fix fix1 substrate setforce 0.0 0.0 0.0
min_style       cg
minimize        1.0e-6 1.0e-8 100000 10000000
unfix fix1


###Soft pushoff
pair_style		soft 1.122
pair_coeff              * * 0.0 1.122
variable		prefactor1 equal 10.0*elapsed/1000
fix			fixs polymer adapt 1 pair soft a * * v_prefactor1
fix			fix1 polymer nve
fix			temp polymer temp/rescale 1 0.0 $T 0.05 1.0
dump			dump1 all custom 1000 pushoff.lammpstrj id mol type q x y z ix iy iz

run			2000

unfix			fixs
unfix	  	        fix1
unfix		        temp
undump		        dump1


fix fix1 substrate setforce 0.0 0.0 0.0
min_style       cg
minimize        1.0e-6 1.0e-8 100000 10000000
unfix fix1


pair_style              lj/cut 2.5 
pair_modify             shift yes
pair_coeff              1 1  1.0   1.0    2.5
pair_coeff              1 2  {epsab}  1.0    2.5
pair_coeff              1 3  1.0   1.0    2.5
pair_coeff              1 4  {epsab}   1.0    2.5

pair_coeff              2 2  1.0   1.0    2.5
pair_coeff              2 3  {epsab}   1.0    2.5
pair_coeff              2 4  1.0   1.0    2.5

pair_coeff              3 3  1.0   1.0    2.5
pair_coeff              3 4  {epsab}   1.0    2.5

pair_coeff              4 4  1.0   1.0    2.5

fix fix1 substrate setforce 0.0 0.0 0.0
min_style       cg
minimize        1.0e-6 1.0e-8 100000 10000000
unfix fix1


restart			        2000000 diblock.restart
fix                     fix6 polymer nvt temp $T $T 0.5 
fix                     fixg polymer gravity 1 vector 0 0 -1
compute                 temp_polymer polymer temp
thermo_style            custom step temp c_temp_polymer vol press ke pe ebond eimp evdwl ecoul elong
dump                    dump3 all custom 20000 3_ht.lammpstrj id mol type q xs ys zs ix iy iz
dump_modify             dump3  first yes
run                     20000
unfix                   fixg
run                     $t
"""

    with open(filename, 'w') as f:
        f.write(content)
    print(f"LAMMPS input file '{filename}' generated, temperature = {temperature}, timestep = {timestep}.")

def generate_slurm_script(lmp_path, input_file, output, domainsize=True, hours=1, run_dir="."):
    """generate work script"""
    output_filename = os.path.join(run_dir, output)
    processors = 96  # 4x4x4 processors
    
    content = f"""#!/bin/bash
#SBATCH --job-name=1_singleion
#SBATCH --account=PAA0027
#SBATCH --nodes=1
#SBATCH --ntasks-per-node={processors}
#SBATCH --time={hours}:00:00
#SBATCH --output=1_jobinfo_%j.out

cd {run_dir}
srun {lmp_path} < {input_file}

"""
    if domainsize:
        content += "python std_domainsize.py\n"
   
    with open(output_filename, 'w', newline='\n') as f:
        f.write(content)
    
    print(f"SLURM script '{output_filename}' generated")
    print(f"LAMMPS direc: {lmp_path}")
    print(f"input file: {input_file}")
