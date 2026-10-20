import os
import subprocess
import time
import numpy as np
from utils import download_domainsize_script, read_klong_data, find_next_run_number
from nosub_polymerbuilder import consbcp
from lammps_generator import generate_bcp_input_file, generate_slurm_script
import sys
sys.dont_write_bytecode = True

def generate_parallel_slurm_script(lmp_path, input_file, output_filename, hours, run_dirs, cores):

    with open(output_filename, 'w') as f:

        f.write(f"#!/bin/bash\n")
        f.write(f"#SBATCH --nodes=1\n")
        f.write(f"#SBATCH --ntasks-per-node=96\n")
        f.write(f"#SBATCH --cluster=cardinal\n")
        f.write(f"#SBATCH --time={hours}:00:00\n\n")
        f.write("module purge\n")
        f.write("module load gcc/12.3.0\n")
        f.write("module load mvapich/4.1\n\n")
                

        for r_dir in run_dirs:

            f.write(f"srun --ntasks={cores} --exclusive --chdir={r_dir} {lmp_path} -in {input_file} > {r_dir}/out.lmp &\n")
            f.write(f"python std_domainsize.py \n")

        f.write("\nwait\n")




def run_simulation(output_filename, cluster, lmp_path, input_file, run_dir="."):
    """doing simulations"""
    original_dir = os.getcwd()
    os.chdir(run_dir)
    
    try:
        if cluster:
            try:
                result = subprocess.run(["sbatch", output_filename], check=True, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                print(f"job submitted:\n{result.stdout}")
            except subprocess.CalledProcessError as e:
                print(f"job failed:\n{e.stderr}")
            except FileNotFoundError:
                print("error: did not find 'sbatch', make sure run the simulation on Cardinal")    
        else:
            subprocess.run(f"{lmp_path} < {input_file}", shell=True, check=True)
    finally:
        os.chdir(original_dir)

def main_multiple_runs(parameter_matrix,cores_per_job=96, start_run_number=None, 
                      temperature=1.0, timestep=200000, hours=3,
                      lmp_path="/users/PAA0027/yuanhaozhang/lammpstest/lammps_latest/lammps_src/build/lmp_mpi",
                      input_file="bcp.inp", output_filename="submit_job.sh",
                      cluster=True, domainsize=True):
    """run multiple simulations"""
    n_runs = parameter_matrix.shape[0]
    cores_per_job = 96 // n_runs
    results = []
    run_directories = []
    
    if start_run_number is None:
        start_run_number = find_next_run_number()
    
    print(f"=== simulation start, initial index: {start_run_number} ===")
    download_domainsize_script(".")        
    print(f"=== submitting {n_runs}th job ===")
    
    # for i in range(n_runs):
    #     params = parameter_matrix[i]
    #     a, b, c, d, A, B, C, D = params
    #     run_number = start_run_number + i
    #     run_dir = f"run_{run_number:03d}"
    #     os.makedirs(run_dir, exist_ok=True)
    #     run_directories.append(run_dir)
        
    #     print(f"\n--- submitting {i+1}/{n_runs}th job (index: {run_number:03d}) ---")
    #     print(f"input variables: M1A={a}, M1B={b}, M2A={c}, M2B={d}, L1A={A}, L1B={B}, L2A={C}, L2B={D}")
    #     print(f"directory: {run_dir}")
        
    #     if os.path.exists("std_domainsize.py"):
    #         import shutil
    #         shutil.copy("std_domainsize.py", run_dir)

    #     target_data_path = os.path.join(run_dir, "data.bcp")
    #     consbcp(a, b, c, d, A, B, C, D, save_path=target_data_path)
    #     generate_bcp_input_file(temperature, timestep, filename=input_file, run_dir=run_dir)
        
    #     generate_slurm_script(lmp_path, input_file, output_filename, domainsize, hours, run_dir)
    #     run_simulation(output_filename, cluster, lmp_path, input_file, run_dir)
    
    # print(f"\n=== all {n_runs} jobs have been submitted ===")
    for i in range(n_runs):
        params = parameter_matrix[i]
        a, b, c, d, A, B, C, D = params
        run_number = start_run_number + i
        run_dir = f"run_{run_number:03d}"
        os.makedirs(run_dir, exist_ok=True)
        run_directories.append(run_dir)
    
        download_domainsize_script(run_dir)
        consbcp(a, b, c, d, A, B, C, D, run_dir=run_dir)
        generate_bcp_input_file(temperature, timestep, filename=input_file, run_dir=run_dir)
        log_file = os.path.join(run_dir, "lammps.out")
        cmd = f"srun --ntasks={cores_per_job} --exclusive --chdir={run_dir} {lmp_path} -in {input_file} > {log_file} &"
        subprocess.Popen(cmd, shell=True)
        print(f"Parallel job {i+1} launched in {run_dir} with {cores_per_job} cores.")

    print(f"\n=== All {n_runs} jobs are running in background. Waiting for klong.data... ===")   

    wait_timeout = hours*3600  
    start_wait_time = time.time()

    print(f"starting time: {time.strftime('%Y-%m-%d %H:%M:%S')}")
    print(f"expected finished at: {time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(time.time() + wait_timeout))}")

    all_completed = False
    for minute in range(wait_timeout // 60):
        elapsed_minutes = minute + 1
        elapsed_seconds = elapsed_minutes * 60
        remaining_seconds = wait_timeout - elapsed_seconds
        progress_percent = (elapsed_seconds / wait_timeout) * 100
        #calcuating remaining time
        remaining_hours = remaining_seconds // 3600
        remaining_minutes = (remaining_seconds % 3600) // 60
        
        completed_count = 0
        for run_dir in run_directories:
            klong_file = os.path.join(run_dir, "klong.data")
            if os.path.exists(klong_file):
                try:
                    with open(klong_file, 'r') as f:
                        content = f.read().strip()
                        if content:
                            completed_count += 1
                except:
                    pass
        
        print(f"\rprogess: {progress_percent:.1f}% | "
            f"have run for: {elapsed_minutes} mins | "
            f"remaining: {int(remaining_hours)}  hours {int(remaining_minutes)}  mins | "
            f"finished: {completed_count}/{n_runs}", 
            end="", flush=True)
        
        # if finished, quit
        if completed_count == n_runs:
            print(f"\nall job done! real waiting time: {elapsed_minutes} mins")
            all_completed = True
            break
        
        # wait for 60 s
        time.sleep(60)
    
    if not all_completed:
        print(f"\nJob done! total waiting time: {wait_timeout//60}min")
        # check whether the jobs are done
        completed_count = 0
        for run_dir in run_directories:
            klong_file = os.path.join(run_dir, "klong.data")
            if os.path.exists(klong_file):
                try:
                    with open(klong_file, 'r') as f:
                        content = f.read().strip()
                        if content:
                            completed_count += 1
                except:
                    pass
        print(f"job completed: {completed_count}/{n_runs} overall run")
    
    print(f"\n=== collecting data ===")
    #collecting the data fron klong.data
    for i, run_dir in enumerate(run_directories):
        print(f"\n--- analysising {i+1}/{n_runs}th job results ---")
        print(f"directory: {run_dir}")            
        klong_values = read_klong_data(run_dir)             
        if klong_values is not None:
            results.append(klong_values)
            print(f"{i+1} th result: {klong_values}")
        else:
            print(f"{i+1} th result: failed (filled with NaN)")
    
    results_array = np.array(results)
    print(f"\n=== all analysis done ===")
    print(f"result matrix shape is: {results_array.shape}")
    
    return results_array
