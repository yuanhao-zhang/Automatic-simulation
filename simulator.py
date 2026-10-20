import torch
import numpy as np
import pandas as pd
import os
import sys
from job_manager import main_multiple_runs
from BOExperimentRec import printNewPoint

# def md(x):
#     parameter_matrix = x.astype(int)     
#     results = main_multiple_runs(parameter_matrix)
    
#     file_name = "training_data.csv"
    
#     if os.path.exists(file_name):
#         with open(file_name, 'a', encoding='utf-8', newline='') as f:
#             for i in range(len(parameter_matrix)):
#                 params = [int(p) for p in parameter_matrix[i]]
#                 sim_res = [float(r) for r in results[i]]
                
#                 param_str = ",".join(map(str, params))
#                 res_str = ",".join(["{:.9f}".format(r) for r in sim_res])
                
#                 full_line = f"{param_str},{res_str}\n"
#                 f.write(full_line)
              
#     return results

def md(x):

    parameter_matrix = x.astype(int) 
    num_jobs = len(parameter_matrix)
    cores_per_job = 96 // num_jobs
    
    print(f" detect {num_jobs} jobs in parallel，assigning {cores_per_job} nodes per job")   
    results = main_multiple_runs(parameter_matrix, cores_per_job=cores_per_job)
    
    return results


def run_bo_and_wait(timeout=10):
    """
    Step 5: test if there exist suggested file and process BO code
    """
    print("\n Excuting BOExperimentRec to get the suggested values for next trail...")
    printNewPoint()
    print("BO activated, waiting for the next recommendation.txt...")

    rec_file = "recommendation.txt"
    # 2. test if there are suggested files
    while not os.path.exists(rec_file):
        current_time = time.strftime("%H:%M:%S", time.localtime())
        print(f"[{current_time}] no detected, finding suggested file after {timeout} s...")
        time.sleep(timeout)
    
    print(f"--- find the next suggested value ---")

if __name__ == "__main__":

    total_iterations = 1# number of iterations
    rec_file = "recommendation.txt"

    for i in range(total_iterations):
        print(f"\n{'#'*50}")
        print(f"### doing  {i+1} / {total_iterations}th iteration ###")
        print(f"{'#'*50}")

        # 1. Reading the input parameters from "recommendation.txt"
        if os.path.exists(rec_file):
            try:
                best_x = np.loadtxt(rec_file, delimiter=',')
                best_x = np.atleast_2d(best_x)
                num_jobs = len(best_x)

                
                print(f" {num_jobs} jobs detected, running the jobs in parallel")
                
                # 2. Remove the previous recommendation.txt file
                os.remove(rec_file)
                print(f"deleting {rec_file}")

                # 3 & 4. Doing simulation
                md(best_x)

                # 5. Do BO code and Wait for the next suggested value 
                # Scaning the BO code will no be excuted if this is the final run
                if i < total_iterations - 1:
                    run_bo_and_wait(timeout=30)
                if i == total_iterations - 1:
                    printNewPoint()


            except Exception as e:
                print(f"Error during processing: {e}")
                break 
        else:
            print(f"did not find suggesting file {rec_file}")
            break

    print("\n[Finish] loop finished or terminated")
