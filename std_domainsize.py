import numpy as np
import os
import time
import shutil
import requests
import re
import sys
from collections import deque
sys.dont_write_bytecode = True




def find_next_run_number(base_dir="."):
    """find next one"""
    pattern = re.compile(r'^run_(\d+)$')
    existing_numbers = []
    
    for item in os.listdir(base_dir):
        if os.path.isdir(item):
            match = pattern.match(item)
            if match:
                existing_numbers.append(int(match.group(1)))
    
    if existing_numbers:
        return max(existing_numbers) + 1
    else:
        return 1

#def download_domainsize_script(run_dir="."):

#    dsurl = "https://raw.githubusercontent.com/yuanhao-zhang/PPPMD/main/domainsize"
#    save_path = os.path.join(run_dir, "domainsize.py")
#
#    try:            
#        response = requests.get(dsurl, stream=True)
#        response.raise_for_status()            
#        with open(save_path, 'wb') as f: 
#            for chunk in response.iter_content(chunk_size=8192):
#                f.write(chunk)           
#        print(f"file downloaded {os.path.abspath(save_path)}")             
#    except Exception as e:
#        print(f"warning: failed to download: {e}")


def download_domainsize_script(run_dir="."):
    source_file = "std_domainsize.py"
    target_path = os.path.join(run_dir, source_file)
    if not os.path.exists(source_file):
        print(f"warning: cannot find {source_file}")
        return
    if os.path.abspath(source_file) == os.path.abspath(target_path):
        return
    try:
        shutil.copy2(source_file, target_path)
        print(f"file have been copied to: {os.path.abspath(target_path)}")
    except Exception as e:
        print(f"failed to copy the file: {e}")

def read_klong_data(run_dir):
    """read the last four numbers in klong.data"""
    klong_file = os.path.join(run_dir, "klong.data")
    try:
        with open(klong_file, 'r') as f:
            last_line = deque(f, 1)
            if not last_line:
                return None
            
            line_content = last_line[0].strip()
            parts = line_content.split()
            if len(parts) >= 4:
                values = [float(x) for x in parts[-4:]]
                return values
            else:
                print(f"Warning: {run_dir} final results less than 4 elements")
                return None               
    except FileNotFoundError:
        print(f"Warning: did not find klong.data in {run_dir}")
        return None