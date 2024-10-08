"""
@Shihab
Weak analysis: Submit jobs with varying node size and 
mesh size. It copies config files and slurm script, modifies
them as necessary, and submits the jobs.

Requires: Path to config file. And slurm script `basic_job.sh` 
to be in same repo.
"""
import yaml
import os
import subprocess

file_path = '../configs/rk-32.yaml'

NODES = [1,2,4,8,16]
MESH  = [[9, 9, 9], [10, 9, 9], [10, 10, 9], [10, 10, 10],
        [11, 10, 10]]

rank_job = []
for node, mesh in zip(NODES, MESH):
    ## Modify YAML
    with open(file_path, 'r') as file:
        data = yaml.safe_load(file)
        data['levels'] = mesh
    
    with open(f"xcf_{node}.yaml", 'w') as out:
        yaml.dump(data, out, default_flow_style=True)

    # modify submit_script
    s = None
    with open('basic_job.sh', 'r') as file:
        lines = file.readlines()
        lines[4] = f"#SBATCH --nodes={node}\n" #BE CRAEFULL abt INDEX
        s = ''.join(lines)
    
    with open(f"xslurm_multi_{node}.sh", 'w') as out:
        out.write(s)

    # Submit job and remember job id
    try:
        result = subprocess.run(['sbatch', f"xslurm_multi_{node}.sh"], 
            stdout=subprocess.PIPE, stderr=subprocess.PIPE,
            universal_newlines = True)        
        if result.returncode == 0:
            output = result.stdout.strip()
            rank_job.append(f"{node}: " + output)
        else:
            print("Error submitting job:", result.stderr)
    
    except Exception as e:
        print("An error occurred while submitting the job:", node, str(e))

with open("weak_exp.txt", 'w') as fh:
    fh.write('\n'.join(rank_job))
