#!/bin/bash 
#SBATCH --job-name=pd
#SBATCH -A default
#SBATCH -t 60:00:00
#SBATCH --nodes=8
#SBATCH --ntasks-per-node=16
#SBATCH -o /home/salvadord/pd/data/pd_scale-1.0_DC-0_TH-0_Balanced-1_1sec.run
#SBATCH -e /home/salvadord/pd/data/pd_scale-1.0_DC-0_TH-0_Balanced-1_1sec.err
#SBATCH --mail-user=salvadordura@gmail.com
#SBATCH --mail-type=end
#SBATCH --exclude=compute[17-64000]

source ~/.bashrc
cd /home/salvadord/pd
mpirun -np 128 nrniv -python -mpi init.py
wait
                            
