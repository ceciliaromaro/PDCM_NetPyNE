#!/bin/bash 
#SBATCH --job-name=pd
#SBATCH -A default
#SBATCH -t 60:00:00
#SBATCH --nodes=16
#SBATCH --ntasks-per-node=64
#SBATCH -o /home/salvadord/pd/data/pd_scale-1.0_DC-0_TH-0_Balanced-1_60sec_1024.run
#SBATCH -e /home/salvadord/pd/data/pd_scale-1.0_DC-0_TH-0_Balanced-1_60sec_1024.err
#SBATCH --mail-user=salvadordura@gmail.com
#SBATCH --mail-type=end
#SBATCH --exclude=compute[17-64000]

source ~/.bashrc
cd /home/salvadord/pd
mpirun -np 1024 nrniv -python -mpi init.py
wait
                            
