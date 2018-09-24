#!/bin/bash 
#SBATCH --job-name=pd
#SBATCH -A default
#SBATCH -t 60:00:00
#SBATCH --nodes=16
#SBATCH --ntasks-per-node=8
#SBATCH -o /home/salvadord/potjans/pd.run
#SBATCH -e /home/salvadord/potjans/pd.err
#SBATCH --mail-user=salvadordura@gmail.com
#SBATCH --mail-type=end
#SBATCH --exclude=compute[17-64000]

source ~/.bashrc
cd /home/salvadord/potjans
mpirun -np 128 nrniv -python -mpi PD_NEURON6.py
wait
                            
