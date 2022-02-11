# NetPyNE implementation of the Potjans and Diesmann Cortical Microcircuit model with scaling option

This code reproduces the simulations for the following paper:

Cecilia Romaro, Fernando Araujo Najman, William W. Lytton, Antonio C. Roque, Salvador Dura-Bernal. **NetPyNE Implementation and Scaling of the Potjans-Diesmann Cortical Microcircuit Model**. Neural Comput 2021; 33 (7): 1993–2032. doi: https://doi.org/10.1162/neco_a_01400

The Potjans-Diesmann cortical microcircuit model is a widely used model originally implemented in NEST. Here, we reimplemented the model using NetPyNE, a high-level Python interface to the NEURON simulator, and reproduced the findings of the original publication. We also implemented a method for scaling the network size that preserves first- and second-order statistics, building on existing work on network theory. Our new implementation enabled the use of more detailed neuron models with multicompartmental morphologies and multiple biophysically realistic ion channels. This opens the model to new research, including the study of dendritic processing, the influence of individual channel parameters, the relation to local field potentials, and other multiscale interactions. The scaling method we used provides flexibility to increase or decrease the network size as needed when running these CPU-intensive detailed simulations. Finally, NetPyNE facilitates modifying or extending the model using its declarative language; optimizing model parameters; running efficient, large-scale parallelized simulations; and analyzing the model through built-in methods, including local field potential calculation and information flow measures.

The code requires NEURON (> v7.8) and Python (>v3.5) and has been tested on Linux machines (CentOS and Ubuntu) and Mac OS X. 
To run the model follow these steps:

1. Unzip all files.
2. Change to model folder and compile the IntFire_PD.mod file via "nrnivmodl .". This should create a directory called either i686 or x86_64, depending on your computer's architecture.
3. Run the init.py file. This can be done directly via `python -i init.py` or using MPI, e.g. `mpiexec -np 4 nrniv -mpi init.py`. Alternatively, if running in an HPC via SLURM you can adapt and run `sbatch pd.sbatch`   

List of files:

- `netParams.py`: High-level specifications of the network parameters. Includes the definition of the populations (`popParams`) and connectivity (`connParams`). Note the external inputs are also implemented as populations of NetStims. This file also contains the function to scale the connections and inputs (`Reescale`) based on the `cfg.ScaleFactor`

- `cfg.py`: Simulation configuration options, including the simulation duration (`cfg.duration`) and the network scale factor (`cfg.ScaleFactor`). This file also contains the flags to reproduce the different conditions shown in the paper: `cfg.DC`, `cfg.TH` and `cfg.Balanced`. Additionally, the user can select different recording, saving and plotting options.

- `init.py`: File to run the simulation; also includes the function calls for creating the network, saving data and analyzing results.

- `pd.sbatch` and `pd_128.sbatch` etc: Script files to run on HPCs via SLURM workload manager using different number of cores.

- `IntFire_PD.mod`: NMODL file for he integrate and fire artificial cell used in this network.

To run the multicompartment neuron version of the network use the files in the `/multicompartment` folder. The running instructions and file description is the same as above, except that you need to compile the `mod` folder via `nrnivmodl mod`. 

The model is also available in Github: https://github.com/suny-downstate-medical-center/PDCM_NetPyNE ; and the multicompartment version is available in the "multicompartment" branch: https://github.com/suny-downstate-medical-center/PDCM_NetPyNE/tree/multicompartment 


For any questions or further assistance please contact: salvadordura at gmail.com

 