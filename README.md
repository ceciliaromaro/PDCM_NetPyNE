# PD_in_NetPyNE
Implementation of the Potjans-Diesmann cortical microcircuit model in NetPyNE/NEURON with rescaling option


To reproduce the PD model, a new NEURON LIF neuron NMODL implementation was required since the default one did not allow setting the membrane time constant higher than the synaptic decay time constant, as in thePD LIF model. The network Python code and the LIF neuron NMODL (.mod) code are made publicly available in GitHub and will be presented in CNS 2018.

Working in progress. New detailed description will be available soon.
