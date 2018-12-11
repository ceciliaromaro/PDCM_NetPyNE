'''
NetPyNE version of Potjans and Diesmann thalamocortical network

init.py -- code to run the simulation

to compile mod files: 
	nrnivmodl  
to run on single core: 
	python init.py
to run on multiple cores: 
	mpiexec -n 2 nrniv -python -mpi init.py 
'''


from netpyne import specs, sim
from neuron import h
from cfg import cfg
from netParams import netParams


############################################################
#               Create network and run simulation
############################################################

sim.initialize(
    simConfig = cfg,  
    netParams = netParams)          # create network object and set cfg and net params
sim.net.createPops()                    # instantiate network populations
sim.net.createCells()                   # instantiate network cells based on defined populations

# randomize m parameter of cells
# rand=h.Random()
# for c in sim.net.cells:
# 	if c.tags['cellModel'] == 'IntFire_PD':
# 		rand.Random123(c.gid, cfg.seeds['m'])
# 		c.hPointp.m = rand.normal(-58,10)

sim.net.addStims()              # add network stimulation
sim.net.connectCells()                  # create connections between cells based on params
sim.setupRecording()                    # setup variables to record for each cell (spikes, V traces, etc)
sim.runSim()                            # run parallel Neuron simulation  
sim.gatherData()                        # gather spiking data and cell info from each node
sim.saveData()                          # save params, cell info and sim output to file (pickle,mat,txt,etc)#
#sim.analysis.plotData()               # plot spike raster etc




