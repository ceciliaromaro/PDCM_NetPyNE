import matplotlib
matplotlib.use('Agg')
from netpyne import specs, sim
from neuron import h
from numpy import *


simLabel = 'pd_scale-1.0_DC-0_TH-0_Balanced-1'

###nr n iv -python PD_NEURON.py
#nrnivmodl #Compilar .mod
#mpiexec -n 2 nrniv -python -mpi PD_NEURON.py 

#from xxx import Reescale
def Reescale(ScaleFactor, C, N_Full, w_p, f_ext, tau_syn, Inp, InpDC):
	if ScaleFactor<1.0 : 
		F_out=array([0.971, 2.868, 4.746, 5.396, 8.142, 9.078, 0.991, 7.523])
		#Para a nao balanceada o F_out deve ser outro. O nao balanceado. Apos corrigir conexoes testar para caso so de aleracao no InpPoiss
		Ncon=vstack(column_stack(0 for i in range(0,8)) for i in range(0,8))
		for r in range(0,8): 
			for c in range(0,8): 
				Ncon[r][c]=(log(1.-C[r][c])/log(1. -1./(N_Full[r]*N_Full[c])) ) /N_Full[r]

		w=w_p*column_stack(vstack( [1.0, -4.0] for i in range(0,8)) for i in range(0,4))
		w[0][2]=2.0*w[0][2]

		x1_all = w * Ncon * F_out
		x1_sum = sum(x1_all, axis=1)
		x1_ext = w_p * Inp * f_ext
		I_ext=column_stack(0 for i in range(0,8))
		I_ext = 0.001 * tau_syn * (
		        (1. - sqrt(ScaleFactor)) * x1_sum + 
		        (1. - sqrt(ScaleFactor)) * x1_ext)
		#if ScaleFactor<0.09:  #Aqui adaptado
		#	I_ext=1.2*I_ext #1.05 da ativo 1.4
				        
		#Inp=ScaleFactor*Inp*w_p*f_ext*tau_syn*0.001
		InpDC=sqrt(ScaleFactor)*InpDC*w_p*f_ext*tau_syn*0.001 #pA
		w_p=w_p/sqrt(ScaleFactor) #pA
		InpDC=InpDC+I_ext
		N_=int(ScaleFactor*N_Full)
	else:
		InpDC=InpDC*w_p*f_ext*tau_syn*0.001
		N_=N_Full	
	
	return InpDC, N_, w_p


###########################################################
#                      Network Options
###########################################################

# Size of Network. Adjust this constants, please!
ScaleFactor=1.0  # 1.0 = 80.000 

# External input DC or Poisson
DC=False #True = DC // False = Poisson

# Thalamic input in 4th and 6th layer on or off
TH=False #True = on // False = off

# Balanced and Unbalanced external input as PD article
Balanced=True #True=Balanced // False=Unbalanced

# DC=True ;  TH=False; Balanced=True   => Reproduce Figure 7 A1 and A2
# DC=False;  TH=False; Balanced=False  => Reproduce Figure 7 B1 and B2
# DC=False ; TH=False; Balanced=True   => Reproduce Figure 8 A, B, C and D
# DC=False ; TH=False; Balanced=True   and run to 60 s to => Table 6 
# DC=False ; TH=True;  Balanced=True   => Figure 10A. But I want a partial reproduce so I guess Figure 10C is not necessary

###########################################################
#                      Network Constants
###########################################################
# Frequency of external input
f_ext=8.0 # (Hz)
# Postsynaptic current time constant
tau_syn=0.5 # (ms)
# Membrane time constant
tau_m=10 # (ms)
# Other constants defined inside .mod
'''
#Absolute refractory period
tauref =2 (ms) 
#Reset potential
Vreset = -65 (mV) : -49 (mV) :
#Fixed firing threshold
Vteta  = -50 (mV)'''
# Membrane capacity
C_m=250 #pF
# Mean amplitude of the postsynaptic potential (in mV).
w_v=0.15
# Mean amplitude of the postsynaptic potential (in pA).
w_p= (((C_m) ** (-1) * tau_m * tau_syn / (
        tau_syn - tau_m) * ((tau_m / tau_syn) ** (
            - tau_m / (tau_m - tau_syn)) - (tau_m / tau_syn) ** (
            - tau_syn / (tau_m - tau_syn)))) ** (-1))
w_p=w_v*w_p #(pA)

#C probability of connection
C=array([[0.1009, 0.1689, 0.0437, 0.0818, 0.0323, 0.,     0.0076, 0.],
         [0.1346, 0.1371, 0.0316, 0.0515, 0.0755, 0.,     0.0042, 0.],
         [0.0077, 0.0059, 0.0497, 0.135,  0.0067, 0.0003, 0.0453, 0.],
         [0.0691, 0.0029, 0.0794, 0.1597, 0.0033, 0.,     0.1057, 0.],
         [0.1004, 0.0622, 0.0505, 0.0057, 0.0831, 0.3726, 0.0204, 0.],
         [0.0548, 0.0269, 0.0257, 0.0022, 0.06,   0.3158, 0.0086, 0.],
         [0.0156, 0.0066, 0.0211, 0.0166, 0.0572, 0.0197, 0.0396, 0.2252],
         [0.0364, 0.001,  0.0034, 0.0005, 0.0277, 0.008,  0.0658, 0.1443]])	
         
#Population size N
#L=['L2e', 'L2i', 'L4e', 'L4i', 'L5e', 'L5i', 'L6e', 'L6i', 'Th']  
#L=['L2e', 'L2i', 'L4e', 'L4i', 'L5e', 'L5i', 'L6e', 'L6i', 'bkg_IF']
L=['L2e', 'L2i', 'L4e', 'L4i', 'L5e', 'L5i', 'L6e', 'L6i']
N_Full=array([20683, 5834, 21915, 5479, 4850, 1065, 14395, 2948, 902])

#Number of Input per Layer
Inp=array([1600, 1500, 2100, 1900, 2000, 1900, 2900, 2100])
#InpUnbalanced=array([2000, 1850, 2000, 1850, 2000, 1850, 2000, 1850])
if Balanced == False:
	InpUnb=array([2000, 1850, 2000, 1850, 2000, 1850, 2000, 1850])



###########################################################
#                      Reescale calculation
###########################################################
# Aquii
if DC == True:
	InpDC=Inp
	#InpPoiss=zeros(8)
	if Balanced== False:
		InpDC=InpUnb
else: 
	InpDC=zeros(8)
	InpPoiss=Inp*ScaleFactor
	if Balanced== False:
		InpPoiss=InpUnb*ScaleFactor
	
	
InpDC, N_, w_p = Reescale(ScaleFactor, C, N_Full, w_p, f_ext, tau_syn, Inp, InpDC)


############################################################
#                    Network Construction
############################################################
# Network parameters
netParams = specs.NetParams()  # object of class NetParams to store the network parameters

# Population parameters
#v_in=array([-58, -58, -58, -58, -58, -58, -58, -58])
#v_in=array([-49.9, -58, -49.9, -58, -49.9, -58, -49.9, -58])

netParams.delayMin_e = 1.5
netParams.ddelay = 0.5
netParams.delayMin_i = 0.75
netParams.weightMin = w_p
netParams.dweight = 0.1


for i in range(0,8):
	netParams.popParams[L[i]] = {'cellType': str(L[i]), 'numCells': int(N_[i]), 'cellModel': 'IntFire_PD', 'm':0, 'Iext':float(InpDC[i])}

# To atualization of Point Neurons
netParams.popParams['bkg_IF'] = {'numCells': 1, 'cellModel': 'NetStim','rate': 40000,  'start':0.0, 'noise': 0.0, 'delay':0}


if DC == False: # External Input as Poisson
	for r in range(0,8):
		netParams.popParams['poiss'+str(L[r])] = {
						'numCells': N_[r], 'cellModel': 'NetStim',
						'rate': InpPoiss[r]*f_ext,   
						'start':0.0, 'noise': 1.0, 'delay':0}
		
		auxConn=array([range(0,N_[r],1),range(0,N_[r],1)])
		netParams.connParams['poiss->'+str(L[r])] = {
			'preConds': {'pop': 'poiss'+str(L[r])},  'postConds': {'pop': L[r]},
			'connList': auxConn.T,   
			'weight':'max(0, weightMin +normal(0,dweight*weightMin))',  'delay': netParams.delayMin_e} 
			
# Thalamus Input: increased of 15Hz that lasts 10 ms
# 0.15 fires in 10 ms each 902 cells -> number of spikes = T*f*N_ = 0.15*902 -> 1 spike each N_*0.15
if TH == True:
	fth=15 #Hz
	Tth=10 #ms
	InTH=[0, 0, 93, 84, 0, 0, 47, 34]
	for r in [2,3,6,7]:
		nTH=int(sqrt(ScaleFactor)*InTH[r]*fth*Tth/1000)
		netParams.popParams['bkg_TH'+str(L[r])] = {'numCells': N_[r], 'cellModel': 'NetStim','rate': 2*(1000*nTH)/Tth ,  'start': 200.0, 'noise': 1.0, 'number': nTH, 'delay':0}
		auxConn=array([range(0,N_[r],1),range(0,N_[r],1)])
		netParams.connParams['bkg_TH->'+str(L[r])] = {
			'preConds': {'pop': 'bkg_TH'+str(L[r])},  'postConds': {'pop': L[r]},
			'connList': auxConn.T,   
			'weight':'max(0, weightMin +normal(0,dweight*weightMin))',  'delay': netParams.delayMin_e} 
		


## Cell connectivity rules
for r in range(0,8):
	for c in range(0,8):
		if (c % 2) == 0:
			if c == 2 and r == 0:
				netParams.connParams[str(L[c])+'->'+str(L[r])] = { 
        			'preConds': {'pop': L[c]},                         # conditions of presyn cells
        			'postConds': {'pop': L[r]},                        # conditions of postsyn cells
        			#'probability'  C[r][c],                            # probability of connection
        			#'synsPerConn': 2,
        			#'divergence': ScaleFactor*log(1. - C[r][c])/( (log(1. - 1./(N_Full[r]*N_Full[c])) )*(N_Full[c]*ScaleFactor)),
					'divergence': ScaleFactor*(log(1.-C[r][c])/log(1. -1./(N_Full[r]*N_Full[c])) ) /N_Full[c],
					                          #log(1.-C[r][c])/log(1. -1./(N_Full[r]*N_Full[c])) ) /N_Full[r
        			'weight':'2*max(0, weightMin +normal(0,dweight*weightMin))', # synaptic weight
        			'delay':'max(0.1, delayMin_e +normal(0,ddelay*delayMin_e))',  # transmission delay (ms)
        			}
			else:
				netParams.connParams[str(L[c])+'->'+str(L[r])] = { 
        			'preConds': {'pop': L[c]},                         # conditions of presyn cells
        			'postConds': {'pop': L[r]},                        # conditions of postsyn cells
        			'divergence': ScaleFactor*(log(1.-C[r][c])/log(1. -1./(N_Full[r]*N_Full[c])) ) /N_Full[c],
        			'weight':'max(0, weightMin +normal(0,dweight*weightMin))', # synaptic weight
        			'delay':'max(0.1, delayMin_e +normal(0,ddelay*delayMin_e))',  # transmission delay (ms)
        			}                                                # synaptic mechanism
		else:
			netParams.connParams[str(L[c])+'->'+str(L[r])] = { 
        		'preConds': {'pop': L[c]},                         # conditions of presyn cells
        		'postConds': {'pop': L[r]},                        # conditions of postsyn cells
        		'divergence': ScaleFactor*(log(1.-C[r][c])/log(1. -1./(N_Full[r]*N_Full[c])) ) /N_Full[c],
        		'weight':'-4*max(0, weightMin +normal(0,dweight*weightMin))', # synaptic weight
        		'delay':'max(0.1, delayMin_i +normal(0,ddelay*delayMin_i))',  # transmission delay (ms)
        		}                                                  # synaptic mechanism
        
netParams.connParams['S2>M'] = {
	'preConds': {'pop': 'bkg_IF'}, 'postConds': {'cellModel': 'IntFire_PD'},
	'probability': 1, 'weight':0,	'delay':netParams.delayMin_e}   #w=80 		

############################################################
#                    Simulation options and plots
############################################################
simConfig = specs.SimConfig() # object of class SimConfig to store simulation configuration
simConfig.seeds['stim']=3
simConfig.duration = 60*1e3 #6*1e2   # Duration of the simulation, in ms
simConfig.dt = 0.025          # Internal integration timestep to use
simConfig.verbose = 0     # Show detailed messages

#Aqui para plotar
#simConfig.recordTraces = {'m': {'var': 'm', 'conds':{'pop': ['L2e', 'L2i']}}}
#simConfig.analysis['plotTraces'] = {'include':[('L2e', [0, 1, 2, 3]),('L2i', [0, 1])], 'timeRange': [0,100],'overlay': True,'oneFigPer': 'trace', 'showFig':False, 'saveFig': 'traceEscala3'+str(ScaleFactor)+'.png'}

simConfig.recordStep = 0.1         # Step size in ms to save data (e.g. V traces, LFP, etc)
simConfig.filename = simLabel  # Set file output name
simConfig.savePickle = False         # Save params, network and sim output to pickle file
simConfig.saveJson = True
simConfig.printSynsAfterRule = True

### Testar reducao de consumo memoria
simConfig.gatherOnlySimData = False #Original
simConfig.saveCellSecs=False
simConfig.saveCellConns=False

#simConfig.recordStim = False
#simConfig.addSynMechs=False
#compactConnFormat=True
simConfig.seeds['m'] = 123

#'include': ['all']
scale = 10
include = [(pop, list(range(0, netParams.popParams[pop]['numCells'], scale))) for pop in L]

#simConfig.analysis['plotRaster']={'include': include, 'timeRange': [100,600], 'popRates' : True , 'orderInverse':True, 'showFig':False, 'saveFig': 'rasterEscala8'+str(ScaleFactor)+'.png'}
#simConfig.analysis['plotRaster']={'include': include, 'timeRange': [100,600], 'popRates' : False , 'orderInverse':True, 'showFig':False, 'saveFig': 'rasterEscala8'+str(ScaleFactor)+'.png'}
simConfig.analysis['plotRaster']={'include': include, 'timeRange': [100,600], 'popRates' : False, 'figSize' : (5,10),  'labels':'overlay', 'orderInverse':True, 'showFig':False, 'saveFig': simLabel+'.png'}

###########################################AQUI######################################################################
simConfig.analysis['plotSpikeStats'] = {'include' : L, 'stats' : ['rate','isicv', 'sync'], 'timeRange' : [100,600],'showFig':False, 'saveFig': 'rasterEscala8'+str(ScaleFactor)+'.png'}
simConfig.printPopAvgRates = True
#####################################################################################################################



############################################################
#                    Create network and run simulation
############################################################

sim.initialize(
    simConfig = simConfig,  
    netParams = netParams)          # create network object and set cfg and net params
sim.net.createPops()                    # instantiate network populations
sim.net.createCells()                   # instantiate network cells based on defined populations

# randomize m parameter of cells
rand=h.Random()
for c in sim.net.cells:
	rand.Random123(c.gid, simConfig.seeds['m'])
	if c.tags['cellModel'] == 'IntFire_PD':
		c.hPointp.m = rand.normal(-58,10)

sim.net.addStims()              # add network stimulation

sim.cfg.createPyStruct = 0      # save memory by not saving py data structure for connections

sim.net.connectCells()                  # create connections between cells based on params
sim.setupRecording()                    # setup variables to record for each cell (spikes, V traces, etc)
sim.runSim()                            # run parallel Neuron simulation  
sim.gatherData()                        # gather spiking data and cell info from each node
sim.saveData()                          # save params, cell info and sim output to file (pickle,mat,txt,etc)#
sim.analysis.plotData()               # plot spike raster etc




