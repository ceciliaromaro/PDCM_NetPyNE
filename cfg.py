'''
NetPyNE version of Potjans and Diesmann thalamocortical network

cfg.py -- contains the simulation configuration (cfg object)

'''

from netpyne import specs

############################################################
#
#                    SIMULATION CONFIGURATION
#
############################################################

cfg = specs.SimConfig() # object of class SimConfig to store simulation configuration

############################################################
# Run options
############################################################

#cfg.seeds['stim']=6
cfg.seeds['stim'] = 0
cfg.duration = 60*1e3 #6*1e2   # Duration of the simulation, in ms
cfg.dt = 0.025          # Internal integration timestep to use
cfg.verbose = 0     # Show detailed messages
#cfg.seeds['m'] = 123
cfg.seeds['m'] = 0
cfg.printPopAvgRates = True
cfg.printRunTime = 1

### Options to save memory in large-scale ismulations
cfg.gatherOnlySimData = False #Original
cfg.saveCellSecs=False
cfg.saveCellConns=False
cfg.createPyStruct = False      # save memory by not saving py data structure for connections


###########################################################
# Network Options
###########################################################

# DC=True ;  TH=False; Balanced=True   => Reproduce Figure 7 A1 and A2
# DC=False;  TH=False; Balanced=False  => Reproduce Figure 7 B1 and B2
# DC=False ; TH=False; Balanced=True   => Reproduce Figure 8 A, B, C and D
# DC=False ; TH=False; Balanced=True   and run to 60 s to => Table 6 
# DC=False ; TH=True;  Balanced=True   => Figure 10A. But I want a partial reproduce so I guess Figure 10C is not necessary

# Size of Network. Adjust this constants, please!
cfg.ScaleFactor = 1.0 # 1.0 = 80.000 

# External input DC or Poisson
cfg.DC = False #True = DC // False = Poisson

# Thalamic input in 4th and 6th layer on or off
cfg.TH = False #True = on // False = off

# Balanced and Unbalanced external input as PD article
cfg.Balanced = True #True=Balanced // False=Unbalanced

cfg.simLabel = 'pd_scale-%s_DC-%d_TH-%d_Balanced-%d_dur-%d_Seed-%d'%(str(cfg.ScaleFactor), int(cfg.DC), int(cfg.TH), int(cfg.Balanced), int(cfg.duration/1e3), int(cfg.seeds['stim']))

###########################################################
# Recording and plotting options
###########################################################

cfg.recordStep = 0.1         # Step size in ms to save data (e.g. V traces, LFP, etc)
cfg.filename = cfg.simLabel  # Set file output name
cfg.saveFolder = 'data/'
cfg.savePickle = True         # Save params, network and sim output to pickle file
cfg.saveJson = False
cfg.recordStim = False
cfg.printSynsAfterRule = False
cfg.recordCellsSpikes = ['L2e', 'L2i', 'L4e', 'L4i', 'L5e', 'L5i','L6e', 'L6i'] # record only spikes of cells (not ext stims)

# raster plot (include update in netParams.py)
'''
cfg.analysis['plotRaster']={'include': [], 'timeRange': [100,600], 'popRates' : False, 'figSize' : (6,7),  
	'labels':'overlay', 'orderInverse': True, 'fontSize':16, 'showFig':False, 'saveFig': True, 'saveData':True}
# raster plot (include update in netParams.py)
'''
# cfg.analysis['plotRaster']={'include': [], 'timeRange': [100,cfg.duration], 'popRates' : False, 'figSize' : (6,7),  
# 	'labels':'overlay', 'orderInverse': True, 'fontSize':16, 'showFig':False, 'saveFig': True, 'saveData':True}

# statistics plot (include update in netParams.py)
cfg.analysis['plotSpikeStats'] = {'include' : [], 'stats' : ['rate','isicv'], 'legendLabels':['L6i', 'L6e', 'L5i', 'L5e', 'L4i', 'L4e','L2i', 'L2e'],
	'timeRange' : [100,cfg.duration], 'fontSize': 16, 'figSize': (6,9),'showFig':False, 'saveFig': True, 'saveData':True}

# plot traces
#cfg.recordTraces = {'m': {'var': 'm', 'conds':{'pop': ['L2e', 'L2i']}}}
#cfg.analysis['plotTraces'] = {'include':[('L2e', [0, 1, 2, 3]),('L2i', [0, 1])], 'timeRange': [0,100],'overlay': True,'oneFigPer': 'trace', 'showFig':False, 'saveFig': 'traceEscala3'+str(ScaleFactor)+'.png'}




