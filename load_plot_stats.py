from netpyne import sim

filenames = ['data/trials/trials_0_0.pkl',
            'data/trials/trials_1_1.pkl',
            'data/trials/trials_2_2.pkl',
            'data/trials/trials_3_3.pkl',
            'data/trials/trials_4_4.pkl',
            'data/trials/trials_5_5.pkl',
            'data/trials/trials_6_6.pkl',
            'data/trials/trials_7_7.pkl',
            'data/trials/trials_8_8.pkl',
            'data/trials/trials_9_9.pkl']

allpops = ['L2e', 'L2i', 'L4e', 'L4i', 'L5e', 'L5i','L6e', 'L6i']

for filename in filenames:
    try:
        sim.load(filename, instantiate=0)
    except:
        pass
    sim.cfg.createNEURONObj=0
    sim.net.createCells()
    sim.analysis.plotSpikeStats(include=allpops, legendLabels=allpops, timeRange=[100, 60000], stats=['rate', 'isicv'], 
                                fontSize= 16, figSize= (6,9), saveFig=1, saveData=1, showFig=0)


