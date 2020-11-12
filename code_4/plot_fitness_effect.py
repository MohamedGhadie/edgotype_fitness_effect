#----------------------------------------------------------------------------------------
# Plot the fitness effect for an edgotype among multiple interactomes.
#----------------------------------------------------------------------------------------

import os
import pickle
import numpy as np
from pathlib import Path
from plot_tools import multi_bar_plot

def main():
    
    # mutation edgotype for which fitness effect is calculated
    # options: quasi-null, edgetic, quasi-wild-type
    edgotype = 'edgetic'
    
    # reference interactome names
    interactome_names = ['HuRI', 'IntAct', 'experiment']
    
    # structural interactome names for plot labels
    legend = ['Y2H-SI', 'Lit-SI', 'Experiment']
    
    # fitness effects
    fitness_effects = ['Effectively neutral', 'Mildly deleterious', 'Strongly detrimental']
    
    # bar colors for different fitness effects
    barColors = {'Effectively neutral':'limegreen',
                 'Mildly deleterious':'orange',
                 'Strongly detrimental':'red'}
    
    # bar hatches for different interactomes
    barHatches = {'HuRI':'o', 'IntAct':'//', 'experiment':''}
    
    # plot confidence interval for the fraction of dispensable PPIs
    plotConfidenceIntervals = True
    
    # show figure legend
    showLeg = False
    
    # show figures
    showFigs = False
    
    # parent directory of all data files
    dataDir = Path('../data')
    dataDir = Path('/Volumes/MG_Samsung/edgotype_fitness_effect/data')
    
    # parent directory of processed data files
    procDir = dataDir / 'processed'
        
    # figure directory
    figDir = Path('../figures') / 'combined'
    
    # input data files
    inFile = '%s_mut_fitness_effect.pkl' % edgotype
    
    # create output directories if not existing
    if not figDir.exists():
        os.makedirs(figDir)
    
    allresults = {}
    for name in interactome_names:
        inPath = procDir / name / inFile
        with open(inPath, 'rb') as f:
            allresults[name] = pickle.load(f)
    
    allFE, allconf = [], []
    for name in interactome_names:
        results = allresults[name]
        if 'FE' in results:
            fe, conf = [], []
            for e in fitness_effects:
                fe.append(100 * results['FE'][e] if e in results['FE'] else np.nan)
                conf.append([100 * c for c in results['CI'][e]] if e in results['CI'] else (np.nan, np.nan))
        else:
            fe = [np.nan for e in fitness_effects]
            conf = [[np.nan, np.nan] for e in fitness_effects]
        allFE.append(fe)
        allconf.append(conf)
    
    multi_bar_plot (allFE,
                    errors = allconf if plotConfidenceIntervals else None,
                    xlabels = [f.replace(' ', '\n') for f in fitness_effects],
                    ylabel = 'Fraction (%)',
                    colors = [[barColors[e] for e in fitness_effects]] * len(interactome_names),
                    hatches = [barHatches[e] for e in interactome_names],
                    barwidth = 0.2,
                    linewidth = 1.5,
                    bargap = 0.05,
                    capsize = 8 if plotConfidenceIntervals else 0,
                    fmt = '.k',
                    msize = 16,
                    ewidth = 2,
                    edgecolor = 'black',
                    ecolors = 'k',
                    fontsize = 20,
                    ylim = [0, 100],
                    yMinorTicks = 4,
                    leg = legend if showLeg else None,
                    show = showFigs,
                    figdir = figDir,
                    figname = '%s_mut_fitness_effect' % edgotype)

if __name__ == "__main__":
    main()
