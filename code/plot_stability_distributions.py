#----------------------------------------------------------------------------------------
# Plot distributions for folding ∆∆G upon mutation.
#----------------------------------------------------------------------------------------

import numpy as np
import scipy.stats as stats
from pathlib import Path
from energy_tools import read_protein_mutation_ddg
from stat_tools import normality_test
from plot_tools import multi_histogram_plot

def main():
    
    # reference interactome name: HI-II-14, HuRI, IntAct
    interactome_name = 'HuRI'
    
    # homology modelling method used to create structural models
    # options: template_based, model_based
    model_method = 'model_based'
    
    # method used to perform calculations; 'geometry' or 'physics'
    edgetic_method = 'physics'
    
    # method that was used to calculate edgetic mutation ∆∆G
    # options: bindprofx, foldx
    edgetic_ddg = 'foldx'
    
    # show figures
    showFigs = False
    
    # parent directory of all data files
    dataDir = Path('../data')
    
    # parent directory of all processed data files
    procDir = dataDir / 'processed'
    
    # directory of processed data files specific to interactome
    interactomeDir = procDir / interactome_name
    
    # directory of processed model-related data files specific to interactome
    modellingDir = interactomeDir / model_method
    
    # directory of edgetic mutation calculation method
    edgeticDir = modellingDir / edgetic_method
    
    # figure directory
    figDir = Path('../figures') / interactome_name / model_method/ edgetic_method
    
    if edgetic_method is 'physics':
        edgeticDir = edgeticDir / (edgetic_ddg + '_edgetics')
        figDir = figDir / ('%s_edgetics' % edgetic_ddg)
    
    # input data files
    natMutDdgFile = edgeticDir / 'nondis_mut_folding_ddg_foldx.txt'
    disMutDdgFile = edgeticDir / 'dis_mut_folding_ddg_foldx.txt'
    
    # create output directories if not existing
    if not figDir.exists():
        os.makedirs(figDir)
    
    #------------------------------------------------------------------------------------
    # load mutation folding ∆∆G
    #------------------------------------------------------------------------------------
    
    natMutDDG = read_protein_mutation_ddg (natMutDdgFile, type = 'folding')
    disMutDDG = read_protein_mutation_ddg (disMutDdgFile, type = 'folding')
    
    natMutDDG = [ddg for _, _, _, ddg in natMutDDG.values() if not np.isnan(ddg)]
    disMutDDG = [ddg for _, _, _, ddg in disMutDDG.values() if not np.isnan(ddg)]
    
    #------------------------------------------------------------------------------------
    # Change in protein stability in relation to mutation distance to center
    #------------------------------------------------------------------------------------
    
    print('\n' + 'Mutations with calculated stability effect:')
    print('non-disease mutations: %d' % len(natMutDDG))
    print('disease mutations: %d' % len(disMutDDG))
    
    print('\n' + 'Non-disease mutations')
    normality_test (natMutDDG, 0.05)
    print('skew = %f' % stats.skew(natMutDDG))
    print('\n' + 'Disease mutations')
    normality_test (disMutDDG, 0.05)
    print('skew = %f' % stats.skew(disMutDDG))
    
    numbins = 30
    multi_histogram_plot (natMutDDG,
                          'c',
                          xlabel = 'Change in protein stability (kcal/mol)',
                          ylabel = 'Number of mutations',
                          edgecolor = 'k',
                          bins = numbins,
                          alpha = 1,
                          xlim = [-10, 50],
                          show = showFigs,
                          figdir = figDir,
                          figname = 'nondisease_mut_ddg_distribution')
    multi_histogram_plot (disMutDDG,
                          'm',
                          xlabel = 'Change in protein stability (kcal/mol)',
                          ylabel = 'Number of mutations',
                          edgecolor = 'k',
                          bins = numbins,
                          alpha = 1,
                          xlim = [-10, 50],
                          show = showFigs,
                          figdir = figDir,
                          figname = 'disease_mut_ddg_distribution')
    multi_histogram_plot ([disMutDDG, natMutDDG],
                          ['m', 'c'],
                          xlabel = 'Change in protein stability (kcal/mol)',
                          ylabel = 'Number of mutations',
                          leg = ('Disease mutations', 'Non-disease mutations'),
                          edgecolor = 'k',
                          bins = numbins,
                          alpha = 0.5,
                          xlim = [-10, 50],
                          show = showFigs,
                          figdir = figDir,
                          figname = 'ddg_distribution')

if __name__ == "__main__":
    main()
