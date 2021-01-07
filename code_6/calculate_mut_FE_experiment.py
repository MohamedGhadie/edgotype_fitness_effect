#----------------------------------------------------------------------------------------
# Calculate the fitness effect for different mutation edgotypes in the experimental dataset 
# of Sahni et al. (2015), i.e., probabilities for edgetic (E), quasi-null (Q) and 
# quasi-wild-type (W) mutations to be effectively neutral (N), mildly deleterious (M) or 
# strongly detrimental (S).
#----------------------------------------------------------------------------------------

import os
import pickle
import pandas as pd
import numpy as np
from pathlib import Path
from math_tools import fitness_effect
from stat_tools import fisher_test, sderror_on_fraction
from plot_tools import pie_plot, curve_plot, bar_plot

def main():
        
    # reference interactome name
    interactome_name = 'experiment'
    
    # mutation edgotype for which fitness effect is calculated
    # options: quasi-null, edgetic, quasi-wild-type
    edgotype = 'edgetic'
    
    # possible fitness effects
    fitness_effects = ['Effectively neutral', 'Mildly deleterious', 'Strongly detrimental']
    
    # bar colors for different fitness effects
    barColors = {'Effectively neutral':'limegreen',
                 'Mildly deleterious':'orange',
                 'Strongly detrimental':'red'}
    
    # assume edgotype probabilities of strongly detrimental (S) mutations to be similar to 
    # those of mildly deleterious (M) mutations. If set to False, strongly detrimental 
    # mutations are assumed to be all quasi-null
    assume_S_as_M = True
    
    # calculate confidence intervals
    computeConfidenceIntervals = True
    
    # confidence interval (%)
    CI = 95
    
    # show figures
    showFigs = False
    
    # parent directory of all data files
    dataDir = Path('../data')
    
    # directory of data files from external sources
    extDir = dataDir / 'external'
    
    # parent directory of all processed data files
    procDir = dataDir / 'processed'
    
    # directory of processed data files specific to interactome
    interactomeDir = procDir / interactome_name
        
    # figure directory
    figDir = Path('../figures') / interactome_name / 'mutation_fitness_effect'
    
    if assume_S_as_M:
        figDir = figDir / 'assume_S_as_M'
    else:
        figDir = figDir / 'assume_S_quasi_null'
    
    # input data files
    mutationsFile = extDir / "Sahni_2015_Table_S3.xlsx"
    
    # output data files
    outputFile = interactomeDir / ('%s_mut_fitness_effect.pkl' % edgotype)
    
    # create output directories if not existing
    if not interactomeDir.exists():
        os.makedirs(interactomeDir)
    if not figDir.exists():
        os.makedirs(figDir)
    
    #------------------------------------------------------------------------------------
    # Load experimental edgotype data
    #------------------------------------------------------------------------------------
    
    expMutations = pd.read_excel (mutationsFile, sheet_name = 'Table S3C')
    expMutations_ppi = pd.read_excel (mutationsFile, sheet_name = 'Table S3A')
    
    expMutations = expMutations [expMutations["Edgotype_class"].apply(
                                 lambda x: x in {'Edgetic', 'Quasi-null', 'Quasi-wild-type'})
                                ].reset_index(drop = True)
    
    expMutations["partners"] = expMutations["Allele_ID"].apply(
        lambda x: expMutations_ppi.loc[expMutations_ppi["Allele_ID"] == x, "Interactor_Gene_ID"].values)
    expMutations["Mut_interaction"] = expMutations["Allele_ID"].apply(
        lambda x: expMutations_ppi.loc[expMutations_ppi["Allele_ID"] == x, "Y2H_score" ].values)
    
    expMutations["WT_interaction"] = expMutations.apply(
        lambda x: np.array([expMutations_ppi.loc[(expMutations_ppi["Category"] == 'Wild-type') & 
                                                 (expMutations_ppi["Entrez_Gene_ID"] == x["Entrez_Gene_ID"]) & 
                                                 (expMutations_ppi["Interactor_Gene_ID"] == p),
                                                 "Y2H_score"].item()
                            for p in x["partners"]]), axis=1)
    
    expMutations["perturbations"] = expMutations.apply(
        lambda x:  0 + ((x["WT_interaction"] == 1) & (x["Mut_interaction"] == 0)), axis=1)
    
    naturalMutations = expMutations [expMutations["Category"] == 'Non-disease variant'].reset_index(drop = True)
    diseaseMutations = expMutations [expMutations["Category"] == 'Disease mutation'].reset_index(drop = True)
    
    #------------------------------------------------------------------------------------
    # Fraction of quasi-null, edgetic and quasi-wild-type mutations
    #------------------------------------------------------------------------------------
    
    numNaturalMut_type = {'Q': sum(naturalMutations["Edgotype_class"] == 'Quasi-null'),
                          'E': sum(naturalMutations["Edgotype_class"] == 'Edgetic'),
                          'W': sum(naturalMutations["Edgotype_class"] == 'Quasi-wild-type')}
    
    numDiseaseMut_type = {'Q': sum(diseaseMutations["Edgotype_class"] == 'Quasi-null'),
                          'E': sum(diseaseMutations["Edgotype_class"] == 'Edgetic'),
                          'W': sum(diseaseMutations["Edgotype_class"] == 'Quasi-wild-type')}
    
    numNaturalMut_considered = sum(numNaturalMut_type.values())
    numDiseaseMut_considered = sum(numDiseaseMut_type.values())
    
    fracNaturalMut_type = {'Q': numNaturalMut_type['Q'] / numNaturalMut_considered,
                           'E': numNaturalMut_type['E'] / numNaturalMut_considered,
                           'W': numNaturalMut_type['W'] / numNaturalMut_considered}
    fracDiseaseMut_type = {'Q': numDiseaseMut_type['Q'] / numDiseaseMut_considered,
                           'E': numDiseaseMut_type['E'] / numDiseaseMut_considered,
                           'W': numDiseaseMut_type['W'] / numDiseaseMut_considered}
    
    print()
    print('Fraction of edgotypes among non-disease mutations:')
    print('quasi-null: %f%% (SE = %g, %d out of %d)' % (100 * fracNaturalMut_type['Q'],
                                                        sderror_on_fraction (numNaturalMut_type['Q'], numNaturalMut_considered),
                                                        numNaturalMut_type['Q'],
                                                        numNaturalMut_considered))
    print('edgetic: %f%% (SE = %g, %d out of %d)' % (100 * fracNaturalMut_type['E'],
                                                     sderror_on_fraction (numNaturalMut_type['E'], numNaturalMut_considered),
                                                     numNaturalMut_type['E'],
                                                     numNaturalMut_considered))
    print('quasi-wild-type: %f%% (SE = %g, %d out of %d)' % (100 * fracNaturalMut_type['W'],
                                                             sderror_on_fraction (numNaturalMut_type['W'], numNaturalMut_considered),
                                                             numNaturalMut_type['W'],
                                                             numNaturalMut_considered))
    
    print()
    print('Fraction of edgotypes among disease mutations:')
    print('quasi-null: %f%% (SE = %g, %d out of %d)' % (100 * fracDiseaseMut_type['Q'],
                                                        sderror_on_fraction (numDiseaseMut_type['Q'], numDiseaseMut_considered),
                                                        numDiseaseMut_type['Q'],
                                                        numDiseaseMut_considered))
    print('edgetic: %f%% (SE = %g, %d out of %d)' % (100 * fracDiseaseMut_type['E'],
                                                     sderror_on_fraction (numDiseaseMut_type['E'], numDiseaseMut_considered),
                                                     numDiseaseMut_type['E'],
                                                     numDiseaseMut_considered))
    print('quasi-wild-type: %f%% (SE = %g, %d out of %d)' % (100 * fracDiseaseMut_type['W'],
                                                             sderror_on_fraction (numDiseaseMut_type['W'], numDiseaseMut_considered),
                                                             numDiseaseMut_type['W'],
                                                             numDiseaseMut_considered))
    
    fisher_test ([numNaturalMut_type['Q'] + numNaturalMut_type['E'], numNaturalMut_type['W']],
                 [numDiseaseMut_type['Q'] + numDiseaseMut_type['E'], numDiseaseMut_type['W']])
        
    pie_plot([numNaturalMut_type['W'], numNaturalMut_type['E'], numNaturalMut_type['Q']],
             angle = 90,
             colors = ['mediumslateblue', 'purple', 'red'],
             edgewidth = 2,
             show = showFigs,
             figdir = figDir,
             figname = 'non_disease_mut_edgotype_experiment')
    pie_plot([numDiseaseMut_type['W'], numDiseaseMut_type['E'], numDiseaseMut_type['Q']],
             angle = 90,
             colors = ['mediumslateblue', 'purple', 'red'],
             edgewidth = 2,
             show = showFigs,
             figdir = figDir,
             figname = 'disease_mut_edgotype_experiment')
        
    #------------------------------------------------------------------------------------
    # Apply Bayes' theorem to calculate the probabilities for a mutation of selected
    # edgotype to be effectively neutral, mildly deleterious or strongly detrimental
    #------------------------------------------------------------------------------------
    
    # edgotype abbreviations
    etype_symbol = {'quasi-null':'Q', 'edgetic':'E', 'quasi-wild-type':'W'}
    T = etype_symbol [edgotype]
    
    # Probability for new missense mutations to be neutral (N)
    pN = 0.27

    # Probability for new missense mutations to be mildly deleterious (M)
    pM = 0.53

    # Probability for new missense mutations to be strongly detrimental (S)
    pS = 0.20
    
    # Probability for strongly detrimental mutations (S) to be of selected edgotype (T)
    if assume_S_as_M:
        pT_S = numDiseaseMut_type[T] / numDiseaseMut_considered
    else:
        pT_S = 1 if edgotype is 'quasi-null' else 0
    
    allresults = fitness_effect (pN,
                                 pM,
                                 pS,
                                 numNaturalMut_type[T],
                                 numNaturalMut_considered,
                                 numDiseaseMut_type[T],
                                 numDiseaseMut_considered,
                                 pT_S = pT_S,
                                 edgotype = edgotype,
                                 CI = 95,
                                 output = True)
    with open(outputFile, 'wb') as fOut:
        pickle.dump(allresults, fOut)
    
    # plot fitness effect
#     posteriors = ['P(%s|%s)' % (p, T) for p in ('N','M','S')]
#     prob = [100 * allresults[p] for p in posteriors]
#     conf = []
#     if computeConfidenceIntervals:
#         for p in posteriors:
#             if p + '_CI' in allresults:
#                 lower, upper = allresults[p + '_CI']
#                 conf.append((100 * lower, 100 * upper))
#             else:
#                 conf.append((0, 0))
    
    prob = [100 * allresults['FE'][f] for f in fitness_effects]
    if computeConfidenceIntervals and ('CI' in allresults):
        conf = []
        for f in fitness_effects:
            lower, upper = allresults['CI'][f]
            conf.append((100 * lower, 100 * upper))
    
    curve_plot (prob,
                error = conf if computeConfidenceIntervals else None,
                xlim = [0.8, 3.1],
                ylim = [0, 100],
                styles = '.k',
                capsize = 10 if computeConfidenceIntervals else 0,
                msize = 28,
                ewidth = 3,
                ecolors = 'k',
                ylabel = 'Fraction (%)',
                yMinorTicks = 4,
                xticks = [1, 2, 3],
                xticklabels = [f.replace(' ', '\n') for f in fitness_effects],
                fontsize = 24,
                adjustBottom = 0.2,
                shiftBottomAxis = -0.1,
                xbounds = (1, 3),
                show = showFigs,
                figdir = figDir,
                figname = '%s_mut_fitness_effect_dot' % edgotype)
    
    bar_plot (prob,
              error = conf if computeConfidenceIntervals else None,
              xlabels = [f.replace(' ', '\n') for f in fitness_effects],
              ylabels = None,
              ylabel = 'Fraction (%)',
              barwidth = 0.5,
              colors = [barColors[f] for f in fitness_effects],
              capsize = 10 if computeConfidenceIntervals else 0,
              fmt = '.k',
              msize = 24,
              ewidth = 3,
              edgecolor = 'black',
              ecolors = 'k',
              fontsize = 21,
              xlim = None,
              ylim = [0, 100],
              xticks = None,
              yMinorTicks = 4,
              adjustBottom = False,
              shiftBottomAxis = None,
              xbounds = None,
              show = showFigs,
              figdir = figDir,
              figname = '%s_mut_fitness_effect_bar' % edgotype)

if __name__ == "__main__":
    main()
