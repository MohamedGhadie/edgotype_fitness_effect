#----------------------------------------------------------------------------------------
# Calculate the the fitness effect for different mutation edgotypes, i.e., probabilities 
# for edgetic (E), quasi-null (Q) or quasi-wild-type (W) mutations to be effectively 
# neutral (N), mildly deleterious (M) or strongly detrimental (S).
#----------------------------------------------------------------------------------------

import os
import pickle
import pandas as pd
from pathlib import Path
from math_tools import fitness_effect
from plot_tools import curve_plot, bar_plot

def main():
    
    # reference interactome name: HI-II-14, HuRI, IntAct
    interactome_name = 'IntAct'
    
    # mutation edgotype for which fitness effect is calculated
    # options: quasi-null, edgetic, quasi-wild-type
    edgotype = 'quasi-wild-type'
    
    # homology modelling method used to create structural models
    # options: template_based, model_based
    model_method = 'model_based_4'
    
    # method used to predict edgetic perturbations
    # options: geometry, physics
    edgetic_method = 'physics'
    
    # method of calculating edgetic mutation ∆∆G
    # options: bindprofx, foldx
    edgetic_ddg = 'foldx'
    
    # possible fitness effects
    fitness_effects = ['Effectively neutral', 'Mildly deleterious', 'Strongly detrimental']
    
    # bar colors for different fitness effects
    barColors = {'Effectively neutral':'limegreen',
                 'Mildly deleterious':'orange',
                 'Strongly detrimental':'red'}
    
    # assume edgotype probabilities of strongly detrimental (S) mutations to be similar to 
    # those of mildly deleterious (M) mutations. If False, strongly detrimental 
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
    
    # parent directory of all processed data files
    procDir = dataDir / 'processed'
    
    # directory of processed data files specific to interactome
    interactomeDir = procDir / interactome_name
    
    # directory of processed model-related data files specific to interactome
    modellingDir = interactomeDir / model_method
    
    # directory of calculation method
    edgeticDir = modellingDir / edgetic_method
    
    # figure directory
    figDir = Path('../figures') / interactome_name / model_method / edgetic_method
    
    if edgetic_method is 'physics':
        edgeticDir = edgeticDir / (edgetic_ddg + '_edgetics')
        figDir = figDir / (edgetic_ddg + '_edgetics') / 'mutation_fitness_effect'
    elif edgetic_method is 'geometry':
        figDir = figDir / 'mutation_fitness_effect'
    
    if assume_S_as_M:
        figDir = figDir / 'assume_S_as_M'
    else:
        figDir = figDir / 'assume_S_quasi_null'
    
    # input data files
    natMutEdgotypeFile = edgeticDir / 'nondisease_mutation_edgotype.txt'
    disMutEdgotypeFile = edgeticDir / 'disease_mutation_edgotype.txt'
    
    # output data files
    outputFile = edgeticDir / ('%s_mut_fitness_effect.pkl' % edgotype)
    
    # create output directories if not existing
    if not edgeticDir.exists():
        os.makedirs(edgeticDir)
    if not figDir.exists():
        os.makedirs(figDir)
    
    #------------------------------------------------------------------------------------
    # load mutations and count edgotypes
    #------------------------------------------------------------------------------------
    
    naturalMutations = pd.read_table (natMutEdgotypeFile, sep='\t')
    diseaseMutations = pd.read_table (disMutEdgotypeFile, sep='\t')
    
    naturalMutations = naturalMutations [naturalMutations["edgotype"] != '-'].reset_index(drop=True)
    diseaseMutations = diseaseMutations [diseaseMutations["edgotype"] != '-'].reset_index(drop=True)
    
    numNaturalMut_type = {'quasi-null': sum(naturalMutations["edgotype"] == 'quasi-null'),
                          'edgetic': sum(naturalMutations["edgotype"] == 'edgetic'),
                          'quasi-wild-type': sum(naturalMutations["edgotype"] == 'quasi-wild-type')}
    
    numDiseaseMut_type = {'quasi-null': sum(diseaseMutations["edgotype"] == 'quasi-null'),
                          'edgetic': sum(diseaseMutations["edgotype"] == 'edgetic'),
                          'quasi-wild-type': sum(diseaseMutations["edgotype"] == 'quasi-wild-type')}
    
    numNaturalMut_considered = len(naturalMutations)
    numDiseaseMut_considered = len(diseaseMutations)
    
    #------------------------------------------------------------------------------------
    # Apply Bayes' theorem to calculate the probabilities for a mutation of selected
    # edgotype to be effectively neutral, mildly deleterious or strongly detrimental
    #------------------------------------------------------------------------------------
    
    # Probability for new missense mutations to be neutral (N)
    pN = 0.27

    # Probability for new missense mutations to be mildly deleterious (M)
    pM = 0.53

    # Probability for new missense mutations to be strongly detrimental (S)
    pS = 0.20
    
    # Probability for strongly detrimental mutations (S) to be of selected edgotype (T)
    if assume_S_as_M:
        pT_S = numDiseaseMut_type[edgotype] / numDiseaseMut_considered
    else:
        pT_S = 1 if edgotype is 'quasi-null' else 0
    
    allresults = fitness_effect (pN,
                                 pM,
                                 pS,
                                 numNaturalMut_type[edgotype],
                                 numNaturalMut_considered,
                                 numDiseaseMut_type[edgotype],
                                 numDiseaseMut_considered,
                                 pT_S = pT_S,
                                 edgotype = edgotype,
                                 CI = 95,
                                 output = True)
    
    with open(outputFile, 'wb') as fOut:
        pickle.dump(allresults, fOut)
    
    # plot fitness effect
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
