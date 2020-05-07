import os
import pickle
import pandas as pd
from pathlib import Path
from math_tools import fitness_effect
from plot_tools import curve_plot, bar_plot

def main():
    
    # reference interactome name
    # options: HI-II-14, HuRI, IntAct
    interactome_name = 'IntAct'
    
    # mutation edgotype for which fitness effect is calculated
    # options: quasi-null, edgetic, quasi-wild-type
    edgotype = 'quasi-null'
    
    # homology modelling method used to create structural models
    # options: template_based, model_based
    model_method = 'model_based'
    
    # method used to predict edgetic perturbations
    # options: geometry, physics
    edgetic_method = 'physics'
    
    # method of calculating edgetic mutation ∆∆G
    # options: bindprofx, foldx
    edgetic_ddg = 'foldx'
    
    # assume edgotype probabilities of strongly detrimental (S) mutations to be similar to 
    # those of mildly deleterious (M) mutations. If set to False, strongly detrimental 
    # mutations are assumed to be all quasi-null
    assume_S_as_M = False
    
    # calculate confidence intervals
    computeConfidenceIntervals = True
    
    # confidence interval (%)
    CI = 95
    
    # calculate dispensable PPI content using fraction of mono-edgetic mutations 
    # instead of fraction of edgetic mutations
    mono_edgetic = False
    
    # show figures
    showFigs = False
    
    # parent directory of all data files
    #dataDir = Path('../data')
    dataDir = Path('/Volumes/MG_Samsung/edgotype_fitness_effect_full_model/data')
    
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
    
    numNaturalMut_type = {'Q': sum(naturalMutations["edgotype"] == 'quasi-null'),
                          'E': sum(naturalMutations["edgotype"] == 'edgetic'),
                          'W': sum(naturalMutations["edgotype"] == 'quasi-wild-type')}
    
    numDiseaseMut_type = {'Q': sum(diseaseMutations["edgotype"] == 'quasi-null'),
                          'E': sum(diseaseMutations["edgotype"] == 'edgetic'),
                          'W': sum(diseaseMutations["edgotype"] == 'quasi-wild-type')}
    
    numNaturalMut_considered = len(naturalMutations)
    numDiseaseMut_considered = len(diseaseMutations)
    
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
    posteriors = ['P(%s|%s)' % (p, T) for p in ('N','M','S')]
    prob = [100 * allresults[p] for p in posteriors]
    conf = []
    if computeConfidenceIntervals:
        for p in posteriors:
            if p + '_CI' in allresults:
                lower, upper = allresults[p + '_CI']
                conf.append((100 * lower, 100 * upper))
            else:
                conf.append((0, 0))
    
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
                xticklabels = ('Effectively\nneutral', 'Mildly\ndeleterious', 'Strongly\ndetrimental'),
                fontsize = 24,
                adjustBottom = 0.2,
                shiftBottomAxis = -0.1,
                xbounds = (1, 3),
                show = showFigs,
                figdir = figDir,
                figname = '%s_mut_fitness_effect_dot' % T)
    
    bar_plot (prob,
              error = conf if computeConfidenceIntervals else None,
              xlabels = ('Effectively\nneutral', 'Mildly\ndeleterious', 'Strongly\ndetrimental'),
              ylabels = None,
              ylabel = 'Fraction (%)',
              barwidth = 0.5,
              colors = ('limegreen', 'orange', 'red'),
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
              figname = '%s_mut_fitness_effect_bar' % T)

if __name__ == "__main__":
    main()
