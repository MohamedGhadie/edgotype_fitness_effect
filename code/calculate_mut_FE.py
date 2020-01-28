import os
import pickle
import pandas as pd
from pathlib import Path
from math_tools import fitness_effect
from plot_tools import curve_plot

def main():
    
    # reference interactome name
    # options: HI-II-14, IntAct
    interactome_name = 'HI-II-14'
    
    # mutation edgotype for which fitness effect is calculated
    # options: quasi-null, edgetic, quasi-wild-type
    edgotype = 'edgetic'
    
    # homology modelling method used to create structural models
    # options: template_based, model_based
    model_method = 'model_based'
    
    # method used to predict edgetic perturbations
    # options: geometry, physics
    edgetic_method = 'geometry'
    
    # assume edgotype probabilities of strongly detrimental (S) mutations to be similar to 
    # those of mildly deleterious (M) mutations. If set to False, strongly detrimental 
    # mutations are assumed to be all quasi-null
    assume_S_as_M = False
    
    # minimum change in protein free energy required for quasi-null mutations
    qnMinDDG = 5
    
    # calculate confidence intervals
    computeConfidenceIntervals = True
    
    # confidence interval (%)
    CI = 95
    
    # set to True to calculate dispensable PPI content using fraction of mono-edgetic mutations 
    # instead of edgetic mutations
    mono_edgetic = False
    
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
    methodDir = modellingDir / edgetic_method
    
    # figure directory
    figDir = Path('../figures') / interactome_name / model_method / edgetic_method / 'mutation_fitness_effect'
    
    if assume_S_as_M:
        figDir = figDir / 'assume_S_as_M'
    else:
        figDir = figDir / 'assume_S_quasi_null'
    
    # input data files
#     natMutLocFile = methodDir / 'nondisease_mutation_struc_loc.txt'
#     disMutLocFile = methodDir / 'disease_mutation_struc_loc.txt'
#     mutationPerturbsFile = interactomeDir / 'unique_mutation_perturbs_geometry.pkl' # new
#     natMutDdgFile = interactomeDir / 'nondisease_mutations_ddg.txt'
#     disMutDdgFile = interactomeDir / 'disease_mutations_ddg.txt'
    natMutEdgotypeFile = methodDir / 'nondisease_mutation_edgotype.txt'
    disMutEdgotypeFile = methodDir / 'disease_mutation_edgotype.txt'
    
    # output data files
    outputFile = methodDir / ('%s_mut_fitness_effect.pkl' % edgotype)
    
    # create output directories if not existing
    if not methodDir.exists():
        os.makedirs(methodDir)
    if not figDir.exists():
        os.makedirs(figDir)
    
    #------------------------------------------------------------------------------------
    # load mutations and count edgotypes
    #------------------------------------------------------------------------------------
    
    naturalMutations = pd.read_table (natMutLocFile, sep='\t')
    diseaseMutations = pd.read_table (disMutLocFile, sep='\t')
    
    naturalMutations = naturalMutations [naturalMutations["edgotype"] != '-'].reset_index(drop=True)
    diseaseMutations = diseaseMutations [diseaseMutations["edgotype"] != '-'].reset_index(drop=True)
    
    numNaturalMut_type = {'QN': sum(naturalMutations["edgotype"] == 'quasi-null'),
                          'E': sum(naturalMutations["edgotype"] == 'edgetic'),
                          'QW': sum(naturalMutations["edgotype"] == 'quasi-wild-type')}
    
    numDiseaseMut_type = {'QN': sum(diseaseMutations["edgotype"] == 'quasi-null'),
                          'E': sum(diseaseMutations["edgotype"] == 'edgetic'),
                          'QW': sum(diseaseMutations["edgotype"] == 'quasi-wild-type')}
    
    numNaturalMut_considered = sum(numNaturalMut_type.values())
    numDiseaseMut_considered = sum(numDiseaseMut_type.values())
    
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
        pT_S = pT_M
    else:
        pT_S = 1 if edgotype is 'quasi-null' else 0
    
    # edgotype abbreviations
    etype_symbol = {'quasi-null':'QN', 'edgetic':'E', 'quasi-wild-type':'QW'}
    T = etype_symbol [edgotype]
    
    allresults = fitness_effect (pN,
                                 pM,
                                 pS,
                                 numNaturalMut_type[T],
                                 numNaturalMut_considered,
                                 numDiseaseMut_type[T],
                                 numDiseaseMut_considered,
                                 pT_S = pT_S,
                                 edgotype = 'edgetic',
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
                msize = 16,
                ewidth = 1.25,
                ecolors = 'k',
                ylabel = 'Fraction of %s mutations (%%)' % edgotype,
                yMinorTicks = 4,
                xticks = [1, 2, 3],
                xticklabels = ('Effectively\nneutral', 'Mildly\ndeleterious', 'Strongly\ndetrimental'),
                fontsize = 20,
                adjustBottom = 0.2,
                shiftBottomAxis = -0.1,
                xbounds = (1, 3),
                show = showFigs,
                figdir = figDir,
                figname = '%s_mut_fitness_effect' % T)

if __name__ == "__main__":
    main()
