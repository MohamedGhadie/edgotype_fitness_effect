#----------------------------------------------------------------------------------------
# Calculate the the fitness effect for different mutation edgotypes, specifically among 
# essential genes.
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
    edgotype = 'quasi-null'
    
    # possible fitness effects
    fitness_effects = ['Effectively neutral', 'Mildly deleterious', 'Strongly detrimental']
    
    # bar colors for different fitness effects
    barColors = {'Effectively neutral':'limegreen',
                 'Mildly deleterious':'orange',
                 'Strongly detrimental':'red'}
    
    # assume edgotype probabilities of strongly detrimental (S) mutations to be similar to 
    # those of mildly deleterious (M) mutations. If False, strongly detrimental 
    # mutations are assumed to be all quasi-null
    assume_S_as_M = False
    
    Pe = {'HuRI':0.21, 'IntAct':0.22, 'experiment':0.11}
    
    # calculate confidence intervals
    computeConfidenceIntervals = True
    
    # confidence interval (%)
    CI = 95
    
    # show figures
    showFigs = False
    
    # parent directory of all data files
    dataDir = Path('../data')
    
    interactomeDir = dataDir / interactome_name
    
    # figure directory
    figDir = Path('../figures') / interactome_name / 'essential_genes'
    
    if assume_S_as_M:
        figDir = figDir / 'assume_S_as_M'
    else:
        figDir = figDir / 'assume_S_quasi_null'
    
    # input data files
    uniprotIDmapFile = dataDir / 'to_human_uniprotID_map.pkl'
    essentialsFile = dataDir / 'Achilles_common_essentials.csv'
    natMutEdgotypeFile = interactomeDir / 'nondisease_mutation_edgotype.txt'
    disMutEdgotypeFile = interactomeDir / 'disease_mutation_edgotype.txt'
    
    # output data files
    outputFile = interactomeDir / ('%s_mut_fitness_effect.pkl' % edgotype)
    
    # create output directories if not existing
    if not interactomeDir.exists():
        os.makedirs(interactomeDir)
    if not figDir.exists():
        os.makedirs(figDir)
    
    with open(uniprotIDmapFile, 'rb') as f:
        uniprotID = pickle.load(f)
    
    with open(essentialsFile, 'r') as fin:
        essentialGenes = list(fin.read().split('\n'))[1:]
    
    essentials = set()
    for id in essentialGenes:
        idsplit = id.split()
        if len(idsplit) == 2:
            name, entrezID = id.split()
            entrezID = entrezID[1:-1]
            if name in uniprotID:
                essentials.add(uniprotID[name])
            elif entrezID in uniprotID:
                essentials.add(uniprotID[entrezID])
            else:
                essentials.add(name)
    
    #------------------------------------------------------------------------------------
    # load mutations and count edgotypes
    #------------------------------------------------------------------------------------
    
    naturalMutations = pd.read_table (natMutEdgotypeFile, sep='\t')
    diseaseMutations = pd.read_table (disMutEdgotypeFile, sep='\t')
    
    if interactome_name is 'experiment':
        naturalMutations["edgotype"] = naturalMutations["Edgotype_class"].apply(lambda x: x.lower())
        diseaseMutations["edgotype"] = diseaseMutations["Edgotype_class"].apply(lambda x: x.lower())
        naturalMutations["Entrez_Gene_ID"] = naturalMutations["Entrez_Gene_ID"].apply(str)
        diseaseMutations["Entrez_Gene_ID"] = diseaseMutations["Entrez_Gene_ID"].apply(str)
        naturalMutations["protein"] = naturalMutations["Entrez_Gene_ID"].apply(lambda x: uniprotID[x]
                                                                                if x in uniprotID
                                                                                else x)
        diseaseMutations["protein"] = diseaseMutations["Entrez_Gene_ID"].apply(lambda x: uniprotID[x]
                                                                                if x in uniprotID
                                                                                else x)
    
    naturalMutations = naturalMutations [naturalMutations["edgotype"] != '-'].reset_index(drop=True)
    diseaseMutations = diseaseMutations [diseaseMutations["edgotype"] != '-'].reset_index(drop=True)
    
    natEss = naturalMutations["protein"].apply(lambda x: x in essentials)
    disEss = diseaseMutations["protein"].apply(lambda x: x in essentials)
    
    natQ = naturalMutations["edgotype"] == 'quasi-null'
    disQ = diseaseMutations["edgotype"] == 'quasi-null'
    natE = naturalMutations["edgotype"] == 'edgetic'
    disE = diseaseMutations["edgotype"] == 'edgetic'
    natW = naturalMutations["edgotype"] == 'quasi-wild-type'
    disW = diseaseMutations["edgotype"] == 'quasi-wild-type'
    
    numNaturalMut_type = {'quasi-null': sum(natQ & natEss),
                          'edgetic': sum(natE & natEss),
                          'quasi-wild-type': sum(natW | (natEss == False))}
    
    numDiseaseMut_type = {'quasi-null': sum(disQ & disEss),
                          'edgetic': sum(disE & disEss),
                          'quasi-wild-type': sum(disW | (disEss == False))}
    
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
        pT_S = Pe[interactome_name] if edgotype is 'quasi-null' else 0
    
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