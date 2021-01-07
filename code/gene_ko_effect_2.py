#----------------------------------------------------------------------------------------
# Calculate gene knockout effect for proteins carrying quasi-null mutations.
#----------------------------------------------------------------------------------------

import os
import io
import pickle
import pandas as pd
import numpy as np
from pathlib import Path
from stat_tools import t_test, sderror, sderror_on_fraction
from plot_tools import multi_bar_plot

def main():
    
    # reference interactome names
    interactome_names = ['HuRI', 'IntAct', 'Experiment']
    
    # structural interactome names for plot labels
    xlabels = ['Y2H-SI', 'Lit-SI', 'Experiment']
    
    # fitness effects
    legend = ['Genes with non-pathogenic quasi-null mutations',
              'Genes with pathogenic quasi-null mutations',
              'All interactome genes']
    
    # bar colors for different fitness effects
    barColors = ['limegreen', 'orange', 'wheat']
    
    # bar hatches for different interactomes
    barHatches = ['o', '//', '']
    
    # show figure legend
    showLeg = True
    
    # show figures
    showFigs = False
    
    # parent directory of all data files
    dataDir = Path('../data')
    
    # figure directory
    figDir = Path('../figures') / 'combined'
    
    # input data files
    uniprotIDmapFile = dataDir / 'to_human_uniprotID_map.pkl'
    geneEffectFile = dataDir / 'Achilles_gene_effect.csv'
    essentialsFile = dataDir / 'Achilles_common_essentials.csv'
    naturalMutationsFile = 'nondisease_mutation_edgotype.txt'
    diseaseMutationsFile = 'disease_mutation_edgotype.txt'
    
    # output data files
    effectDictFile = dataDir / 'gene_effect.pkl'
    
    # create output directories if not existing
    if not figDir.exists():
        os.makedirs(figDir)
    
    #------------------------------------------------------------------------------------
    # Calculate gene knockout effect
    #------------------------------------------------------------------------------------
    
    with open(uniprotIDmapFile, 'rb') as f:
        uniprotID = pickle.load(f)
    
    if not effectDictFile.is_file():
        print('Producing gene effect dictionary')
        geneEffect = pd.read_table (geneEffectFile, sep='\t')
    
        with io.open(geneEffectFile, "r") as f:
            proteins = []
            headerLine = f.readline().split(',')
            for id in headerLine[1:]:
                name, entrezID = id.split()
                entrezID = entrezID[1:-1]
                if name in uniprotID:
                    proteins.append(uniprotID[name])
                elif entrezID in uniprotID:
                    proteins.append(uniprotID[entrezID])
                else:
                    proteins.append(name)
        
            cellLines = []
            effect = {p:[] for p in proteins}
            for line in f:
                linesplit = line.split(',')
                cellLines.append(linesplit[0])
                for p, e in zip(proteins, linesplit[1:]):
                    try:
                        effect[p].append(float(e))
                    except:
                        effect[p].append(np.nan)
        
        with open(effectDictFile, 'wb') as fOut:
            pickle.dump([cellLines, effect], fOut)
    
    with open(effectDictFile, 'rb') as f:
        cellLines, effect = pickle.load(f)
    
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
    # process mutations
    #------------------------------------------------------------------------------------
    
    natEffect, disEffect, geneEffect = [], [], []
    natSE, disSE, geneSE = [], [], []
    
    natEss, disEss, allEss = [], [], []
    natEssSE, disEssSE, allEssSE = [], [], []
    
    for name in interactome_names:
        if name is 'Experiment':
            interactome = pd.read_table (dataDir / name / 'reference_interactome.txt', sep='\t')
            interactome["Protein_1"] = interactome["Protein_1"].apply(lambda x: uniprotID[x[1:]]
                                                                                if x[1:] in uniprotID
                                                                                else x[1:])
            interactome["Protein_2"] = interactome["Protein_2"].apply(lambda x: uniprotID[x[1:]]
                                                                                if x[1:] in uniprotID
                                                                                else x[1:])
        else:
            interactome = pd.read_table (dataDir / name / 'structural_interactome.txt', sep='\t')
        
        naturalMutations = pd.read_table (dataDir / name / naturalMutationsFile, sep='\t')
        diseaseMutations = pd.read_table (dataDir / name / diseaseMutationsFile, sep='\t')
        
        if name is 'Experiment':
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
        
        naturalMutations = naturalMutations [naturalMutations["edgotype"] == 'quasi-null'].reset_index(drop=True)
        diseaseMutations = diseaseMutations [diseaseMutations["edgotype"] == 'quasi-null'].reset_index(drop=True)
        
        print('\n------------------------------------------------------------------------')
        print('Interactome: %s' % name)
        print('------------------------------------------------------------------------\n')
        
        print('non-disease mutations: %d' % len(naturalMutations))
        print('disease mutations: %d' % len(diseaseMutations))
        
        interactomeProteins = set(interactome[["Protein_1", "Protein_2"]].values.flatten())
        natMutationProteins = set(naturalMutations["protein"].tolist())
        disMutationProteins = set(diseaseMutations["protein"].tolist())
        overlapProteins = natMutationProteins & disMutationProteins
        
        print()
        print('all interactome proteins: %d' % len(interactomeProteins))
        print('proteins carrying non-disease mutations: %d' % len(natMutationProteins))
        print('proteins carrying disease mutations: %d' % len(disMutationProteins))
        print('overlapping proteins: %d' % len(overlapProteins))
    
        natE = []
        for p in natMutationProteins:
            if p in effect:
                natE.extend([-e for e in effect[p] if not np.isnan(e)])
    
        disE = []
        for p in disMutationProteins:
            if p in effect:
                disE.extend([-e for e in effect[p] if not np.isnan(e)])
    
        geneE = []
        for p in interactomeProteins:
            if p in effect:
                geneE.extend([-e for e in effect[p] if not np.isnan(e)])
        
        natEffect.append(np.mean(natE))
        disEffect.append(np.mean(disE))
        geneEffect.append(np.mean(geneE))
        
        natSE.append(sderror(natE))
        disSE.append(sderror(disE))
        geneSE.append(sderror(geneE))
        
        # print results
        print()
        print('Average gene knockout effect:')
        print('non-disease QN mutations: %.2f (SE = %g, n = %d)' % (natEffect[-1], natSE[-1], len(natE)))
        print('disease QN mutations: %.2f (SE = %g, n = %d)' % (disEffect[-1], disSE[-1], len(disE)))
        print('all interactome genes: %.2f (SE = %g, n = %d)' % (geneEffect[-1], geneSE[-1], len(geneE)))
        
        print()
        print('Statistical significance, disease QN mutations and non-disease QN mutations')
        t_test(natE, disE)
        
        print()
        print('Statistical significance, non-disease QN mutations and all interactome genes')
        t_test(natE, geneE)
        
        print()
        print('Statistical significance, disease QN mutations and all interactome genes')
        t_test(disE, geneE)
        
        natEssProteins = natMutationProteins & essentials
        disEssProteins = disMutationProteins & essentials
        allEssProteins = interactomeProteins & essentials
        
        natEss.append (len(natEssProteins) / len(natMutationProteins))
        disEss.append (len(disEssProteins) / len(disMutationProteins))
        allEss.append (len(allEssProteins) / len(interactomeProteins))
        
        natEssSE.append (sderror_on_fraction (len(natEssProteins), len(natMutationProteins)))
        disEssSE.append (sderror_on_fraction (len(disEssProteins), len(disMutationProteins)))
        allEssSE.append (sderror_on_fraction (len(allEssProteins), len(interactomeProteins)))
    return
#     genomeE = []
#     for v in effect.values():
#         genomeE.extend([-e for e in v if not np.isnan(e)])
#     
#     print()
#     print('Average knockout effect for all human genes: %.2f (SE = %g, n = %d)' % (np.mean(genomeE),
#                                                                                    sderror(genomeE),
#                                                                                    len(genomeE)))
    return
    multi_bar_plot ([natEffect, disEffect, geneEffect],
                    errors = [natSE, disSE, geneSE],
                    xlabels = xlabels,
                    ylabels = [0, 0.1, 0.2, 0.3],
                    ylabel = 'Gene knockout effect',
                    colors = barColors,
                    hatches = barHatches,
                    barwidth = 0.2,
                    linewidth = 1.5,
                    bargap = 0.05,
                    capsize = 8,
                    fmt = '.k',
                    msize = 14,
                    ewidth = 1,
                    edgecolor = 'black',
                    ecolors = 'k',
                    fontsize = 20,
                    bottomPos = 'zero',
                    ybounds = [0, 0.3],
                    leg = legend if showLeg else None,
                    show = showFigs,
                    figdir = figDir,
                    figname = 'qn_mut_gene_knockout_effect')
    
    multi_bar_plot ([natEss, disEss, allEss],
                    errors = [natEssSE, disEssSE, allEssSE],
                    xlabels = xlabels,
                    ylabels = [0, 0.1, 0.2, 0.3],
                    ylabel = 'Fraction of essential genes',
                    colors = barColors,
                    hatches = barHatches,
                    barwidth = 0.2,
                    linewidth = 1.5,
                    bargap = 0.05,
                    capsize = 8,
                    fmt = '.k',
                    msize = 14,
                    ewidth = 1,
                    edgecolor = 'black',
                    ecolors = 'k',
                    fontsize = 20,
                    leg = legend if showLeg else None,
                    show = showFigs,
                    figdir = figDir,
                    figname = 'qn_mut_gene_essential')

if __name__ == "__main__":
    main()
