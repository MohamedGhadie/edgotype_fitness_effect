#----------------------------------------------------------------------------------------
# Calculate gene knockout effect for proteins carrying quasi-null mutations.
#----------------------------------------------------------------------------------------

import os
import io
import pickle
import pandas as pd
import numpy as np
from pathlib import Path
from stat_tools import t_test, sderror

def main():
    
    # reference interactome name: HI-II-14, HuRI, IntAct
    interactome_name = 'IntAct'
        
    # parent directory of all data files
    dataDir = Path('../data')
    
    # directory of processed data files specific to interactome
    interactomeDir = dataDir / interactome_name
    
    # input data files
    uniprotIDmapFile = dataDir / 'to_human_uniprotID_map.pkl'
    geneEffectFile = dataDir / 'Achilles_gene_effect.csv'
    naturalMutationsFile = interactomeDir / 'nondisease_mutation_edgotype.txt'
    diseaseMutationsFile = interactomeDir / 'disease_mutation_edgotype.txt'
    
    effectDictFile = dataDir / 'gene_effect.pkl'
    
    #------------------------------------------------------------------------------------
    # load mutations
    #------------------------------------------------------------------------------------
    
    print('Reading mutations')
    naturalMutations = pd.read_table (naturalMutationsFile, sep='\t')
    diseaseMutations = pd.read_table (diseaseMutationsFile, sep='\t')
    
    naturalMutations = naturalMutations [naturalMutations["edgotype"] == 'quasi-null'].reset_index(drop=True)
    diseaseMutations = diseaseMutations [diseaseMutations["edgotype"] == 'quasi-null'].reset_index(drop=True)
    
    print('Number of mutations:')
    print('non-disease mutations: %d' % len(naturalMutations))
    print('disease mutations: %d' % len(diseaseMutations))
    print()
    
    natMutationProteins = set(naturalMutations["protein"].tolist())
    disMutationProteins = set(diseaseMutations["protein"].tolist())
    overlapProteins = natMutationProteins & disMutationProteins
    
    print('Number of proteins carrying mutations:')
    print('non-disease mutations: %d' % len(natMutationProteins))
    print('disease mutations: %d' % len(disMutationProteins))
    print('overlapping proteins: %d' % len(overlapProteins))
    print()
    
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
            pickle.dump(effect, fOut)
    
    with open(effectDictFile, 'rb') as f:
        effect = pickle.load(f)
        
    natEffect = []
    for p in natMutationProteins:
        if p in effect:
            natEffect.extend([e for e in effect[p] if not np.isnan(e)])
    
    disEffect = []
    for p in disMutationProteins:
        if p in effect:
            disEffect.extend([e for e in effect[p] if not np.isnan(e)])
    
    allEffect = []
    for v in effect.values():
        allEffect.extend([e for e in v if not np.isnan(e)])
        
    # print results
    print('Average gene knockout effect for non-disease QN mutations: %.2f (SE = %g, n = %d)' 
                % (np.mean(natEffect), sderror(natEffect), len(natEffect)))
    
    print('Average gene knockout effect for disease QN mutations: %.2f (SE = %g, n = %d)' 
                % (np.mean(disEffect), sderror(disEffect), len(disEffect)))
    
    print('Average gene knockout effect for all genes: %.2f (SE = %g, n = %d)' 
                % (np.mean(allEffect), sderror(allEffect), len(allEffect)))
    
    print('Statistical significance, natural QN mutations and all genes')
    t_test(natEffect, allEffect)
    
    print('Statistical significance, disease QN mutations and all genes')
    t_test(disEffect, allEffect)

if __name__ == "__main__":
    main()
