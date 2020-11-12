#----------------------------------------------------------------------------------------
# Produce mutation overlap in two interactomes.
#----------------------------------------------------------------------------------------

import os
import pickle
import pandas as pd
from pathlib import Path

def main():
    
    # names of reference interactomes
    # options: any pair from HI-II-14, HuRI, IntAct, experiment
    interactome_names = ['IntAct', 'experiment']
    
    # homology modelling method used to create structural models
    # options: template_based, model_based
    model_method = 'model_based'
    
    # method used to perform calculations; 'geometry' or 'physics'
    edgetic_method = 'physics'
    
    # method that was used to calculate edgetic mutation ∆∆G
    # options: bindprofx, foldx
    edgetic_ddg = 'foldx'
        
    # name of gene name column to be used if protein column not found
    gene_col = 'Symbol'
    
    # name of mutation column to be used if processed mut_position column not found
    mut_RefSeq_col = 'Mutation_RefSeq_AA'
    
    # show figures
    showFigs = False
    
    # parent directory of all data files
    dataDir = Path('../data')
    
    # parent directory of all processed data files
    procDir = dataDir / 'processed'
    
    # figure directory
    figDir = Path('../figures') / 'combined'
    
    # input data files
    uniprotIDmapFile = procDir / 'to_human_uniprotID_map.pkl'
    
    natMutFile = 'nondisease_mutation_edgotype.txt'
    disMutFile = 'disease_mutation_edgotype.txt'
    
    # output data files
    natMutOutFile = procDir / ('natural_mutation_overlap_%s' % '_'.join(interactome_names))
    disMutOutFile = procDir / ('disease_mutation_overlap_%s' % '_'.join(interactome_names))
    
    # create output directories if not existing
    if not procDir.exists():
        os.makedirs(procDir)
    
    with open(uniprotIDmapFile, 'rb') as f:
        uniprotID = pickle.load(f)
    
    #----------------------------------------------------------------------------------------
    # process mutations in both interactomes
    #----------------------------------------------------------------------------------------
    
    for name in interactome_names:
        if name is 'experiment':
            fileDir = procDir / name
        else:
            fileDir = procDir / name / model_method / edgetic_method / (edgetic_ddg + '_edgetics')
        
        naturalMutations = pd.read_table (fileDir / natMutFile, sep='\t')
        diseaseMutations = pd.read_table (fileDir / disMutFile, sep='\t')
    
        if 'protein' not in naturalMutations.columns:
            naturalMutations["protein"] = naturalMutations[gene_col].apply(lambda x:
                                                                            uniprotID[x] if x in uniprotID 
                                                                            else x + 'unknown')
        if 'protein' not in diseaseMutations.columns:
            diseaseMutations["protein"] = diseaseMutations[gene_col].apply(lambda x:
                                                                            uniprotID[x] if x in uniprotID 
                                                                            else x + 'unknown')
    
        if 'mut_position' not in naturalMutations.columns:
            naturalMutations["mutation"] = naturalMutations[mut_RefSeq_col].apply(lambda x: x[x.find('p.')+2:])
            naturalMutations["mutation"] = naturalMutations["mutation"].apply(lambda x: (x[0], int(x[1:-1]), x[-1]))
            (naturalMutations["wt_res"],
             naturalMutations["mut_position"],
             naturalMutations["mut_res"]) = zip(* naturalMutations["mutation"].values)
    
        if 'mut_position' not in diseaseMutations.columns:
            diseaseMutations["mutation"] = diseaseMutations[mut_RefSeq_col].apply(lambda x: x[x.find('p.')+2:])
            diseaseMutations["mutation"] = diseaseMutations["mutation"].apply(lambda x: (x[0], int(x[1:-1]), x[-1]))
            (diseaseMutations["wt_res"],
             diseaseMutations["mut_position"],
             diseaseMutations["mut_res"]) = zip(* diseaseMutations["mutation"].values)
        
        naturalMutations.rename (columns = {'edgotype': 'edgotype_' + name,
                                            'partners': 'partners_' + name,
                                            'perturbations': 'perturbations_' + name}, inplace=True)
        diseaseMutations.rename (columns = {'edgotype': 'edgotype_' + name,
                                            'partners': 'partners_' + name,
                                            'perturbations': 'perturbations_' + name}, inplace=True)
        
        print()
        print('%s mutations:' % name)
        print('Non-disease: %d' % len(diseaseMutations))
        print('Disease: %d' % len(diseaseMutations))
        
        if name is interactome_names[0]:
            naturalMutations1, diseaseMutations1 = naturalMutations.copy(), diseaseMutations.copy()
        else:
            naturalMutations2, diseaseMutations2 = naturalMutations.copy(), diseaseMutations.copy()
    
    #----------------------------------------------------------------------------------------
    # calculate overlap in protein space
    #----------------------------------------------------------------------------------------
    
    natMut_overlap = pd.merge (naturalMutations1,
                               naturalMutations2,
                               how = 'inner',
                               on = ['protein', 'mut_position', 'wt_res', 'mut_res'])
    disMut_overlap = pd.merge (diseaseMutations1,
                               diseaseMutations2,
                               how = 'inner',
                               on = ['protein', 'mut_position', 'wt_res', 'mut_res'])
    
    print()
    print('Overlapping mutations:')
    print('Non-disease: %d' % len(natMut_overlap))
    print('Disease: %d' % len(disMut_overlap))
    
    cols = ["protein", 'mut_position', 'wt_res', 'mut_res']
    cols = cols + ['edgotype_' + n for n in interactome_names]
    cols = cols + ['partners_' + n for n in interactome_names]
    cols = cols + ['perturbations_' + n for n in interactome_names]
    
    natMut_overlap[cols].to_csv (natMutOutFile, sep='\t', index=False)
    disMut_overlap[cols].to_csv (disMutOutFile, sep='\t', index=False)    

if __name__ == '__main__':
    main()
