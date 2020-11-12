#----------------------------------------------------------------------------------------
# Predict edgetic mutations based on physics.
#----------------------------------------------------------------------------------------

import os
import numpy as np
import pandas as pd
from pathlib import Path

def main():
    
    # reference interactome name: HI-II-14, HuRI, IntAct
    interactome_name = 'IntAct'
    
    # sample size factor
    factor = 20
    
    # homology modelling method used to create structural models
    # options: template_based, model_based
    model_method = 'model_based'
    
    # method of calculating mutation ∆∆G for which results will be used
    # options: bindprofx, foldx
    ddg_method = 'foldx'
    
    # parent directory of all data files
    dataDir = Path('../data')
    
    # parent directory of all processed data files
    procDir = dataDir / 'processed'
    
    # directory of processed data files specific to interactome
    interactomeDir = procDir / interactome_name
    
    # directory of processed model-related data files specific to interactome
    modellingDir = interactomeDir / model_method
    
    # directory of edgetic mutation calculation method
    edgeticDir = modellingDir / 'physics' / (ddg_method + '_edgetics')
    
    # input data files
    naturalMutationsFile = edgeticDir / 'nondisease_mutation_struc_loc.txt'
    diseaseMutationsFile = edgeticDir / 'disease_mutation_struc_loc.txt'
    
    # output data files
    natMutOutFile = modellingDir / 'geometry' / 'nondisease_mutation_struc_loc_sample.txt'
    disMutOutFile = modellingDir / 'geometry' / 'disease_mutation_struc_loc_sample.txt'
    
    # create output directories if not existing
    if not modellingDir.exists():
        os.makedirs(modellingDir)
    
    #------------------------------------------------------------------------------------
    # Read mutation geometry-based edgetic perturbations and mutation ∆∆G
    #------------------------------------------------------------------------------------
    
    naturalMutations = pd.read_table (naturalMutationsFile, sep='\t')
    diseaseMutations = pd.read_table (diseaseMutationsFile, sep='\t')
    
    #------------------------------------------------------------------------------------
    # predict PPI perturbations based on physics
    #------------------------------------------------------------------------------------
       
    natMutSample = pd.DataFrame()
    for r in ['interface', 'buried', 'exposed-noninterface']:
        subset = naturalMutations [naturalMutations["structural_location"] == r].reset_index(drop=True)
        numMut = len(subset)
        s = int(np.round(numMut/factor))
        ind = np.random.choice(numMut, size = s, replace=False)
        natMutSample = natMutSample.append(subset.iloc[ind,:])
        print('Non-disease %s mutations selected: %d' % (r, s))
    
    print()
    disMutSample = pd.DataFrame()
    for r in ['interface', 'buried', 'exposed-noninterface']:
        subset = diseaseMutations [diseaseMutations["structural_location"] == r].reset_index(drop=True)
        numMut = len(subset)
        s = int(np.round(numMut/factor))
        ind = np.random.choice(numMut, size = s, replace=False)
        disMutSample = disMutSample.append(subset.iloc[ind,:])
        print('Disease %s mutations selected: %d' % (r, s))
    
    dropCol = [c for c in ["edgotype", "partners", "perturbations"] if c in natMutSample.columns]
    natMutSample.drop (dropCol, axis=1, inplace=True)
    dropCol = [c for c in ["edgotype", "partners", "perturbations"] if c in disMutSample.columns]
    disMutSample.drop (dropCol, axis=1, inplace=True)
    
    print()
    print('Overall mutations selected')
    print('Non-disease: %d' % len(natMutSample))
    print('Disease: %d' % len(disMutSample))
    
    natMutSample.to_csv (natMutOutFile, index=False, sep='\t')
    disMutSample.to_csv (disMutOutFile, index=False, sep='\t')

if __name__ == "__main__":
    main()
