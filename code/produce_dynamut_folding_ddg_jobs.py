#----------------------------------------------------------------------------------------
# Produce jobs for mutations mapped onto protein structures to be submitted 
# to DynaMut2 for folding ∆∆G calculations. Jobs may contain multiple mutations. 
# Mutations with existing ∆∆G values in the input file are skipped.
#----------------------------------------------------------------------------------------

import os
from pathlib import Path
from energy_tools import (append_mutation_ddg_files,
                          read_unprocessed_ddg_mutations,
                          produce_dynamut_jobs)

def main():
    
    # reference interactome name: HI-II-14, HuRI, IntAct
    interactome_name = 'IntAct'
    
    # homology modelling method used to create structural models
    # options: template_based, model_based
    model_method = 'model_based'
    
    # method used to perform calculations; 'geometry' or 'physics'
    edgetic_method = 'physics'
    
    # method that was used to calculate edgetic mutation binding ∆∆G
    # options: bindprofx, foldx, mCSM
    edgetic_ddg = 'mCSM'
    
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
    
    if edgetic_method is 'physics':
        edgeticDir = edgeticDir / (edgetic_ddg + '_edgetics')
    
    # directory of dynamut output jobs
    outDir = edgeticDir / 'dynamut'
    
    # directory of PDB structure files
    pdbDir = Path('../../pdb_files')
    
    if model_method is 'model_based':
        modelDir = modellingDir / 'protein_models'
    else:
        modelDir = pdbDir
    
    # input file containing mutations to submit to dynamut
    nondiseaseMutFile = edgeticDir / 'nondis_mut_folding_ddg_dynamut.txt'
    diseaseMutFile = edgeticDir / 'dis_mut_folding_ddg_dynamut.txt'
    
    # temporary files
    allMutFile = edgeticDir / 'all_mut_folding_ddg_dynamut.txt'
    
    # create output directories if not existing
    if not outDir.exists():
        os.makedirs(outDir)
    
    append_mutation_ddg_files (nondiseaseMutFile, diseaseMutFile, allMutFile)
    mutations = read_unprocessed_ddg_mutations (allMutFile, type = 'folding')
    produce_dynamut_jobs (mutations, modelDir, outDir)
    os.remove(allMutFile)

if __name__ == "__main__":
    main()
