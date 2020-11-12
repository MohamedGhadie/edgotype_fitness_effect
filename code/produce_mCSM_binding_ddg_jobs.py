#----------------------------------------------------------------------------------------
# Produce jobs for mutations mapped onto PPI structures to be submitted 
# to mCSM_PPI2 for binding ∆∆G calculations. Produced jobs may contain multiple mutations. 
# Mutations with existing ∆∆G values in the input file are skipped.
#----------------------------------------------------------------------------------------

import os
from pathlib import Path
from energy_tools import (append_mutation_ddg_files,
                          read_unprocessed_ddg_mutations,
                          produce_mCSM_jobs)

def main():
    
    # reference interactome name: HI-II-14, HuRI, IntAct
    interactome_name = 'IntAct'
    
    # homology modelling method used to create structural models
    # options: template_based, model_based
    model_method = 'model_based'
    
    # parent directory of all data files
    dataDir = Path('../data')
    
    # parent directory of all processed data files
    procDir = dataDir / 'processed'
    
    # directory of processed data files specific to interactome
    interactomeDir = procDir / interactome_name
    
    # directory of processed model-related data files specific to interactome
    modellingDir = interactomeDir / model_method
    
    # directory of foldx output jobs
    outDir = modellingDir / 'mCSM'
    
    # directory of PDB structure files
    pdbDir = Path('../../pdb_files')
    
    if model_method is 'model_based':
        modelDir = modellingDir / 'ppi_models'
    else:
        modelDir = pdbDir
    
    # input file containing mutations to submit to foldx
    nondiseaseMutFile = modellingDir / 'nondis_mut_binding_ddg_mCSM.txt'
    diseaseMutFile = modellingDir / 'dis_mut_binding_ddg_mCSM.txt'
    
    # temporary files
    allMutFile = modellingDir / 'all_mut_binding_ddg_mCSM.txt'
    
    # create output directories if not existing
    if not outDir.exists():
        os.makedirs(outDir)
    
    append_mutation_ddg_files (nondiseaseMutFile, diseaseMutFile, allMutFile)
    mutations = read_unprocessed_ddg_mutations (allMutFile, 'binding')
    produce_mCSM_jobs (mutations, modelDir, outDir)
    os.remove(allMutFile)

if __name__ == "__main__":
    main()
