#----------------------------------------------------------------------------------------
# Produce jobs for mutations mapped onto PPI structures to be submitted 
# to FoldX for binding ∆∆G calculations. Each produced job contains only one mutation. 
# Mutations with existing ∆∆G values in the input file are skipped.
#----------------------------------------------------------------------------------------

import os
from pathlib import Path
from energy_tools import (append_mutation_ddg_files,
                          read_unprocessed_ddg_mutations,
                          produce_foldx_and_beluga_jobs)

def main():
    
    # reference interactome name: HI-II-14, HuRI, IntAct
    interactome_name = 'HuRI'
    
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
    outDir = modellingDir / 'foldx'
    
    # directory of PDB structure files
    pdbDir = Path('../../pdb_files')
    
    if model_method is 'model_based':
        modelDir = modellingDir / 'ppi_models'
    else:
        modelDir = pdbDir
    
    # input file containing mutations to submit to foldx
    nondiseaseMutFile = modellingDir / 'nondis_mut_binding_ddg_foldx.txt'
    diseaseMutFile = modellingDir / 'dis_mut_binding_ddg_foldx.txt'
    
    # temporary files
    allMutFile = modellingDir / 'all_mut_binding_ddg_foldx.txt'
    
    # create output directories if not existing
    if not outDir.exists():
        os.makedirs(outDir)
    
    append_mutation_ddg_files (nondiseaseMutFile, diseaseMutFile, allMutFile)
    mutations = read_unprocessed_ddg_mutations (allMutFile, 'binding')
    
    produce_foldx_and_beluga_jobs (mutations,
                                   modelDir,
                                   outDir,
                                   'binding',
                                   account = 'ctb-yxia',
                                   walltime = '1-00',
                                   ntasks = 1,
                                   nodes = 1,
                                   ntasks_per_node = 1,
                                   cpus_per_task = 1,
                                   mem_per_cpu = '4G',
                                   outputfile = '/project/ctb-yxia/ghadie84/foldx/data/%x-%j.out',
                                   username = 'ghadie84',
                                   serverDataDir = '/project/ctb-yxia/ghadie84/foldx/data')
    os.remove(allMutFile)

if __name__ == "__main__":
    main()
