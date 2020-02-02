import os
from pathlib import Path
from ddg_tools import (append_mutation_ddg_files,
                       read_unprocessed_ddg_mutations,
                       produce_foldx_and_beluga_jobs)

def main():
    
    # reference interactome name
    # options: HI-II-14, IntAct
    interactome_name = 'HI-II-14'
    
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
    pdbDir = Path('/Volumes/MG_Samsung/pdb_files')
    
    if model_method is 'model_based':
        modelDir = Path('../models')
    else
        modelDir = pdbDir
    
    # input file containing mutations to submit to bindprofx
    nondiseaseMutFile = modellingDir / 'nondis_mut_folding_ddg_foldx.txt'
    diseaseMutFile = modellingDir / 'dis_mut_folding_ddg_foldx.txt'
    
    # temporary files
    allMutFile = modellingDir / 'all_mut_folding_ddg_foldx.txt'
    
    # create output directories if not existing
    if not outDir.exists():
        os.makedirs(outDir)
    
    append_mutation_ddg_files (nondiseaseMutFile, diseaseMutFile, allMutFile)
    mutations = read_unprocessed_ddg_mutations (allMutFile, type = 'folding')
    
    produce_foldx_and_beluga_jobs (mutations,
                                   modelDir,
                                   outDir,
                                   'folding',
                                   account = 'ctb-yxia',
                                   walltime = '1-00',
                                   ntasks = 1,
                                   nodes = 1,
                                   ntasks_per_node = 1,
                                   cpus_per_task = 1,
                                   mem_per_cpu = '8000M',
                                   username = 'ghadie84',
                                   outputfile = '/project/def-yxia/ghadie84/foldx/data/%x-%j.out',
                                   serverDataDir = '/project/def-yxia/ghadie84/foldx/data')

if __name__ == "__main__":
    main()
