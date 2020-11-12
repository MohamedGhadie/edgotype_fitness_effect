#----------------------------------------------------------------------------------------
# Process results from folding ∆∆G calculations by FoldX. Failed jobs with multiple 
# mutations will be split again into multiple single-mutation jobs.
#----------------------------------------------------------------------------------------

import os
from pathlib import Path
from energy_tools import (read_foldx_results,
                          write_mutation_ddg_tofile,
                          produce_foldx_and_beluga_jobs)

def main():
    
    # reference interactome name: HI-II-14, HuRI, IntAct
    interactome_name = 'HuRI'
    
    # homology modelling method used to create structural models
    # options: template_based, model_based
    model_method = 'model_based'
    
    # method used to perform calculations; 'geometry' or 'physics'
    edgetic_method = 'physics'
    
    # method that was used to calculate edgetic mutation binding ∆∆G
    # options: bindprofx, foldx
    edgetic_ddg = 'foldx'
    
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
    
    # directory of foldx results
    inDir = edgeticDir / 'foldx' / 'results'
    
    # directory of foldx output jobs
    outDir = edgeticDir / 'foldx'
    
    # directory of PDB structure files
    pdbDir = Path('../../pdb_files')
    
    if model_method is 'model_based':
        modelDir = modellingDir / 'protein_models'
    else:
        modelDir = pdbDir
    
    # create output directories if not existing
    if not outDir.exists():
        os.makedirs(outDir)
    
    processed, unprocessed = read_foldx_results (inDir, type = 'folding')
    
    write_mutation_ddg_tofile (processed,
                               edgeticDir / 'nondis_mut_folding_ddg_foldx.txt',
                               edgeticDir / 'nondis_mut_folding_ddg_foldx_2.txt',
                               type = 'folding')
    write_mutation_ddg_tofile (processed,
                               edgeticDir / 'dis_mut_folding_ddg_foldx.txt',
                               edgeticDir / 'dis_mut_folding_ddg_foldx_2.txt',
                               type = 'folding')
    
    os.remove (edgeticDir / 'nondis_mut_folding_ddg_foldx.txt')
    os.remove (edgeticDir / 'dis_mut_folding_ddg_foldx.txt')
    os.rename (edgeticDir / 'nondis_mut_folding_ddg_foldx_2.txt',
               edgeticDir / 'nondis_mut_folding_ddg_foldx.txt')
    os.rename (edgeticDir / 'dis_mut_folding_ddg_foldx_2.txt',
               edgeticDir / 'dis_mut_folding_ddg_foldx.txt')
    
    produce_foldx_and_beluga_jobs (unprocessed,
                                   modelDir,
                                   outDir,
                                   'folding',
                                   account = 'ctb-yxia',
                                   walltime = '1-00',
                                   ntasks = 1,
                                   nodes = 1,
                                   ntasks_per_node = 1,
                                   cpus_per_task = 1,
                                   mem_per_cpu = '4G',
                                   username = 'ghadie84',
                                   outputfile = '/project/ctb-yxia/ghadie84/foldx/data/%x-%j.out',
                                   serverDataDir = '/project/ctb-yxia/ghadie84/foldx/data')

if __name__ == "__main__":
    main()
