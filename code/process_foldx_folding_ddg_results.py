import os
from pathlib import Path
from ddg_tools import (read_foldx_results,
                       write_mutation_ddg_tofile,
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
    
    # directory of foldx results
    inDir = modellingDir / 'foldx' / 'results'
    
    # directory of foldx output jobs
    outDir = modellingDir / 'foldx'
    
    # directory of PDB structure files
    pdbDir = Path('/Volumes/MG_Samsung/pdb_files')
    
    if model_method is 'model_based':
        modelDir = modellingDir / 'protein_models'
    else:
        modelDir = pdbDir
    
    # create output directories if not existing
    if not outDir.exists():
        os.makedirs(outDir)
    
    processed, unprocessed = read_foldx_results (inDir, type = 'folding')
    
    write_mutation_ddg_tofile (processed,
                               modellingDir / 'nondis_mut_folding_ddg_foldx.txt',
                               modellingDir / 'nondis_mut_folding_ddg_foldx_2.txt',
                               type = 'folding')
    write_mutation_ddg_tofile (processed,
                               modellingDir / 'dis_mut_folding_ddg_foldx.txt',
                               modellingDir / 'dis_mut_folding_ddg_foldx_2.txt',
                               type = 'folding')
    
    os.remove (modellingDir / 'nondis_mut_folding_ddg_foldx.txt')
    os.remove (modellingDir / 'dis_mut_folding_ddg_foldx.txt')
    os.rename (modellingDir / 'nondis_mut_folding_ddg_foldx_2.txt',
               modellingDir / 'nondis_mut_folding_ddg_foldx.txt')
    os.rename (modellingDir / 'dis_mut_folding_ddg_foldx_2.txt',
               modellingDir / 'dis_mut_folding_ddg_foldx.txt')
    
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
                                   mem_per_cpu = '4000M',
                                   username = 'ghadie84',
                                   outputfile = '/project/ctb-yxia/ghadie84/foldx/data/%x-%j.out',
                                   serverDataDir = '/project/ctb-yxia/ghadie84/foldx/data')

if __name__ == "__main__":
    main()
