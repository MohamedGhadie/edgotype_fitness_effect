#----------------------------------------------------------------------------------------
# Process results of folding ∆∆G calculations by DynaMut2.
#----------------------------------------------------------------------------------------

import os
from pathlib import Path
from energy_tools import read_dynamut_results, write_mutation_ddg_tofile

def main():
    
    # reference interactome name: HI-II-14, HuRI, IntAct
    interactome_name = 'IntAct'
    
    # homology modelling method used to create structural models
    # options: template_based, model_based
    model_method = 'model_based'
    
    # method used to perform calculations; 'geometry' or 'physics'
    edgetic_method = 'physics'
    
    # method that was used to calculate edgetic mutation binding ∆∆G
    # options: bindprofx, foldx
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
    
    # directory of DynaMut2 results
    inDir = edgeticDir / 'dynamut'
    
    processed = read_dynamut_results (inDir)
    
    write_mutation_ddg_tofile (processed,
                               edgeticDir / 'nondis_mut_folding_ddg_dynamut.txt',
                               edgeticDir / 'nondis_mut_folding_ddg_dynamut_2.txt',
                               type = 'folding')
    write_mutation_ddg_tofile (processed,
                               edgeticDir / 'dis_mut_folding_ddg_dynamut.txt',
                               edgeticDir / 'dis_mut_folding_ddg_dynamut_2.txt',
                               type = 'folding')
    
    os.remove (edgeticDir / 'nondis_mut_folding_ddg_dynamut.txt')
    os.remove (edgeticDir / 'dis_mut_folding_ddg_dynamut.txt')
    os.rename (edgeticDir / 'nondis_mut_folding_ddg_dynamut_2.txt',
               edgeticDir / 'nondis_mut_folding_ddg_dynamut.txt')
    os.rename (edgeticDir / 'dis_mut_folding_ddg_dynamut_2.txt',
               edgeticDir / 'dis_mut_folding_ddg_dynamut.txt')

if __name__ == "__main__":
    main()
