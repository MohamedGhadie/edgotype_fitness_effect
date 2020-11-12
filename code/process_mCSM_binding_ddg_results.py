#----------------------------------------------------------------------------------------
# Process binding ∆∆G results from mCSM.
#----------------------------------------------------------------------------------------

import os
from pathlib import Path
from energy_tools import read_mCSM_results, write_mutation_ddg_tofile

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
    
    # results file
    inPath = modellingDir / 'mCSM' / 'results.txt'
    
    processed = read_mCSM_results (inPath)
    
    write_mutation_ddg_tofile (processed,
                               modellingDir / 'nondis_mut_binding_ddg_mCSM.txt',
                               modellingDir / 'nondis_mut_binding_ddg_mCSM_2.txt',
                               type = 'binding')
    write_mutation_ddg_tofile (processed,
                               modellingDir / 'dis_mut_binding_ddg_mCSM.txt',
                               modellingDir / 'dis_mut_binding_ddg_mCSM_2.txt',
                               type = 'binding')
    
    os.remove (modellingDir / 'nondis_mut_binding_ddg_mCSM.txt')
    os.remove (modellingDir / 'dis_mut_binding_ddg_mCSM.txt')
    os.rename (modellingDir / 'nondis_mut_binding_ddg_mCSM_2.txt',
               modellingDir / 'nondis_mut_binding_ddg_mCSM.txt')
    os.rename (modellingDir / 'dis_mut_binding_ddg_mCSM_2.txt',
               modellingDir / 'dis_mut_binding_ddg_mCSM.txt')

if __name__ == "__main__":
    main()
