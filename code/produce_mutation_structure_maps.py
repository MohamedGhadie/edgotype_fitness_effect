#----------------------------------------------------------------------------------------
# Map mutations in the structural interactome onto structural models to be submitted for 
# protein folding ∆∆G calculations.
#
# Run the following scripts before running this script:
# - produce_data_mappings.py
# - produce_structural_interactome.py
# - calculate_dist_to_center.py
#----------------------------------------------------------------------------------------

import os
import pandas as pd
from pathlib import Path
from threeD_structure_tools import write_mutation_structure_maps

def main():
    
    # reference interactome name
    # options: HI-II-14, IntAct
    interactome_name = 'HI-II-14'
    
    # homology modelling method used to create structural models
    # options: template_based, model_based
    model_method = 'model_based'
    
    # method of calculating mutation ∆∆G for which results will be used
    # options: foldx
    ddg_method = 'foldx'
    
    # parent directory of all data files
    dataDir = Path('../data')
    
    # parent directory of all processed data files
    procDir = dataDir / 'processed'
    
    # directory of processed data files specific to interactome
    interactomeDir = procDir / interactome_name
    
    # directory of processed model-related data files specific to interactome
    modellingDir = interactomeDir / model_method
    
    # directory of calculation method
    methodDir = modellingDir / 'physics'
    
    # directory for PDB structure files
    pdbDir = Path('/Volumes/MG_Samsung/pdb_files')
    #pdbDir = Path('../pdb_files')
    
    if model_method is 'model_based':
        modelDir = Path('../models')
    else:
        modelDir = pdbDir
    
    # input data files
    chainSeqFile = modellingDir / 'protein_chain_sequences.pkl'
    chainStrucResFile = modellingDir / 'protein_chain_strucRes.pkl'
    naturalMutationsFile = methodDir / 'nondisease_mutation_dist_to_center.txt'
    diseaseMutationsFile = methodDir / 'disease_mutation_dist_to_center.txt'
    
    # output data files
    natural_mutations_ddg_file = modellingDir / ('nondis_mut_folding_ddg_%s.txt' % ddg_method)
    disease_mutations_ddg_file = modellingDir / ('dis_mut_folding_ddg.txt' % ddg_method)
    
    # create output directories if not existing
    if not methodDir.exists():
        os.makedirs(methodDir)
    
    #------------------------------------------------------------------------------------
    # write mutations mapped onto structural models to file
    #------------------------------------------------------------------------------------
    
    naturalMutations = pd.read_table (naturalMutationsFile, sep='\t')
    diseaseMutations = pd.read_table (diseaseMutationsFile, sep='\t')
        
    # write non-disease mutations
    if not natural_mutations_ddg_file.is_file():
        print('Writing non-disease mutations to file %s' % str(natural_mutations_ddg_file))
        naturalMutationMaps = list (zip(naturalMutations["protein"].tolist(),
                                        naturalMutations["mut_position"].tolist(),
                                        naturalMutations["mut_res"].tolist(),
                                        naturalMutations["model_chain"].tolist(),
                                        naturalMutations["chain_pos"].tolist()))
        write_mutation_structure_maps (naturalMutationMaps,
                                       chainSeqFile,
                                       chainStrucResFile,
                                       pdbDir,
                                       natural_mutations_ddg_file)
    
    # write disease mutations
    if not disease_mutations_ddg_file.is_file():
        print('Writing disease mutations to file %s' % str(disease_mutations_ddg_file))
        diseaseMutationMaps = list (zip(diseaseMutations["protein"].tolist(),
                                        diseaseMutations["mut_position"].tolist(),
                                        diseaseMutations["mut_res"].tolist(),
                                        diseaseMutations["model_chain"].tolist(),
                                        diseaseMutations["chain_pos"].tolist()))
        write_mutation_structure_maps (diseaseMutationMaps,
                                       chainSeqFile,
                                       chainStrucResFile,
                                       pdbDir,
                                       disease_mutations_ddg_file)
    
if __name__ == "__main__":
    main()
