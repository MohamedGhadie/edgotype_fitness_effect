#----------------------------------------------------------------------------------------
# Produce mutation mappings onto structural models of single proteins to be submitted to 
# FoldX for protein folding ∆∆G calculations.
#----------------------------------------------------------------------------------------

import os
import pandas as pd
from pathlib import Path
from threeD_structure_tools import write_mutation_structure_maps

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
    
    # directory for PDB structure files
    pdbDir = Path('../../pdb_files')
    
    if model_method is 'model_based':
        modelDir = modellingDir / 'protein_models'
    else:
        modelDir = pdbDir
    
    # input data files
    chainSeqFile = modellingDir / 'protein_chain_sequences.pkl'
    chainStrucResFile = modellingDir / 'protein_chain_strucRes.pkl'
    naturalMutationsFile = edgeticDir / 'nondisease_mutation_struc_loc.txt'
    diseaseMutationsFile = edgeticDir / 'disease_mutation_struc_loc.txt'
    
    # output data files
    natural_mutations_ddg_file = edgeticDir / 'nondis_mut_folding_ddg_dynamut.txt'
    disease_mutations_ddg_file = edgeticDir / 'dis_mut_folding_ddg_dynamut.txt'
    
    # create output directories if not existing
    if not edgeticDir.exists():
        os.makedirs(edgeticDir)
    
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
                                       modelDir,
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
                                       modelDir,
                                       disease_mutations_ddg_file)
    
if __name__ == "__main__":
    main()
