#----------------------------------------------------------------------------------------
# This script constructs a structural interactome from a reference interactome by mapping 
# interaction binding interfaces at amino acid resolution from experimentally determined 
# three-dimensional structural models in PDB onto interactions in the reference interactome.
#
# Run the following scripts before running this script:
# - produce_data_mappings.py
# - process_interactome.py
# - BLAST human protein sequences against PDB sequences and save results into path
#   ../data/external/human_pdb_e-5
#----------------------------------------------------------------------------------------

import os
from pathlib import Path
from interactome_tools import read_single_interface_annotated_interactome
from modelling_tools import (set_model_dir,
                             produce_model_annotated_interactome,
                             produce_ppi_model_chainSeq_dict,
                             produce_ppi_chain_pos_mapping,
                             produce_model_chain_strucRes_dict)
from structural_annotation import (produce_interface_annotated_interactome,
                                   merge_interactome_interface_annotations,
                                   remove_duplicate_interface_annotations)

def main():
    
    # reference interactome name
    # options: HI-II-14, IntAct
    interactome_name = 'HI-II-14'
    
    # max binding distance for interface residues in PDB structure
    bindingDist = 5
    
    # show figures
    showFigs = False
    
    # parent directory of all data files
    dataDir = Path('../data')
    
    # parent directory of all processed data files
    procDir = dataDir / 'processed'
    
    # directory of processed data files specific to interactome
    interactomeDir = procDir / interactome_name
    
    # directory of processed model-related data files specific to interactome
    modellingDir = interactomeDir / 'model_based'
    
    # figure directory
    figDir = Path('../figures') / interactome_name / 'model_based'
    
    # directory for model structure files
    modelDir = Path('../models')
    
    # input data files
    templateAnnotatedInteractomeFile = modellingDir / 'human_template_annotated_interactome.txt'
    #pdbBlastFile = extDir / 'human_pdb_e-5'
    proteinSeqFile = procDir / 'human_reference_sequences.pkl'
    #chainSeqFile = procDir / 'chain_sequences.pkl'
    #pdbChainsFile = procDir / 'pdb_seqres_chains.pkl'
    #chainListFile = procDir / 'pdb_seqres_chains.list'
    #chainStrucResFile = procDir / 'chain_strucRes.pkl'
    #interactomeFile = interactomeDir / 'human_interactome.txt'
    
    # output data files
    modelAnnotatedInteractomeFile = modellingDir / 'human_model_annotated_interactome.txt'
    chainSeqFile = modellingDir / 'interactome_model_chain_sequences.pkl'
    chainMapFile = modellingDir / 'interactome_model_chain_map.txt'
    modelInterfaceFile = procDir / 'model_interfaces.txt'
    chainStrucResFile = modellingDir / 'interactome_model_chain_strucRes.pkl'

#     proteinChainsFile = procDir / 'protein_chains.pkl'
#     alignmentEvalueFile = procDir / 'human_protein_chain_min_alignment_evalues.pkl'
#     chainInterfaceFile = procDir / 'pdb_interfaces.txt'
#     chainAnnotatedInteractomeFile = interactomeDir / 'human_chain_annotated_interactome.txt'
#     chainIDFile = interactomeDir / 'interactome_chainIDs.txt'
#     pdbIDFile = interactomeDir / 'interactome_pdbIDs.txt'
    structuralInteractomeFile1 = modellingDir / 'human_structural_interactome_withDuplicates.txt'
    structuralInteractomeFile = modellingDir / 'human_structural_interactome.txt'
#     refInteractomeChainMapFile = interactomeDir / 'ref_interactome_pdb_chain_map.txt'
#     strucInteractomeChainMapFile = interactomeDir / 'struc_interactome_pdb_chain_map.txt'
    
    if not modellingDir.exists():
        os.makedirs(modellingDir)
    if not figDir.exists():
        os.makedirs(figDir)
    
    set_model_dir (modelDir)
    
    if not modelAnnotatedInteractomeFile.is_file():
        print('producing model annotated interactome')
        produce_model_annotated_interactome (templateAnnotatedInteractomeFile,
                                             modelAnnotatedInteractomeFile)
    
    if not chainSeqFile.is_file():
        print('producing model chain sequence dictionary')
        produce_ppi_model_chainSeq_dict (modelAnnotatedInteractomeFile,
                                         proteinSeqFile,
                                         chainSeqFile)
    
    if not chainMapFile.is_file():
        print('producing model chain position mapping file')
        produce_ppi_chain_pos_mapping (modelAnnotatedInteractomeFile,
                                       chainSeqFile,
                                       chainMapFile)
    
    if not chainStrucResFile.is_file():
        print('producing model chain structured residue label file')
        produce_model_chain_strucRes_dict (chainSeqFile, chainStrucResFile)
    
    if not structuralInteractomeFile1.is_file():
        print('mapping model interfaces onto model-annotated interactome')
        produce_interface_annotated_interactome(modelAnnotatedInteractomeFile,
                                                modelDir,
                                                chainSeqFile,
                                                chainMapFile,
                                                modelInterfaceFile,
                                                chainStrucResFile,
                                                1,
                                                1,
                                                False,
                                                0,
                                                bindingDist,
                                                structuralInteractomeFile1,
                                                downloadPDB = False,
                                                suppressWarnings = False)
        
    print('merging interface annotations for each PPI')
    merge_interactome_interface_annotations (structuralInteractomeFile1,
                                             structuralInteractomeFile)
    
    structuralInteractome = read_single_interface_annotated_interactome( structuralInteractomeFile )
    interactomeProteins = list(set(structuralInteractome[["Protein_1", "Protein_2"]].values.flatten()))
    print( '\n' + 'Structural interactome:' )
    print( '%d PPIs' % len(structuralInteractome) )
    print( '%d proteins' % len(interactomeProteins) )
    print()

if __name__ == "__main__":
    main()
