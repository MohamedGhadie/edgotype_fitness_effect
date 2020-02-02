import os
from pathlib import Path
from modelling_tools import (set_model_dir,
                             produce_protein_model_chainSeq_dict,
                             produce_model_chain_strucRes_dict,
                             produce_protein_chain_pos_mapping)

def main():
    
    # reference interactome name
    # options: HI-II-14, IntAct
    interactome_name = 'HI-II-14'
    
    # parent directory of all data files
    dataDir = Path('../data')
    
    # parent directory of all processed data files
    procDir = dataDir / 'processed'
    
    # directory of processed data files specific to interactome
    interactomeDir = procDir / interactome_name
    
    # directory of processed model-related data files specific to interactome
    modelBasedDir = interactomeDir / 'model_based'
    
    # input data files
    templateMapFile = modelBasedDir / 'single_template_map_per_protein.txt'
    proteinSeqFile = procDir / 'human_reference_sequences.pkl'
    
    # output data files
    chainSeqFile = modelBasedDir / 'protein_chain_sequences.pkl'
    chainStrucResFile = modelBasedDir / 'protein_chain_strucRes.pkl'
    chainMapFile = modelBasedDir / 'single_chain_map_per_protein.txt'
    
    if not modelBasedDir.exists():
        os.makedirs(modelBasedDir)
    
    if not chainSeqFile.is_file():
        print('producing protein model chain sequence dictionary')
        produce_protein_model_chainSeq_dict (templateMapFile, proteinSeqFile, chainSeqFile)
    
    if not chainStrucResFile.is_file():
        print('producing protein model chain structured residue label file')
        produce_model_chain_strucRes_dict (chainSeqFile, chainStrucResFile)
    
    if not chainMapFile.is_file():
        print('producing protein model chain position mapping file')
        produce_protein_chain_pos_mapping (templateMapFile, chainSeqFile, chainMapFile)

if __name__ == "__main__":
    main()
