import os
from pathlib import Path
from modelling_tools import (set_model_dir,
                             produce_protein_fullmodel_chainSeq_dict,
                             produce_fullmodel_chain_strucRes_dict,
                             produce_protein_fullmodel_pos_mapping)

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
    proteinSeqFile = procDir / 'human_reference_sequences.pkl'
    templateMapFile = modelBasedDir / 'single_template_map_per_protein.txt'
    
    # output data files
    chainSeqFile = modelBasedDir / 'protein_chain_sequences.pkl'
    chainStrucResFile = modelBasedDir / 'protein_chain_strucRes.pkl'
    chainMapFile = modelBasedDir / 'single_chain_map_per_protein.txt'
    
    if not modelBasedDir.exists():
        os.makedirs(modelBasedDir)
    
    print('producing protein model chain sequence dictionary')
    produce_protein_fullmodel_chainSeq_dict (templateMapFile, proteinSeqFile, chainSeqFile)
    
    print('producing protein model chain structured residue label file')
    produce_fullmodel_chain_strucRes_dict (chainSeqFile, chainStrucResFile)
    
    print('producing protein model chain position mapping file')
    produce_protein_fullmodel_pos_mapping (templateMapFile, chainSeqFile, chainMapFile)

if __name__ == "__main__":
    main()
