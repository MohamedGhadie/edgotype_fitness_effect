import os
import pandas as pd
from pathlib import Path
from id_mapping import produce_chainSeq_dict
from modelling_tools import (set_pdb_dir,
                             set_template_dir,
                             enable_pdb_downloads,
                             disable_pdb_warnings,
                             write_template_sequences,
                             produce_template_files)

def main():
    
    # reference interactome name
    # options: HI-II-14, IntAct
    interactome_name = 'HI-II-14'
    
    # homology modelling method used to create structural models
    # options: template_based, model_based
    #model_method = 'model_based'
    
    # parent directory of all data files
    dataDir = Path('../data')
    
    # directory of data files from external sources
    extDir = dataDir / 'external'
    
    # parent directory of all processed data files
    procDir = dataDir / 'processed'
    
    # directory of processed data files specific to interactome
    interactomeDir = procDir / interactome_name
    
    # directory of processed template-related data files specific to interactome
    templateBasedDir = interactomeDir / 'template_based'
    
    # directory of processed model-related data files specific to interactome
    modelBasedDir = interactomeDir / 'model_based'
    
    # directory for PDB structure files
    pdbDir = Path('/Volumes/MG_Samsung/pdb_files')
    
    # directory for template structure files
    templateDir = Path('../templates')
    
    # input data files
    templateMapFile = templateBasedDir / 'single_chain_map_per_protein.txt'
    chainSeqresFile = templateBasedDir / 'protein_chain_sequences.pkl'
    chainStrucResFile = templateBasedDir / 'protein_chain_strucRes.pkl'
    
    # output data files
    chainStrucSeqFastaFile = modelBasedDir / 'protein_template_sequences.fasta'
    chainStrucSeqFile = modelBasedDir / 'protein_template_sequences.pkl'
    
    # create output directories if not existing
    if not modellingDir.exists():
        os.makedirs(modellingDir)
    if not pdbDir.exists():
        os.makedirs(pdbDir)
    if not templateDir.exists():
        os.makedirs(templateDir)
    
    # set directory of raw PDB coordinate files for modelling tools
    set_pdb_dir (pdbDir)
    
    # set directory of template coordinate files for modelling tools
    set_template_dir (templateDir)
    
    # enable or disable PDB downloads
    enable_pdb_downloads (allow_pdb_downloads)
    
    # suppress or allow PDB warnings
    disable_pdb_warnings (suppress_pdb_warnings)
    
    templateMap = pd.read_table(templateMapFile, sep='\t')
    templateIDs = templateMap["Subject"].tolist()
    
    if not chainStrucSeqFastaFile.is_file():
        print('writing protein template sequences to Fasta file')
        write_template_sequences (templateIDs,
                                  chainSeqresFile,
                                  chainStrucResFile,
                                  chainStrucSeqFastaFile)
    
    if not chainStrucSeqFile.is_file():
        print('producing protein template sequence dictionary')
        produce_chainSeq_dict (chainStrucSeqFastaFile, chainStrucSeqFile)
    
    print('Extracting coordinate files for protein templates')
    produce_template_files (templateIDs, chainSeqresFile, chainStrucResFile)

if __name__ == "__main__":
    main()
