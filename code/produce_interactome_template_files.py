import os
from pathlib import Path
from id_mapping import produce_chainSeq_dict
from modelling_tools import (set_pdb_dir,
                             set_template_dir,
                             enable_pdb_downloads,
                             disable_pdb_warnings,
                             write_interactome_template_sequences,
                             produce_interactome_template_files)

def main():
    
    # reference interactome name
    # options: HI-II-14, IntAct
    interactome_name = 'HI-II-14'
    
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
    interactomeFile = templateBasedDir / 'human_structural_interactome.txt'
    chainSeqresFile = procDir / 'chain_seqres.pkl'
    chainStrucResFile = procDir / 'chain_strucRes.pkl'
    
    # output data files
    chainStrucSeqFastaFile = modelBasedDir / 'chain_struc_sequences.fasta'
    chainStrucSeqFile = modelBasedDir / 'chain_struc_sequences.pkl'
    
    # create output directories if not existing
    if not modelBasedDir.exists():
        os.makedirs(mdoelBasedDir)
    if not templateDir.exists():
        os.makedirs(templateDir)
    if not pdbDir.exists():
        os.makedirs(pdbDir)
    
    # set directory of raw PDB coordinate files for modelling tools
    set_pdb_dir (pdbDir)
    
    # set directory of template coordinate files for modelling tools
    set_template_dir (templateDir)
    
    # enable or disable PDB downloads
    enable_pdb_downloads (allow_pdb_downloads)
    
    # suppress or allow PDB warnings
    disable_pdb_warnings (suppress_pdb_warnings)
    
    if not chainStrucSeqFastaFile.is_file():
        print('writing interactome template sequences to Fasta file')
        write_interactome_template_sequences (interactomeFile,
                                              chainSeqresFile,
                                              chainStrucResFile,
                                              chainStrucSeqFastaFile)
    
    if not chainStrucSeqFile.is_file():
        print('producing template sequence dictionary')
        produce_chainSeq_dict (chainStrucSeqFastaFile, chainStrucSeqFile)
    
    print('Extracting coordinate files for PPI templates')
    produce_interactome_template_files (interactomeFile, chainSeqresFile, chainStrucResFile)

if __name__ == "__main__":
    main()
