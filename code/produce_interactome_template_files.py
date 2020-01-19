import os
from pathlib import Path
from text_tools import parse_blast_file
from structural_annotation import filter_chain_annotations
from modelling_tools import (set_pdb_dir,
                             set_template_dir,
                             set_alignment_dir,
                             enable_pdb_downloads,
                             disable_pdb_warnings,
                             extend_alignments,
                             produce_interactome_template_files,
                             produce_interactome_alignment_files)

def main():
    
    # reference interactome name
    # options: HI-II-14, IntAct
    interactome_name = 'HI-II-14'
    
    # Maximum e-value cutoff to filter out protein-chain annotations
    evalue = 1e-10
    
    # Minimum protein coverage fraction required for protein-chain annotation
    proteinCov = 0
    
    # Minimum chain coverage fraction required for protein-chain annotation
    chainCov = 0
    
    # allow downloading of PDB structures while constructing the structural interactome
    allow_pdb_downloads = True
    
    # suppress PDB warnings when constructing the structural interactome
    suppress_pdb_warnings = True
    
    # verbosity for Modeller
    verbose = True
    
    # parent directory of all data files
    dataDir = Path('../data')
    
    # parent directory of all processed data files
    procDir = dataDir / 'processed'
    
    # directory of processed data files specific to interactome
    interactomeDir = procDir / interactome_name
        
    # directory for PDB structure files
    pdbDir = Path('/Volumes/MG_Samsung/pdb_files')
    
    # directory for template structure files
    templateDir = Path('../templates')
    
    # directory for alignment files
    alignmentDir = Path('../alignments')
    
    # input data files
    blastFile = interactomeDir / 'interactome_template_blast_alignments_e-5'
    proteinSeqFile = procDir / 'human_reference_sequences.pkl'
    chainStrucSeqFile = interactomeDir / 'chain_struc_sequences.pkl'
    chainSeqresFile = procDir / 'chain_seqres.pkl'
    chainStrucResFile = procDir / 'chain_strucRes.pkl'
    interactomeFile = interactomeDir / 'human_template_annotated_interactome.txt'
    
    # output data files
    alignmentFile1 = interactomeDir / 'interactome_template_alignments.txt'
    alignmentFile2 = interactomeDir / 'interactome_template_filtered_alignments.txt'
    alignmentFile3 = interactomeDir / 'interactome_template_extended_alignments.txt'
    annotatedInteractomeFile = interactomeDir / 'human_model_annotated_interactome.txt'    
    
    # create output directories if not existing
    if not dataDir.exists():
        os.makedirs(dataDir)
    if not procDir.exists():
        os.makedirs(procDir)
    if not interactomeDir.exists():
        os.makedirs(interactomeDir)
    
    # set directory of raw PDB coordinate files for modelling tools
    set_pdb_dir (pdbDir)
    
    # set directory of template coordinate files for modelling tools
    set_template_dir (templateDir)
    
    # set directory of alignment files for modelling tools
    set_alignment_dir (alignmentDir)
    
    # enable or disable PDB downloads
    enable_pdb_downloads (allow_pdb_downloads)
    
    # suppress or allow PDB warnings
    disable_pdb_warnings (suppress_pdb_warnings)
    
    if not alignmentFile1.is_file():
        print( 'Parsing BLAST protein-chain alignment file' )
        parse_blast_file (blastFile, alignmentFile1)
    
    if not alignmentFile2.is_file():
        print('Filtering alignments')
        filter_chain_annotations (alignmentFile1,
                                  alignmentFile2,
                                  evalue = evalue,
                                  prCov = proteinCov,
                                  chCov = chainCov)
    
    if not alignmentFile3.is_file():
        print('Extending alignments to full sequences')
        extend_alignments (alignmentFile2,
                           proteinSeqFile,
                           chainStrucSeqFile,
                           alignmentFile3)
    
    print('Extracting coordinate files for PPI templates')
    produce_interactome_template_files (interactomeFile,
                                        chainSeqresFile,
                                        chainStrucResFile,
                                        templateDir)
    
    print('Producing PPI alignment files')
    produce_interactome_alignment_files (interactomeFile,
                                         alignmentFile3,
                                         chainSeqresFile,
                                         chainStrucResFile,
                                         annotatedInteractomeFile,
                                         verbose = verbose)
    
if __name__ == "__main__":
    main()