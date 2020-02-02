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
                             produce_interactome_alignment_files)

def main():
    
    # reference interactome name
    # options: HI-II-14, IntAct
    interactome_name = 'HI-II-14'
    
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
    
    # directory of processed template-related data files specific to interactome
    templateBasedDir = interactomeDir / 'template_based'
    
    # directory of processed model-related data files specific to interactome
    modelBasedDir = interactomeDir / 'model_based'
    
    # directory for PDB structure files
    pdbDir = Path('/Volumes/MG_Samsung/pdb_files')
    
    # directory for template structure files
    templateDir = Path('../templates')
    
    # directory for alignment files
    alignmentDir = Path('../alignments')
    
    # input data files
    proteinSeqFile = procDir / 'human_reference_sequences.pkl'
    chainStrucSeqFile = modelBasedDir / 'interactome_template_sequences.pkl'
    chainSeqresFile = templateBasedDir / 'protein_chain_sequences.pkl'
    chainStrucResFile = templateBasedDir / 'protein_chain_strucRes.pkl'
    interactomeFile = templateBasedDir / 'human_structural_interactome.txt'
    blastFile = modelBasedDir / 'interactome_template_blast_alignments_e100'
    
    # output data files
    alignmentFile1 = modelBasedDir / 'ppi_template_alignments.txt'
    alignmentFile2 = modelBasedDir / 'ppi_template_filtered_alignments.txt'
    alignmentFile3 = modelBasedDir / 'ppi_template_extended_alignments.txt'
    annotatedInteractomeFile = modelBasedDir / 'human_template_annotated_interactome.txt'    
    
    # create output directories if not existing
    if not modelBasedDir.exists():
        os.makedirs(modelBasedDir)
    if not templateDir.exists():
        os.makedirs(templateDir)
    if not alignmentDir.exists():
        os.makedirs(alignmentDir)
    
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
        # This is only to calculate alignment coverage, and remove duplicate alignments
        # for each protein-chain pair. Filtering by evalue is not needed at this point.
        filter_chain_annotations (alignmentFile1,
                                  alignmentFile2,
                                  evalue = 1000,
                                  prCov = 0,
                                  chCov = 0)
    
    if not alignmentFile3.is_file():
        print('Extending alignments to full sequences')
        extend_alignments (alignmentFile2,
                           proteinSeqFile,
                           chainStrucSeqFile,
                           alignmentFile3)
    
    print('Producing PPI alignment files')
    produce_interactome_alignment_files (interactomeFile,
                                         alignmentFile3,
                                         chainSeqresFile,
                                         chainStrucResFile,
                                         annotatedInteractomeFile,
                                         verbose = verbose)
    
if __name__ == "__main__":
    main()