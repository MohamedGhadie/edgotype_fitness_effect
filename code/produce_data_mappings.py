import os
from pathlib import Path
from text_tools import parse_fasta_file, reduce_fasta_headers
from id_mapping import (produce_geneName_dict,
                        produce_uniqueGene_swissProtIDs,
                        produce_uniqueGene_sequences,
                        produce_proteinSeq_dict,
                        produce_uniprotID_dict,
                        produce_chainSeq_dict,
                        produce_chain_strucRes_dict)
from interactome_tools import write_interactome_sequences
from structural_annotation import write_interactome_template_sequences

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
        
    # directory for PDB structure files
    pdbDir = Path('/Volumes/MG_Samsung/pdb_files')
    
    # input data files
    uniprotRefSeqFile = extDir / 'UP000005640_9606.fasta'
    idmapFile = extDir / 'HUMAN_9606_idmapping.dat'
    proteomeListFile = extDir / 'uniprot_reviewed_human_proteome.list'
    pdbSeqresFile = extDir / 'pdb_seqres.txt'
    chainResAnnotFile = extDir / 'ss_dis.txt'
    interactomeFile = interactomeDir / 'human_template_annotated_interactome.txt'
    
    # output data files
    refSeqFile = procDir / 'human_reference_sequences.fasta'
    sequenceFile = procDir / 'human_reference_sequences.txt'
    geneMapFile = procDir / 'to_human_geneName_map.pkl'
    uniqueGeneSwissProtIDFile = procDir / 'uniprot_unique_gene_reviewed_human_proteome.list'
    uniqueGeneSequenceFile = procDir / 'human_unique_gene_reference_sequences.txt'
    proteinSeqFile = procDir / 'human_reference_sequences.pkl'
    uniprotIDmapFile = procDir / 'to_human_uniprotID_map.pkl'
    seqresFile = procDir / 'pdb_seqres_reduced.fasta'
    chainSeqresFile = procDir / 'chain_seqres.pkl'
    chainStrucResFile = procDir / 'chain_strucRes.pkl'
    interactomeSequenceFile = interactomeDir / 'interactome_sequences.fasta'
    chainStrucSeqFastaFile = interactomeDir / 'chain_struc_sequences.fasta'
    chainStrucSeqFile = interactomeDir / 'chain_struc_sequences.pkl'
    
    # create output directories if not existing
    if not dataDir.exists():
        os.makedirs(dataDir)
    if not extDir.exists():
        os.makedirs(extDir)
    if not procDir.exists():
        os.makedirs(procDir)
    if not interactomeDir.exists():
        os.makedirs(interactomeDir)
    if not pdbDir.exists():
        os.makedirs(pdbDir)
    
    if not refSeqFile.is_file():
        print('Reducing headers in protein sequence fasta file')
        reduce_fasta_headers (uniprotRefSeqFile, '|', 2, 2, refSeqFile)
    
    if not sequenceFile.is_file():
        print('reading protein sequence fasta file')
        parse_fasta_file (refSeqFile, sequenceFile)
    
    if not geneMapFile.is_file():
        print('producing UniProtID-to-geneName dictionary')
        produce_geneName_dict (idmapFile, proteomeListFile, geneMapFile)
    
    if not uniqueGeneSwissProtIDFile.is_file():
        print('producing list of unique-gene UniProt IDs')
        produce_uniqueGene_swissProtIDs (proteomeListFile, geneMapFile, uniqueGeneSwissProtIDFile)
    
    if not uniqueGeneSequenceFile.is_file():
        print('producing sequence file for unique-gene UniProt IDs')
        produce_uniqueGene_sequences (sequenceFile, uniqueGeneSwissProtIDFile, geneMapFile, uniqueGeneSequenceFile)
    
    if not proteinSeqFile.is_file():
        print('producing protein sequence dictionary')
        produce_proteinSeq_dict (refSeqFile, proteinSeqFile)
    
    if not uniprotIDmapFile.is_file():
        print('producing to-UniProt-ID dictionary')
        produce_uniprotID_dict (idmapFile, uniqueGeneSwissProtIDFile, uniprotIDmapFile)
    
    if not seqresFile.is_file():
        print('reducing headers in PDB chain sequence file')
        reduce_fasta_headers (pdbSeqresFile, ' ', 1, 1, seqresFile)
    
    if not chainSeqresFile.is_file():
        print('producing PDB chain sequence dictionary from fasta records')
        produce_chainSeq_dict (seqresFile, chainSeqresFile)
    
    if not chainStrucResFile.is_file():
        print('parsing PDB chain structured-residue order file')
        produce_chain_strucRes_dict (chainResAnnotFile, chainStrucResFile)
    
    if not interactomeSequenceFile.is_file():
        print('writing interactome protein sequences')
        write_interactome_sequences (interactomeFile,
                                     uniqueGeneSequenceFile,
                                     interactomeSequenceFile)
    
    if not chainStrucSeqFastaFile.is_file():
        print('writing interactome template sequences to Fasta file')
        write_interactome_template_sequences (interactomeFile,
                                              chainSeqresFile,
                                              chainStrucResFile,
                                              pdbDir,
                                              chainStrucSeqFastaFile)
    
    if not chainStrucSeqFile.is_file():
        print('producing template sequence dictionary')
        produce_chainSeq_dict (chainStrucSeqFastaFile, chainStrucSeqFile)

if __name__ == "__main__":
    main()
