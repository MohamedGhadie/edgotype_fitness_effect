import os
import sys
import io
import warnings
import pickle
import pandas as pd
from pathlib import Path
from interactome_tools import read_single_interface_annotated_interactome
from pdb_tools import (allow_pdb_downloads,
                       suppress_pdb_warnings,
                       load_pdbtools_chain_sequences,
                       load_pdbtools_chain_strucRes_labels,
                       write_partial_structure,
                       get_chain_IDs,
                       ordered_residue_IDs,
                       ordered_residue_ID)

# directory for PDB structure files
pdbDir = Path('../pdb_files')

# directory for template structure files
templateDir = Path('../templates')

# directory for alignment files
alignmentDir = Path('../alignments')

def set_pdb_dir (dir):
    
    global pdbDir
    pdbDir = dir
    if not pdbDir.exists():
        os.makedirs(pdbDir)

def set_template_dir (dir):
    
    global templateDir
    templateDir = dir
    if not templateDir.exists():
        os.makedirs(templateDir)

def set_alignment_dir (dir):
    
    global alignmentDir
    alignmentDir = dir
    if not alignmentDir.exists():
        os.makedirs(alignmentDir)

def enable_pdb_downloads (download):
    """Set global variable in pdb_tools module to allow PDB downloads.

    Args:
        download (bool): if true, allow PDB downloads.

    """
    allow_pdb_downloads (download)

def disable_pdb_warnings (suppress):
    """Set global variable in pdb_tools module to suppress PDB warnings.

    Args:
        suppress (bool): if true, suppress PDB warnings.

    """
    suppress_pdb_warnings (suppress)

def extend_alignments (inPath, querySeqFile, subjectSeqFile, outPath):
    
    alignments = pd.read_table (inPath, sep="\t")
    with open(querySeqFile, 'rb') as f:
        querySeq = pickle.load(f)
    with open(subjectSeqFile, 'rb') as f:
        subjectSeq = pickle.load(f)
    
    query_alignments, subject_alignments = [], []
    for _, row in alignments.iterrows():
        Qalign, Salign = extend_alignment (querySeq[row.Query],
                                           subjectSeq[row.Subject],
                                           row.Qseq,
                                           row.Sseq,
                                           row.Qstart,
                                           row.Qend,
                                           row.Sstart,
                                           row.Send)
        query_alignments.append(Qalign)
        subject_alignments.append(Salign)
    
    alignments["Qseq"] = query_alignments
    alignments["Sseq"] = subject_alignments
    alignments.to_csv (outPath, index=False, sep='\t')

def extend_alignment (FullQseq,
                      FullSseq,
                      Qseq,
                      Sseq,
                      Qstart,
                      Qend,
                      Sstart,
                      Send):
    
    leftExt, rightExt = Qstart - 1, len(FullQseq) - Qend
    Qalign = FullQseq[:leftExt] + Qseq
    Salign = ('-' * leftExt) + Sseq
    if rightExt > 0:
        Qalign = Qalign + FullQseq[-rightExt:]
        Salign = Salign + ('-' * rightExt)
    
    leftExt, rightExt = Sstart - 1, len(FullSseq) - Send
    Salign = FullSseq[:leftExt] + Salign
    Qalign = ('-' * leftExt) + Qalign
    if rightExt > 0:
        Salign = Salign + FullSseq[-rightExt:]
        Qalign = Qalign + ('-' * rightExt)
    
    return Qalign, Salign

def produce_interactome_template_files (inPath,
                                        chainSeqresFile,
                                        chainStrucResFile,
                                        outDir):
    
    load_pdbtools_chain_sequences (chainSeqresFile)
    load_pdbtools_chain_strucRes_labels (chainStrucResFile)
    interactome = read_single_interface_annotated_interactome (inPath)
    n = len(interactome)
    
    for i, row in interactome.iterrows():
        sys.stdout.write('  PPI %d out of %d (%.2f%%) \r' % (i+1, n, 100*(i+1)/n))
        sys.stdout.flush()
        for chain1, chain2 in row.Chain_pairs:
            pdbid, chainID1 = chain1.split('_')
            _, chainID2 = chain2.split('_')
            selectChains = sorted([chainID1, chainID2])
            templateID = '-'.join([pdbid] + selectChains)
            outFile = outDir / ('pdb' + templateID + '.ent')
            if not outFile.is_file():
                resIDs = {c:ordered_residue_IDs (pdbid, c, pdbDir) for c in selectChains}
                write_partial_structure (pdbid,
                                         selectChains,
                                         pdbDir,
                                         outFile,
                                         resIDs = resIDs)
    print()

def produce_interactome_alignment_files (interactomeFile,
                                         alignmentFile,
                                         chainSeqresFile,
                                         chainStrucResFile,
                                         outPath,
                                         numModels = 1,
                                         verbose = True):
    
    load_pdbtools_chain_sequences (chainSeqresFile)
    load_pdbtools_chain_strucRes_labels (chainStrucResFile)
    interactome = read_single_interface_annotated_interactome (interactomeFile)
    alignments = pd.read_table (alignmentFile, sep="\t")
    n = len(interactome)
    
    allIDs = []
    for i, row in interactome.iterrows():
        sys.stdout.write('  PPI %d out of %d (%.2f%%) \r' % (i+1, n, 100*(i+1)/n))
        sys.stdout.flush()
        IDs = create_complex_alignment ([row.Protein_1, row.Protein_2],
                                        row.Chain_pairs[0],
                                        alignments,
                                        numModels = numModels,
                                        verbose = verbose)
        allIDs.append(IDs)
    print()
    interactome["Complex_ID"], interactome["Template_ID"], interactome["Alignment_ID"] = zip(* allIDs)
    interactome = interactome [interactome["Alignment_ID"] != '-']
    interactome.drop(["Interfaces",	"Chain_pairs"], axis=1, inplace=True)
    interactome.to_csv(outPath, index=False, sep='\t')

def create_complex_alignment (proteins,
                              templates,
                              alignments,
                              numModels = 1,
                              verbose = True):
    
    pdbid, _ = templates[0].split('_')
    chainIDs = [id.split('_')[1] for id in templates]
    templateMap = {c:p for c, p in zip(chainIDs, proteins)}
    templateID = '-'.join([pdbid] + sorted(chainIDs))
    
    chainIDs = get_chain_IDs (templateID, templateDir)
    orderedProteins = [templateMap[c] for c in chainIDs]
    
    complexID = '-'.join(orderedProteins)
    alignmentID = '_'.join([complexID, '-'.join([pdbid] + chainIDs)])
    
    pr_alignments, ch_alignments = [], []
    for protein, chainID in zip(orderedProteins, chainIDs):
        align = alignments [(alignments["Query"] == protein) & 
                            (alignments["Subject"] == '_'.join([pdbid, chainID]))]
        if not align.empty:
            prAlign, chAlign = align[["Qseq", "Sseq"]].values[0]
            pr_alignments.append(prAlign)
            ch_alignments.append(chAlign)
    
    if len(pr_alignments) == len(orderedProteins):
        write_alignment (complexID,
                         templateID,
                         chainIDs,
                         [pr_alignments, ch_alignments],
                         alignmentDir / alignmentID)
        return complexID, templateID, alignmentID
    else:
        warnings.warn("Full alignment for complex %s not found. Alignment file not created" % complexID)
        return '-', '-', '-'

def write_alignment (complexID,
                     templateID,
                     chainIDs,
                     alignments,
                     outPath):

    pr_alignments, ch_alignments = alignments
    het, resNum, icode = ordered_residue_ID ('first', templateID, chainIDs[0], templateDir)
    firstResID = str(resNum) if icode == ' ' else str(resNum) + icode
    het, resNum, icode = ordered_residue_ID ('last', templateID, chainIDs[-1], templateDir)
    lastResID = str(resNum) if icode == ' ' else str(resNum) + icode
    with io.open(outPath, "w") as fout:
        fout.write ('>P1;%s\n' % templateID)
        fout.write ('structure:%s:%s:%s:%s:%s::::\n' % (templateID,
                                                        firstResID,
                                                        chainIDs[0],
                                                        lastResID,
                                                        chainIDs[-1]))
        fout.write ('/'.join(ch_alignments) + '*\n\n')
        fout.write ('>P1;%s\n' % complexID)
        fout.write ('sequence:%s:.:.:.:.::::\n' % complexID)
        fout.write ('/'.join(pr_alignments) + '*\n')
