import os
import io
import pickle
import warnings
import pandas as pd
from pathlib import Path
from modeller import *
from modeller.automodel import *
from interactome_tools import read_single_interface_annotated_interactome
from pdb_tools import (allow_pdb_downloads,
                       suppress_pdb_warnings,
                       load_pdbtools_chain_sequences,
                       load_pdbtools_chain_strucRes_labels,
                       write_partial_structure,
                       get_chain_IDs,
                       ordered_residue_ID_at_index)

# directory for PDB structure files
pdbDir = Path('../pdb_files')

# directory for template structure files
templateDir = Path('../templates')

# directory for alignment files
alignmentDir = Path('../alignments')

# directory for output models
modelDir = Path('../models')

def set_pdb_dir (dir):
    
    global pdbDir
    pdbDir = dir
    if not pdbDir.exists():
        os.makedirs(str(pdbDir))

def set_template_dir (dir):
    
    global templateDir
    templateDir = dir
    if not templateDir.exists():
        os.makedirs(str(templateDir))

def set_alignment_dir (dir):
    
    global alignmentDir
    alignmentDir = dir
    if not alignmentDir.exists():
        os.makedirs(str(alignmentDir))

def set_model_dir (dir):
    
    global modelDir
    modelDir = dir
    if not modelDir.exists():
        os.makedirs(str(modelDir))

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

def produce_interactome_template_files (inPath, outDir):
    
    interactome = read_single_interface_annotated_interactome (inPath)
    n = len(interactome)
    
    for i, row in interactome.iterrows():
        sys.stdout.write('  PPI %d out of %d (%.2f%%) \r' % (i+1, n, 100*(i+1)/n))
        sys.stdout.flush()
        for chain1, chain2 in row.Chain_pairs:
            pdbid, chainID1 = chain1.split('_')
            _, chainID2 = chain2.split('_')
            selectChains = sorted([chainID1, chainID2])
            templateID = '_'.join([pdbid] + selectChains)
            write_partial_structure (pdbid,
                                     selectChains,
                                     pdbDir,
                                     outDir / (templateID + '.ent'))

def create_interactome_models (interactomeFile,
                               alignmentFile,
                               chainSeqresFile,
                               chainStrucResFile,
                               numModels = 1,
                               verbose = True):
    
    load_pdbtools_chain_sequences (chainSeqresFile)
    load_pdbtools_chain_strucRes_labels (chainStrucResFile)
    interactome = read_single_interface_annotated_interactome (interactomeFile)
    alignments = pd.read_table (alignmentFile, sep="\t")
    n = len(interactome)
    
    for i, row in interactome.iterrows():
        sys.stdout.write('  PPI %d out of %d (%.2f%%) \r' % (i+1, n, 100*(i+1)/n))
        sys.stdout.flush()
        create_complex_model ([row.Protein_1, row.Protein_2],
                              row.Chain_pairs[0],
                              alignments,
                              numModels = numModels,
                              verbose = verbose)

def create_complex_model (proteins,
                          templates,
                          alignments,
                          numModels = 1,
                          verbose = True):
    
    pdbid, _ = templates[0].split('_')
    chainIDs = [id.split('_')[1] for id in templates]
    templateMap = {c:p for c, p in zip(chainIDs, proteins)}
    templateID = '_'.join([pdbid] + sorted(chainIDs))
    
    chainIDs = get_chain_IDs (templateID, templateDir)
    orderedProteins = [templateMap[c] for c in chainIDs]
    
    complexID = '-'.join(orderedProteins)
    alignmentID = '_'.join([complexID, '-'.join([pdbid] + chainIDs)])
    
#     create_complex_alignment (orderedProteins,
#                               templateID,
#                               orderedChainIDs,
#                               alignments,
#                               alignmentDir / alignmentID)
    
#     create_protein_model (complexID,
#                           templateID,
#                           alignmentDir / alignmentID,
#                           starting_model = 1,
#                           ending_model = numModels,
#                           verbose = True)

# def create_complex_alignment (proteins,
#                               templateID,
#                               chainIDs,
#                               alignments,
#                               outPath):
    
#     complexID = '-'.join(orderedProteins)
#     pdbid = templateID.split('_')[0] if '_' in templateID else templateID
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
    else:
        warnings.warn("Full alignment for complex %s not found. Alignment file not created" % complexID)

def write_alignment (complexID,
                     templateID,
                     chainIDs,
                     alignments,
                     outPath):

    pr_alignments, ch_alignments = alignments
    het, resNum, icode = ordered_residue_ID_at_index ('first', templateID, chainIDs[0], templateDir)
    firstResID = str(resNum) if icode == ' ' else str(resNum) + icode
    het, resNum, icode = ordered_residue_ID_at_index ('last', templateID, chainIDs[-1], templateDir)
    lastResID = str(resNum) if icode == ' ' else str(resNum) + icode
    with io.open(outPath, "w") as fout:
        fout.write ('>P1;%s\n' % pdbid)
        fout.write ('structure:%s:%s:%s:%s:%s::::\n' % (templateID,
                                                        firstResID,
                                                        chainIDs[0],
                                                        lastResID,
                                                        chainIDs[-1]))
        fout.write ('/'.join(ch_alignments) + '*\n\n')
        fout.write ('>P1;%s\n' % complexID)
        fout.write ('sequence:%s:.:.:.:.::::\n' % complexID)
        fout.write ('/'.join(pr_alignments) + '*\n')
