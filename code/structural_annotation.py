#----------------------------------------------------------------------------------------
# Modules for protein structural annotation.
#----------------------------------------------------------------------------------------

import io
import time
import pickle
import pandas as pd
import numpy as np
from pathlib import Path
from interactome_tools import read_single_interface_annotated_interactome
from pdb_tools import (load_pdbtools_chain_sequences,
                       load_pdbtools_chain_strucRes_labels,
                       produce_chain_struc_sequences)

def filter_chain_annotations (inPath,
                              outPath,
                              evalue = 1e-5,
                              prCov = 0,
                              chCov = 0):
    """Filter protein chain annotations by coverage and e-value.

    Args:
        inPath (Path): path to tab-deleimited file containing protein-chain alignments.
        outPath (Path): file path to save filtered alignments to.
        evalue (numeric): maximum alignment e-value cutoff.
        prCov (numeric): minimum alignment coverage fraction required on protein sequence.
        chCov (numeric): minimum alignment coverage fraction required on chain sequence.

    """
    with io.open(inPath, "r", encoding="utf8", errors='ignore') as f, io.open(outPath, "w") as fout:
        headers = f.readline().strip()
        fout.write(headers + '\t' + 'Qcov' + '\t' + 'Scov' + '\n')
        headerSplit = headers.split('\t')
        numCol = len(headerSplit)
        evalueCol = headerSplit.index('Expect')
        QseqCol = headerSplit.index('Qseq')
        SseqCol = headerSplit.index('Sseq')
        QlenCol = headerSplit.index('Qlen')
        SlenCol = headerSplit.index('Slen')
    with io.open(inPath, "r", encoding="utf8", errors='ignore') as f, io.open(outPath, "a") as fout:
        next(f)
        for line in f:
            line = line.strip()
            linesplit = line.split('\t')
            if len(linesplit) == numCol:
                if float(linesplit[evalueCol]) < evalue:
                    Qseq = linesplit[QseqCol]
                    Sseq = linesplit[SseqCol]
                    Qcov = (len(Qseq) - Qseq.count('-')) / int(linesplit[QlenCol])
                    Scov = (len(Sseq) - Sseq.count('-')) / int(linesplit[SlenCol])
                    if (Qcov >= prCov) and (Scov >= chCov):
                        fout.write('\t'.join([line, str(Qcov), str(Scov)]) + '\n')
    remove_duplicate_chain_annotations(outPath, outPath)

def remove_duplicate_chain_annotations (inPath, outPath):
    """Keep only one alignment for each protein-chain pair, the one with the smallest 
        e-value.

    Args:
        inPath (Path): path to tab-deleimited file containing protein-chain alignments.
        outPath (Path): file path to save filtered alignments to.

    """
    chainMap = pd.read_table(inPath, sep='\t')
    chainMap = chainMap.sort_values("Expect", axis=0, ascending=True)
    chainMap = chainMap.drop_duplicates(subset = ["Query", "Subject"], keep='first')
    chainMap = chainMap.sort_values(['Query','Subject'], axis=0, ascending=True)
    chainMap.to_csv(outPath, index=False, sep='\t')

def filter_chain_annotations_by_protein (inPath, proteins, outPath):
    """Filter protein chain annotations by proteins.

    Args:
        inPath (Path): path to tab-deleimited file containing protein-chain alignments.
        proteins (list): proteins to be selected.
        outPath (Path): file path to save filtered alignments to.

    """
    with io.open(inPath, "r", errors='ignore') as f, io.open(outPath, "w") as fout:
        headers = f.readline()
        fout.write(headers)
        headerSplit = headers.strip().split('\t')
        numCol = len(headerSplit)
        queryPos = headerSplit.index('Query')
    with io.open(inPath, "r", errors='ignore') as f, io.open(outPath, "a") as fout:
        next(f)
        for line in f:
            linesplit = line.strip().split('\t')
            if len(linesplit) == numCol:
                if linesplit[queryPos] in proteins:
                    fout.write(line)

def single_chain_per_protien (inPath, outPath):
    """Keep only one chain alignment for each protein, the one with the smallest e-value.

    Args:
        inPath (Path): path to tab-deleimited file containing protein-chain alignments.
        outPath (Path): file path to save filtered alignments to.

    """
    chainMap = pd.read_table (inPath, sep='\t')
    chainMap = chainMap.sort_values ("Expect", axis=0, ascending=True)
    chainMap = chainMap.drop_duplicates (subset = "Query", keep='first')
    chainMap = chainMap.sort_values ("Query", axis=0, ascending=True)
    chainMap.to_csv (outPath, index=False, sep='\t')

def write_interactome_template_sequences (interactomeFile,
                                          chainSeqresFile,
                                          chainStrucResFile,
                                          pdbDir,
                                          outPath):

    load_pdbtools_chain_sequences (chainSeqresFile)
    load_pdbtools_chain_strucRes_labels (chainStrucResFile)
    interactome = read_single_interface_annotated_interactome (interactomeFile)
    chainIDs = []
    for templates in interactome["Chain_pairs"].values:
        for t in templates:
            chainIDs.extend(t)
    produce_chain_struc_sequences (chainIDs, pdbDir, outPath)
