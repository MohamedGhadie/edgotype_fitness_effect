import os
import io
import sys
import warnings
import pickle
import pandas as pd
import numpy as np
from pathlib import Path
from Bio.PDB import DSSP
from simple_tools import segment_overlaps
from text_tools import read_list_table
from interactome_tools import (read_single_interface_annotated_interactome,
                               get_protein_interface)
from pdb_tools import (allow_pdb_downloads,
                       suppress_pdb_warnings,
                       load_pdbtools_chain_sequences,
                       load_pdbtools_chain_strucRes_labels,
                       return_structure,
                       count_neighbors_by_chainIDs,
                       return_chain_sequence,
                       return_multichain_res_posToIDs,
                       return_chain_res_IDToPos,
                       return_chain_res_posToID,
                       return_chain_res_posToIDs,
                       get_distance)
from rsa_tools import (load_empirical_maxAcc,
                       read_dssp_file,
                       acc_to_rsa,
                       calculate_structure_rsa)

pdbDir = Path('../pdb_files')
accDir = Path('../res_acc')
downloadStructures = True
suppressWarnings = True
mutcount = proteincount = 0
#centerDist = {}
distFile = Path('../data/processed/tmpDist.pkl')

def set_pdb_dir (dir):
    
    global pdbDir
    pdbDir = dir
    if not pdbDir.exists():
        os.makedirs(pdbDir)

def set_res_accDir (dir):
    
    global accDir
    accDir = dir
    if not accDir.exists():
        os.makedirs(accDir)

def load_dictionaries (chainSequencePath = None, empMaxAccPath = None):
    
    if chainSequencePath:
        load_pdbtools_chain_sequences(chainSequencePath)
    if empMaxAccPath:
        load_empirical_maxAcc (empMaxAccPath)

def produce_empirical_maxAcc (pdbIDs, dsspDir, outPath, perc = 99.99):
    
    residues = ['A','R','N','D','C','E','Q','G','H','I',
                'L','K','M','F','P','S','T','W','Y','V']
    resAcc = { res : [] for res in residues }
    l = len(pdbIDs)
    for pdbid in pdbIDs:
        dsspFile = dsspDir / (pdbid + '.dssp')
        if dsspFile.is_file():
            acc = read_dssp_file(dsspFile)
            for saTup in acc.values():
                res, sa = saTup[ 1 ], saTup[ 3 ]
                if res in resAcc:
                    resAcc[ res ].append(sa)
    maxAcc = {}
    for res, sa in resAcc.items():
        maxAcc[ res ] = np.percentile(sa, perc) if sa else np.nan
    print('Maximum residue accessibility values:')
    print(maxAcc)
    with open(outPath, 'wb') as fOut:
        pickle.dump(maxAcc, fOut)

def protein_residue_RSA (positions,
                         pdbChainsFile,
                         chainSeqFile,
                         chainMapFile,
                         chainStrucResFile,
                         pdbDir,
                         accDir,
                         dsspDir = None,
                         maxAccFile = None,
                         pdbMapFile = None,
                         prCov = 0,
                         chCov = 0,
                         allChainMap = False,
                         singleChain = True,
                         mapToProtein = True,
                         downloadPDB = True,
                         suppressWarnings = True):
    
    set_pdb_dir (pdbDir)
    set_res_accDir (accDir)
    allow_pdb_downloads (downloadPDB)
    suppress_pdb_warnings (suppressWarnings)
    with open(pdbChainsFile, 'rb') as f:
        pdbChains = pickle.load(f)
    if pdbMapFile:
        with open(pdbMapFile, 'rb') as f:
            toPDB = pickle.load(f)
    else:
        toPDB = {}
    if maxAccFile:
        load_empirical_maxAcc (maxAccFile)
    chainMap = read_list_table (chainMapFile, ["Qpos", "Spos"], [int, int], '\t')
    load_pdbtools_chain_sequences (chainSeqFile)
    load_pdbtools_chain_strucRes_labels (chainStrucResFile)
    
    resRSA = {}
    n = len(positions)
    for i, (protein, pos) in enumerate(positions):
        sys.stdout.write('  protein %d out of %d (%.2f%%) \r' % (i+1, n, 100*(i+1)/n))
        sys.stdout.flush()
        rsa =  res_rsa (protein,
                        [pos],
                        pdbChains,
                        chainMap,
                        dsspDir = dsspDir,
                        pdbIDs = toPDB[protein] if protein in toPDB else None,
                        prCov = prCov,
                        chCov = chCov,
                        allChainMap = allChainMap,
                        singleChain = singleChain,
                        mapToProtein = mapToProtein)
        if rsa:
            resRSA[(protein, pos)] = rsa[pos]
    print()
    return resRSA

def produce_protein_model_RSA (prLen,
                               pdbChainsFile,
                               chainSeqFile,
                               chainMapFile,
                               chainStrucResFile,
                               pdbDir,
                               accDir,
                               outPath,
                               dsspDir = None,
                               maxAccFile = None,
                               pdbMapFile = None,
                               prCov = 0,
                               chCov = 0,
                               allChainMap = False,
                               singleChain = True,
                               mapToProtein = True,
                               downloadPDB = True,
                               suppressWarnings = True):
    
    set_pdb_dir (pdbDir)
    set_res_accDir (accDir)
    allow_pdb_downloads (downloadPDB)
    suppress_pdb_warnings (suppressWarnings)
    with open(pdbChainsFile, 'rb') as f:
        pdbChains = pickle.load(f)
    if pdbMapFile:
        with open(pdbMapFile, 'rb') as f:
            toPDB = pickle.load(f)
    else:
        toPDB = {}
    if maxAccFile:
        load_empirical_maxAcc (maxAccFile)
    chainMap = read_list_table (chainMapFile, ["Qpos", "Spos"], [int, int], '\t')
    load_pdbtools_chain_sequences (chainSeqFile)
    load_pdbtools_chain_strucRes_labels (chainStrucResFile)
    
    rsa = {}
    n = len(prLen)
    for i, (protein, ln) in enumerate(prLen.items()):
        sys.stdout.write('  protein %d out of %d (%.2f%%) \r' % (i+1, n, 100*(i+1)/n))
        sys.stdout.flush()
        rsa[protein] =  res_rsa (protein,
                                 np.arange(1, ln + 1),
                                 pdbChains,
                                 chainMap,
                                 dsspDir = dsspDir,
                                 pdbIDs = toPDB[protein] if protein in toPDB else None,
                                 prCov = prCov,
                                 chCov = chCov,
                                 allChainMap = allChainMap,
                                 singleChain = singleChain,
                                 mapToProtein = mapToProtein)
    with open(outPath, 'wb') as fOut:
        pickle.dump(rsa, fOut)
    print()

def produce_protein_surface_residues (rsaFile,
                                      outPath,
                                      minRSA = 0.25,
                                      downloadPDB = True,
                                      suppressWarnings = True):
    
    allow_pdb_downloads (downloadPDB)
    suppress_pdb_warnings (suppressWarnings)
    with open(rsaFile, 'rb') as f:
        rsaDict = pickle.load(f)
    
    surfaceRes = {}
    for p, rsa in rsaDict.items():
        surfaceRes[p] = list(np.where(rsa >= minRSA)[0] + 1)
    with open(outPath, 'wb') as fOut:
        pickle.dump(surfaceRes, fOut)
    
def mutation_dist_to_surface (mutations,
                              surfaceResFile,
                              pdbChainsFile,
                              chainMapFile,
                              chainStrucResFile,
                              pdbDir,
                              prCov = 0,
                              chCov = 0,
                              pdbMapFile = None,
                              allChainMap = False,
                              singleChain = True,
                              downloadPDB = True,
                              suppressWarnings = True):
    
    set_pdb_dir (pdbDir)
    allow_pdb_downloads (downloadPDB)
    suppress_pdb_warnings (suppressWarnings)
    with open(surfaceResFile, 'rb') as f:
        surfaceRes = pickle.load(f)
    with open(pdbChainsFile, 'rb') as f:
        pdbChains = pickle.load(f)
    if pdbMapFile:
        with open(pdbMapFile, 'rb') as f:
            toPDB = pickle.load(f)
    else:
        toPDB = {}
    chainMap = read_list_table(chainMapFile, ["Qpos", "Spos"], [int, int], '\t')
    load_pdbtools_chain_strucRes_labels (chainStrucResFile)
    
    distances = []
    n = len(mutations)
    for i, (protein, mut_pos) in enumerate(mutations[["protein", "mut_position"]].values):
        sys.stdout.write('  mutation %d out of %d (%.2f%%) \r' % (i+1, n, 100*(i+1)/n))
        sys.stdout.flush()
        dist = res_dist_to_surface (protein,
                                    mut_pos,
                                    surfaceRes[protein] if protein in surfaceRes else [],
                                    pdbChains,
                                    chainMap,
                                    pdbIDs = toPDB[protein] if protein in toPDB else None,
                                    prCov = prCov,
                                    chCov = chCov,
                                    allChainMap = allChainMap,
                                    singleChain = singleChain)
        distances.append(dist)
    print()
    return distances

def multiprotein_res_dist_to_surface (proteins,
                                      surfaceResFile,
                                      pdbChainsFile,
                                      chainMapFile,
                                      chainStrucResFile,
                                      pdbDir,
                                      prCov = 0,
                                      chCov = 0,
                                      pdbMapFile = None,
                                      allChainMap = False,
                                      singleChain = True,
                                      downloadPDB = True,
                                      suppressWarnings = True):
    
    set_pdb_dir (pdbDir)
    allow_pdb_downloads (downloadPDB)
    suppress_pdb_warnings (suppressWarnings)
    with open(surfaceResFile, 'rb') as f:
        surfaceRes = pickle.load(f)
    with open(pdbChainsFile, 'rb') as f:
        pdbChains = pickle.load(f)
    if pdbMapFile:
        with open(pdbMapFile, 'rb') as f:
            toPDB = pickle.load(f)
    else:
        toPDB = {}
    chainMap = read_list_table (chainMapFile, ["Qpos", "Spos"], [int, int], '\t')
    load_pdbtools_chain_strucRes_labels (chainStrucResFile)
    
    distances = []
    n = len(proteins)
    for i, (protein, pr_len) in enumerate(proteins[["protein", "len"]].values):
        print('protein %d out of %d (%f%%)' % (i+1, n, 100*(i+1)/n))
        dist = multires_dist_to_surface (protein,
                                         np.arange(1, pr_len + 1),
                                         surfaceRes[protein] if protein in surfaceRes else [],
                                         pdbChains,
                                         chainMap,
                                         #pdbIDs = toPDB[protein] if protein in toPDB else None,
                                         prCov = prCov,
                                         chCov = chCov,
                                         allChainMap = allChainMap,
                                         singleChain = singleChain)
        distances.append(dist)
    return distances

def multires_dist_to_surface (protein,
                              resPos,
                              surface,
                              pdbChains,
                              chainMap,
                              pdbIDs = None,
                              prCov = 0,
                              chCov = 0,
                              allChainMap = False,
                              singleChain = True):
    
    distances = multires_dist_to_residues (protein,
                                           resPos,
                                           surface,
                                           pdbChains,
                                           chainMap,
                                           pdbIDs = pdbIDs,
                                           prCov = prCov,
                                           chCov = chCov,
                                           allChainMap = allChainMap,
                                           singleChain = singleChain)
    minDist = []
    if distances:
        for chainID, chainPos, dist in distances:
            dist = dist[np.isnan(dist) == False]
            if dist.size > 0:
                minDist.append( (chainID, chainPos, min(dist)) )
            else:
                minDist.append( (chainID, chainPos, np.nan) )
    return minDist
#     return np.array( list( map( min, distances ) ) )

def res_dist_to_surface (protein,
                         resPos,
                         surface,
                         pdbchains,
                         chainMap,
                         pdbIDs = None,
                         prCov = 0,
                         chCov = 0,
                         allChainMap = False,
                         singleChain = True):
    
    distances = dist_to_residues (protein,
                                  resPos,
                                  surface,
                                  pdbchains,
                                  chainMap,
                                  pdbIDs = pdbIDs,
                                  prCov = prCov,
                                  chCov = chCov,
                                  allChainMap = allChainMap,
                                  singleChain = singleChain)
    if distances:
        chainID, chainPos, dist = distances
        dist = dist[np.isnan(dist) == False]
        if dist.size > 0:
            return chainID, chainPos, min(dist)
        else:
            return chainID, chainPos, np.nan
    else:
        return ()

def mutation_dist_to_center (mutations,
                             pdbChainsFile,
                             chainSeqFile,
                             chainMapFile,
                             chainStrucResFile,
                             pdbDir,
                             prCov = 0,
                             chCov = 0,
                             pdbMapFile = None,
                             allChainMap = False,
                             singleChain = True,
                             downloadPDB = True,
                             suppressWarnings = True,
                             tmpFile = None):
    
    set_pdb_dir (pdbDir)
    allow_pdb_downloads (downloadPDB)
    suppress_pdb_warnings (suppressWarnings)
    with open(pdbChainsFile, 'rb') as f:
        pdbChains = pickle.load(f)
    if pdbMapFile:
        with open(pdbMapFile, 'rb') as f:
            toPDB = pickle.load(f)
    else:
        toPDB = {}
    load_pdbtools_chain_sequences (chainSeqFile)
    load_pdbtools_chain_strucRes_labels (chainStrucResFile)
    chainMap = read_list_table (chainMapFile, ["Qpos", "Spos"], [int, int], '\t')
    
#     global centerDist, distFile, mutcount
#     if tmpFile:
#         mutcount = 0
#         distFile = tmpFile
#         if distFile.is_file():
#             with open(distFile, 'rb') as f:
#                 centerDist = pickle.load(f)
    
    distances = []
    n = len(mutations)
    for i, (protein, mut_pos) in enumerate(mutations[["protein", "mut_position"]].values):
        sys.stdout.write('  mutation %d out of %d (%.2f%%) \r' % (i+1, n, 100*(i+1)/n))
        sys.stdout.flush()
        dist = res_dist_to_center (protein,
                                   mut_pos,
                                   pdbChains,
                                   chainMap,
                                   pdbIDs = toPDB[protein] if protein in toPDB else None,
                                   prCov = prCov,
                                   chCov = chCov,
                                   allChainMap = allChainMap,
                                   singleChain = singleChain)
        distances.append(dist)
    print()
    return distances

def multiprotein_model_res_dist_to_center (proteins,
                                           pdbChainsFile,
                                           chainSeqFile,
                                           chainMapFile,
                                           chainStrucResFile,
                                           pdbDir,
                                           prCov = 0,
                                           chCov = 0,
                                           pdbMapFile = None,
                                           allChainMap = False,
                                           singleChain = True,
                                           downloadPDB = True,
                                           suppressWarnings = True):
    
    set_pdb_dir (pdbDir)
    allow_pdb_downloads (downloadPDB)
    suppress_pdb_warnings (suppressWarnings)
    with open(pdbChainsFile, 'rb') as f:
        pdbChains = pickle.load(f)
    if pdbMapFile:
        with open(pdbMapFile, 'rb') as f:
            toPDB = pickle.load(f)
    else:
        toPDB = {}
    load_pdbtools_chain_sequences (chainSeqFile)
    load_pdbtools_chain_strucRes_labels (chainStrucResFile)
    chainMap = read_list_table (chainMapFile, ["Qpos", "Spos"], [int, int], '\t')
    
    distances = []
    n = len(proteins)
    for i, protein in enumerate(proteins):
        sys.stdout.write('  protein %d out of %d (%.2f%%) \r' % (i+1, n, 100*(i+1)/n))
        sys.stdout.flush()
        dist = protein_model_res_dist_to_center (protein,
                                                 pdbChains,
                                                 chainMap,
                                                 pdbIDs = toPDB[protein] if protein in toPDB else None,
                                                 prCov = prCov,
                                                 chCov = chCov,
                                                 allChainMap = allChainMap,
                                                 singleChain = singleChain)
        distances.append(dist)
    print()
    return distances

def protein_model_res_dist_to_center (protein,
                                      pdbChains,
                                      chainMap,
                                      pdbIDs = None,
                                      prCov = 0,
                                      chCov = 0,
                                      allChainMap = False,
                                      singleChain = True):
    
    if pdbIDs:
        pdbIDs = sorted(set(pdbIDs))
    else:
        mappingChains = chainMap.loc[chainMap["Query"] == protein, "Subject"].tolist()
        pdbIDs = sorted({x.split('_')[0] for x in mappingChains})
    
    strucMap, cov = best_structure_map (protein,
                                        pdbIDs,
                                        pdbChains,
                                        chainMap,
                                        prCov = prCov,
                                        chCov = chCov,
                                        allChainMap = allChainMap,
                                        singleChain = singleChain)
    
    if strucMap is not None:
        chainIDs = strucMap["Subject"].values
        pdbid = chainIDs[0].split('_')[0]
        chainIDs = {id.split('_')[1] for id in chainIDs}
        distances = multichain_allres_dist_to_center (pdbid, chainIDs)
        distList = []
        for chainID, posdist in distances.items():
            distList.extend([(chainID, pos, dist) for pos, dist in posdist.items()])
        return distList
    else:
        return []
    
def mutation_dist_to_interface (mutations,
                                interactomeFile,
                                pdbChainsFile,
                                chainSeqFile,
                                chainMapFile,
                                chainStrucResFile,
                                pdbDir,
                                prCov = 0,
                                chCov = 0,
                                pdbMapFile = None,
                                allChainMap = False,
                                singleChain = True,
                                downloadPDB = True,
                                suppressWarnings = True):
    
    set_pdb_dir (pdbDir)
    allow_pdb_downloads (downloadPDB)
    suppress_pdb_warnings (suppressWarnings)
    interactome = read_single_interface_annotated_interactome (interactomeFile)
    with open(pdbChainsFile, 'rb') as f:
        pdbChains = pickle.load(f)
    if pdbMapFile:
        with open(pdbMapFile, 'rb') as f:
            toPDB = pickle.load(f)
    else:
        toPDB = {}
    load_pdbtools_chain_sequences (chainSeqFile)
    load_pdbtools_chain_strucRes_labels (chainStrucResFile)
    chainMap = read_list_table (chainMapFile, ["Qpos", "Spos"], [int, int], '\t')
    
    distances = []
    n = len(mutations)
    for i, (protein, mut_pos) in enumerate(mutations[["protein", "mut_position"]].values):
        sys.stdout.write('  mutation %d out of %d (%.2f%%) \r' % (i+1, n, 100*(i+1)/n))
        sys.stdout.flush()
        dist = res_dist_to_interface (protein,
                                      mut_pos,
                                      get_protein_interface (interactome, protein),
                                      pdbChains,
                                      chainMap,
                                      pdbIDs = toPDB[protein] if protein in toPDB else None,
                                      prCov = prCov,
                                      chCov = chCov,
                                      allChainMap = allChainMap,
                                      singleChain = singleChain)
        distances.append(dist)
    print()
    return distances

def multiprotein_model_res_dist_to_interface (proteins,
                                              interactomeFile,
                                              pdbChainsFile,
                                              chainSeqFile,
                                              chainMapFile,
                                              chainStrucResFile,
                                              pdbDir,
                                              prCov = 0,
                                              chCov = 0,
                                              pdbMapFile = None,
                                              allChainMap = False,
                                              singleChain = True,
                                              downloadPDB = True,
                                              suppressWarnings = True):
    
    set_pdb_dir (pdbDir)
    allow_pdb_downloads (downloadPDB)
    suppress_pdb_warnings (suppressWarnings)
    interactome = read_single_interface_annotated_interactome (interactomeFile)
    with open(pdbChainsFile, 'rb') as f:
        pdbChains = pickle.load(f)
    if pdbMapFile:
        with open(pdbMapFile, 'rb') as f:
            toPDB = pickle.load(f)
    else:
        toPDB = {}
    load_pdbtools_chain_strucRes_labels (chainStrucResFile)
    load_pdbtools_chain_sequences (chainSeqFile)
    chainMap = read_list_table (chainMapFile, ["Qpos", "Spos"], [int, int], '\t')
    
    distances = []
    n = len(proteins)
    for i, protein in enumerate(proteins):
        sys.stdout.write('  protein %d out of %d (%.2f%%) \r' % (i+1, n, 100*(i+1)/n))
        sys.stdout.flush()
        dist = protein_model_res_dist_to_interface (protein,
                                                    get_protein_interface(interactome, protein),
                                                    pdbChains,
                                                    chainMap,
                                                    pdbIDs = toPDB[protein] if protein in toPDB else None,
                                                    prCov = prCov,
                                                    chCov = chCov,
                                                    allChainMap = allChainMap,
                                                    singleChain = singleChain)
        distances.append(dist)
    print()
    return distances

def protein_model_res_dist_to_interface (protein,
                                         interface,
                                         pdbChains,
                                         chainMap,
                                         pdbIDs = None,
                                         prCov = 0,
                                         chCov = 0,
                                         allChainMap = False,
                                         singleChain = True):
    
    if pdbIDs:
        pdbIDs = sorted(set(pdbIDs))
    else:
        mappingChains = chainMap.loc[chainMap["Query"] == protein, "Subject"].tolist()
        pdbIDs = sorted({x.split('_')[0] for x in mappingChains})
    
    strucMap, cov = best_structure_map (protein,
                                        pdbIDs,
                                        pdbChains,
                                        chainMap,
                                        prCov = prCov,
                                        chCov = chCov,
                                        allChainMap = allChainMap,
                                        singleChain = singleChain)
    
    if strucMap is not None:
        posMapping = multires_map_to_structure (interface, strucMap)
        chainIDs = strucMap["Subject"].values
        pdbID = chainIDs[0].split('_')[0]
        chainIDs = {id.split('_')[1] for id in chainIDs}
        interfaceMap = [(chainID, chainPos) for pos, pdbid, chainID, chainPos in posMapping]
        return multichain_res_min_dist_to_residues (pdbID, chainIDs, interfaceMap)
    else:
        return []
    
def multires_map_to_structure (resPos, strucMap):
    
    mapping = []
    for pos in resPos:
        posMap = position_map (pos, strucMap)
        try:
            ind = list(np.isnan(posMap)).index(False)
            chainID = strucMap.iloc[ind]["Subject"]
            pdbID, chainID = chainID.split('_')
            mapPos = posMap[ind]
            mapping.append( (pos, pdbID, chainID, mapPos) )
        except ValueError:
            pass
    return mapping

def multichain_res_min_dist_to_residues (pdbID, chainIDs, targetPositions):
    
    allChainDist = []
    for chainID in chainIDs:
        fullID = pdbID + '_' + chainID
        distances = chain_allres_min_dist_to_residues (pdbID, chainID, targetPositions)
        distances = [(fullID, pos, dist) for chainID, pos, dist in distances]
        allChainDist.extend(distances)
    return allChainDist

def chain_allres_min_dist_to_residues (pdbID, chainID, targetPositions):
    
    chainSeq = return_chain_sequence (pdbID + '_' + chainID)
    if chainSeq:
        sourcePositions = [(chainID, pos) for pos in np.arange(1, len(chainSeq) + 1)]
        return chain_multires_min_dist_to_residues (pdbID, sourcePositions, targetPositions)
    else:
        return []

def chain_multires_min_dist_to_residues (pdbID, sourcePositions, targetPositions):
    
    distances = []
    n = len(sourcePositions)
    for i, sourcePos in enumerate(sourcePositions):
        dist = chain_res_min_dist_to_residues (pdbID, sourcePos, targetPositions)
        if not np.isnan(dist):
            distances.append(sourcePos + (dist,))
    return distances

def chain_res_min_dist_to_residues (pdbID, sourcePos, targetPositions):
    
    distances = chain_res_dist_to_residues (pdbID, sourcePos, targetPositions)
    distances = distances [np.isnan(distances) == False]
    if distances.size > 0:
        return min(distances)
    else:
        return np.nan

def chain_res_dist_to_residues (pdbID, sourcePos, targetPositions):
    
    distances = []
    chain, pos = sourcePos
    for targetChain, targetPos in targetPositions:
        dist = res_distance_by_chainIDs (pdbID,
                                         chain,
                                         pos,
                                         targetChain,
                                         targetPos,
                                         pdbDir)
        distances.append(dist)
    return np.array(distances)

def multires_dist_to_interface (protein,
                                resPos,
                                interface,
                                pdbChains,
                                chainMap,
                                pdbIDs = None,
                                prCov = 0,
                                chCov = 0,
                                allChainMap = False,
                                singleChain = True):
    
    distances = multires_dist_to_residues (protein,
                                           resPos,
                                           interface,
                                           pdbChains,
                                           chainMap,
                                           pdbIDs = pdbIDs,
                                           prCov = prCov,
                                           chCov = chCov,
                                           allChainMap = allChainMap,
                                           singleChain = singleChain)
    minDist = []
    if distances:
        for chainID, chainPos, dist in distances:
            dist = dist[np.isnan(dist) == False]
            if dist.size > 0:
                minDist.append( (chainID, chainPos, min(dist)) )
            else:
                minDist.append( () )
    return minDist
    #return np.array( list( map( min, distances ) ) )

def res_dist_to_interface (protein,
                           resPos,
                           interface,
                           pdbchains,
                           chainMap,
                           pdbIDs = None,
                           prCov = 0,
                           chCov = 0,
                           allChainMap = False,
                           singleChain = True):
    
    distances = dist_to_residues (protein,
                                  resPos,
                                  interface,
                                  pdbchains,
                                  chainMap,
                                  pdbIDs = pdbIDs,
                                  prCov = prCov,
                                  chCov = chCov,
                                  allChainMap = allChainMap,
                                  singleChain = singleChain)
    if distances:
        chainID, chainPos, dist = distances
        dist = dist[np.isnan(dist) == False]
        if dist.size > 0:
            return chainID, chainPos, min(dist)
        else:
            return ()
    else:
        return ()

def multires_dist_to_residues (protein,
                               resPos,
                               otherPos,
                               pdbChains,
                               chainMap,
                               pdbIDs = None,
                               prCov = 0,
                               chCov = 0,
                               allChainMap = False,
                               singleChain = True):
    
    if pdbIDs:
        pdbIDs = sorted(set(pdbIDs))
    else:
        mappingChains = chainMap.loc[chainMap["Query"] == protein, "Subject"].tolist()
        pdbIDs = sorted({x.split('_')[0] for x in mappingChains})
    
    strucMap, cov = best_structure_map (protein,
                                        pdbIDs,
                                        pdbChains,
                                        chainMap,
                                        prCov = prCov,
                                        chCov = chCov,
                                        allChainMap = allChainMap,
                                        singleChain = singleChain)
    
    multiResDistances = []
    if strucMap is not None:
        n = len(resPos)
        for i, pos in enumerate(resPos):
            distances = dist_to_residues_mapped_struc (pos, otherPos, strucMap)
            if distances:
                multiResDistances.append(distances)
    return multiResDistances
#     return [ np.array( [ np.nan ] * len( otherPos ) ) ] * len( resPos )

def dist_to_residues (protein,
                      resPos,
                      otherPos,
                      pdbChains,
                      chainMap,
                      pdbIDs = None,
                      prCov = 0,
                      chCov = 0,
                      allChainMap = False,
                      singleChain = True):
    
    if pdbIDs:
        pdbIDs = sorted(set(pdbIDs))
    else:
        mappingChains = chainMap.loc[chainMap["Query"] == protein, "Subject"].tolist()
        pdbIDs = sorted({x.split('_')[0] for x in mappingChains})
    
    if len(otherPos) > 0:
        strucMap, cov = best_structure_map (protein,
                                            pdbIDs,
                                            pdbChains,
                                            chainMap,
                                            prCov = prCov,
                                            chCov = chCov,
                                            allChainMap = allChainMap,
                                            singleChain = singleChain)
        if strucMap is not None:
            return dist_to_residues_mapped_struc (resPos, otherPos, strucMap)
        #return np.array( [ np.nan ] * len( otherPos ) )
    else:
        warnings.warn("Requesting distance to empty list of residues")
        #return np.array( [ np.nan ] )
    return ()
    
def dist_to_residues_mapped_struc (resPos, otherPos, strucMap):
    
    mappedPos = position_map (resPos, strucMap)
    try:
        ind = list(np.isnan(mappedPos)).index(False)
        mapPos = mappedPos[ind]
        pdbid, hostChain = strucMap.iloc[ind]["Subject"].split('_')
        distances = []
        for pos in otherPos:
            mappedPos = position_map (pos, strucMap)
            try:
                ind = list(np.isnan(mappedPos)).index(False)
                _, otherChain = strucMap.iloc[ind]["Subject"].split('_')
                mapPos2 = mappedPos[ind]
                dist = res_distance_by_chainIDs (pdbid,
                                                 hostChain,
                                                 mapPos,
                                                 otherChain,
                                                 mapPos2,
                                                 pdbDir)
                distances.append(dist)
            except ValueError:
                distances.append(np.nan)
        return pdbid + '_' + hostChain, mapPos, np.array(distances)
    except ValueError:
        #return np.array( [ np.nan ] * len( otherPos ) )
        return ()

def multires_dist_to_center (protein,
                             resPos,
                             pdbChains,
                             chainMap,
                             pdbIDs = None,
                             prCov = 0,
                             chCov = 0,
                             allChainMap = False,
                             singleChain = True):
    
    if pdbIDs:
        pdbIDs = sorted(set(pdbIDs))
    else:
        mappingChains = chainMap.loc[chainMap["Query"] == protein, "Subject"].tolist()
        pdbIDs = sorted({x.split('_')[0] for x in mappingChains})
    
    strucMap, cov = best_structure_map (protein,
                                        pdbIDs,
                                        pdbChains,
                                        chainMap,
                                        prCov = prCov,
                                        chCov = chCov,
                                        allChainMap = allChainMap,
                                        singleChain = singleChain)
    distances = {}
    if strucMap is not None:
        n = len(resPos)
        for i, pos in enumerate(resPos):
            dist = dist_to_center_mapped_struc (pos, strucMap)
            if dist:
                distances[pos] = dist
    return distances

def res_dist_to_center (protein,
                        resPos,
                        pdbChains,
                        chainMap,
                        pdbIDs = None,
                        prCov = 0,
                        chCov = 0,
                        allChainMap = False,
                        singleChain = True):
    
    # global centerDist, distFile
#     k = protein + '_' + str(resPos)
#     if k in centerDist:
#         return centerDist[k]
#     else:
    if pdbIDs:
        pdbIDs = sorted(set(pdbIDs))
    else:
        mappingChains = chainMap.loc[chainMap["Query"] == protein, "Subject"].tolist()
        pdbIDs = sorted({x.split('_')[0] for x in mappingChains})

    strucMap, cov = best_structure_map (protein,
                                        pdbIDs,
                                        pdbChains,
                                        chainMap,
                                        prCov = prCov,
                                        chCov = chCov,
                                        allChainMap = allChainMap,
                                        singleChain = singleChain)
    if strucMap is not None:
        dist = dist_to_center_mapped_struc (resPos, strucMap)
    else:
        dist = ()
#         centerDist[k] = dist
#         if not (mutcount % 50): 
#             with open(distFile, 'wb') as fOut:
#                 pickle.dump(centerDist, fOut)
    return dist

def dist_to_center_mapped_struc (resPos, strucMap):
    
    mappedPos = position_map (resPos, strucMap)
    try:
        ind = list(np.isnan(mappedPos)).index(False)
        mapPos = mappedPos[ind]
        fullChainID = strucMap.iloc[ind]["Subject"]
        pdbid, hostChainID = fullChainID.split('_')
        allChainIDs = {id.split('_')[1] for id in strucMap["Subject"].values}
        center = multichain_geometric_center (pdbid, allChainIDs)
        dist = resPos_dist_to_point (pdbid, hostChainID, mapPos, center)
        if not np.isnan(dist):
            return fullChainID, mapPos, dist
        else:
            return ()
    except ValueError:
        return ()

def multichain_geometric_center (pdbid, chainIDs):
    
    struc = return_structure (pdbid, pdbDir)
    if struc:
        points = []
        for chainID in chainIDs:
            chain = struc[0][chainID]
            for residue in chain:
                for atom in residue:
                    points.append(atom.get_coord())
        return geometric_center (points)
    else:
        return ()

def geometric_center (points):
    
    if points:
        x, y, z = list(zip(* points))
        return np.mean(x), np.mean(y), np.mean(z)
    else:
        return ()
    
def resPos_dist_to_point (pdbid, chainID, resPos, point):
    
    struc = return_structure (pdbid, pdbDir)
    if struc:
        chain = struc[0][chainID]
        resID = return_chain_res_posToID (pdbid, chainID, resPos, pdbDir)
        if resID:
            return res_dist_to_point (chain[resID], point)
    return np.nan

def multichain_allres_dist_to_center (pdbid, chainIDs):
    
    center = multichain_geometric_center (pdbid, chainIDs)
    return multichain_allres_dist_to_point (pdbid, chainIDs, center)

def multichain_allres_dist_to_point (pdbid, chainIDs, point):
    
    distances = {}
    for chainID in chainIDs:
        fullID = pdbid + '_' + chainID
        distances[fullID] = chain_allres_dist_to_point (pdbid, chainID, point)
    return distances

def chain_allres_dist_to_point (pdbid, chainID, point):
    
    chainSeq = return_chain_sequence (pdbid + '_' + chainID)
    if chainSeq:
        chainSeqPositions = np.arange(1, len(chainSeq) + 1)
        return multiResPos_dist_to_point (pdbid, chainID, chainSeqPositions, point)
    else:
        return {}
    
def multiResPos_dist_to_point (pdbid, chainID, resPos, point):
    
    distances = {}
    struc = return_structure (pdbid, pdbDir)
    if struc:
        chain = struc[0][chainID]
        idmaps = return_chain_res_posToIDs (pdbid, chainID, pdbDir)
        n = len(resPos)
        for i, pos in enumerate(resPos):
            if pos in idmaps:
                resID = idmaps[pos]
                residue = chain[resID]
                distances[pos] = res_dist_to_point (residue, point)
    return distances

def multiRes_dist_to_point (residues, point):
    
    distance = []
    n = len(residues)
    for i, residue in enumerate(residues):
        distances.append(res_dist_to_point (residue, point))
    return distances
    
def res_dist_to_point (residue, point):
    
    mindist = np.inf
    px, py, pz = point
    for atom in residue:
        x, y, z = atom.get_coord()
        dist = (x - px)**2 + (y - py)**2 + (z - pz)**2
        if dist < mindist:
            mindist = dist
    return np.sqrt(mindist)

def res_distance_by_chainIDs (pdbid,
                              chainID1,
                              pos1,
                              chainID2,
                              pos2,
                              pdbDir):
    
    struc = return_structure (pdbid, pdbDir)
    if struc:
        resID1 = return_chain_res_posToID (pdbid, chainID1, pos1, pdbDir)
        resID2 = return_chain_res_posToID (pdbid, chainID2, pos2, pdbDir)
        if resID1 and resID2:
            res1 = struc[0][chainID1][resID1]
            res2 = struc[0][chainID2][resID2]
            return get_distance (res1, res2)
    return np.nan

def res_rsa (protein,
             resPos,
             pdbChains,
             chainMap,
             dsspDir = None,
             pdbIDs = None,
             prCov = 0,
             chCov = 0,
             allChainMap = False,
             singleChain = True,
             mapToProtein = True):
    
    if pdbIDs:
        pdbIDs = sorted(set(pdbIDs))
    else:
        mappingChains = chainMap.loc[chainMap["Query"] == protein, "Subject"].tolist()
        pdbIDs = sorted({x.split('_')[0] for x in mappingChains})
    
    strucMap, cov = best_structure_map (protein,
                                        pdbIDs,
                                        pdbChains,
                                        chainMap,
                                        prCov = prCov,
                                        chCov = chCov,
                                        allChainMap = allChainMap,
                                        singleChain = singleChain)
    
    if strucMap is not None:
        pdbid = strucMap.iloc[0]["Subject"].split('_')[0]
        rsa = res_rsa_mapped_struc (resPos,
                                    strucMap,
                                    allChains = pdbChains[pdbid],
                                    dsspDir = dsspDir,
                                    mapToProtein = mapToProtein)
        rsaReduced = {}
        if mapToProtein:
            for pos, val in rsa.items():
                chainID, chainPos, rsaInfo = val
                rsaReduced[pos] = chainID, chainPos, rsaInfo[3]
        else:
            for pos, rsaInfo in rsa.items():
                rsaReduced[pos] = rsaInfo[3]
        return rsaReduced
    else:
        return {}

def res_rsa_mapped_struc (resPos,
                          strucMap,
                          allChains = None,
                          dsspDir = None,
                          wholeStructureRSA = False,
                          mapToProtein = True):
    
    if wholeStructureRSA:
        return precalculated_res_rsa (resPos,
                                      strucMap,
                                      dsspDir,
                                      mapToProtein = mapToProtein)
    else:
        return calculate_res_rsa (resPos,
                                  strucMap,
                                  allChains = allChains,
                                  mapToProtein = mapToProtein)

def calculate_res_rsa (resPos,
                       strucMap,
                       allChains = None,
                       mapToProtein = True):
    
    pdbIDs = sorted({x.split('_')[0] for x in strucMap["Subject"].values})
    if len(pdbIDs) == 0:
        warnings.warn("Requesting residue RSA using empty structure map")
        return {}
    else:
        if len(pdbIDs) > 1:
            warnings.warn("Requesting residue RSA using structure maps of multiple PDB IDs. " +
                          "First encountered PDB ID will be used ...")
        pdbid = pdbIDs.pop()
        selMap = strucMap [strucMap["Subject"].apply(lambda x:
                            x.startswith(pdbid + '_'))].reset_index(drop=True)
        selChains = {id.split('_')[1] for id in selMap["Subject"].values}
        if selChains == allChains:
            rsa = calculate_structure_rsa (accDir, pdbDir, pdbid)
        else:
            rsa = calculate_structure_rsa (accDir, pdbDir, pdbid, selChains = selChains)
        if mapToProtein:
            rsa = rsa_map (rsa, resPos, selMap)
        else:
            rsa = chain_rsa_map (rsa, pdbid)
        return rsa

def precalculated_res_rsa (resPos,
                           strucMap,
                           dsspDir,
                           mapToProtein = True):
    
    rsa = {}
    pdbIDs = sorted({x.split('_')[0] for x in strucMap["Subject"].values})
    if len(pdbIDs) == 0:
        warnings.warn("Requesting residue RSA using empty structure map")
    else:
        if len(pdbIDs) > 1:
            warnings.warn("Requesting residue RSA using structure maps of multiple PDB IDs" +
                          "Using first encountered PDB ID ...")
        pdbid = pdbIDs.pop()
        dsspFile = dsspDir / (pdbid + '.dssp')
        if dsspFile.is_file():
            acc = read_dssp_file (dsspFile)
            rsa = acc_to_rsa (acc)
            if mapToProtein:
                rsa = rsa_map (rsa, resPos, strucMap)
            else:
                rsa = chain_rsa_map (rsa, pdbid)
    return rsa
    
def rsa_map (rsa, resPos, strucMap):
    
    rsaMap = {}
    chainResPos, mappedPositions = [], []
    for pos in resPos:
        mappedPos = position_map (pos, strucMap)
        try:
            ind = list(np.isnan(mappedPos)).index(False)
            pdbid, mapChain = strucMap.iloc[ind]["Subject"].split('_')
            mapPos = mappedPos[ind]
            chainResPos.append( (mapChain, mapPos) )
            mappedPositions.append(pos)
        except ValueError:
            pass
    if len(chainResPos) > 0:
        chainResIDs = return_multichain_res_posToIDs (pdbid, chainResPos, pdbDir)
        for (chainID, chainPos, resID), pos in zip(chainResIDs, mappedPositions):
            k = chainID, resID
            if k in rsa:
                rsaMap[pos] = pdbid + '_' + chainID, chainPos, rsa[k]
    return rsaMap

def chain_rsa_map (rsa, pdbid):
    
    rsaMap = {}
    for k in rsa.keys():
        chainID, resID = k
        pos = return_chain_res_IDToPos (pdbid, chainID, resID, pdbDir)
        rsaMap[ (pdbid + '_' + chainID, pos) ] = rsa[k]
    return rsaMap

def count_mutation_neighbors (mutations,
                              pdbChainsFile,
                              chainMapFile,
                              chainStrucResFile,
                              pdbDir,
                              nbDist = 5,
                              prCov = 0,
                              chCov = 0,
                              pdbMapFile = None,
                              allChainMap = False,
                              singleChain = True,
                              downloadPDB = True,
                              suppressWarnings = True):
    
    allow_pdb_downloads (downloadPDB)
    suppress_pdb_warnings (suppressWarnings)
    with open(pdbChainsFile, 'rb') as f:
        pdbChains = pickle.load(f)
    if pdbMapFile:
        with open(pdbMapFile, 'rb') as f:
            toPDB = pickle.load(f)
    else:
        toPDB = {}
    chainMap = read_list_table (chainMapFile, ["Qpos", "Spos"], [int, int], '\t')
    load_pdbtools_chain_strucRes_labels (chainStrucResFile)
    
    neighbors = []
    n = len(mutations)
    for i, (protein, mut_pos) in enumerate(mutations[["protein", "mut_position"]].values):
        sys.stdout.write('  mutation %d out of %d (%.2f%%) \r' % (i+1, n, 100*(i+1)/n))
        sys.stdout.flush()
        nb = count_res_neighbors (protein,
                                  mut_pos,
                                  pdbChains,
                                  chainMap,
                                  pdbDir,
                                  nbDist = nbDist,
                                  prCov = prCov,
                                  chCov = chCov,
                                  pdbIDs = toPDB[protein] if protein in toPDB else None,
                                  allChainMap = allChainMap,
                                  singleChain = singleChain)
        neighbors.append(nb)
    print()
    return neighbors

def count_multiprotein_res_neighbors (proteins,
                                      pdbChainsFile,
                                      chainMapFile,
                                      chainStrucResFile,
                                      pdbDir,
                                      nbDist = 5,
                                      prCov = 0,
                                      chCov = 0,
                                      pdbMapFile = None,
                                      allChainMap = False,
                                      singleChain = True,
                                      downloadPDB = True,
                                      suppressWarnings = True):
    
    allow_pdb_downloads (downloadPDB)
    suppress_pdb_warnings (suppressWarnings)
    with open(pdbChainsFile, 'rb') as f:
        pdbChains = pickle.load(f)
    if pdbMapFile:
        with open(pdbMapFile, 'rb') as f:
            toPDB = pickle.load(f)
    else:
        toPDB = {}
    chainMap = read_list_table (chainMapFile, ["Qpos", "Spos"], [int, int], '\t')
    load_pdbtools_chain_strucRes_labels (chainStrucResFile)
    
    neighbors = []
    n = len(proteins)
    for i, (protein, positions) in enumerate(proteins[["protein", "res_positions"]].values):
        sys.stdout.write('  protein %d out of %d (%.2f%%) \r' % (i+1, n, 100*(i+1)/n))
        sys.stdout.flush()
        nb = count_protein_res_neighbors (protein,
                                          positions,
                                          pdbChains,
                                          chainMap,
                                          pdbDir,
                                          nbDist = nbDist,
                                          prCov = prCov,
                                          chCov = chCov,
                                          pdbIDs = toPDB[protein] if protein in toPDB else None,
                                          allChainMap = allChainMap,
                                          singleChain = singleChain)
        neighbors.append(nb)
    print()
    return neighbors

def count_protein_res_neighbors (protein,
                                 resPos,
                                 pdbChains,
                                 chainMap,
                                 pdbDir,
                                 nbDist = 5,
                                 prCov = 0,
                                 chCov = 0,
                                 pdbIDs = None,
                                 allChainMap = False,
                                 singleChain = True):
    
    neighbors = []
    n = len(resPos)
    for i, pos in enumerate(resPos):
        nb = count_res_neighbors (protein,
                                  pos,
                                  pdbChains,
                                  chainMap,
                                  pdbDir,
                                  nbDist = nbDist,
                                  prCov = prCov,
                                  chCov = chCov,
                                  pdbIDs = pdbIDs,
                                  allChainMap = allChainMap,
                                  singleChain = singleChain)
        neighbors.append(nb)
    return np.array(neighbors)

def count_res_neighbors (protein,
                         resPos,
                         pdbChains,
                         chainMap,
                         pdbDir,
                         nbDist = 5,
                         prCov = 0,
                         chCov = 0,
                         pdbIDs = None,
                         allChainMap = False,
                         singleChain = True):
    
    if pdbIDs:
        pdbIDs = sorted(set(pdbIDs))
    else:
        mappingChains = chainMap.loc[chainMap["Query"] == protein, "Subject"].tolist()
        pdbIDs = sorted({x.split('_')[0] for x in mappingChains})
    
    strucMap, cov = best_structure_map (protein,
                                        pdbIDs,
                                        pdbChains,
                                        chainMap,
                                        prCov = prCov,
                                        chCov = chCov,
                                        allChainMap = allChainMap,
                                        singleChain = singleChain)
    if strucMap is not None:
        return count_res_neighbors_mapped_struc (resPos,
                                                 strucMap,
                                                 pdbDir,
                                                 nbDist = nbDist)
    return np.nan

def count_res_neighbors_mapped_struc (resPos, strucMap, pdbDir, nbDist = 5):
    
    mappedPos = position_map (resPos, strucMap)
    try:
        ind = list(np.isnan(mappedPos)).index(False)
        mapPos = mappedPos[ind]
        hostChain = strucMap.iloc[ind]["Subject"]
        allChains = set(strucMap["Subject"].values)
        pdbID = hostChain.split('_')[0]
        return count_neighbors_by_chainIDs (pdbID,
                                            mapPos,
                                            hostChain,
                                            allChains,
                                            pdbDir,
                                            nbDist = nbDist)
    except ValueError:
        return np.nan

def best_structure_map (protein,
                        pdbIDs,
                        pdbChains,
                        chainMap,
                        prCov = 0,
                        chCov = 0,
                        allChainMap = False,
                        singleChain = True):
    
    bestMap, bestCov = None, 0
    for id in pdbIDs:
        strucMap, mapped, cov = structure_map (protein,
                                               id,
                                               pdbChains,
                                               chainMap,
                                               prCov = prCov,
                                               chCov = chCov,
                                               allChainMap = allChainMap,
                                               singleChain = singleChain)
        if mapped and (cov > bestCov):
            bestMap, bestCov = strucMap, cov
    return bestMap, bestCov

def structure_map (protein,
                   pdbid,
                   pdbChains,
                   chainMap,
                   prCov = 0,
                   chCov = 0,
                   allChainMap = False,
                   singleChain = True):
    
    mapInd = ((chainMap["Query"] == protein) & 
                chainMap["Subject"].apply(lambda x: x.startswith(pdbid + '_')))
    strucMap = chainMap[mapInd].sort_values("Expect", axis=0, ascending=True)
    if (pdbid in pdbChains) and not strucMap.empty:
        if singleChain:
            strucMap = strucMap [(strucMap["Qcov"] >= prCov) & (strucMap["Scov"] >= chCov)]
            if not strucMap.empty:
                strucMap = strucMap.sort_values("Qcov", axis=0, ascending=False)
                return strucMap.iloc[:1], True, strucMap.iloc[0]["Qcov"]
        elif (not allChainMap) or (set(strucMap["Subject"].values) == pdbChains[pdbid]):
            strucMap = strucMap.drop_duplicates(subset = "Subject")
            if all(strucMap["Scov"] >= chCov):
                map_coord = list(zip(strucMap["Qstart"], strucMap["Qend"]))
                overlap_frac = []
                for coord in map_coord:
                    overlap_frac.append(sum(segment_overlaps(coord, map_coord) > 0))
                if max(overlap_frac) <= 1:
                    totalPrCov = strucMap["Qcov"].sum()
                    if totalPrCov >= prCov:
                        return strucMap, True, totalPrCov
    return strucMap, False, 0

def pos_structure_map (strucMap,
                       protein,
                       chainID,
                       pos,
                       firstMap = False):
    """Return protein-chain sequence alignments that include a specific protein 
        sequence position.

    Args:
        strucMap (DataFrame): table of protein-chain sequence alignments.
        protein (str): protein UniProt ID.
        chainID (str): PDB chain ID.
        pos (int): protein sequence position to be passed through alignments.
        firstMap (bool): if True, only first successful alignment returned.

    Returns:
        DataFrame

    """
    mappings = strucMap [(strucMap["Query"] == protein) & 
                         (strucMap["Subject"] == chainID)].reset_index(drop=True)
    mappings["posMaps"] = position_map (pos, mappings)
    mappings = mappings [np.isnan(mappings["posMaps"]) == False]
    if firstMap and not mappings.empty:
        return mappings.iloc[0]
    else:
        return mappings if not mappings.empty else None

def position_map (resPos, strucMap):
    """Map a sequence position to positions through multiple sequence alignments.

    Args:
        resPos (int): reference sequence position to be mapped.
        strucMap (DataFrame): table of protein-chain sequence alignments.

    Returns:
        array

    """
    posMap = []
    for Qpos, Spos in strucMap[["Qpos", "Spos"]].values:
        if resPos in Qpos:
            posMap.append(Spos[Qpos.index(resPos)])
        else:
            posMap.append(np.nan)
    return np.array(posMap)

def write_mutation_structure_maps (mutations,
                                   chainSeqFile,
                                   chainStrucResFile,
                                   pdbDir,
                                   outPath,
                                   downloadPDB = True,
                                   suppressWarnings = True):
    
    set_pdb_dir (pdbDir)
    allow_pdb_downloads (downloadPDB)
    suppress_pdb_warnings (suppressWarnings)
    load_pdbtools_chain_sequences (chainSeqFile)
    load_pdbtools_chain_strucRes_labels (chainStrucResFile)
    with io.open(outPath, "w") as fout:
        fout.write('\t'.join(['protein',
                              'partner',
                              'protein_pos',
                              'pdb_id',
                              'chain_id',
                              'chain_pos',
                              'chain_mutation',
                              'partner_chain']) + '\n')
        n = len(mutations)
        for i, (protein, proteinPos, mutRes, chainID, chainPos) in enumerate(mutations):
            sys.stdout.write('  mutation %d out of %d (%.2f%%) \r' % (i+1, n, 100*(i+1)/n))
            sys.stdout.flush()
            seq = return_chain_sequence (chainID)
            if seq:
                chainWT = seq[chainPos - 1]
                if chainWT != mutRes:
                    pdbid, chainID = chainID.split('_')
                    resID = return_chain_res_posToID (pdbid, chainID, chainPos, pdbDir)
                    if resID:
                        _, resNum, _ = resID
                        fout.write('\t'.join([protein,
                                              '-',
                                              str(proteinPos),
                                              pdbid,
                                              chainID,
                                              str(chainPos),
                                              ''.join([chainWT, chainID, str(resNum), mutRes]),
                                              '-']) + '\n')
        print()

def remove_mutation_RSA_bias (mutations):
    
    edgetic = mutations["edgotype"] == 'edgetic'
    nonEdgetic = mutations["edgotype"] == 'non-edgetic'
    numNonEdgetic = sum(nonEdgetic)
    validRSA = np.isnan(mutations["RSA"]) == False
    numNonEdgeticRSA = sum(nonEdgetic & validRSA)
    mappingFactor = numNonEdgeticRSA / numNonEdgetic
    
    mutations = mutations [edgetic | validRSA].reset_index(drop=True)
    
    edgetic = mutations["edgotype"] == 'edgetic'
    numEdgetic = sum(edgetic)
    numEdgeticsToRemove = numEdgetic - int(np.round(numEdgetic * mappingFactor))
    
    nans = np.isnan(mutations["RSA"])
    numEdgeticNaN = sum(edgetic & nans)
    
    if numEdgeticsToRemove <= numEdgeticNaN:
        edgeticInd = [i for i, v in enumerate(edgetic & nans) if v]
    else:
        mutations = mutations[nans == False].reset_index(drop=True)
        numEdgeticsToRemove -= numEdgeticNaN
        edgeticInd = [i for i, v in enumerate(mutations["edgotype"] == 'edgetic') if v]
    
    sampleIndex = np.random.choice (edgeticInd, size = numEdgeticsToRemove, replace = False)
    toKeep = [i not in sampleIndex for i in range(len(mutations))]
    mutations = mutations[toKeep].reset_index(drop=True)
    return mutations
