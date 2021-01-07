#----------------------------------------------------------------------------------------
# Calculate mutation residue RSA.
#----------------------------------------------------------------------------------------

import os
import pickle
import pandas as pd
import numpy as np
from pathlib import Path
from stat_tools import t_test, sderror
from threeD_structure_tools import produce_protein_model_RSA

def main():
    
    # reference interactome name: HI-II-14, HuRI, IntAct
    interactome_name = 'IntAct'
    
    # homology modelling method used to create structural models
    # options: template_based, model_based
    model_method = 'model_based_6'
    
    # method used to perform calculations; 'geometry' or 'physics'
    edgetic_method = 'physics'
    
    # method that was used to calculate edgetic mutation ∆∆G
    # options: bindprofx, foldx
    edgetic_ddg = 'foldx'
    
    # allow downloading of PDB structures
    allow_pdb_downloads = False
    
    # suppress PDB warnings when loading structures
    suppress_pdb_warnings = True
    
    # parent directory of all data files
    dataDir = Path('../data')
    
    # parent directory of all processed data files
    procDir = dataDir / 'processed'
    
    # directory of processed data files specific to interactome
    interactomeDir = procDir / interactome_name
    
    # directory of processed model-related data files specific to interactome
    modellingDir = interactomeDir / model_method
    
    # directory of edgetic mutation calculation method
    edgeticDir = modellingDir / edgetic_method
    
    if edgetic_method is 'physics':
        edgeticDir = edgeticDir / (edgetic_ddg + '_edgetics')
        
    # directory of PDB structures
    pdbDir = Path('../../pdb_files')
    
    if model_method is 'model_based_6':
        modelDir = interactomeDir / 'model_based' / 'protein_models'
    else:
        modelDir = pdbDir
    
    # directory for calculated solvent accessibility files
    accDir = modellingDir / 'res_acc'
    
    # input data files
    proteinSeqFile = procDir / 'human_reference_sequences.pkl'
    maxAccFile = procDir / 'empirical_maxAcc_99_99.pkl'
    modelChainsFile = modellingDir / 'protein_model_chains.pkl'
    chainSeqFile = modellingDir / 'protein_chain_sequences.pkl'
    chainStrucResFile = modellingDir / 'protein_chain_strucRes.pkl'
    proteinModelFile = modellingDir / 'single_chain_map_per_protein.txt'
    naturalMutationsFile = edgeticDir / 'nondisease_mutation_edgetics.txt'
    diseaseMutationsFile = edgeticDir / 'disease_mutation_edgetics.txt'
    
    # output data files
    proteinRSAFile = edgeticDir / 'protein_RSA_99_99.pkl'
    natMutRSAFile = edgeticDir / 'nondisease_mutation_RSA.txt'
    disMutRSAFile = edgeticDir / 'disease_mutation_RSA.txt'
    
    # create output directories if not existing
    if not edgeticDir.exists():
        os.makedirs(edgeticDir)
    
    #------------------------------------------------------------------------------------
    # load mutations
    #------------------------------------------------------------------------------------
    
    print('Reading mutations')
    naturalMutations = pd.read_table (naturalMutationsFile, sep='\t')
    diseaseMutations = pd.read_table (diseaseMutationsFile, sep='\t')
    
    naturalMutations = naturalMutations [(naturalMutations["edgotype"] == 'edgetic') |
                                         (naturalMutations["edgotype"] == 'non-edgetic')].reset_index(drop=True)
    diseaseMutations = diseaseMutations [(diseaseMutations["edgotype"] == 'edgetic') |
                                         (diseaseMutations["edgotype"] == 'non-edgetic')].reset_index(drop=True)
    
    print('Number of mutations:')
    print('non-disease mutations: %d' % len(naturalMutations))
    print('disease mutations: %d' % len(diseaseMutations))
    print()
    
    natMutationProteins = set(naturalMutations["protein"].tolist())
    disMutationProteins = set(diseaseMutations["protein"].tolist())
    mutationProteins = natMutationProteins | disMutationProteins
    
    print('Number of proteins carrying mutations:')
    print('non-disease mutations: %d' % len(natMutationProteins))
    print('disease mutations: %d' % len(disMutationProteins))
    print('all mutations: %d' % len(mutationProteins))
    print()
    
    #------------------------------------------------------------------------------------
    # Calculate mutation RSA
    #------------------------------------------------------------------------------------    
    
    if not proteinRSAFile.is_file():
        with open(proteinSeqFile, 'rb') as f:
            prSeq = pickle.load(f)
        prLen = {}
        for p in mutationProteins:
            if p in prSeq:
                prLen[p] = len(prSeq[p])
        
        print('producing protein RSA dictionary')
        produce_protein_model_RSA (prLen,
                                   modelChainsFile,
                                   chainSeqFile,
                                   proteinModelFile,
                                   chainStrucResFile,
                                   modelDir,
                                   accDir,
                                   proteinRSAFile,
                                   maxAccFile = maxAccFile,
                                   mapToProtein = True,
                                   downloadPDB = allow_pdb_downloads,
                                   suppressWarnings = suppress_pdb_warnings)
        
    with open(proteinRSAFile, 'rb') as f:
        proteinRSA = pickle.load(f)
    
    naturalMutations["RSA"] = [()] * len(naturalMutations)
    for i, row in naturalMutations.iterrows():
        if row.protein in proteinRSA:
            if row.mut_position in proteinRSA[row.protein]:
                naturalMutations.at[i, "RSA"] = proteinRSA[row.protein][row.mut_position]
            else:
                naturalMutations.at[i, "RSA"] = ('-', -1, np.nan)
    
    diseaseMutations["RSA"] = [()] * len(diseaseMutations)
    for i, row in diseaseMutations.iterrows():
        if row.protein in proteinRSA:
            if row.mut_position in proteinRSA[row.protein]:
                diseaseMutations.at[i, "RSA"] = proteinRSA[row.protein][row.mut_position]
            else:
                diseaseMutations.at[i, "RSA"] = ('-', -1, np.nan)
    
    naturalMutations = naturalMutations [naturalMutations["RSA"].apply(bool)]
    natMutModelChain, natMutChainPos, natMutRSA = list(zip(* naturalMutations["RSA"].values))
    diseaseMutations = diseaseMutations [diseaseMutations["RSA"].apply(bool)]
    disMutModelChain, disMutChainPos, disMutRSA = list(zip(* diseaseMutations["RSA"].values))
    
    naturalMutations = naturalMutations.assign (model_chain = natMutModelChain,
                                                chain_pos = natMutChainPos,
                                                RSA = natMutRSA)
    
    diseaseMutations = diseaseMutations.assign (model_chain = disMutModelChain,
                                                chain_pos = disMutChainPos,
                                                RSA = disMutRSA)
    
    naturalMutations.to_csv (natMutRSAFile, index=False, sep='\t')
    diseaseMutations.to_csv (disMutRSAFile, index=False, sep='\t')
    
    natMutRSA = [rsa for rsa in natMutRSA if not np.isnan(rsa)]
    disMutRSA = [rsa for rsa in disMutRSA if not np.isnan(rsa)]
    
    # total number of mutation residues with calculated distance to center
    numNatural, numDisease = len(natMutRSA), len(disMutRSA)
    
    # print results
    if numNatural > 0:
        print('Average RSA for non-disease mutations: %f Å (SE = %g, n = %d)' 
                % (np.mean(natMutRSA), sderror(natMutRSA), numNatural))
    else:
        print('non-disease mutations: NA')
    
    if numDisease > 0:
        print('Average RSA for disease mutations: %f Å (SE = %g, n = %d)' 
                % (np.mean(disMutRSA), sderror(disMutRSA), numDisease))
    else:
        print('disease mutations: NA')
    
    if (numNatural > 1) and (numDisease > 1):
        t_test(natMutRSA, disMutRSA)
    print()

if __name__ == "__main__":
    main()
