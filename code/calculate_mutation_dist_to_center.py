import os
import numpy as np
import pandas as pd
from pathlib import Path
from threeD_structure_tools import mutation_dist_to_center
from stat_tools import t_test, sderror

def main():
    
    # reference interactome name
    # options: HI-II-14, IntAct
    interactome_name = 'HI-II-14'
    
    # homology modelling method used to create structural models
    # options: template_based, model_based
    model_method = 'model_based'
    
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
    pdbDir = Path('/Volumes/MG_Samsung/pdb_files')
    #pdbDir = Path('../../pdb_files')
    
    if model_method is 'model_based':
        modelDir = modellingDir / 'protein_models'
    else:
        modelDir = pdbDir
        
    # input data files
    modelChainsFile = modellingDir / 'protein_model_chains.pkl'
    chainSeqFile = modellingDir / 'protein_chain_sequences.pkl'
    proteinModelFile = modellingDir / 'single_chain_map_per_protein.txt'
    chainStrucResFile = modellingDir / 'protein_chain_strucRes.pkl'
    naturalMutationsFile = edgeticDir / 'nondisease_mutation_edgetics.txt'
    diseaseMutationsFile = edgeticDir / 'disease_mutation_edgetics.txt'

    # output data files
    natMutCenterDistFile = edgeticDir / 'nondisease_mutation_dist_to_center.txt'
    disMutCenterDistFile = edgeticDir / 'disease_mutation_dist_to_center.txt'
    
    # create output directories if not existing
    if not edgeticDir.exists():
        os.makedirs(edgeticDir)
        
    #------------------------------------------------------------------------------------
    # load mutations
    #------------------------------------------------------------------------------------
    
    print('Reading mutations')
    naturalMutations = pd.read_table (naturalMutationsFile, sep='\t')
    diseaseMutations = pd.read_table (diseaseMutationsFile, sep='\t')
    
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
    # calculate mutation distance to structure center
    #------------------------------------------------------------------------------------
    
    print('Calculating distance to protein center for non-disease mutations')
    naturalMutations["dist_to_center"] = mutation_dist_to_center (naturalMutations,
                                                                  modelChainsFile,
                                                                  chainSeqFile,
                                                                  proteinModelFile,
                                                                  chainStrucResFile,
                                                                  modelDir,
                                                                  downloadPDB = allow_pdb_downloads,
                                                                  suppressWarnings = suppress_pdb_warnings)
    print('Calculating distance to protein center for disease mutations')
    diseaseMutations["dist_to_center"] = mutation_dist_to_center (diseaseMutations,
                                                                  modelChainsFile,
                                                                  chainSeqFile,
                                                                  proteinModelFile,
                                                                  chainStrucResFile,
                                                                  modelDir,
                                                                  downloadPDB = allow_pdb_downloads,
                                                                  suppressWarnings = suppress_pdb_warnings)
    
    naturalMutations = naturalMutations [naturalMutations["dist_to_center"].apply(bool)]
    natMutModelChain, natMutChainPos, natMutDist = list(zip(* naturalMutations["dist_to_center"].values))
    diseaseMutations = diseaseMutations [diseaseMutations["dist_to_center"].apply(bool)]
    disMutModelChain, disMutChainPos, disMutDist = list(zip(* diseaseMutations["dist_to_center"].values))
    
    naturalMutations = naturalMutations.assign (model_chain = natMutModelChain,
                                                chain_pos = natMutChainPos,
                                                dist_to_center = natMutDist)
    
    diseaseMutations = diseaseMutations.assign (model_chain = disMutModelChain,
                                                chain_pos = disMutChainPos,
                                                dist_to_center = disMutDist)
    
    naturalMutations.to_csv (natMutCenterDistFile, index=False, sep='\t')
    diseaseMutations.to_csv (disMutCenterDistFile, index=False, sep='\t')
    
    # total number of mutation residues with calculated distance to center
    numNatural, numDisease = len(natMutDist), len(disMutDist)
    
    # print results
    print()
    if numNatural > 0:
        print('Average distance to protein center for non-disease mutations: %.1f Å (SE = %g, n = %d)' 
                % (np.mean(natMutDist), sderror(natMutDist), numNatural))
    else:
        print('non-disease mutations: NA')
    
    if numDisease > 0:
        print('Average distance to protein center for disease mutations: %.1f Å (SE = %g, n = %d)' 
                % (np.mean(disMutDist), sderror(disMutDist), numDisease))
    else:
        print('disease mutations: NA')
    
    if (numNatural > 1) and (numDisease > 1):
        t_test(natMutDist, disMutDist)
    print()

if __name__ == "__main__":
    main()
