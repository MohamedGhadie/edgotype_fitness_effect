#----------------------------------------------------------------------------------------
# Compare mutation residue RSA with other residues.
#----------------------------------------------------------------------------------------

import os
import pickle
import pandas as pd
import numpy as np
from pathlib import Path
from threeD_structure_tools import produce_protein_model_RSA
from stat_tools import t_test, sderror, round_data, cont_pdf
from plot_tools import bar_plot, multi_bar_plot

def main():
    
    # reference interactome name: HI-II-14, HuRI, IntAct
    interactome_name = 'IntAct'
    
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
    
    # show figures
    showFigs = False
    
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
    
    # directory of PDB structures
    pdbDir = Path('../../pdb_files')
    
    if model_method is 'model_based':
        modelDir = modellingDir / 'protein_models'
    else:
        modelDir = pdbDir
    
    # directory for calculated solvent accessibility files
    accDir = modellingDir / 'res_acc'
    
    # figure directory
    figDir = Path('../figures') / interactome_name / model_method / edgetic_method
    
    if edgetic_method is 'physics':
        edgeticDir = edgeticDir / (edgetic_ddg + '_edgetics')
        figDir = figDir / (edgetic_ddg + '_edgetics')
    
    # input data files
    proteinSeqFile = procDir / 'human_reference_sequences.pkl'
    maxAccFile = procDir / 'empirical_maxAcc_99_99.pkl'
    modelChainsFile = modellingDir / 'protein_model_chains.pkl'
    chainSeqFile = modellingDir / 'protein_chain_sequences.pkl'
    chainStrucResFile = modellingDir / 'protein_chain_strucRes.pkl'
    proteinModelFile = modellingDir / 'single_chain_map_per_protein.txt'
    naturalMutationsFile = edgeticDir / 'nondisease_mutation_RSA.txt'
    diseaseMutationsFile = edgeticDir / 'disease_mutation_RSA.txt'
    
    # output data files
    proteinModelRSAFile = edgeticDir / 'protein_model_RSA_99_99.pkl'
    
    # create output directories if not existing
    if not edgeticDir.exists():
        os.makedirs(edgeticDir)
    if not figDir.exists():
        os.makedirs(figDir)
    
    #------------------------------------------------------------------------------------
    # Load mutations and RSA values
    #------------------------------------------------------------------------------------    
    
    print('Reading all mutations')
    naturalMutations = pd.read_table (naturalMutationsFile, sep='\t')
    diseaseMutations = pd.read_table (diseaseMutationsFile, sep='\t')
    
    print('Total number of mutations:')
    print('non-disease mutations: %d' % len(naturalMutations))
    print('disease mutations: %d' % len(diseaseMutations))
    print()
    
    natMutationProteins = set(naturalMutations["protein"].tolist())
    disMutationProteins = set(diseaseMutations["protein"].tolist())
    mutationProteins = natMutationProteins | disMutationProteins
    
    print('Total number of proteins carrying mutations:')
    print('non-disease mutations: %d' % len(natMutationProteins))
    print('disease mutations: %d' % len(disMutationProteins))
    print('all mutations: %d' % len(mutationProteins))
    print()
    
    #------------------------------------------------------------------------------------
    # Calculate RSA for all residues in structural models carrying mutations
    #------------------------------------------------------------------------------------ 

    if not proteinModelRSAFile.is_file():
        with open(proteinSeqFile, 'rb') as f:
            prSeq = pickle.load(f)
        prLen = {}
        for p in mutationProteins:
            if p in prSeq:
                prLen[p] = len(prSeq[p])
        
        print('Calculating RSA for all model residues')
        produce_protein_model_RSA (prLen,
                                   modelChainsFile,
                                   chainSeqFile,
                                   proteinModelFile,
                                   chainStrucResFile,
                                   modelDir,
                                   accDir,
                                   proteinModelRSAFile,
                                   maxAccFile = maxAccFile,
                                   mapToProtein = False,
                                   downloadPDB = allow_pdb_downloads,
                                   suppressWarnings = suppress_pdb_warnings)
    
    with open(proteinModelRSAFile, 'rb') as f:
        modelRSA = pickle.load(f)
    
    allResidues = pd.DataFrame({"protein" : sorted(mutationProteins)})
    allResidues["RSA"] = [{}] * len(allResidues)
    allResidues["RSA"] = allResidues["protein"].apply(lambda x: modelRSA[x]
                                                                if x in modelRSA
                                                                else {})
    allResidues = allResidues [allResidues["RSA"].apply(bool)]
    
    allResRSA = []
    for rsa in allResidues["RSA"].values:
        allResRSA.extend(rsa.values())
    allResRSA = [rsa for rsa in allResRSA if not np.isnan(rsa)]
    numAllRes = len(allResRSA)
    
    if numAllRes > 0:
        print('Average RSA for all residues: %f Å (SE = %g, n = %d)' 
                % (np.mean(allResRSA), sderror(allResRSA), numAllRes))
    else:
        print('all residues: NA')
    print()
    
    #------------------------------------------------------------------------------------
    # print results
    #------------------------------------------------------------------------------------
    
    natMutRSA = [rsa for rsa in naturalMutations["RSA"].values if not np.isnan(rsa)]
    disMutRSA = [rsa for rsa in diseaseMutations["RSA"].values if not np.isnan(rsa)]
    
    # total number of mutations with calculated RSA
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
        print('Difference between non-disease and disease mutations:')
        t_test (natMutRSA, disMutRSA)
    if (numNatural > 1) and (numAllRes > 1):
        print('Difference between non-disease mutations and all residues:')
        t_test (natMutRSA, allResRSA)
    if (numDisease > 1) and (numAllRes > 1):
        print('Difference between disease mutations and all residues:')
        t_test (disMutRSA, allResRSA)
    print() 
    
    data = [np.mean(disMutRSA), np.mean(natMutRSA), np.mean(allResRSA)]
    maxY = 0.1 * np.ceil(max(data) / 0.1)
    bar_plot (data,
              error = [sderror(disMutRSA), sderror(natMutRSA), sderror(allResRSA)],
              xlabels = ['Pathogenic\nmutations', 'Non-pathogenic\nmutations', 'All\nresidues'],
              ylabel = 'RSA',
              ylabels = list(np.around(np.linspace(0, maxY, int(maxY*10+1)),1)),
              colors = ['magenta', 'turquoise', 'dodgerblue'],
              edgecolor = 'black',
              ewidth = 2,
              barwidth = 0.6,
              fontsize = 16,
              capsize = 10,
              msize = 20,
              ylim = [0, maxY],
              show = showFigs,
              figdir = figDir,
              figname = 'avg_rsa')
    
    alldata = []
    for rsa in [disMutRSA, natMutRSA, allResRSA]:
        rounded = round_data (rsa, w = 0.1)
        dens, x = cont_pdf (rounded, minVal = 0, maxVal = 1, w = 0.1)
        alldata.append(100 * dens)
    
    barwidth = 0.4
    tickshift = barwidth / 2.
    xticks = [i + tickshift for i in np.arange(1, len(alldata[0]) + 1, 2)]
    xticklabels = [round(f, 1) for f in np.arange(0, 1.1, 0.2)]
    diseaseData, nondiseaseData = alldata[0] - alldata[2], alldata[1] - alldata[2]
    minY = min([i for ls in (diseaseData, nondiseaseData) for i in ls])
    minY = 5 * np.floor(minY / 5) if minY < 0 else 0
    maxY = max([i for ls in (diseaseData, nondiseaseData) for i in ls])
    maxY = 5 * np.ceil(maxY / 5) if maxY > 0 else 0
    multi_bar_plot ([diseaseData, nondiseaseData],
                    xlabel = 'RSA',
                    ylabel = 'Difference in fraction\nfrom all residues (%)',
                    xlabels = xticklabels,
                    colors = ['magenta', 'turquoise'],
                    hatches = ['/', '..'],
                    barwidth = barwidth,
                    fontsize = 22,
                    ylim = [minY, maxY],
                    xticks = xticks,
                    opacity = None,
                    leg = ('Pathogenic mutations', 'Non-pathogenic mutations'),
                    overlap = False,
                    show = showFigs,
                    figdir = figDir,
                    figname = 'rsa_bar_distribution_difference')

if __name__ == "__main__":
    main()
