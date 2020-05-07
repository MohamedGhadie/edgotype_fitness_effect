import os
import pickle
import numpy as np
import pandas as pd
from pathlib import Path
from threeD_structure_tools import multiprotein_model_res_dist_to_center
from stat_tools import t_test, sderror, round_data, cont_pdf
from plot_tools import bar_plot, multi_bar_plot

def main():
    
    # reference interactome name
    # options: HI-II-14, HuRI, IntAct
    interactome_name = 'HuRI'
    
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
    #pdbDir = Path('/Volumes/MG_Samsung/pdb_files')
    pdbDir = Path('../../pdb_files')
    
    if model_method is 'model_based':
        modelDir = modellingDir / 'protein_models'
    else:
        modelDir = pdbDir
    
    # figure directory
    figDir = Path('../figures') / interactome_name / model_method / edgetic_method
    
    if edgetic_method is 'physics':
        edgeticDir = edgeticDir / (edgetic_ddg + '_edgetics')
        figDir = figDir / (edgetic_ddg + '_edgetics')
        
    # input data files
    modelChainsFile = modellingDir / 'protein_model_chains.pkl'
    chainSeqFile = modellingDir / 'protein_chain_sequences.pkl'
    proteinModelFile = modellingDir / 'single_chain_map_per_protein.txt'
    chainStrucResFile = modellingDir / 'protein_chain_strucRes.pkl'
    naturalMutationsFile = edgeticDir / 'nondisease_mutation_dist_to_center.txt'
    diseaseMutationsFile = edgeticDir / 'disease_mutation_dist_to_center.txt'
    
    # output data files
    allResDistFile = edgeticDir / 'protein_model_res_dist_to_center.pkl'
    
    # create output directories if not existing
    if not edgeticDir.exists():
        os.makedirs(edgeticDir)
    if not figDir.exists():
        os.makedirs(figDir)
        
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
    
    natMutDist = naturalMutations["dist_to_center"].tolist()
    disMutDist = diseaseMutations["dist_to_center"].tolist()
    
    # total number of mutation residues with calculated distance to center
    numNatural, numDisease = len(natMutDist), len(disMutDist)
    
    # print results
    print()
    if numNatural > 0:
        print('Average distance to protein center for non-disease mutations: %f Å (SE = %g, n = %d)' 
                % (np.mean(natMutDist), sderror(natMutDist), numNatural))
    else:
        print('non-disease mutations: NA')
    
    if numDisease > 0:
        print('Average distance to protein center for disease mutations: %f Å (SE = %g, n = %d)'
                % (np.mean(disMutDist), sderror(disMutDist), numDisease))
    else:
        print('disease mutations: NA')
    
    if (numNatural > 1) and (numDisease > 1):
        t_test (natMutDist, disMutDist)
    print()
    
    #------------------------------------------------------------------------------------
    # calculate distance to structure center for all residues
    #------------------------------------------------------------------------------------
    
    if allResDistFile.is_file():
        with open(allResDistFile, 'rb') as f:
            allResidues = pickle.load(f)
    else:
        print('Calculating distance to center for all residues')
        mutationProteins = sorted(set(naturalMutations["protein"]) | set(diseaseMutations["protein"]))
        allResidues = pd.DataFrame({"protein" : mutationProteins})
        allResidues["dist_to_center"] = multiprotein_model_res_dist_to_center(
                                                            allResidues["protein"].tolist(),
                                                            modelChainsFile,
                                                            chainSeqFile,
                                                            proteinModelFile,
                                                            chainStrucResFile,
                                                            modelDir,
                                                            downloadPDB = allow_pdb_downloads,
                                                            suppressWarnings = suppress_pdb_warnings)
    
        with open(allResDistFile, 'wb') as fOut:
            pickle.dump(allResidues, fOut)
        print('Distances to center for all protein model residues written to file %s' % str(allResDistFile))
    
    allResidues = allResidues [allResidues["dist_to_center"].apply(bool)]
    
    allResDist = []
    for d in allResidues["dist_to_center"].values:
        allResDist.extend([dist for chainID, chainPos, dist in d if not np.isnan(dist)])
    numAllRes = len(allResDist)
    
    # print results  
    if numAllRes > 0:
        print()
        print('Average distance to protein center for all protein residues: %.1f Å (SE = %g, n = %d)' 
                % (np.mean(allResDist), sderror(allResDist), numAllRes))
    else:
        print('all protein residues: NA')
    
    # statistical significance
    if (numNatural > 1) and (numAllRes > 1):
        print('Difference between non-disease mutations and all residues:')
        t_test (natMutDist, allResDist)
    if (numDisease > 1) and (numAllRes > 1):
        print('Difference between disease mutations and all residues:')
        t_test (disMutDist, allResDist)
    print()
   
    #------------------------------------------------------------------------------------
    # plot results
    #------------------------------------------------------------------------------------ 
    
    means = [np.mean(disMutDist), np.mean(natMutDist), np.mean(allResDist)]
    maxY = 5 * np.ceil(max(means) / 5)
    bar_plot (means,
              error = [sderror(disMutDist), sderror(natMutDist), sderror(allResDist)],
              ylabel = 'Distance to protein center (Å)',
              xlabels = ('Disease\nmutations', 'Non-disease\nmutations', 'All\nresidues'),
              #ylabels = [0, 0.1, 0.2, 0.3],
              colors = ['magenta', 'dodgerblue', 'turquoise'],
              edgecolor = 'black',
              ewidth = 2,
              barwidth = 0.6,
              fontsize = 20,
              capsize = 10,
              msize = 20,
              ylim = [0, maxY],
              show = showFigs,
              figdir = figDir,
              figname = 'Avg_distance_to_center')
    
    maxDist = 100
    w = 10
    alldata = []
    for dist in [disMutDist, natMutDist, allResDist]:
        rounded = round_data (dist, w = w, maxVal = maxDist)
        dens, x = cont_pdf (rounded, minVal = 0, maxVal = maxDist, w = w)
        alldata.append(100 * dens)
    
    maxY = max([i for ls in alldata for i in ls])
    maxY = 5 * np.ceil(maxY / 5)
    multi_bar_plot (alldata,
                    xlabel = 'Distance to protein center (Å)',
                    ylabel = 'Fraction (%)',
                    xlabels = list(map(str, x[:-1])) + [ '≥' + str(x[-1]) ],
                    colors = ['r', 'b', 'g'],
                    barwidth = 0.2,
                    ylim = [0, maxY],
                    leg = ('Disease mutations', 'Non-disease mutations', 'All residues'),
                    show = showFigs,
                    figdir = figDir,
                    figname = 'dist_to_center_bar_distribution')
    
    barwidth = 0.4
    tickshift = barwidth / 2.
    xticks = [i + tickshift for i in np.arange(1, len(alldata[0]) + 1, 2)]
    xticklabels = list(map(str, np.arange(0, maxDist, 20))) + ['≥' + str(maxDist)]
    diseaseData, nondiseaseData = alldata[0] - alldata[2], alldata[1] - alldata[2]
    minY = min([i for ls in (diseaseData, nondiseaseData) for i in ls])
    minY = 5 * np.floor(minY / 5) if minY < 0 else 0
    maxY = max([i for ls in (diseaseData, nondiseaseData) for i in ls])
    maxY = 5 * np.ceil(maxY / 5) if maxY > 0 else 0
    multi_bar_plot ([diseaseData, nondiseaseData],
                    xlabel = 'Distance to protein center (Å)',
                    ylabel = 'Difference in fraction from average residue (%)',
                    xlabels = xticklabels,
                    colors = ['magenta', 'dodgerblue'],
                    barwidth = barwidth,
                    #ylim = [-12, 40],
                    xticks = xticks,
                    opacity = None,
                    #leg = ('Disease mutations', 'Non-disease mutations'),
                    overlap = False,
                    show = showFigs,
                    figdir = figDir,
                    figname = 'dist_to_center_bar_distribution_difference')

if __name__ == "__main__":
    main()
