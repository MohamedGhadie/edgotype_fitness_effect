import pandas as pd
import numpy as np
from pathlib import Path
from energy_tools import read_protein_mutation_ddg
from stat_tools import sderror
from plot_tools import multi_bar_plot, curve_plot

def main():
    
    # reference interactome name
    # options: HI-II-14, HuRI, IntAct
    interactome_name = 'IntAct'
    
    # homology modelling method used to create structural models
    # options: template_based, model_based
    model_method = 'model_based'
    
    # method used to perform calculations; 'geometry' or 'physics'
    edgetic_method = 'physics'
    
    # method that was used to calculate edgetic mutation ∆∆G
    # options: bindprofx, foldx
    edgetic_ddg = 'foldx'
    
    # minimum RSA required for surface residues
    buriedMaxRSA = 0.25
    
    # show figures
    showFigs = False
    
    # parent directory of all data files
    #dataDir = Path('../data')
    dataDir = Path('/Volumes/MG_Samsung/edgotype_fitness_effect_full_model/data')
    
    # parent directory of all processed data files
    procDir = dataDir / 'processed'
    
    # directory of processed data files specific to interactome
    interactomeDir = procDir / interactome_name
    
    # directory of processed model-related data files specific to interactome
    modellingDir = interactomeDir / model_method
    
    # directory of edgetic mutation calculation method
    edgeticDir = modellingDir / edgetic_method
    
    # figure directory
    figDir = Path('../figures') / interactome_name / model_method / edgetic_method
    
    if edgetic_method is 'physics':
        edgeticDir = edgeticDir / (edgetic_ddg + '_edgetics')
        figDir = figDir / (edgetic_ddg + '_edgetics')
    
    # input data files
    naturalMutationsFile = edgeticDir / 'nondisease_mutation_RSA.txt'
    diseaseMutationsFile = edgeticDir / 'disease_mutation_RSA.txt'
    natMutDdgFile = edgeticDir / 'nondis_mut_folding_ddg_foldx.txt'
    disMutDdgFile = edgeticDir / 'dis_mut_folding_ddg_foldx.txt'
    
    # output data files
    natMutOutFile = edgeticDir / 'nondisease_mutation_RSA_ddg.txt'
    disMutOutFile = edgeticDir / 'disease_mutation_RSA_ddg.txt'
    
    #------------------------------------------------------------------------------------
    # load mutations
    #------------------------------------------------------------------------------------
    
    print('Reading mutations')
    naturalMutations = pd.read_table (naturalMutationsFile, sep='\t')
    diseaseMutations = pd.read_table (diseaseMutationsFile, sep='\t')
    
    print()
    print('Number of mutations:')
    print('non-disease mutations: %d' % len(naturalMutations))
    print('disease mutations: %d' % len(diseaseMutations))
    
    natMutationProteins = set(naturalMutations["protein"].tolist())
    disMutationProteins = set(diseaseMutations["protein"].tolist())
    mutationProteins = natMutationProteins | disMutationProteins
    
    print()
    print('Number of proteins carrying mutations:')
    print('non-disease mutations: %d' % len(natMutationProteins))
    print('disease mutations: %d' % len(disMutationProteins))
    print('all mutations: %d' % len(mutationProteins))
    
    #------------------------------------------------------------------------------------
    # load mutation effect on protein stability
    #------------------------------------------------------------------------------------
    
    natMutDdg = read_protein_mutation_ddg (natMutDdgFile, type = 'folding')
    disMutDdg = read_protein_mutation_ddg (disMutDdgFile, type = 'folding')
    
    energy_change = []
    for _, row in naturalMutations.iterrows():
        k = row.protein, row.mut_position, row.mut_res
        if k in natMutDdg:
            pdbid, chainID, ch_mut, ddg = natMutDdg[k]
            energy_change.append(ddg)
        else:
            energy_change.append(np.nan)
    naturalMutations["ddg"] = energy_change
    naturalMutations = naturalMutations [np.isnan(naturalMutations["ddg"]) == False].reset_index(drop=True)
    
    energy_change = []
    for _, row in diseaseMutations.iterrows():
        k = row.protein, row.mut_position, row.mut_res
        if k in disMutDdg:
            pdbid, chainID, ch_mut, ddg = disMutDdg[k]
            energy_change.append(ddg)
        else:
            energy_change.append(np.nan)
    diseaseMutations["ddg"] = energy_change
    diseaseMutations = diseaseMutations [np.isnan(diseaseMutations["ddg"]) == False].reset_index(drop=True)
    
    naturalMutations.to_csv (natMutOutFile, index=False, sep='\t')
    diseaseMutations.to_csv (disMutOutFile, index=False, sep='\t')
    
    natMutRSA = naturalMutations["RSA"].tolist()
    disMutRSA = diseaseMutations["RSA"].tolist()
    natMutDDG = naturalMutations["ddg"].tolist()
    disMutDDG = diseaseMutations["ddg"].tolist()
    
    print()
    print('Mutations with calculated change in stability and RSA:')
    print('non-disease mutations: %d' % len(natMutRSA))
    print('disease mutations: %d' % len(disMutRSA))
    
    xticklabels = [round(i, 1) for i in np.arange(0, 1.1, 0.2)]
    curve_plot ([disMutDDG, natMutDDG],
                xdata = [disMutRSA, natMutRSA],
                styles = ['^m', 'og'],
                fitstyles = ['--m', '--g'],
                capsize = 12,
                msize = 20,
                ewidth = 3,
                ecolors = ['m', 'g'],
                fontsize = 24,
                xticklabels = xticklabels,
                xlabel = 'RSA',
                ylabel = 'Change in protein stability\n∆∆G (kcal / mol)',
                #leg = ('Disease mutations', 'Non-disease mutations'),
                xlim = [-0.05, 1],
                ylim = [-0.25, 6],
                compress = True,
                linefit = False,
                xstart = -0.05,
                perbin = 2,
                show = showFigs,
                figdir = figDir,
                figname = 'Stability_vs_RSA')
    
    surfaceNatMutDDG = [d for d, rsa in zip(natMutDDG, natMutRSA) if rsa > buriedMaxRSA]
    buriedNatMutDDG = [d for d, rsa in zip(natMutDDG, natMutRSA) if rsa <= buriedMaxRSA]
    surfaceDisMutDDG = [d for d, rsa in zip(disMutDDG, disMutRSA) if rsa > buriedMaxRSA]
    buriedDisMutDDG = [d for d, rsa in zip(disMutDDG, disMutRSA) if rsa <= buriedMaxRSA]
    
    multi_bar_plot ([(np.mean(surfaceNatMutDDG), np.mean(buriedNatMutDDG)),
                     (np.mean(surfaceDisMutDDG), np.mean(buriedDisMutDDG))],
                    [(sderror(surfaceNatMutDDG), sderror(buriedNatMutDDG)),
                     (sderror(surfaceDisMutDDG), sderror(buriedDisMutDDG))],
                  ylabel = 'Change in protein stability ∆∆G (kcal / mol)',
                  xlabels = ('Exposed', 'Buried'),
                  leg = ('Non-disease mutations', 'Disease mutations'),
                  colors = ('c', 'm'),
                  edgecolor = 'k',
                  ewidth = 2,
                  bargap = 0.02,
                  fontsize = 16,
                  ylim = [0, 4],
                  show = showFigs,
                  figdir = figDir,
                  figname = 'Stability_vs_surface_buried')

if __name__ == "__main__":
    main()

