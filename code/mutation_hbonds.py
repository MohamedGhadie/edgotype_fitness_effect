#----------------------------------------------------------------------------------------
# Calculate change in number of potential H-bonds upon mutation.
#----------------------------------------------------------------------------------------

import os
import pandas as pd
import numpy as np
from pathlib import Path
from stat_tools import t_test, sderror
from plot_tools import multi_bar_plot

def main():
    
    edgotype = 'edgetic'
    
    # reference interactome names
    interactome_names = ['HuRI', 'IntAct', 'experiment']
    
    # structural interactome names for plot labels
    xlabels = ['Y2H-SI', 'Lit-SI', 'Experiment']
    
    # fitness effects
    legend = ['Non-disease %s mutations' % edgotype,
              'Disease-causing %s mutations' % edgotype]
    
    # bar colors for different fitness effects
    barColors = ['limegreen', 'orange']
    
    # bar hatches for different interactomes
    barHatches = ['.', '//']
    
    # show figure legend
    showLeg = False
    
    # show figures
    showFigs = False
    
    # parent directory of all data files
    dataDir = Path('../data')
    
    # figure directory
    figDir = Path('../figures') / 'combined'
    
    # input data files
    residuePropertyFile = dataDir / 'aa_properties.xlsx'
    naturalMutationsFile = 'nondisease_mutation_edgotype.txt'
    diseaseMutationsFile = 'disease_mutation_edgotype.txt'
    
    # create output directories if not existing
    if not figDir.exists():
        os.makedirs(figDir)
    
    #------------------------------------------------------------------------------------
    # Calculate gene knockout effect
    #------------------------------------------------------------------------------------
    
    properties = pd.read_excel (residuePropertyFile, sheet_name = 'Sheet1')
    
    prop = {}
    for _, row in properties.iterrows():
        prop[row.Residue] = {'charge':row.Charge, 
                             'group': row.Side_chain_group,
                             'class':row.Class,
                             'abundance':row.Abundance,
                             'flexibility':row.Side_chain_flexibility,
                             'interaction_modes':row.Interaction_modes,
                             'H-bonds':row.Potential_H_bonds,
                             'weight':row.Molecular_weight,
                             'isoelectric_point':row.Isoelectric_point,
                             'hydrophobicity':row.Hydrophobicity,
                             'codons':row.Standard_codons}
        
    #------------------------------------------------------------------------------------
    # process mutations
    #------------------------------------------------------------------------------------
    
    natH, disH = [], []
    natSE, disSE = [], []
    
    for name in interactome_names:
        naturalMutations = pd.read_table (dataDir / name / naturalMutationsFile, sep='\t')
        diseaseMutations = pd.read_table (dataDir / name / diseaseMutationsFile, sep='\t')
        
        if name is 'experiment':
            naturalMutations["edgotype"] = naturalMutations["Edgotype_class"].apply(lambda x: x.lower())
            diseaseMutations["edgotype"] = diseaseMutations["Edgotype_class"].apply(lambda x: x.lower())
        
        print('\n------------------------------------------------------------------------')
        print('Interactome: %s' % name)
        print('------------------------------------------------------------------------\n')
        
        print('non-disease mutations: %d' % len(naturalMutations))
        print('disease mutations: %d' % len(diseaseMutations))
        
        naturalMutations = naturalMutations [naturalMutations["edgotype"] == edgotype].reset_index(drop=True)
        diseaseMutations = diseaseMutations [diseaseMutations["edgotype"] == edgotype].reset_index(drop=True)
        
        nat_wt_hbonds = naturalMutations['wt_res'].apply(lambda x: prop[x]['H-bonds'])
        nat_mt_hbonds = naturalMutations['mut_res'].apply(lambda x: prop[x]['H-bonds'])
        nat_hbonds_change = nat_mt_hbonds - nat_wt_hbonds
        
        dis_wt_hbonds = diseaseMutations['wt_res'].apply(lambda x: prop[x]['H-bonds'])
        dis_mt_hbonds = diseaseMutations['mut_res'].apply(lambda x: prop[x]['H-bonds'])
        dis_hbonds_change = dis_mt_hbonds - dis_wt_hbonds
        
        natH.append(np.mean(nat_hbonds_change))
        disH.append(np.mean(dis_hbonds_change))
        
        natSE.append(sderror(nat_hbonds_change))
        disSE.append(sderror(dis_hbonds_change))
                
        # print results
        print()
        print('Average change in number of potential H-bonds:')
        print('non-disease %s mutations: %.2f (SE = %g, n = %d)' % (edgotype, natH[-1], natSE[-1], len(nat_hbonds_change)))
        print('disease %s mutations: %.2f (SE = %g, n = %d)' % (edgotype, disH[-1], disSE[-1], len(dis_hbonds_change)))
        print()
        t_test(nat_hbonds_change.values, dis_hbonds_change.values)
    
    multi_bar_plot ([natH, disH],
                    errors = [natSE, disSE],
                    xlabels = xlabels,
                    #ylabels = [0, 1, 2],
                    #ylim = [0, 2],
                    ylabel = 'Change in number of potential H-bonds',
                    colors = barColors,
                    hatches = barHatches,
                    barwidth = 0.2,
                    linewidth = 1.5,
                    bargap = 0.05,
                    capsize = 8,
                    fmt = '.k',
                    msize = 14,
                    ewidth = 1.5,
                    edgecolor = 'black',
                    ecolors = 'k',
                    fontsize = 20,
                    yMinorTicks = 1,
                    leg = legend if showLeg else None,
                    show = showFigs,
                    figdir = figDir,
                    figname = '%s_mut_hbonds' % edgotype)

if __name__ == "__main__":
    main()