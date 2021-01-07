#----------------------------------------------------------------------------------------
# Calculate change in hydrophobicity upon mutation.
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
    
    natHB, disHB = [], []
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
        
        nat_wt_hb = naturalMutations['wt_res'].apply(lambda x: prop[x]['hydrophobicity'])
        nat_mt_hb = naturalMutations['mut_res'].apply(lambda x: prop[x]['hydrophobicity'])
        nat_wt_hb = nat_wt_hb.apply(lambda x: 0.01 if x == 0 else x)
        nat_mt_hb = nat_mt_hb.apply(lambda x: 0.01 if x == 0 else x)
        nat_hb_change = nat_mt_hb / nat_wt_hb
        
        dis_wt_hb = diseaseMutations['wt_res'].apply(lambda x: prop[x]['hydrophobicity'])
        dis_mt_hb = diseaseMutations['mut_res'].apply(lambda x: prop[x]['hydrophobicity'])
        dis_wt_hb = dis_wt_hb.apply(lambda x: 0.01 if x == 0 else x)
        dis_mt_hb = dis_mt_hb.apply(lambda x: 0.01 if x == 0 else x)
        dis_hb_change = dis_mt_hb / dis_wt_hb
        
        natHB.append(np.mean(nat_hb_change))
        disHB.append(np.mean(dis_hb_change))
        
        natSE.append(sderror(nat_hb_change))
        disSE.append(sderror(dis_hb_change))
                
        # print results
        print()
        print('Average hydrophobicity change:')
        print('non-disease %s mutations: %f (SE = %g, n = %d)' % (edgotype, natHB[-1], natSE[-1], len(nat_hb_change)))
        print('disease %s mutations: %f (SE = %g, n = %d)' % (edgotype, disHB[-1], disSE[-1], len(dis_hb_change)))
        print()
        t_test(nat_hb_change.values, dis_hb_change.values)
       
    multi_bar_plot ([natHB, disHB],
                    errors = [natSE, disSE],
                    xlabels = xlabels,
                    ylabels = [0, 10, 20, 30],
                    ylim = [0, 30],
                    ylabel = 'Fold change in hydrophobicity',
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
                    #topAxis = True,
                    leg = legend if showLeg else None,
                    show = showFigs,
                    figdir = figDir,
                    figname = '%s_mut_hydrophobicity' % edgotype)

if __name__ == "__main__":
    main()
