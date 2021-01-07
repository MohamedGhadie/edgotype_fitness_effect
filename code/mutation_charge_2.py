#----------------------------------------------------------------------------------------
# Calculate change in molecular weight upon mutation.
#----------------------------------------------------------------------------------------

import os
import pandas as pd
import numpy as np
from pathlib import Path
from stat_tools import t_test, sderror, fisher_test
from plot_tools import multi_bar_plot

def main():
    
    # reference interactome names
    interactome_names = ['HuRI', 'IntAct']
    
    # structural interactome names for plot labels
    xlabels = ['Y2H-SI', 'Lit-SI']
    
    # fitness effects
#     legend = ['Non-disease %s mutations' % edgotype,
#               'Disease-causing %s mutations' % edgotype]
    
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
    
    diseaseMutations = pd.read_table (dataDir / 'HuRI' / diseaseMutationsFile, sep='\t')
    diseaseMutations2 = pd.read_table (dataDir / 'IntAct' / diseaseMutationsFile, sep='\t')
    
    diseaseMutations = diseaseMutations.append(diseaseMutations2)
    diseaseMutations = diseaseMutations.drop_duplicates(subset=["protein",
                                                                "mut_position",
                                                                "mut_res"]).reset_index(drop=True)
    
    dis_wt_hb = diseaseMutations['wt_res'].apply(lambda x: prop[x]['charge'])
    dis_mt_hb = diseaseMutations['mut_res'].apply(lambda x: prop[x]['charge'])
#     dis_wt_hb = dis_wt_hb.apply(lambda x: 0.01 if x == 0 else x)
#     dis_mt_hb = dis_mt_hb.apply(lambda x: 0.01 if x == 0 else x)
    #dis_hb_change = dis_mt_hb / dis_wt_hb
    dis_hb_change = dis_mt_hb != dis_wt_hb
    
    dis_hb_change_w = dis_hb_change[diseaseMutations["edgotype"] == 'quasi-wild-type'].values
    dis_hb_change_e = dis_hb_change[diseaseMutations["edgotype"] == 'edgetic'].values
    dis_hb_change_q = dis_hb_change[diseaseMutations["edgotype"] == 'quasi-null'].values
                
    # print results
#     print()
#     print('Average hydrophobicity change:')
#     print('disease quasi-wild-type mutations: %f (SE = %g, n = %d)' % (np.mean(dis_hb_change_w),
#                                                                        sderror(dis_hb_change_w),
#                                                                        len(dis_hb_change_w)))
#     print('disease edgetic mutations: %f (SE = %g, n = %d)' % (np.mean(dis_hb_change_e),
#                                                                sderror(dis_hb_change_e),
#                                                                len(dis_hb_change_e)))
#     print('disease quasi-null mutations: %f (SE = %g, n = %d)' % (np.mean(dis_hb_change_q),
#                                                                   sderror(dis_hb_change_q),
#                                                                   len(dis_hb_change_q)))
#     print()
#     t_test(dis_hb_change_q, dis_hb_change_e)
#     t_test(dis_hb_change_q, dis_hb_change_w)
    
    print()
    print('Fraction of mutations with charge change:')
    print('disease quasi-wild-type mutations: %f (n = %d)' 
            % (sum(dis_hb_change_w) / len(dis_hb_change_w), len(dis_hb_change_w)))
    print('disease edgetic mutations: %f (n = %d)' 
            % (sum(dis_hb_change_e) / len(dis_hb_change_e), len(dis_hb_change_e)))
    print('disease quasi-null mutations: %f (n = %d)' 
            % (sum(dis_hb_change_q) / len(dis_hb_change_q), len(dis_hb_change_q)))
    
    print()
    print('Statistical significance, edgetic and quasi-wild-type')
    fisher_test([sum(dis_hb_change_e), sum(dis_hb_change_e == False)],
                [sum(dis_hb_change_w), sum(dis_hb_change_w == False)])
    
    print()
    print('Statistical significance, quasi-null and edgetic')
    fisher_test([sum(dis_hb_change_q), sum(dis_hb_change_q == False)],
                [sum(dis_hb_change_e), sum(dis_hb_change_e == False)])
    
    print()
    print('Statistical significance, quasi-null and quasi-wild-type')
    fisher_test([sum(dis_hb_change_q), sum(dis_hb_change_q == False)],
                [sum(dis_hb_change_w), sum(dis_hb_change_w == False)])

if __name__ == "__main__":
    main()