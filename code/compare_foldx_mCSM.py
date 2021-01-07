#----------------------------------------------------------------------------------------
# This script compares mutation binding ∆∆G values calculated by FoldX and mCSM..
#----------------------------------------------------------------------------------------

import numpy as np
import pandas as pd
from pathlib import Path
from scipy.stats import pearsonr
from energy_tools import read_protein_mutation_ddg

def main():
    
    # reference interactome name: HI-II-14, HuRI, IntAct
    interactome_name = 'IntAct'
    
    # homology modelling method used to create structural models
    # options: template_based, model_based
    model_method = 'model_based'
    
    # minimum reduction in binding free energy ∆∆G required for PPI disruption
    ddgCutoff = 0.5
    
    # parent directory of all data files
    dataDir = Path('../data')
    
    # parent directory of all processed data files
    procDir = dataDir / 'processed'
    
    # directory of processed data files specific to interactome
    interactomeDir = procDir / interactome_name
    
    # directory of processed model-related data files specific to interactome
    modellingDir = interactomeDir / model_method
    
    # input data files
    foldx_natMutFile = modellingDir / 'nondis_mut_binding_ddg_foldx.txt'
    foldx_disMutFile = modellingDir / 'dis_mut_binding_ddg_foldx.txt'
    mcsm_natMutFile = modellingDir / 'nondis_mut_binding_ddg_mCSM.txt'
    mcsm_disMutFile = modellingDir / 'dis_mut_binding_ddg_mCSM.txt'
    
    foldx = read_protein_mutation_ddg (foldx_natMutFile, type = 'binding')
    foldx_disMut = read_protein_mutation_ddg (foldx_disMutFile, type = 'binding')
    mcsm = read_protein_mutation_ddg (mcsm_natMutFile, type = 'binding')
    mcsm_disMut = read_protein_mutation_ddg (mcsm_disMutFile, type = 'binding')
    
    for k, v in foldx_disMut.items():
        foldx[k] = v
    for k, v in mcsm_disMut.items():
        mcsm[k] = v
    
    mut = pd.DataFrame()
    mcsm_data = [k + (v[-1],) for k, v in mcsm.items()]
    mut["protein"], mut["partner"], mut["mut_position"], mut["mut_res"], mut["mCSM_ddg"] = zip(*mcsm_data)
    
    foldx_ddg = []
    for _, row in mut.iterrows():
        k = row.protein, row.partner, row.mut_position, row.mut_res
        if k in foldx:
            foldx_ddg.append(foldx[k][-1])
        else:
            foldx_ddg.append(np.nan)
    mut["foldx_ddg"] = foldx_ddg
    mut = mut[(np.isnan(mut["foldx_ddg"]) == False) & (np.isnan(mut["mCSM_ddg"]) == False)]
    
    print('Number of mutations = %d' % len(mut))
    print('Pearson correlation = %f, p-value = %f' % pearsonr(mut["foldx_ddg"].values, mut["mCSM_ddg"].values))
    
    mut["edgetic_foldx"] = mut["foldx_ddg"] > ddgCutoff
    mut["edgetic_mCSM"] = mut["mCSM_ddg"] > ddgCutoff
    
    print(mut)

    TP = sum(mut["edgetic_foldx"] & mut["edgetic_mCSM"])
    FP = sum(mut["edgetic_foldx"] & (mut["edgetic_mCSM"] == False))
    TN = sum((mut["edgetic_foldx"] == False) & (mut["edgetic_mCSM"] == False))
    FN = sum((mut["edgetic_foldx"] == False) & mut["edgetic_mCSM"])
    
    TPR = TP / (TP + FN)
    FPR = FP / (FP + TN)
    
    print('TPR = %.2f' % TPR)
    print('FPR = %.2f' % FPR)
    
if __name__ == '__main__':
    main()
