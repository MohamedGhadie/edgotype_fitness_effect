import os
import pandas as pd
import numpy as np
from pathlib import Path
from ddg_tools import read_protein_mutation_ddg
from stat_tools import fisher_test, sderror_on_fraction
from plot_tools import pie_plot

def main():
    
    # reference interactome name
    # options: HI-II-14, IntAct
    interactome_name = 'HI-II-14'
    
    # homology modelling method used to create structural models
    # options: template_based, model_based
    model_method = 'model_based'
    
    # method used to predict edgetic perturbations
    # options: geometry, physics
    edgetic_method = 'geometry'
    
    # minimum change in protein free energy required for quasi-null mutations
    qnMinDDG = 5
        
    # calculate confidence intervals
    computeConfidenceIntervals = True
    
    # confidence interval (%)
    CI = 95
    
    # set to True to calculate dispensable PPI content using fraction of mono-edgetic mutations 
    # instead of edgetic mutations
    mono_edgetic = False
    
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
    
    # directory of calculation method
    methodDir = modellingDir / edgetic_method
    
    # figure directory
    figDir = Path('../figures') / interactome_name / model_method / edgetic_method
    
    # input data files
    natMutLocFile = methodDir / 'nondisease_mutation_struc_loc.txt'
    disMutLocFile = methodDir / 'disease_mutation_struc_loc.txt'
    natMutDdgFile = modellingDir / 'nondis_mut_folding_ddg_foldx.txt'
    disMutDdgFile = modellingDir / 'dis_mut_folding_ddg_foldx.txt'
    
    # output data files
    natMutEdgotypeFile = methodDir / 'nondisease_mutation_edgotype.txt'
    disMutEdgotypeFile = methodDir / 'disease_mutation_edgotype.txt'
    
    # create output directories if not existing
    if not methodDir.exists():
        os.makedirs(methodDir)
    if not figDir.exists():
        os.makedirs(figDir)
    
    #------------------------------------------------------------------------------------
    # load mutations
    #------------------------------------------------------------------------------------
    
    print('Reading mutation locations on protein structure')
    naturalMutations = pd.read_table (natMutLocFile, sep='\t')
    diseaseMutations = pd.read_table (disMutLocFile, sep='\t')
    
    #------------------------------------------------------------------------------------
    # load mutation effect on protein stability
    #------------------------------------------------------------------------------------
    
    natMutDdg = read_protein_mutation_ddg (natMutDdgFile, type = 'folding')
    disMutDdg = read_protein_mutation_ddg (disMutDdgFile, type = 'folding')
    
    ddgValues = {'nat_mut':natMutDdg, 'dis_mut':disMutDdg}
    mutations = {'nat_mut':naturalMutations, 'dis_mut':diseaseMutations}
    for k, mut in mutations.items():
        energy_change = []
        for _, row in mut.iterrows():
            p = row.protein, row.mut_position, row.mut_res
            if p in ddgValues[k]:
                pdbid, chainID, ch_mut, ddg = ddgValues[k][p]
                energy_change.append(ddg)
            else:
                energy_change.append(np.nan)
        mutations[k]["ddg"] = energy_change
    
    #------------------------------------------------------------------------------------
    # Assign mutation edgotypes
    #------------------------------------------------------------------------------------    
    
    mutations = {'nat_mut':naturalMutations, 'dis_mut':diseaseMutations}
    for k, mut in mutations.items():
        edgotypes = []
        for _, row in mut.iterrows():
            if row.edgotype == 'edgetic':
                edgotypes.append('edgetic')
            elif row.edgotype == 'non-edgetic':
                if row.structural_location == 'exposed-noninterface':
                    edgotypes.append('quasi-wild-type')
                elif row.structural_location == 'buried':
                    if np.isnan(row.ddg):
                        edgotypes.append('-')
                    elif row.ddg >= qnMinDDG:
                        edgotypes.append('quasi-null')
                    else:
                        edgotypes.append('quasi-wild-type')
            else:
                edgotypes.append('-')
        mutations[k]['edgotype'] = edgotypes
    
    naturalMutations.to_csv(natMutEdgotypeFile, index=False, sep='\t')
    diseaseMutations.to_csv(disMutEdgotypeFile, index=False, sep='\t')
    
    #------------------------------------------------------------------------------------
    # Display mutation edgotype fractions
    #------------------------------------------------------------------------------------ 
    
    naturalMutations = naturalMutations [naturalMutations["edgotype"] != '-'].reset_index(drop=True)
    diseaseMutations = diseaseMutations [diseaseMutations["edgotype"] != '-'].reset_index(drop=True)
    
    numNaturalMut_type = {'QN': sum(naturalMutations["edgotype"] == 'quasi-null'),
                          'E': sum(naturalMutations["edgotype"] == 'edgetic'),
                          'QW': sum(naturalMutations["edgotype"] == 'quasi-wild-type')}
    
    numDiseaseMut_type = {'QN': sum(diseaseMutations["edgotype"] == 'quasi-null'),
                          'E': sum(diseaseMutations["edgotype"] == 'edgetic'),
                          'QW': sum(diseaseMutations["edgotype"] == 'quasi-wild-type')}
    
    numNaturalMut_considered = sum(numNaturalMut_type.values())
    numDiseaseMut_considered = sum(numDiseaseMut_type.values())
    
    fracNaturalMut_type = {'QN': numNaturalMut_type['QN'] / numNaturalMut_considered,
                           'E': numNaturalMut_type['E'] / numNaturalMut_considered,
                           'QW': numNaturalMut_type['QW'] / numNaturalMut_considered}
    fracDiseaseMut_type = {'QN': numDiseaseMut_type['QN'] / numDiseaseMut_considered,
                           'E': numDiseaseMut_type['E'] / numDiseaseMut_considered,
                           'QW': numDiseaseMut_type['QW'] / numDiseaseMut_considered}
    
    print()
    print('Fraction of edgotypes among non-disease mutations:')
    print('quasi-null: %f%% (SE = %g, %d out of %d)' % (100 * fracNaturalMut_type['QN'],
                                                        sderror_on_fraction (numNaturalMut_type['QN'], numNaturalMut_considered),
                                                        numNaturalMut_type['QN'],
                                                        numNaturalMut_considered))
    print('edgetic: %f%% (SE = %g, %d out of %d)' % (100 * fracNaturalMut_type['E'],
                                                     sderror_on_fraction (numNaturalMut_type['E'], numNaturalMut_considered),
                                                     numNaturalMut_type['E'],
                                                     numNaturalMut_considered))
    print('quasi-wild-type: %f%% (SE = %g, %d out of %d)' % (100 * fracNaturalMut_type['QW'],
                                                             sderror_on_fraction (numNaturalMut_type['QW'], numNaturalMut_considered),
                                                             numNaturalMut_type['QW'],
                                                             numNaturalMut_considered))
    
    print()
    print('Fraction of edgotypes among disease mutations:')
    print('quasi-null: %f%% (SE = %g, %d out of %d)' % (100 * fracDiseaseMut_type['QN'],
                                                        sderror_on_fraction (numDiseaseMut_type['QN'], numDiseaseMut_considered),
                                                        numDiseaseMut_type['QN'],
                                                        numDiseaseMut_considered))
    print('edgetic: %f%% (SE = %g, %d out of %d)' % (100 * fracDiseaseMut_type['E'],
                                                     sderror_on_fraction (numDiseaseMut_type['E'], numDiseaseMut_considered),
                                                     numDiseaseMut_type['E'],
                                                     numDiseaseMut_considered))
    print('quasi-wild-type: %f%% (SE = %g, %d out of %d)' % (100 * fracDiseaseMut_type['QW'],
                                                             sderror_on_fraction (numDiseaseMut_type['QW'], numDiseaseMut_considered),
                                                             numDiseaseMut_type['QW'],
                                                             numDiseaseMut_considered))
    
    fisher_test ([numNaturalMut_type['QN'] + numNaturalMut_type['E'], numNaturalMut_type['QW']],
                 [numDiseaseMut_type['QN'] + numDiseaseMut_type['E'], numDiseaseMut_type['QW']])
        
    pie_plot([numNaturalMut_type['QW'], numNaturalMut_type['E'], numNaturalMut_type['QN']],
             angle = 90,
             colors = ['mediumslateblue', 'purple', 'red'],
             edgewidth = 2,
             show = showFigs,
             figdir = figDir,
             figname = 'non_disease_mut_edgotype')
    pie_plot([numDiseaseMut_type['QW'], numDiseaseMut_type['E'], numDiseaseMut_type['QN']],
             angle = 90,
             colors = ['mediumslateblue', 'purple', 'red'],
             edgewidth = 2,
             show = showFigs,
             figdir = figDir,
             figname = 'disease_mut_edgotype')

if __name__ == "__main__":
    main()