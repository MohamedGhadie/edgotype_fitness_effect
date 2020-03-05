#----------------------------------------------------------------------------------------
# Predict interactome perturbations based on physics. For mutations predicted to perturb 
# interaction based on geometry, this script re-predicts perturbations for those mutations 
# based on change in interaction binding energy caused by mutations at the interface.
#
# Requirements before running this script:
# 
# First run the following scripts:
# - predict_perturbations_geometry.py
# - produce_mutation_structure_maps.py
#
# For physics_based perturbation prediction using ∆∆G values calculated by bindprofx:
# 1- run script produce_bindprofx_jobs.py on the two data files nondisease_mutations_onchains.txt 
#    and disease_mutations_onchains.txt to produce bindprofx jobs for mutations grouped per structure
# 2- submit jobs to bindprofx to calculate mutation ∆∆G
# 3- run script process_bindprofx_results.py to process bindprofx results 
#    and produce second-round jobs for single mutations of failed jobs
# 4- submit second-round jobs to bindprofx to calculate mutation ∆∆G
# 5- run script process_bindprofx_results.py to process bindprofx second-round results
# 6- repeat steps 1-5 until no new jobs are produced
# 7- Files containing final ∆∆G results should be named:
#    'nondisease_mutations_bindprofx_ddg.txt' and 'disease_mutations_bindprofx_ddg.txt'
#
# For physics_based perturbation prediction using ∆∆G values calculated by foldx:
# 1- run script produce_foldx_jobs.py on the two data files nondisease_mutations_onchains.txt 
#    and disease_mutations_onchains.txt to produce foldx jobs for single mutations
# 2- submit jobs to foldx to calculate mutation ∆∆G
# 3- run script process_foldx_jobs.py to process foldx results
# 4- repeat steps 1-3 until no new jobs are produced
# 5- Files containing final ∆∆G results should be named:
#    'nondisease_mutations_foldx_ddg.txt' and 'disease_mutations_foldx_ddg.txt'
#----------------------------------------------------------------------------------------

import os
import pandas as pd
import numpy as np
from pathlib import Path
from ddg_tools import read_protein_mutation_ddg
from stat_tools import t_test, fisher_test, sderror, sderror_on_fraction
from plot_tools import bar_plot, multi_histogram_plot

def main():
    
    # reference interactome name
    # options: HI-II-14, IntAct
    interactome_name = 'IntAct'
    
    # homology modelling method used to create structural models
    # options: template_based, model_based
    model_method = 'model_based'
    
    # method of calculating mutation ∆∆G for which results will be used
    # options: bindprofx, foldx
    ddg_method = 'foldx'
    
    # Minimum reduction in binding free energy DDG required for interaction perturbation
    ddgCutoff = 0.5
    
    # structural interactome names
    struc_name = {'HI-II-14':'Y2H-SI', 'IntAct':'IntAct-SI'}
    
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
    methodDir = modellingDir / 'physics'
    
    # figure directory
    figDir = Path('../figures') / interactome_name / model_method / 'physics' / ('%s_edgetics' % ddg_method)
    
    # input data files
    #geometryPerturbsFile = interactomeDir / 'unique_mutation_perturbs_geometry.pkl'
    natMutDDGinFile = modellingDir / ('nondis_mut_binding_ddg_%s.txt' % ddg_method)
    disMutDDGinFile = modellingDir / ('dis_mut_binding_ddg_%s.txt' % ddg_method)
    structuralInteractomeFile = modellingDir / 'human_structural_interactome.txt'
    
    # output data files
    natMutDDGoutFile = methodDir / ('nondisMut_binding_ddg_table_%s.txt' % ddg_method)
    disMutDDGoutFile = methodDir / ('disMut_binding_ddg_table_%s.txt' % ddg_method)
    
    # create output directories if not existing
    if not methodDir.exists():
        os.makedirs(methodDir)
    if not figDir.exists():
        os.makedirs(figDir)
    
    #------------------------------------------------------------------------------------
    # Fraction of mutation-targeted PPIs with ∆∆G exceeding a specified cutoff
    #------------------------------------------------------------------------------------
    
    # read change in binding free energy for interfacial mutations mapped on PDB chains
    naturalMutationsDDG = read_protein_mutation_ddg (natMutDDGinFile, 'binding')
    diseaseMutationsDDG = read_protein_mutation_ddg (disMutDDGinFile, 'binding')
    
    naturalMutations = pd.DataFrame(columns=["protein", "partner", "protein_pos", "mut_res", 
                                             "pdb_id", "chain_id", "chain_partner", "chain_mut", "ddg"])
    for i, item in enumerate(naturalMutationsDDG.items()):
        naturalMutations.loc[i] = item[0] + item[1]
    
    diseaseMutations = pd.DataFrame(columns=["protein", "partner", "protein_pos", "mut_res", 
                                             "pdb_id", "chain_id", "chain_partner", "chain_mut", "ddg"])
    for i, item in enumerate(diseaseMutationsDDG.items()):
        diseaseMutations.loc[i] = item[0] + item[1]
    
    naturalMutations.to_csv (natMutDDGoutFile, index=False, sep='\t')
    diseaseMutations.to_csv (disMutDDGoutFile, index=False, sep='\t')
    
    natMutDDGs = naturalMutations["ddg"].values
    disMutDDGs = diseaseMutations["ddg"].values
    
    numNatural_ddg_considered = len(natMutDDGs)
    numDisease_ddg_considered = len(disMutDDGs)
    
    print( '\n' + 'Avg. change in binding energy (∆∆G) for mutation-targeted PPIs:' )
    print( 'Non-disease: %.1f (SE = %g, n = %d)' % (np.mean(natMutDDGs),
                                                    sderror(natMutDDGs),
                                                    numNatural_ddg_considered) )
    
    print( 'Disease: %.1f (SE = %g, n = %d)' % (np.mean(disMutDDGs),
                                                sderror(disMutDDGs),
                                                numDisease_ddg_considered) )
    # Statistical significance of difference in means
    t_test(natMutDDGs, disMutDDGs)
    
    multi_histogram_plot ([disMutDDGs, natMutDDGs],
                          ['red', 'green'],
                          xlabel = 'Change in binding free energy (∆∆G)',
                          ylabel = 'Number of mutations',
                          leg = ['Disease interfacial mutations', 'Non-disease interfacial mutations'],
                          bins = 25,
                          alpha = 0.7,
                          fontsize = 24,
                          show = showFigs,
                          figdir = figDir,
                          figname = 'mut_ddg_histogram_%s' % ddg_method)
    
    numNatural_ddg = sum(natMutDDGs > ddgCutoff)
    numDisease_ddg = sum(disMutDDGs > ddgCutoff)
    fracNatural_ddg = numNatural_ddg / numNatural_ddg_considered
    fracDisease_ddg = numDisease_ddg / numDisease_ddg_considered
    fracNatural_ddg_error = sderror_on_fraction( numNatural_ddg, numNatural_ddg_considered )
    fracDisease_ddg_error = sderror_on_fraction( numDisease_ddg, numDisease_ddg_considered )
    
    print( '\n' + 'Fraction of mutation-targeted PPIs with ∆∆G > %.1f kcal/mol:' % ddgCutoff )
    print( 'Non-disease: %.3f (SE = %g, ntot = %d)' % (fracNatural_ddg,
                                                       fracNatural_ddg_error,
                                                       numNatural_ddg_considered) )
    
    print( 'Disease: %.3f (SE = %g, ntot = %d)' % (fracDisease_ddg,
                                                   fracDisease_ddg_error,
                                                   numDisease_ddg_considered) )
    
    # Statistical significance of difference in fractions
    fisher_test([numNatural_ddg, numNatural_ddg_considered - numNatural_ddg],
                [numDisease_ddg, numDisease_ddg_considered - numDisease_ddg])
    
    bar_plot([fracNatural_ddg, fracDisease_ddg],
             error = [fracNatural_ddg_error, fracDisease_ddg_error],
             xlabels = ('%s PPIs\nwith non-disease\nmutations\nat interface' % struc_name[interactome_name],
                        '%s PPIs\nwith disease\nmutations\nat interface' % struc_name[interactome_name]),
             ylabel = ('Fraction with ∆∆G > %.1f kcal / mol' % ddgCutoff),
             ylabels = [0, 0.2, 0.4, 0.6, 0.8],
             ylim = [0, 0.8],
             colors = ['turquoise', 'magenta'],
             edgecolor = 'black',
             ewidth = 2.5,
             barwidth = 0.6,
             fontsize = 24,
             capsize = 10,
             msize = 26,
             show = showFigs,
             figdir = figDir,
             figname = 'mut_ddg_frac_>%.1f_%s' % (ddgCutoff, ddg_method))

if __name__ == "__main__":
    main()
