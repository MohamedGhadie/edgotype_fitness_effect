import os
import pickle
import pandas as pd
import numpy as np
from pathlib import Path
from text_tools import write_list_table
from mutation_interface_edgotype import assign_edgotypes
from ddg_tools import read_protein_mutation_ddg
from stat_tools import fisher_test, sderror_on_fraction, proportion_ratio_CI, proportion_sum_CI
from plot_tools import pie_plot, curve_plot

def main():
    
    # reference interactome name
    # options: HI-II-14, IntAct
    interactome_name = 'HI-II-14'
    
    # mutation edgotype for which fitness effect is calculated
    # options: quasi-null, edgetic, quasi-wild-type
    edgotype = 'edgetic'
    
    # minimum change in protein free energy required for quasi-null mutations
    qnMinDDG = 5
    
    # assume edgotype probabilities of strongly detrimental (S) mutations to be similar to 
    # those of mildly deleterious (M) mutations. If set to False, strongly detrimental 
    # mutations are assumed to be all quasi-null
    assume_S_as_M = False
    
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
    
    # directory of calculation method
    methodDir = interactomeDir / 'geometry'
    
    # figure directory
    figDir = Path('../figures') / interactome_name / 'mutation_fitness_effect_using_geometry'
    
    if assume_S_as_M:
        figDir = figDir / 'assume_S_as_M'
    else:
        figDir = figDir / 'assume_S_quasi_null'
    
    # input data files
#     natMutLocFile = methodDir / 'nondisease_mutation_struc_loc.txt'
#     disMutLocFile = methodDir / 'disease_mutation_struc_loc.txt'
    mutationPerturbsFile = interactomeDir / 'unique_mutation_perturbs_geometry.pkl' # new
    natMutDdgFile = interactomeDir / 'nondisease_mutations_ddg.txt'
    disMutDdgFile = interactomeDir / 'disease_mutations_ddg.txt'
    
    # output data files
    natMutEdgotypeFile = methodDir / 'nondisease_mutation_edgotype.txt'
    disMutEdgotypeFile = methodDir / 'disease_mutation_edgotype.txt'
    outputFile = methodDir / ('%s_mut_fitness_effect.pkl' % edgotype)
    
    # create output directories if not existing
    if not methodDir.exists():
        os.makedirs(methodDir)
    if not figDir.exists():
        os.makedirs(figDir)
    
    #------------------------------------------------------------------------------------
    # terms to be used
    #------------------------------------------------------------------------------------
    
    # abbreviations for different structural regions
    etype_symbol = {'quasi-null':'QN', 'edgetic':'E', 'quasi-wild-type':'QW'}
    
    # abbreviation for selected structural region
    T = etype_symbol [edgotype]
    
    # prior conditional probabilities
    cond_probs = ['P(%s|%s)' % (T, p) for p in ('N','M','S')]
    
    # posterior conditional probabilities
    posteriors = ['P(%s|%s)' % (p, T) for p in ('N','M','S')]
    
    #------------------------------------------------------------------------------------
    # New: load mutations
    #------------------------------------------------------------------------------------
    
    with open(mutationPerturbsFile, 'rb') as f:
        naturalMutations, diseaseMutations = pickle.load(f)
    
    #------------------------------------------------------------------------------------
    # load mutations
    #------------------------------------------------------------------------------------
    
#     print('Reading processed mutations')
#     naturalMutations = pd.read_table (natMutLocFile, sep='\t')
#     diseaseMutations = pd.read_table (disMutLocFile, sep='\t')   
#     
#     naturalMutations = naturalMutations [naturalMutations["edgotype"].apply(
#                         lambda x: x in ['edgetic', 'non-edgetic'])].reset_index(drop=True)
#     diseaseMutations = diseaseMutations [diseaseMutations["edgotype"].apply(
#                         lambda x: x in ['edgetic', 'non-edgetic'])].reset_index(drop=True)
#     
#     print()
#     print('Number of mutations:')
#     print('non-disease mutations: %d' % len(naturalMutations))
#     print('disease mutations: %d' % len(diseaseMutations))
#     
#     natMutationProteins = set(naturalMutations["protein"].tolist())
#     disMutationProteins = set(diseaseMutations["protein"].tolist())
#     mutationProteins = natMutationProteins | disMutationProteins
#     
#     print()
#     print('Number of proteins carrying mutations:')
#     print('non-disease mutations: %d' % len(natMutationProteins))
#     print('disease mutations: %d' % len(disMutationProteins))
#     print('all mutations: %d' % len(mutationProteins))
    
    #------------------------------------------------------------------------------------
    # load mutation effect on protein stability
    #------------------------------------------------------------------------------------
    
#     natMutDdg = read_protein_mutation_ddg (natMutDdgFile, type = 'folding')
#     disMutDdg = read_protein_mutation_ddg (disMutDdgFile, type = 'folding')
#     
#     ddgValues = {'nat_mut':natMutDdg, 'dis_mut':disMutDdg}
#     mutations = {'nat_mut':naturalMutations, 'dis_mut':diseaseMutations}
#     for k, mut in mutations.items():
#         energy_change = []
#         for _, row in mut.iterrows():
#             p = row.protein, row.mut_position, row.mut_res
#             if p in ddgValues[k]:
#                 pdbid, chainID, ch_mut, ddg = ddgValues[k][p]
#                 energy_change.append(ddg)
#             else:
#                 energy_change.append(np.nan)
#         mutations[k]["ddg"] = energy_change
    
    #------------------------------------------------------------------------------------
    # Calculate the fraction of predicted quasi-null mutations
    #------------------------------------------------------------------------------------    
    
#     mutations = {'nat_mut':naturalMutations, 'dis_mut':diseaseMutations}
#     for k, mut in mutations.items():
#         edgotypes = []
#         for _, row in mut.iterrows():
#             if row.edgotype == 'edgetic':
#                 edgotypes.append('edgetic')
#             elif row.structural_location == 'exposed-noninterface':
#                 edgotypes.append('quasi-wild-type')
#             elif row.structural_location == 'buried':
#                 if np.isnan(row.ddg):
#                     edgotypes.append('-')
#                 elif row.ddg >= qnMinDDG:
#                     edgotypes.append('quasi-null')
#                 else:
#                     edgotypes.append('quasi-wild-type')
#             else:
#                 edgotypes.append('-')
#         mutations[k]['edgotype'] = edgotypes
    
    #------------------------------------------------------------------------------------
    # Newly added
    #------------------------------------------------------------------------------------
    # Assign mutation edgotypes
    #------------------------------------------------------------------------------------
    
    print( '\n' + 'Labeling mutation edgotypes:' )
    print( '%d non-disease mutations' % len(naturalMutations) )
    print( '%d disease mutations' % len(diseaseMutations) )
    
    naturalMutations["edgotype"] = assign_edgotypes (naturalMutations["perturbations"].tolist(),
                                                     mono_edgetic = False)
    diseaseMutations["edgotype"] = assign_edgotypes (diseaseMutations["perturbations"].tolist(),
                                                     mono_edgetic = False)
    
    nat_mono_edgotype = assign_edgotypes (naturalMutations["perturbations"].tolist(), mono_edgetic = True)
    dis_mono_edgotype = assign_edgotypes (diseaseMutations["perturbations"].tolist(), mono_edgetic = True)
    
    if mono_edgetic:
        print( '\n' + 'Labeling mono-edgetic mutations' )
        naturalMutations["mono-edgotype"] = nat_mono_edgotype
        diseaseMutations["mono-edgotype"] = dis_mono_edgotype
    
    # write predicted mutation edgotypes to tab-delimited file
    write_list_table (naturalMutations, ["partners", "perturbations"], natMutEdgotypeFile)
    write_list_table (diseaseMutations, ["partners", "perturbations"], disMutEdgotypeFile)
    
    #------------------------------------------------------------------------------------
    
    naturalMutations = naturalMutations [naturalMutations["edgotype"] != '-'].reset_index(drop=True)
    diseaseMutations = diseaseMutations [diseaseMutations["edgotype"] != '-'].reset_index(drop=True)
    
    naturalMutations.to_csv(natMutEdgotypeFile, index=False, sep='\t')
    diseaseMutations.to_csv(disMutEdgotypeFile, index=False, sep='\t')
    
    numNaturalMut_type = {'QN': sum(naturalMutations["edgotype"] == 'quasi-null'),
                          'E': sum(naturalMutations["edgotype"] == 'edgetic'),
                          'QW': sum(naturalMutations["edgotype"] == 'quasi-wild-type')}
    
    numDiseaseMut_type = {'QN': sum(diseaseMutations["edgotype"] == 'quasi-null'),
                          'E': sum(diseaseMutations["edgotype"] == 'edgetic'),
                          'QW': sum(diseaseMutations["edgotype"] == 'quasi-wild-type')}
    
#     numNaturalMut_considered = sum(numNaturalMut_type.values())
#     numDiseaseMut_considered = sum(numDiseaseMut_type.values())
    
    numNaturalMut_considered = len(naturalMutations)
    numDiseaseMut_considered = len(diseaseMutations)
    
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
    
    #------------------------------------------------------------------------------------
    # Apply Bayes' theorem to calculate the probabilities for a mutation of selected
    # edgotype to be effectively neutral, mildly deleterious or strongly detrimental
    #------------------------------------------------------------------------------------
    
    # Probability for new missense mutations to be neutral (N)
    pN = 0.27

    # Probability for new missense mutations to be mildly deleterious (M)
    pM = 0.53

    # Probability for new missense mutations to be strongly detrimental (S)
    pS = 0.20
    
    # Probability for effectively neutral mutations (N) to be of selected edgotype (T)
    pT_N = fracNaturalMut_type[T]
    
    # Probability for mildly deleterious mutations (M) to be of selected edgotype (T)
    pT_M = fracDiseaseMut_type[T]
    
    # Probability for strongly detrimental mutations (S) to be of selected edgotype (T)
    if assume_S_as_M:
        pT_S = pT_M
    else:
        pT_S = 1 if edgotype is 'quasi-null' else 0
    
    # Probability for a new missense mutation to be of selected edgotype
    pT = (pT_N * pN) + (pT_M * pM) + (pT_S * pS)
    
    # Probability for mutations of selected edgotype to be effectively neutral
    pN_T = pT_N * pN / pT
    
    # Probability for mutations of selected edgotype to be mildly deleterious
    pM_T = pT_M * pM / pT
    
    # Probability for mutations of selected edgotype to be strongly detrimental
    pS_T = pT_S * pS / pT
    
    allresults = {prob:p for prob, p in zip(posteriors, (pN_T, pM_T, pS_T))}
    allresults['P(%s)' % T] = pT
    
    print()
    print('Fitness effect calculation for %s (%s) mutations:' % (edgotype, T))
    print('P(N) = %.1f %%' % (100 * pN))
    print('P(M) = %.1f %%' % (100 * pM))
    print('P(S) = %.1f %%' % (100 * pS))
    print('%s = %.1f %%' % (cond_probs[0], 100 * pT_N))
    print('%s = %.1f %%' % (cond_probs[1], 100 * pT_M))
    print('%s = %.1f %%' % (cond_probs[2], 100 * pT_S))
    print('P(%s) = %sP(N) + %sP(M) + %sP(S) = %.1f %%' 
            % (T, cond_probs[0], cond_probs[1], cond_probs[2], 100 * pT))
    print('Probability for %s mutations to be effectively neutral %s = %.1f %%' 
            % (edgotype, posteriors[0], 100 * pN_T))
    print('Probability for %s mutations to be mildly deleterious %s = %.1f %%' 
            % (edgotype, posteriors[1], 100 * pM_T))
    print('Probability for %s mutations to be strongly detrimental %s = %.1f %%' 
            % (edgotype, posteriors[2], 100 * pS_T))
    
    # calculate 95% confidence interval
    if computeConfidenceIntervals:
        n_N, n_M = numNaturalMut_considered, numDiseaseMut_considered
        k_obs_N, k_obs_M = numNaturalMut_type[T], numDiseaseMut_type[T]
        for prob, p in zip(posteriors, (pN_T, pM_T, pS_T)):
#            if p > 0:
            if prob is posteriors[2] and p == 0:
                p_lower, p_upper = 0, 0
            else:
                if prob is posteriors[0]:
                    pr_lower, pr_upper = proportion_ratio_CI (k_obs_M,
                                                              n_M,
                                                              k_obs_N,
                                                              n_N,
                                                              a = pM / pN,
                                                              b = pT_S * pS / pN,
                                                              conf = CI)
                elif prob is posteriors[1]:
                    pr_lower, pr_upper = proportion_ratio_CI (k_obs_N,
                                                              n_N,
                                                              k_obs_M,
                                                              n_M,
                                                              a = pN / pM,
                                                              b = pT_S * pS / pM,
                                                              conf = CI)
                elif prob is posteriors[2]:
                    pr_lower, pr_upper = proportion_sum_CI (k_obs_N,
                                                            n_N,
                                                            k_obs_M,
                                                            n_M,
                                                            a = pN / (pT_S * pS),
                                                            b = pM / (pT_S * pS),
                                                            conf = CI)
                p_lower = 1 / (1 + pr_upper)
                p_upper = 1 / (1 + pr_lower)
            print( '%.1f%% confidence interval for %s = (%f, %f)' 
                    % (CI, prob, 100 * p_lower, 100 * p_upper) )
#             else:
#                 p_lower, p_upper = np.nan, np.nan
            if (not np.isnan(p_lower)) and (not np.isnan(p_upper)):
                allresults[prob + '_CI'] = [p - p_lower, p_upper - p]
    with open(outputFile, 'wb') as fOut:
        pickle.dump(allresults, fOut)
    
    # plot fitness effect
    prob = [100 * allresults[p] for p in posteriors]
    conf = []
    if computeConfidenceIntervals:
        for p in posteriors:
            if p + '_CI' in allresults:
                lower, upper = allresults[p + '_CI']
                conf.append((100 * lower, 100 * upper))
            else:
                conf.append((0, 0))
    
    curve_plot (prob,
                error = conf if computeConfidenceIntervals else None,
                xlim = [0.8, 3.1],
                ylim = [0, 100],
                styles = '.k',
                capsize = 10 if computeConfidenceIntervals else 0,
                msize = 16,
                ewidth = 1.25,
                ecolors = 'k',
                ylabel = 'Fraction of %s mutations (%%)' % edgotype,
                yMinorTicks = 4,
                xticks = [1, 2, 3],
                xticklabels = ('Effectively\nneutral', 'Mildly\ndeleterious', 'Strongly\ndetrimental'),
                fontsize = 20,
                adjustBottom = 0.2,
                shiftBottomAxis = -0.1,
                xbounds = (1, 3),
                show = showFigs,
                figdir = figDir,
                figname = '%s_mut_fitness_effect' % T)

if __name__ == "__main__":
    main()
