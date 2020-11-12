#----------------------------------------------------------------------------------------
# Modules for mathematical computations.
#----------------------------------------------------------------------------------------

import numpy as np
from stat_tools import proportion_ratio_CI, proportion_sum_CI

def fitness_effect (pN,
                    pM,
                    pS,
                    k_N,
                    n_N,
                    k_M,
                    n_M,
                    pT_S = None,
                    edgotype = 'edgetic',
                    CI = 95,
                    output = True):
    """Calculate the fitness effect probabilities with confidence intervals for a given edgotype T.

    Args:
        pN (numeric): prior probability for mutations to be effectively neutral.
        pM (numeric): prior probability for mutations to be mildly deleterious.
        pS (numeric): prior probability for mutations to be strongly detrimental.
        k_N (numeric): number of neutral mutations with edgotype T.
        n_N (numeric): total number of neutral mutations.
        k_M (numeric): number of mildly deleterious mutations with edgotype T.
        n_M (numeric): total number of mildly deleterious mutations.
        pT_S (numeric): probability for strongly detrimental mutations to have edgotype T.
        edgotype (str): edgotype for which to calculate fitness effect, 
                        either 'edgetic', 'quasi-null' or 'quasi-wild-type'.
        CI (numeric): percent confidence interval to be calculated for each fitness effect probability.
        output (bool): print output if True.

    Returns:
        dict

    """
    # abbreviations for different structural regions
    edgotype_symbol = {'quasi-null':'Q', 'edgetic':'E', 'quasi-wild-type':'W'}
    
    # abbreviation for selected structural region
    T = edgotype_symbol [edgotype]
    
    # fitness effects
    fitness_effects = ['Effectively neutral', 'Mildly deleterious', 'Strongly detrimental']
    
    # Probability for effectively neutral mutations (N) to be of selected edgotype (T)
    pT_N = k_N / n_N
    
    # Probability for mildly deleterious mutations (M) to be of selected edgotype (T)
    pT_M = k_M / n_M
    
    # Assume strongly detrimental mutations are similar to mildly deleterious mutations
    if pT_S is None:
        pT_S = pT_M
    
    # Probability for a new missense mutation to be of selected edgotype
    pT = (pT_N * pN) + (pT_M * pM) + (pT_S * pS)
    
    if pT == 0:
        return {}
    else:
        # Probability for mutations of selected edgotype to be effectively neutral
        pN_T = pT_N * pN / pT
    
        # Probability for mutations of selected edgotype to be mildly deleterious
        pM_T = pT_M * pM / pT
    
        # Probability for mutations of selected edgotype to be strongly detrimental
        pS_T = pT_S * pS / pT
    
        #allresults = {prob:p for prob, p in zip(posteriors, (pN_T, pM_T, pS_T))}
        #allresults['P(%s)' % T] = pT
        allresults = {'FE':{f:p for f, p in zip(fitness_effects, (pN_T, pM_T, pS_T))}}
    
        if output:
            # prior conditional probabilities
            cond_probs = ['P(%s|%s)' % (T, p) for p in ('N','M','S')]
    
            # posterior probabilities
            posteriors = ['P(%s|%s)' % (p, T) for p in ('N','M','S')]
            
            print()
            print('Fitness effect calculation for %s (%s) mutations:' % (edgotype, T))
            print('P(N) = %.1f %%' % (100 * pN))
            print('P(M) = %.1f %%' % (100 * pM))
            print('P(S) = %.1f %%' % (100 * pS))
            print('P(%s|N) = %.1f %%' % (T, 100 * pT_N))
            print('P(%s|M) = %.1f %%' % (T, 100 * pT_M))
            print('P(%s|S) = %.1f %%' % (T, 100 * pT_S))
            print('P(%s) = P(%s|N)P(N) + P(%s|M)P(M) + P(%s|S)P(S) = %.2f %%' 
                    % (T, T, T, T, 100 * pT))
            print('Probability for %s mutations (%s) to be effectively neutral (N): P(N|%s) = %.2f %%' 
                    % (edgotype, T, T, 100 * pN_T))
            print('Probability for %s mutations (%s) to be mildly deleterious (M): P(M|%s) = %.2f %%' 
                    % (edgotype, T, T, 100 * pM_T))
            print('Probability for %s mutations (%s) to be strongly detrimental (S): P(S|%s) = %.2f %%' 
                    % (edgotype, T, T, 100 * pS_T))
    
        # calculate 95% confidence interval
        if CI:
            allresults['CI'] = {}
            for f, s, p in zip(fitness_effects, ('N','M','S'), (pN_T, pM_T, pS_T)):
                if f is 'Effectively neutral':
                    p_lower, p_upper = proportion_ratio_CI (k_M,
                                                            n_M,
                                                            k_N,
                                                            n_N,
                                                            a = pM / pN,
                                                            b = pT_S * pS / pN,
                                                            conf = CI)
                elif f is 'Mildly deleterious':
                    p_lower, p_upper = proportion_ratio_CI (k_N,
                                                            n_N,
                                                            k_M,
                                                            n_M,
                                                            a = pN / pM,
                                                            b = pT_S * pS / pM,
                                                            conf = CI)
                elif f is 'Strongly detrimental':
                    if pT_S > 0:
                        p_lower, p_upper = proportion_sum_CI (k_N,
                                                              n_N,
                                                              k_M,
                                                              n_M,
                                                              a = pN / (pT_S * pS),
                                                              b = pM / (pT_S * pS),
                                                              conf = CI)
                    else:
                        p_lower, p_upper = np.nan, np.nan
                        allresults['CI'][f] = [p_lower, p_upper]
                
                if not (np.isnan(p_lower) or np.isnan(p_upper)):
                    p_lower, p_upper = 1 / (1 + p_upper), 1 / (1 + p_lower)
                    allresults['CI'][f] = [p - p_lower, p_upper - p]
                
                if output:
                    print('%.1f%% confidence interval for P(%s|%s) = (%f, %f)' 
                            % (CI, s, T, 100 * p_lower, 100 * p_upper))
        
        return allresults
