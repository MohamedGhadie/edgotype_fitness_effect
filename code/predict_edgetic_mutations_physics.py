#----------------------------------------------------------------------------------------
# Predict edgetic mutations based on physics.
#----------------------------------------------------------------------------------------

import os
from pathlib import Path
from text_tools import read_list_table, write_list_table
from interactome_tools import read_single_interface_annotated_interactome
from energy_tools import read_protein_mutation_ddg
from mutation_interface_edgotype import (energy_based_perturbation,
                                         assign_edgotypes,
                                         create_perturbed_network)
from stat_tools import sderror_on_fraction, fisher_test
from plot_tools import network_plot

def main():
    
    # reference interactome name: HI-II-14, HuRI, IntAct
    interactome_name = 'HuRI'
    
    # homology modelling method used to create structural models
    # options: template_based, model_based
    model_method = 'model_based'
    
    # method of calculating mutation ∆∆G for which results will be used
    # options: bindprofx, foldx
    ddg_method = 'foldx'
    
    # minimum reduction in binding free energy ∆∆G required for PPI disruption
    ddgCutoff = 0.5
    
    # if True predict mono-edgetic mutations instead of edgetic mutations
    mono_edgetic = False
    
    # plot perturbed interactome and produce files for use by Cytoscape
    plot_perturbations = True
    
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
    edgeticDir = modellingDir / 'physics' / (ddg_method + '_edgetics')
    
    # directory of network perturbation output data files for use by Cytoscape
    cytoscapeDir = edgeticDir / 'cytoscape'
    
    # figure directory
    figDir = Path('../figures') / interactome_name / model_method / 'physics' / (ddg_method + '_edgetics')
    
    # input data files
    structuralInteractomeFile = modellingDir / 'structural_interactome.txt'
    natMutGeomEdgotypeFile = modellingDir / 'geometry' / 'nondisease_mutation_edgetics.txt'
    disMutGeomEdgotypeFile = modellingDir / 'geometry' / 'disease_mutation_edgetics.txt'
    natMutDDGFile = modellingDir / ('nondis_mut_binding_ddg_%s.txt' % ddg_method)
    disMutDDGFile = modellingDir / ('dis_mut_binding_ddg_%s.txt' % ddg_method)
    
    # output data files
    natMutEdgotypeFile = edgeticDir / 'nondisease_mutation_edgetics.txt'
    disMutEdgotypeFile = edgeticDir / 'disease_mutation_edgetics.txt'
    naturalMutEdgeFile = cytoscapeDir / 'nondiseaseMut_perturbed_edges'
    naturalMutNodeFile = cytoscapeDir / 'nondiseaseMut_node_colors'
    diseaseMutEdgeFile = cytoscapeDir / 'diseaseMut_perturbed_edges'
    diseaseMutNodeFile = cytoscapeDir / 'diseaseMut_node_colors'
    
    # create output directories if not existing
    if not edgeticDir.exists():
        os.makedirs(edgeticDir)
    if not cytoscapeDir.exists():
        os.makedirs(cytoscapeDir)
    if not figDir.exists():
        os.makedirs(figDir)
    
    #------------------------------------------------------------------------------------
    # Read mutation geometry-based edgetic perturbations and mutation ∆∆G
    #------------------------------------------------------------------------------------
    
    naturalMutations = read_list_table (natMutGeomEdgotypeFile,
                                        ["partners", "perturbations"],
                                        [str, int])
    diseaseMutations = read_list_table (disMutGeomEdgotypeFile,
                                        ["partners", "perturbations"],
                                        [str, int])
    
    # read change in binding free energy for interfacial mutations
    naturalMutationsDDG = read_protein_mutation_ddg (natMutDDGFile, 'binding')
    diseaseMutationsDDG = read_protein_mutation_ddg (disMutDDGFile, 'binding')
    
    #------------------------------------------------------------------------------------
    # predict PPI perturbations based on physics
    #------------------------------------------------------------------------------------
    
    print( '\n' + 'Performing physics-based edgotype prediction for non-disease mutations' )
    naturalMutations["perturbations"], knownDDG, unknownDDG = energy_based_perturbation (naturalMutations,
                                                                                         naturalMutationsDDG,
                                                                                         ddgCutoff)
    print( '\n' + 'Performing physics-based edgotype prediction for disease mutations' )
    diseaseMutations["perturbations"], knownDDG, unknownDDG = energy_based_perturbation (diseaseMutations,
                                                                                         diseaseMutationsDDG,
                                                                                         ddgCutoff)
    
    #------------------------------------------------------------------------------------
    # Assign mutation edgotypes
    #------------------------------------------------------------------------------------
    
    print('\n' + 'Labeling mutation edgotypes')
    print('%d non-disease mutations' % len(naturalMutations))
    print('%d disease mutations' % len(diseaseMutations))
    
    naturalMutations["edgotype"] = assign_edgotypes (naturalMutations["perturbations"].tolist(),
                                                     mono_edgetic = False)
    diseaseMutations["edgotype"] = assign_edgotypes (diseaseMutations["perturbations"].tolist(),
                                                     mono_edgetic = False)
    
    nat_mono_edgotype = assign_edgotypes (naturalMutations["perturbations"].tolist(), mono_edgetic = True)
    dis_mono_edgotype = assign_edgotypes (diseaseMutations["perturbations"].tolist(), mono_edgetic = True)
    
    if mono_edgetic:
        print('\n' + 'Labeling mono-edgetic mutations')
        naturalMutations["mono-edgotype"] = nat_mono_edgotype
        diseaseMutations["mono-edgotype"] = dis_mono_edgotype
    else:
        if "mono-edgotype" in naturalMutations.columns.values:
            naturalMutations = naturalMutations.drop("mono-edgotype", axis=1)
        if "mono-edgotype" in diseaseMutations.columns.values:
            diseaseMutations = diseaseMutations.drop("mono-edgotype", axis=1)
    
    naturalMutations = naturalMutations [naturalMutations["edgotype"] != '-'].reset_index(drop=True)
    diseaseMutations = diseaseMutations [diseaseMutations["edgotype"] != '-'].reset_index(drop=True)
    
    if mono_edgetic:
        numNaturalMut_edgetic = sum(naturalMutations["mono-edgotype"] == 'mono-edgetic')
        numNaturalMut_nonedgetic = sum(naturalMutations["mono-edgotype"].apply(lambda x: 
                                                        x in ('non-edgetic', 'edgetic')))
        numDiseaseMut_edgetic = sum(diseaseMutations["mono-edgotype"] == 'mono-edgetic')
        numDiseaseMut_nonedgetic = sum(diseaseMutations["mono-edgotype"].apply(lambda x: 
                                                        x in ('non-edgetic', 'edgetic')))
    else:
        numNaturalMut_edgetic = sum(naturalMutations["edgotype"] == 'edgetic')
        numNaturalMut_nonedgetic = sum(naturalMutations["edgotype"] == 'non-edgetic')
        numDiseaseMut_edgetic = sum(diseaseMutations["edgotype"] == 'edgetic')
        numDiseaseMut_nonedgetic = sum(diseaseMutations["edgotype"] == 'non-edgetic')
    
    numNaturalMut_considered = numNaturalMut_edgetic + numNaturalMut_nonedgetic
    numDiseaseMut_considered = numDiseaseMut_edgetic + numDiseaseMut_nonedgetic
    
    label = 'monoedgetic' if mono_edgetic else 'edgetic'
    print( '\n' + 'Fraction of predicted %s mutations:' % label )
    print('Non-disease mutations: %.3f (SE = %g, %d out of %d)' 
            % (numNaturalMut_edgetic / numNaturalMut_considered,
               sderror_on_fraction (numNaturalMut_edgetic, numNaturalMut_considered),
               numNaturalMut_edgetic,
               numNaturalMut_considered))
    
    print('Disease mutations: %.3f (SE = %g, %d out of %d)' 
            % (numDiseaseMut_edgetic / numDiseaseMut_considered,
               sderror_on_fraction (numDiseaseMut_edgetic, numDiseaseMut_considered),
               numDiseaseMut_edgetic,
               numDiseaseMut_considered))
    
    fisher_test ([numNaturalMut_edgetic, numNaturalMut_nonedgetic],
                 [numDiseaseMut_edgetic, numDiseaseMut_nonedgetic])
    
    # write predicted mutation edgotypes to tab-delimited file
    write_list_table (naturalMutations, ["partners", "perturbations"], natMutEdgotypeFile)
    write_list_table (diseaseMutations, ["partners", "perturbations"], disMutEdgotypeFile)
    
    #------------------------------------------------------------------------------------
    # plot network perturbations
    #------------------------------------------------------------------------------------
    
    if plot_perturbations:
        structuralInteractome = read_single_interface_annotated_interactome (structuralInteractomeFile)
        
        print( '\n' + 'Creating network perturbed by non-disease mutations' )
        nodes, edges, nodeColors, edgeColors = create_perturbed_network (structuralInteractome,
                                                                         naturalMutations,
                                                                         naturalMutEdgeFile,
                                                                         naturalMutNodeFile)
        network_plot (edges,
                      nodes = nodes,
                      nodeSizes = [20] * len(nodes),
                      edgeWidth = 1,
                      nodeColors = nodeColors,
                      edgeColors = edgeColors,
                      show = showFigs,
                      figdir = figDir,
                      figname = 'nondisease_mut_perturbed_interactome_%s' % ddg_method)
    
        print( '\n' + 'Creating network perturbed by disease mutations' )
        nodes, edges, nodeColors, edgeColors = create_perturbed_network (structuralInteractome,
                                                                         diseaseMutations,
                                                                         diseaseMutEdgeFile,
                                                                         diseaseMutNodeFile)
        network_plot (edges,
                      nodes = nodes,
                      nodeSizes = [20] * len(nodes),
                      edgeWidth = 1,
                      nodeColors = nodeColors,
                      edgeColors = edgeColors,
                      show = showFigs,
                      figdir = figDir,
                      figname = 'disease_mut_perturbed_interactome_%s' % ddg_method)

if __name__ == "__main__":
    main()
