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
import pickle
from pathlib import Path
from text_tools import read_list_table, write_list_table
from interactome_tools import read_single_interface_annotated_interactome
from ddg_tools import read_protein_mutation_ddg
from mutation_interface_edgotype import (energy_based_perturbation,
                                         assign_edgotypes,
                                         create_perturbed_network)
from plot_tools import network_plot

def main():
    
    # reference interactome name
    # options: HI-II-14, IntAct
    interactome_name = 'HI-II-14'
    
    # homology modelling method used to create structural models
    # options: template_based, model_based
    model_method = 'model_based'
    
    # method of calculating mutation ∆∆G for which results will be used
    # options: bindprofx, foldx
    ddg_method = 'foldx'
    
    # Minimum reduction in binding free energy DDG required for interaction perturbation
    ddgCutoff = 0.5
    
    # set to True to calculate dispensable PPI content using fraction of mono-edgetic mutations 
    # instead of edgetic mutations
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
    
    # directory of calculation method
    methodDir = modellingDir / 'physics'
    
    # directory of network perturbation output data files for use by Cytoscape
    cytoscapeDir = methodDir / 'cytoscape'
    
    # figure directory
    figDir = Path('../figures') / interactome_name / model_method / 'physics' / ('%s_edgetics' % ddg_method)
    
    # input data files
    #geometryPerturbsFile = interactomeDir / 'unique_mutation_perturbs_geometry.pkl'
    structuralInteractomeFile = modellingDir / 'human_structural_interactome.txt'
    natMutGeomEdgotypeFile = modellingDir / 'geometry' / 'nondisease_mutation_edgetics.txt'
    disMutGeomEdgotypeFile = modellingDir / 'geometry' / 'disease_mutation_edgetics.txt'
    natMutDDGFile = modellingDir / ('nondis_mut_binding_ddg_%s.txt' % ddg_method)
    disMutDDGFile = modellingDir / ('dis_mut_binding_ddg_%s.txt' % ddg_method)
    
    # output data files
    #physicsPerturbsFile = interactomeDir / ('mutation_perturbs_physics_%s.pkl' % ddg_method)
    natMutEdgotypeFile = methodDir / ('nondisease_mutation_edgetics_%s.txt' % ddg_method)
    disMutEdgotypeFile = methodDir / ('disease_mutation_edgetics_%s.txt' % ddg_method)
    naturalMutEdgeFile = cytoscapeDir / ('nondiseaseMut_perturbed_edges_%s' % ddg_method)
    naturalMutNodeFile = cytoscapeDir / ('nondiseaseMut_node_colors_%s' % ddg_method)
    diseaseMutEdgeFile = cytoscapeDir / ('diseaseMut_perturbed_edges_%s' % ddg_method)
    diseaseMutNodeFile = cytoscapeDir / ('diseaseMut_node_colors_%s' % ddg_method)
    
    # create output directories if not existing
    if not methodDir.exists():
        os.makedirs(methodDir)
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
    # predict PPI perturbations
    #------------------------------------------------------------------------------------

#     if geometryPerturbsFile.is_file():
#         print( '\n' + 'Loading geometry-based PPI perturbation predictions' )
#         with open(geometryPerturbsFile, 'rb') as f:
#             naturalPerturbs, diseasePerturbs = pickle.load(f)
#     else:
#         print( '\n' + 'Geometry-based PPI perturbation prediction file not found' )
#         return

    print( '\n' + 'Performing physics-based edgotype prediction for non-disease mutations' )
    naturalMutations["perturbations"], knownDDG, unknownDDG = energy_based_perturbation (naturalMutations,
                                                                                         naturalMutationsDDG,
                                                                                         ddgCutoff)
    print( '\n' + 'Performing physics-based edgotype prediction for disease mutations' )
    diseaseMutations["perturbations"], knownDDG, unknownDDG = energy_based_perturbation (diseaseMutations,
                                                                                         diseaseMutationsDDG,
                                                                                         ddgCutoff)
#     with open(physicsPerturbsFile, 'wb') as fOut:
#         pickle.dump([naturalMutations, diseaseMutations], fOut)
    
    #------------------------------------------------------------------------------------
    # Assign mutation edgotypes
    #------------------------------------------------------------------------------------
    
    print( '\n' + 'Labeling mutation edgotypes' )
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
    else:
        if "mono-edgotype" in naturalMutations.columns.values:
            naturalMutations = naturalMutations.drop("mono-edgotype", axis=1)
        if "mono-edgotype" in diseaseMutations.columns.values:
            diseaseMutations = diseaseMutations.drop("mono-edgotype", axis=1)
    
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
