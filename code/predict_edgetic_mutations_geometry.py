#----------------------------------------------------------------------------------------
# Predict interactome perturbations based on geometry. Given a structural interactome 
# with PPI binding interface annotations, common neutral mutations not associated 
# with disease as well as Mendelian disease-causing mutations are mapped onto the 
# structural interactome. Then, PPI perturbations caused by mutations are predicted based 
# on mutation location relative to interaction interface.
#
# Run the following scripts before running this script:
# - produce_structural_interactome.py
# - process_dbsnp_mutations.py
# - process_clinvar_mutations.py
#----------------------------------------------------------------------------------------

import os
import pickle
import numpy as np
from pathlib import Path
from text_tools import write_list_table
from interactome_tools import read_interface_annotated_interactome
from mutation_processing_tools import remove_mutation_overlaps
from mutation_interface_edgotype import (mutation_PPI_interface_perturbations, 
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
    methodDir = modellingDir / 'geometry'
    
    # directory of network perturbation output data files for use by Cytoscape
    cytoscapeDir = methodDir / 'cytoscape'
        
    # figure directory
    figDir = Path('../figures') / interactome_name / model_method / 'geometry'
    
    # input data files
    naturalMutationsFile = procDir / 'dbsnp_mutations4.txt'
    diseaseMutationsFile = procDir / 'clinvar_mutations6.txt'
    structuralInteractomeFile = modellingDir / 'structural_interactome.txt'
    
    # output data files
    #uniqueMutationPerturbsFile = interactomeDir / 'unique_mutation_perturbs_geometry.pkl'
    natMutEdgotypeFile = methodDir / 'nondisease_mutation_edgetics.txt'
    disMutEdgotypeFile = methodDir / 'disease_mutation_edgetics.txt'
    naturalMutEdgeFile = cytoscapeDir / 'nondiseaseMut_perturbed_edges'
    naturalMutNodeFile = cytoscapeDir / 'nondiseaseMut_node_colors'
    diseaseMutEdgeFile = cytoscapeDir / 'diseaseMut_perturbed_edges'
    diseaseMutNodeFile = cytoscapeDir / 'diseaseMut_node_colors'
    
    # create output directories if not existing
    if not methodDir.exists():
        os.makedirs(methodDir)
    if not cytoscapeDir.exists():
        os.makedirs(cytoscapeDir)
    if not figDir.exists():
        os.makedirs(figDir)
    
    #------------------------------------------------------------------------------------
    # further process mutations
    #------------------------------------------------------------------------------------
    
    naturalMutations, diseaseMutations = remove_mutation_overlaps (naturalMutationsFile,
                                                                   diseaseMutationsFile)
    
    #------------------------------------------------------------------------------------
    # predict PPI perturbations
    #------------------------------------------------------------------------------------
    
    # Consider PPI perturbations only for PPIs with this maximum number of interfaces 
    # mapped from distinct PDB binding chain pairs.
    # This parameter is irrelevent if the flag "merge_interfaces" is True.
    # Set to inf for unlimited number of interfaces.
    maxInterfaces = np.inf
    
    # Predict PPI perturbation if mutation is this number of residues away in sequence 
    # from an interface residue. Set to 0 if mutation must be exactly at interface residue.
    numResFromInterface = 0
    
    # Consider PPI perturbations only for PPIs with this minimum number of partners
    minPartners = 1
    
    structuralInteractome = read_interface_annotated_interactome (structuralInteractomeFile)
    
    print( '\n' + 'Predicting PPI perturbations by non-disease mutations based on geometry' )
    naturalMutation_perturbs = mutation_PPI_interface_perturbations (naturalMutations,
                                                                     structuralInteractome,
                                                                     maxInterfaces = maxInterfaces,
                                                                     dist = numResFromInterface)
    print( '\n' + 'Predicting PPI perturbations by disease mutations based on geometry' )
    diseaseMutation_perturbs = mutation_PPI_interface_perturbations (diseaseMutations,
                                                                     structuralInteractome,
                                                                     maxInterfaces,
                                                                     dist = numResFromInterface)
    
    naturalMutations["partners"], naturalMutations["perturbations"] = zip(* naturalMutation_perturbs)
    naturalMutations = naturalMutations[naturalMutations["partners"].apply(len) >= minPartners]
    
    diseaseMutations["partners"], diseaseMutations["perturbations"] = zip(* diseaseMutation_perturbs)
    diseaseMutations = diseaseMutations[diseaseMutations["partners"].apply(len) >= minPartners]
    
    # drop duplicate mutations based on location, regardless of residue type 
    naturalMutations = naturalMutations.drop_duplicates(subset=["protein", "mut_position"]).reset_index(drop=True)
    diseaseMutations = diseaseMutations.drop_duplicates(subset=["protein", "mut_position"]).reset_index(drop=True)
    
    print( '\n' + 'Number of mutations with PPI perturbation predictions after removing duplicates by position' )
    print( 'non-disease: %d' % len(naturalMutations) )
    print( 'disease: %d' % len(diseaseMutations) )
    
    #------------------------------------------------------------------------------------
    # Assign mutation edgotypes
    #------------------------------------------------------------------------------------
    
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
    
    # write predicted mutation edgotypes to file
    write_list_table (naturalMutations, ["partners", "perturbations"], natMutEdgotypeFile)
    write_list_table (diseaseMutations, ["partners", "perturbations"], disMutEdgotypeFile)
    
#     with open(uniqueMutationPerturbsFile, 'wb') as fOut:
#         pickle.dump([naturalMutations, diseaseMutations], fOut)
    
    #------------------------------------------------------------------------------------
    # plot network perturbations
    #------------------------------------------------------------------------------------
    
    if plot_perturbations:
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
                      figname = 'nondisease_mut_perturbed_interactome')
    
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
                      figname = 'disease_mut_perturbed_interactome')

if __name__ == "__main__":
    main()
