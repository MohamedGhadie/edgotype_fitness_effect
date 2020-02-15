import os
import pandas as pd
from pathlib import Path
from stat_tools import sderror_on_fraction
from plot_tools import pie_plot

def main():
    
    # reference interactome name
    # options: HI-II-14, IntAct
    interactome_name = 'HI-II-14'
    
    # homology modelling method used to create structural models
    # options: template_based, model_based
    model_method = 'model_based'
    
    # method used to perform calculations; 'geometry' or 'physics'
    edgetic_method = 'physics'
    
    # method that was used to calculate edgetic mutation binding ∆∆G
    # options: bindprofx, foldx
    edgetic_ddg = 'foldx'
    
    # maximum RSA for buried mutations
    burialMaxRSA = 0.25
    
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
    edgeticDir = modellingDir / edgetic_method
    
    # figure directory
    figDir = Path('../figures') / interactome_name / model_method / edgetic_method
    
    if edgetic_method is 'physics':
        edgeticDir = edgeticDir / (edgetic_ddg + '_edgetics')
        figDir = figDir / (edgetic_ddg + '_edgetics')
        
    # input data files
    naturalMutationsFile = edgeticDir / 'nondisease_mutation_RSA.txt'
    diseaseMutationsFile = edgeticDir / 'disease_mutation_RSA.txt'
    natMutGeometryFile = modellingDir / 'geometry' / 'nondisease_mutation_edgetics.txt'
    disMutGeometryFile = modellingDir / 'geometry' / 'disease_mutation_edgetics.txt'
    
    # output data files
    natMutLocFile = edgeticDir / 'nondisease_mutation_struc_loc.txt'
    disMutLocFile = edgeticDir / 'disease_mutation_struc_loc.txt'
    
    # create output directories if not existing
    if not edgeticDir.exists():
        os.makedirs(edgeticDir)
    if not figDir.exists():
        os.makedirs(figDir)
    
    #------------------------------------------------------------------------------------
    # load mutations
    #------------------------------------------------------------------------------------
    
    print('Reading processed mutations')
    naturalMutations = pd.read_table (naturalMutationsFile, sep='\t')
    diseaseMutations = pd.read_table (diseaseMutationsFile, sep='\t')
    
    print()
    print('Number of mutations:')
    print('non-disease mutations: %d' % len(naturalMutations))
    print('disease mutations: %d' % len(diseaseMutations))
    
    natMutationProteins = set(naturalMutations["protein"].tolist())
    disMutationProteins = set(diseaseMutations["protein"].tolist())
    mutationProteins = natMutationProteins | disMutationProteins
    
    print()
    print('Number of proteins carrying mutations:')
    print('non-disease mutations: %d' % len(natMutationProteins))
    print('disease mutations: %d' % len(disMutationProteins))
    print('all mutations: %d' % len(mutationProteins))
    
    #------------------------------------------------------------------------------------
    # Identify mutation locations on protein structure
    #------------------------------------------------------------------------------------
    
    #edgetic_region_label = 'interface' if edgetic_method is 'geometry' else 'edgetic'
    #exposed_region_label = 'exposed-noninterface' if edgetic_method is 'geometry' else 'exposed-nonedgetic'
    natMutGeometry = pd.read_table (natMutGeometryFile, sep='\t')
    disMutGeometry = pd.read_table (disMutGeometryFile, sep='\t')
    natMutGeometryEdgotype = {(row.protein, row.mut_position, row.mut_res):row.edgotype
                                for _, row in natMutGeometry.iterrows()}
    disMutGeometryEdgotype = {(row.protein, row.mut_position, row.mut_res):row.edgotype
                                for _, row in disMutGeometry.iterrows()}
    mutGeomEdgotype = {'nat_mut':natMutGeometryEdgotype, 'dis_mut':disMutGeometryEdgotype}
    
    mutations = {'nat_mut':naturalMutations, 'dis_mut':diseaseMutations}
    for k, mut in mutations.items():
        locations = []
        geomEdg = mutGeomEdgotype[k]
        for _, row in mut.iterrows():
            t = row.protein, row.mut_position, row.mut_res
            if row.edgotype == '-':
                locations.append('-')
            elif (row.edgotype == 'edgetic') or (geomEdg[t] == 'edgetic'):
                locations.append('interface')
            elif row.RSA <= burialMaxRSA:
                locations.append('buried')
            else:
                locations.append('exposed-noninterface')
        mut["structural_location"] = locations
    
    naturalMutations.to_csv (natMutLocFile, index=False, sep='\t')
    diseaseMutations.to_csv (disMutLocFile, index=False, sep='\t')
    
    #------------------------------------------------------------------------------------
    # Calculate fraction of mutations in each region on protein structure
    #------------------------------------------------------------------------------------
    
    numNaturalMut_interface = sum(naturalMutations["structural_location"] == 'interface')
    numNaturalMut_exposed = sum(naturalMutations["structural_location"] == 'exposed-noninterface')
    numNaturalMut_buried = sum(naturalMutations["structural_location"] == 'buried')
    
    numDiseaseMut_interface = sum(diseaseMutations["structural_location"] == 'interface')
    numDiseaseMut_exposed = sum(diseaseMutations["structural_location"] == 'exposed-noninterface')
    numDiseaseMut_buried = sum(diseaseMutations["structural_location"] == 'buried')    
    
    numNaturalMut_considered = numNaturalMut_interface + numNaturalMut_exposed + numNaturalMut_buried
    numDiseaseMut_considered = numDiseaseMut_interface + numDiseaseMut_exposed + numDiseaseMut_buried
    
    #------------------------------------------------------------------------------------
    # Print results
    #------------------------------------------------------------------------------------
    
    print()
    print('Fraction of non-disease mutations in each region:')
    print('exposed-noninterface: %f%% (SE = %g, %d out of %d)' 
            % (100 * numNaturalMut_exposed / numNaturalMut_considered,
               sderror_on_fraction (numNaturalMut_exposed, numNaturalMut_considered),
               numNaturalMut_exposed,
               numNaturalMut_considered))
    print('interface: %f%% (SE = %g, %d out of %d)' 
            % (100 * numNaturalMut_interface / numNaturalMut_considered,
               sderror_on_fraction (numNaturalMut_interface, numNaturalMut_considered),
               numNaturalMut_interface,
               numNaturalMut_considered))
    print('buried: %f%% (SE = %g, %d out of %d)' 
            % (100 * numNaturalMut_buried / numNaturalMut_considered,
               sderror_on_fraction (numNaturalMut_buried, numNaturalMut_considered),
               numNaturalMut_buried,
               numNaturalMut_considered))
    
    print()
    print('Fraction of disease mutations in each region:')
    print('exposed-noninterface: %f%% (SE = %g, %d out of %d)' 
            % (100 * numDiseaseMut_exposed / numDiseaseMut_considered,
               sderror_on_fraction (numDiseaseMut_exposed, numDiseaseMut_considered),
               numDiseaseMut_exposed,
               numDiseaseMut_considered))
    print('interface: %f%% (SE = %g, %d out of %d)' 
            % (100 * numDiseaseMut_interface / numDiseaseMut_considered,
               sderror_on_fraction (numDiseaseMut_interface, numDiseaseMut_considered),
               numDiseaseMut_interface,
               numDiseaseMut_considered))
    print('buried: %f%% (SE = %g, %d out of %d)' 
            % (100 * numDiseaseMut_buried / numDiseaseMut_considered,
               sderror_on_fraction (numDiseaseMut_buried, numDiseaseMut_considered),
               numDiseaseMut_buried,
               numDiseaseMut_considered))
    
    pie_plot ([numNaturalMut_exposed, numNaturalMut_interface, numNaturalMut_buried],
              angle = 90,
              #labels = [exposed_region_label.title(), edgetic_region_label.title(), 'Buried'],
              #labels = ['exposed-noninterface', 'interface', 'Buried'],
              colors = ['mediumslateblue', 'purple', 'red'],
              edgewidth = 2,
              show = showFigs,
              figdir = figDir,
              figname = 'non_disease_mutation_structural_locations')
    pie_plot ([numDiseaseMut_exposed, numDiseaseMut_interface, numDiseaseMut_buried],
              angle = 90,
              #labels = [exposed_region_label.title(), edgetic_region_label.title(), 'Buried'],
              #labels = ['exposed-noninterface', 'interface', 'Buried'],
              colors = ['mediumslateblue', 'purple', 'red'],
              edgewidth = 2,
              show = showFigs,
              figdir = figDir,
              figname = 'disease_mutation_structural_locations')

if __name__ == "__main__":
    main()
