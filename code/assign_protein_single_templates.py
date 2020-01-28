import os
import pandas as pd
from pathlib import Path
from structural_annotation import single_chain_per_protien
from pdb_tools import download_structures

def main():
    
    # maximum e-value cutoff to filter out protein-chain annotations
    #evalue = 1e-10
    
    # minimum protein sequence coverage fraction required for protein-chain annotation
    #prCov = 0.5
    
    # minimum chain sequence coverage fraction required for protein-chain annotation
    #chCov = 0.5
    
    # download missing template structures from PDB
    download_template_structures = True
    
    # parent directory of all data files
    dataDir = Path('../data')
    
    # parent directory of all processed data files
    procDir = dataDir / 'processed'
    
    # directory of PDB structure files
    pdbDir = Path('/Volumes/MG_Samsung/pdb_files')
    
    # directory of processed data files specific to interactome
    interactomeDir = procDir / interactome_name
      
    # input data files
    #chainMapFile = procDir / 'human_pdb_chain_map.txt'
    chainMapFile = interactomeDir / 'struc_interactome_pdb_chain_map.txt'
    
    # output data files
    #filteredChainMapFile = interactomeDir / 'human_pdb_chain_map_filtered.txt'
    singleChainMapFile = interactomeDir / 'human_single_chain_map_per_protein.txt'
    chainIDFile = interactomeDir / 'unique_template_chainIDs.txt'
    pdbIDFile = interactomeDir / 'unique_template_pdbIDs.txt'
    
    # create output directories if not existing
    if not interactomeDir.exists():
        os.makedirs(interactomeDir)
    if not pdbDir.exists():
        os.makedirs(pdbDir)

    #------------------------------------------------------------------------------------
    # select one chain model per protein
    #------------------------------------------------------------------------------------
    
#     if not filteredChainMapFile.is_file():
#         print('filtering chain annotations')
#         filter_chain_annotations (chainMapFile,
#                                   filteredChainMapFile,
#                                   evalue = evalue,
#                                   prCov = prCov,
#                                   chCov = chCov)
    
    if not singleChainMapFile.is_file():
        print('selecting one chain model per protein')
        single_chain_per_protien (chainMapFile, singleChainMapFile)
    
    proteinModels = pd.read_table (singleChainMapFile, sep='\t')
    uniqueChains = set(proteinModels["Subject"].values)
    uniquePDBs = {id.split('_')[0] for id in uniqueChains}
    print('%d unique chains in %d unique PDB structures' % (len(uniqueChains), len(uniquePDBs)))
    
    with open(chainIDFile, 'w') as f:
        for i in sorted(uniqueChains):
            f.write("%s\n" % i)
    with open(pdbIDFile, 'w') as f:
        for i in sorted(uniquePDBs):
            f.write("%s\n" % i)
    
    if download_template_structures:
        print('downloading structures for selected protein models')
        download_structures (pdbIDFile, pdbDir)

if __name__ == "__main__":
    main()
