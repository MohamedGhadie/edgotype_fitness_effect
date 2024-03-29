#----------------------------------------------------------------------------------------
# Produce dictionary of empirically calculated residue maximum solvent accessibility.
#----------------------------------------------------------------------------------------

import pickle
from pathlib import Path
from pdb_tools import produce_chain_list, produce_chain_dict
from threeD_structure_tools import produce_empirical_maxAcc

def main():
    
    # parent directory of all data files
    dataDir = Path('../data')
    
    # parent directory of all processed data files
    procDir = dataDir / 'processed'
    
    # local directory of precalculated RSA files
    dsspDir = Path('../../dssp')
    
    # input data files
    pdbSeqresFile = procDir / 'pdb_seqres_reduced.fasta'
    
    # output data files
    chainListFile = procDir / 'pdb_chains.list'
    pdbChainsFile = procDir / 'pdb_chains.pkl'
    maxAccFile = procDir / 'empirical_maxAcc_99_99.pkl'
    
    # create output directories if not existing
    if not procDir.exists():
        os.makedirs(procDir)
    
    #------------------------------------------------------------------------------------
    # Produce empirical maximum solvent accessibility dictionary
    #------------------------------------------------------------------------------------    
    
    if not chainListFile.is_file():
        print('producing PDB chain ID file from fasta records')
        produce_chain_list (pdbSeqresFile, chainListFile)
    
    if not pdbChainsFile.is_file():
        print('producing PDB chain dictionary from chain list file')
        produce_chain_dict (chainListFile, pdbChainsFile)
    
    if not maxAccFile.is_file():
        print('producing empirical maximum solvent accessibility dictionary')
        with open(pdbChainsFile, 'rb') as f:
            pdbChains = pickle.load(f)
        pdbIDs = sorted(pdbChains.keys())
        produce_empirical_maxAcc (pdbIDs, dsspDir, maxAccFile)

if __name__ == "__main__":
    main()
