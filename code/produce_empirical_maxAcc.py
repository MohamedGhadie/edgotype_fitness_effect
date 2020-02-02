import pickle
from pathlib import Path
from id_mapping import produce_chain_dict
from pdb_tools import produce_chain_list
from threeD_structure_tools import produce_empirical_maxAcc

def main():
    
    # parent directory of all data files
    dataDir = Path('../data')
    
    # parent directory of all processed data files
    procDir = dataDir / 'processed'
    
    # directory of PDB structures
    pdbDir = Path('/Volumes/MG_Samsung/pdb_files')
    #pdbDir = Path('../../pdb_files')
    
    # directory of precalculated RSA files on local computer
    dsspDir = Path('/Volumes/MG_Samsung/dssp')
    #dsspDir = Path('../../dssp')
    
    # input data files
    pdbSeqresFile = extDir / 'pdb_seqres.txt'
    
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
