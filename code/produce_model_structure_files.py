import os
import sys
from pathlib import Path
from interactome_tools import read_single_interface_annotated_interactome
from pdb_tools import write_partial_structure

def main():
    
    # reference interactome name
    # options: HI-II-14, IntAct
    interactome_name = 'HI-II-14'
    
    # parent directory of all data files
    dataDir = Path('../data')
    
    # parent directory of all processed data files
    procDir = dataDir / 'processed'
    
    # directory of processed data files specific to interactome
    interactomeDir = procDir / interactome_name
        
    # directory for PDB structure files
    pdbDir = Path('/Volumes/MG_Samsung/pdb_files')
    
    # directory for template structure files
    templateDir = Path('../templates')
    
    # input data files
    interactomeFile = interactomeDir / 'human_interface_annotated_interactome.txt'
    
    # create output directories if not existing
    if not dataDir.exists():
        os.makedirs(dataDir)
    if not procDir.exists():
        os.makedirs(procDir)
    if not interactomeDir.exists():
        os.makedirs(interactomeDir)
    if not pdbDir.exists():
        os.makedirs(pdbDir)
    if not templateDir.exists():
        os.makedirs(templateDir)
    
    interactome = read_single_interface_annotated_interactome (interactomeFile)
    n = len(interactome)
    
    for i, row in interactome.iterrows():
        sys.stdout.write('  PPI %d out of %d (%.2f%%) \r' % (i+1, n, 100*(i+1)/n))
        sys.stdout.flush()
        for chain1, chain2 in row.Chain_pairs:
            pdbid, chainID1 = chain1.split('_')
            _, chainID2 = chain2.split('_')
            modelID = '_'.join([pdbid, chainID1, chainID2])
            write_partial_structure (pdbid,
                                     [chainID1, chainID2],
                                     pdbDir,
                                     templateDir / (modelID + '.ent'))

if __name__ == "__main__":
    main()
