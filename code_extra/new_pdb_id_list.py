import io
from pathlib import Path

def main():
    
    # parent directory of all data files
    dataDir = Path('../data')
    
    # parent directory of all processed data files
    procDir = dataDir / 'processed'
    
    # directory of PDB structures
    pdbDir = Path('/Volumes/MG_Samsung/pdb_files')
    
    # input data files
    inPath = procDir / 'HI-II-14' / 'template_based' / 'interactome_pdbIDs.txt'
    compareToPath = procDir / 'IntAct' / 'template_based' / 'interactome_pdbIDs.txt'
    
    # output data files
    outPath = procDir / 'HI-II-14' / 'template_based' / 'new_pdbIDs.txt'
    
    with open(inPath, 'r') as f:
        inIDs = set(f.read().split())
    
    with open(compareToPath, 'r') as f:
        compareIDs = set(f.read().split())
    
    newIDs = inIDs - compareIDs
    
    nofile = []
    with io.open(outPath, 'w') as fout:
        for id in newIDs:
            filename = pdbDir / ('pdb' + id + '.ent')
            if filename.is_file():
                fout.write(id + '\n')
            else:
                nofile.append(id)
    print('New PDB IDs with no structure file available:')
    print(', '.join(nofile))

if __name__ == "__main__":
    main()
