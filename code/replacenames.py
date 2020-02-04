import io
from pathlib import Path

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
    
    # directory of processed model-related data files specific to interactome
    modelBasedDir = interactomeDir / 'model_based'
    
    # input data files
    #inPath = modelBasedDir / 'human_structural_interactome_withDuplicates.txt'
    #inPath = modelBasedDir / 'human_structural_interactome.txt'
    inPath = modelBasedDir / 'model_interfaces.txt'
    
    # output data files
    #outPath = modelBasedDir / 'human_structural_interactome_withDuplicates_new.txt'
    #outPath = modelBasedDir / 'human_structural_interactome_new.txt'
    outPath = modelBasedDir / 'model_interfaces_new.txt'
    
    with io.open(inPath, "r", encoding="utf-8") as f, io.open(outPath, "w") as fout:
        for line in f:
            line = line.replace('#','=')
            fout.write(line)

if __name__ == "__main__":
    main()