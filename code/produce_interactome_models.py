import os
from pathlib import Path
from modeller_tools import produce_interactome_models

def main():
    
    # reference interactome name
    # options: HI-II-14, IntAct
    interactome_name = 'HI-II-14'
    
    # verbosity for Modeller
    verbosity = 'none'
    
    # parent directory of all data files
    dataDir = Path('../data')
    
    # parent directory of all processed data files
    procDir = dataDir / 'processed'
    
    # directory of processed data files specific to interactome
    interactomeDir = procDir / interactome_name
    
    # directory of processed model-related data files specific to interactome
    modelBasedDir = interactomeDir / 'model_based'
    
    # directory for template structure files
    templateDir = Path('../templates')
    
    # directory for alignment files
    alignmentDir = Path('../alignments')
    
    # directory for output models
    modelDir = Path('../models')
    
    # input data files
    interactomeFile = modelBasedDir / 'human_template_annotated_interactome.txt'
    
    # create output directories if not existing
    if not modelDir.exists():
        os.makedirs(str(modelDir))
    
    print('Creating PPI models')
    produce_interactome_models (interactomeFile,
                                alignmentDir,
                                templateDir,
                                modelDir,
                                numModels = 1,
                                verbosity = verbosity)

if __name__ == "__main__":
    main()