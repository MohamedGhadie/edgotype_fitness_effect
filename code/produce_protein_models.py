import os
from pathlib import Path
from modeller_tools import produce_protein_models

def main():
    
    # reference interactome name
    # options: HI-II-14, IntAct
    interactome_name = 'HI-II-14'
    
    # verbosity for Modeller
    verbosity = 'none'
    
    # parent directory of all data files
    dataDir = Path('../../../../')
    
    # parent directory of all processed data files
    procDir = dataDir / 'processed'
    
    # directory of processed data files specific to interactome
    interactomeDir = procDir / interactome_name
    
    # directory of processed model-related data files specific to interactome
    modelBasedDir = interactomeDir / 'model_based'
    
    # directory for template structure files
    templateDir = modelBasedDir / 'protein_templates'
    
    # directory for alignment files
    alignmentDir = modelBasedDir / 'protein_alignments'
    
    # directory for output models
    modelDir = modelBasedDir / 'protein_models'
    
    # input data files
    templateMapFile = modelBasedDir / 'single_template_map_per_protein.txt'
    
    # create output directories if not existing
    if not modelDir.exists():
        os.makedirs(str(modelDir))
    
    print('Creating protein models')
    produce_protein_models (templateMapFile,
                            alignmentDir,
                            templateDir,
                            modelDir,
                            numModels = 1,
                            verbosity = verbosity)

if __name__ == "__main__":
    main()