import os
from pathlib import Path
from modeller_tools import produce_protein_models

def main():
    
    # reference interactome name
    # options: HI-II-14, IntAct
    interactome_name = 'IntAct'
    
    # verbosity for Modeller
    verbosity = 'none'
    
    modellerTimeout = 30
    
    # parent directory of all data files
    #dataDir = Path('../data')
    dataDir = Path('../../../../')
    
    # parent directory of all processed data files
    procDir = dataDir / 'processed'
    
    # directory of processed data files specific to interactome
    interactomeDir = procDir / interactome_name
    
    # directory of processed model-related data files specific to interactome
    modelBasedDir = interactomeDir / 'model_based'
    
    # directory for template structure files
    templateDir = modelBasedDir / 'ppi_templates'
    
    # directory for alignment files
    alignmentDir = modelBasedDir / 'ppi_alignments'
    
    # directory for output models
    modelDir = modelBasedDir / 'ppi_models_3'
    
    # input data files
    interactomeFile = modelBasedDir / 'template_annotated_interactome_3.txt'
    
    # create output directories if not existing
    if not modelDir.exists():
        os.makedirs(str(modelDir))
    
    print('Creating PPI models')
    produce_protein_models (interactomeFile,
                            alignmentDir,
                            templateDir,
                            modelDir,
                            numModels = 1,
                            verbosity = verbosity,
                            modellerTimeout = modellerTimeout * 60)
    print('done')

if __name__ == "__main__":
    main()