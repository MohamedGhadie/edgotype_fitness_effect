import os
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
    
    # directory for output models
    modelDir = modelBasedDir / 'protein_models'
    
    modelFiles = os.listdir(modelDir)
    for name in modelFiles:
        if name.endswith('.B99990001.pdb'):
            modelID, _, _ = name.split('.')
            os.rename(modelDir / name, modelDir / ('pdb' + modelID + '.ent'))
    
#     modelFiles = os.listdir(modelDir)
#     for name in modelFiles:
#         if name.endswith('.pdb'):
#             namesplit = name.split('.')
#             if len(namesplit) == 2:
#                 os.rename(modelDir / name, modelDir / ('pdb' + namesplit[0] + '.ent'))
    
#     modelFiles = os.listdir(modelDir)
#     for name in modelFiles:
#         if name.endswith('.ent'):
#             namesplit = name.split('-')
#             if len(namesplit) == 2:
#                 os.rename(modelDir / name, modelDir / '#'.join(namesplit))
    
#     modelFiles = os.listdir(modelDir)
#     for name in modelFiles:
#         if name.endswith('.ent'):
#             namesplit = name.split('+')
#             if len(namesplit) == 2:
#                 os.rename(modelDir / name, modelDir / '#'.join(namesplit))

#     modelFiles = os.listdir(modelDir)
#     for name in modelFiles:
#         if '#' in name:
#             newname = name.replace('#','=')
#             os.rename(modelDir / name, modelDir / newname)

if __name__ == "__main__":
    main()