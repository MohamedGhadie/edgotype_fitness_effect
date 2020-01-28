import os
from pathlib import Path

def main():
    
    # directory for output models
    modelDir = Path('../models')
    
    if not modelDir.exists():
        print("Model directory not found")
        return
    
    modelFiles = os.listdir(modelDir)
    for name in modelFiles:
        if name.endswith('.pdb'):
            namesplit = name.split('.')
            if len(namesplit) == 3:
                os.rename(modelDir / name, modelDir / ('pdb' + namesplit[0] + '.ent'))
    
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
#                 os.rename(modelDir / name, modelDir / '+'.join(namesplit))
    
#     modelFiles = os.listdir(modelDir)
#     for name in modelFiles:
#         if name.endswith('.ent'):
#             namesplit = name.split('+')
#             if len(namesplit) == 2:
#                 os.rename(modelDir / name, modelDir / '='.join(namesplit))

if __name__ == "__main__":
    main()