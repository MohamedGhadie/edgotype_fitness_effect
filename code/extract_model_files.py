import os
import subprocess
from pathlib import Path

def main():
    
    # directory for output models
    modelDir = Path('../data/processed/HI-II-14/model_based/ppi_models')
    
    outDir = Path('../models_only')
    
    if not modelDir.exists():
        print("Model directory not found")
        return
    
    if not outDir.exists():
        os.makedirs(outDir)
    
    modelFiles = os.listdir(modelDir)
    for name in modelFiles:
        if name.endswith('.ent'):
            cmd = ['cp', str(modelDir / name), str(outDir / name)]
            subprocess.run(cmd)

if __name__ == "__main__":
    main()