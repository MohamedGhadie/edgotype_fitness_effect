
import pandas as pd
from pathlib import Path

def main():
    
    interactome_name = 'experiment'
    
    # parent directory of all data files
    dataDir = Path('../data')
    
    interactomeDir = dataDir / interactome_name
    
    # input data files
    mutationFile = interactomeDir / 'processed_sahni_mutations.txt'
    
    # output data files
    naturalMutationsFile = interactomeDir / 'nondisease_mutation_edgotype.txt'
    diseaseMutationsFile = interactomeDir / 'disease_mutation_edgotype.txt'
    
    mutations = pd.read_table (mutationFile, sep='\t')
    naturalMutations = mutations [mutations["Category"] == 'Non-disease variant'].reset_index(drop = True)
    diseaseMutations = mutations [mutations["Category"] == 'Disease mutation'].reset_index(drop = True)
    
    naturalMutations.to_csv (naturalMutationsFile, index=False, sep='\t')
    diseaseMutations.to_csv (diseaseMutationsFile, index=False, sep='\t')

if __name__ == "__main__":
    main()
