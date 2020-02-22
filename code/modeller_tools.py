import pandas as pd
from modeller import *
from modeller.automodel import *

def produce_protein_models (inPath,
                            alignmentDir,
                            templateDir,
                            modelDir,
                            numModels = 1,
                            verbosity = 'minimal'):
    
    templateMap = pd.read_table (inPath, sep='\t')
    n = len(templateMap)
    
    for i, row in templateMap.iterrows():
        print('\n********************************************************************')
        print('Protein complex %d out of %d (%.2f%%)' % (i+1, n, 100.*(i+1)/n))
        print('********************************************************************\n')
        modelFile = modelDir / (row.Complex_ID + '.B99990001.pdb')
        if not modelFile.is_file():
            create_protein_model (row.Complex_ID,
                                  row.Template_file_ID,
                                  str(alignmentDir / row.Alignment_file_ID),
                                  str(templateDir),
                                  str(modelDir),
                                  starting_model = 1,
                                  ending_model = numModels,
                                  verbosity = verbosity)

def create_protein_model (protein,
                          templateIDs,
                          alignmentFile,
                          templateDir,
                          modelDir,
                          starting_model = 1,
                          ending_model = 1,
                          verbosity = 'minimal'):
    
    if verbosity is 'verbose':
        log.verbose()
    elif verbosity is 'minimal':
        log.minimal()
    else:
        log.none()
    
    env = environ()
    
    # read heteroatoms
    env.io.hetatm = True
    
    # directories for input atom files
    env.io.atom_files_directory = [modelDir, templateDir]

    # Be sure to use 'MyModel' rather than 'automodel' here!
    a = MyModel(env,
                alnfile  = alignmentFile,
                knowns   = templateIDs,
                sequence = protein)

    a.starting_model = starting_model
    a.ending_model = ending_model
    a.make()

# Override MyModel methods
class MyModel (automodel):
    
    def special_patches (self, aln):
        # number of chains in model
        numChains = len(self.chains)
        
        # Rename both chains and renumber the residues in each
        self.rename_segments (segment_ids = [c.name for c in self.chains],
                              renumber_residues = [1] * numChains)
        
        # Rename chain for single chain model:
        if numChains == 1:
            self.chains[0].name = 'A'
