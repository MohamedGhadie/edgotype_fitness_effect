#----------------------------------------------------------------------------------------
# Parse a single flat file of mutations from the dbSNP database.
#----------------------------------------------------------------------------------------

import os
from pathlib import Path
from mutation_processing_tools import parse_dbsnp_flatfile_keys, parse_dbsnp_flatfile

def main():
    
    # chromosome number
    chr = '20'
    
    # parent directory of all data files
    dataDir = Path('../data')
    
    # directory of data files from external sources
    extDir = dataDir / 'external'
    
    # parent directory of all processed data files
    procDir = dataDir / 'processed'
    
    # directory for dbSNP raw data files
    dbsnpInDir = extDir / 'dbsnp'
    
    # directory to save dbSNP intermediate processed data files
    dbsnpOutDir = procDir / 'dbsnp_intermediate'
    
    # input data files
    dbsnpInFile = dbsnpInDir / ('ds_flat_ch%s.flat' % chr)
    
    # output data files
    dbsnpOutFile = dbsnpOutDir / ('dbsnp_chr%s.txt' % chr)
    
    # create output directories if not existing
    if not procDir.exists():
        os.makedirs(procDir)
    if not dbsnpOutDir.exists():
        os.makedirs(dbsnpOutDir)
    
    #------------------------------------------------------------------------------------
    # process dbSNP mutations
    #------------------------------------------------------------------------------------
    
    if chr is '1':
        print('collecting dbSNP flatfile keys')
        parse_dbsnp_flatfile_keys (dbsnpInFile, dbsnpOutDir, pausetime = 0)
    
    print('parsing dbSNP flatfile')
    parse_dbsnp_flatfile (dbsnpInFile, dbsnpOutDir, dbsnpOutFile, pausetime = 0)
    
    print('done')
    
if __name__ == '__main__':
    main()
