#----------------------------------------------------------------------------------------
# Modules for text processing.
#----------------------------------------------------------------------------------------

import io
import pandas as pd
import copy as cp
from pathlib import Path

def write_fasta_file (data, idCol, seqCol, outPath):
    """Produce FASTA file from DataFrame columns.

    Args:
        data (DataFrame): table to be written to FASTA file.
        idCol (str): name of column containing IDs to be written as FASTA headers.
        seqCol (str): name of column containing sequences to be written to FASTA file.
        outPath (Path): path to save FASTA file to.
        
    """
    with io.open(outPath, "w") as fout:
        for _, row in data.iterrows():
            fout.write('>' + row[idCol] + '\n')
            fout.write(row[seqCol] + '\n')

def produce_item_list (inPath, cols, outPath):
    """Write unique items from specific table columns to file.

    Args:
        inPath (Path): path to file containing table to be written.
        cols (list, str): names of columns to be saved to file.
        outPath (Path): path to save list of items to.

    """
    df = pd.read_table(inPath, sep='\t')
    items = sorted(set(df[cols].values.flatten()))
    with open(outPath, 'w') as fout:
        for item in items:
            fout.write('%s\n' % item)

def read_list_table (inPath, cols, dtyp, delm = '\t'):
    """Read a table from a file into a DataFrame with specified column entries converted 
        to lists.

    Args:
        inPath (Path): path to file containing table to be read.
        cols (list, str): names of columns whose entries are to be converted to lists.
        dtyp (list, datatype): data types to convert each column list elements to.
        delm (str): delimiter used to read table.

    """
    df = pd.read_table(inPath, sep = delm)
    if not isinstance(cols, (list, tuple)):
        cols = [cols]
        dtyp = [dtyp]
    for col, typ in zip(cols, dtyp):
        strcol = df[col].apply(str)
        df[col] = strcol.apply(lambda x: list(map(typ, map(str.strip, x.split(',')))))
    return df

def write_list_table (df, cols, outPath, delm = '\t'):
    """Write a DataFrame to text file with specified columns with list entries converted 
        to comma-separated strings.

    Args:
        df (DataFrame): table to be writen to file.
        cols (list, str): names of columns whose list entries are to be converted to strings.
        outPath (Path): path to save table to.
        delm (str): delimiter used to write table to file.

    """
    table = cp.deepcopy(df)
    if isinstance(cols, (list, tuple)):
        for col in cols:
            table[col] = table[col].apply(lambda x: ','.join(map(str, x)))
    else:
        table[cols] = table[cols].apply(lambda x: ','.join(map(str, x)))
    table.to_csv(outPath, index=False, sep=delm)

def write_guillimin_job (outPath,
                         nodes = 1,
                         ppn = 1,
                         pmem = 7700,
                         walltime = '1:00:00:00',
                         outputfile = 'outputfile',
                         errorfile = 'errorfile',
                         rapid = None,
                         jobid = None,
                         commands = None):
    """Create a job file to be run by McGill HPC server.
        See http://www.hpc.mcgill.ca

    Args:
        outPath (Path): path to save job file to.
        nodes (numeric): number of server nodes to be allocated.
        ppn (numeric): total number of CPU cores to be allocated.
        pmem (numeric): default random access memory (RAM) in MB to be reserved per core.
        walltime (str: 'days:hr:min:sec'): maximum time allowed for job to run.
        outputfile (Path): path to file where standard output is written.
        errorfile (Path): path to file where runtime error is written.
        rapid (str): resource allocation project identifier (RAPid).
        jobid ('str'): job name.
        commands (list): additional commands in string format to be written to job file.

    """
    with io.open(outPath, "w") as fout:
        fout.write('#!/bin/bash')
        fout.write('\n' + '#PBS -l nodes=%d:ppn=%d,pmem=%dm,walltime=%s' % (nodes, ppn, pmem, walltime))
        if rapid:
            fout.write('\n' + '#PBS -A %s' % rapid)
        fout.write('\n' + '#PBS -o %s' % outputfile)
        fout.write('\n' + '#PBS -e %s' % errorfile)
        if jobid:
            fout.write('\n' + '#PBS -N %s' % jobid)
        fout.write('\n\n' + 'cd $PBS_O_WORKDIR')
        if commands:
            for cmd in commands:
                fout.write('\n' + cmd)

def write_beluga_job (outPath,
                      account = 'ctb-yxia',
                      walltime = '1-00',
                      ntasks = 1,
                      nodes = 1,
                      ntasks_per_node = 1,
                      cpus_per_task = 1,
                      mem = None,
                      mem_per_cpu = '4000M',
                      outputfile = '%x-%j.out',
                      errorfile = None,
                      jobname = None,
                      commands = None):
    """Create a job file to be run on Compute Canada Beluga cluster.
        See https://docs.computecanada.ca/wiki/Beluga/en

    Args:
        outPath (Path): path to save job file to.
        account (str): project account name.
        walltime (str): maximum time allowed for job to run in the format days:hr:min:sec.
        ntasks (numeric): number of processes to be allocated.
        nodes (numeric): number of server nodes to be allocated.
        ntasks_per_node (numeric): number of processes to be allocated per node.
        cpus_per_task (numeric): number of nodes to be allocated per process.
        mem (numeric): memory per node.
        mem_per_cpu (numeric): memory per core.
        outputfile (Path): path to file where standard output is written.
        errorfile (Path): path to file where runtime error is written.
        jobname (str): job name.
        commands (list): commands in string format to be executed by job.

    """
    with io.open(outPath, "w") as fout:
        fout.write('#!/bin/bash')
        fout.write('\n' + '#SBATCH --account=%s' % account)
        fout.write('\n' + '#SBATCH --time=%s' % walltime)
        fout.write('\n' + '#SBATCH --ntasks=%d' % ntasks)
        fout.write('\n' + '#SBATCH --nodes=%d' % nodes)
        fout.write('\n' + '#SBATCH --ntasks-per-node=%d' % ntasks_per_node)
        fout.write('\n' + '#SBATCH --cpus-per-task=%d' %cpus_per_task)
        if mem:
            fout.write('\n' + '#SBATCH --mem=%s' % mem)
        fout.write('\n' + '#SBATCH --mem-per-cpu=%s' % mem_per_cpu)
        fout.write('\n' + '#SBATCH --output=%s' % outputfile)
        if errorfile:
            fout.write('\n' + '#SBATCH --error=%s' % errorfile)
        if jobname:
            fout.write('\n' + '#SBATCH --job-name=%s' % jobname)
        if commands:
            fout.write('\n')
            for cmd in commands:
                fout.write('\n' + cmd)
