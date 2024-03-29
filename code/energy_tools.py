#----------------------------------------------------------------------------------------
# Modules for calculating and processing structural energy.
#----------------------------------------------------------------------------------------

import os
import io
import re
import sys
from pathlib import Path
from simple_tools import toOneLetterAA
from text_tools import write_beluga_job
from pdb_tools import pdbfile_id, solve_pdbfile_id, clear_structures, write_partial_structure

def read_unprocessed_ddg_mutations (inPath, type = 'binding'):
    """Read PDB chain mutations with missing ∆∆G values from file.

    Args:
        inPath (Path): path to file containing mutations.
        type (str): type of ∆∆G values; 'binding' for interface ∆∆G, 'folding' for protein folding ∆∆G.

    Returns:
        dict

    """
    mutations = {}
    done = set()
    with io.open(inPath, "r", encoding="utf-8") as f:
        next(f)
        for line in f:
            strsplit = list( map ( str.strip, line.split('\t') ) )
            if len(strsplit) >= 8:
                protein, partner, pr_pos, pdbid, chainID, ch_pos, mut, ch_partner = strsplit[:8]
                if type is 'binding':
                    mutation = '-'.join( [protein, partner, pr_pos, mut[-1]] )
                elif type is 'folding':
                    mutation = '-'.join( [protein, pr_pos, mut[-1]] )
                if mutation not in done:
                    if len(strsplit) == 8:
                        done.add(mutation)
                        if type is 'binding':
                            struc = (pdbid,) + tuple(sorted([chainID, ch_partner]))
                        elif type is 'folding':
                            struc = pdbid, chainID
                        if struc in mutations:
                            mutations[struc].add(mut)
                        else:
                            mutations[struc] = {mut}
                    elif strsplit[8] is not 'X':
                        done.add(mutation)
    return {k:list(v) for k, v in mutations.items()}

def produce_foldx_and_beluga_jobs (mutations,
                                   pdbDir,
                                   outDir,
                                   type,
                                   foldxParam = None,
                                   account = 'ctb-yxia',
                                   walltime = '1-00',
                                   ntasks = 1,
                                   nodes = 1,
                                   ntasks_per_node = 1,
                                   cpus_per_task = 1,
                                   mem = None,
                                   mem_per_cpu = '4G',
                                   outputfile = '%x-%j.out',
                                   errorfile = None,
                                   username = '',
                                   extraCommands = None,
                                   serverDataDir = '../data'):
    """Produce jobs for the FoldX method for ∆∆G calculations using several commands.

    Args:
        mutations (dict): mutations associated with each structural model.
        pdbDir (Path): file directory containing PDB structures.
        outDir (Path): file directory to save FoldX data files and Beluga jobs to.
        type (str): type of ∆∆G values; 'binding' for interface ∆∆G, 'folding' for protein folding ∆∆G.
        foldxParam (dict): parameters used for each FoldX job, otherwise default.
        account (str): project account name.
        walltime (str): maximum time allowed for job to run.
        ntasks (numeric): number of processes to be allocated.
        nodes (numeric): number of server nodes to be allocated.
        ntasks_per_node (numeric): number of processes to be allocated per node.
        cpus_per_task (numeric): number of nodes to be allocated per process.
        mem (str): memory per node.
        mem_per_cpu (str): memory per core.
        outputfile (Path): path to file where standard output is written.
        errorfile (Path): path to file where runtime error is written.
        username (str): user name to associate with Beluga server job.
        extraCommands (list): additional non-foldx commands to be written to job file.
        serverDataDir (Path): data directory used by Beluga server relative to job directory.

    """
    if type is 'folding':
        produce_foldx_buildmodel_jobs (mutations,
                                       pdbDir,
                                       outDir / 'data',
                                       parameters = foldxParam)
        produce_beluga_foldx_jobs (mutations,
                                   type,
                                   outDir / 'jobs',
                                   account = account,
                                   walltime = walltime,
                                   ntasks = ntasks,
                                   nodes = nodes,
                                   ntasks_per_node = ntasks_per_node,
                                   cpus_per_task = cpus_per_task,
                                   mem = mem,
                                   mem_per_cpu = mem_per_cpu,
                                   outputfile = outputfile,
                                   errorfile = errorfile,
                                   username = username,
                                   extraCommands = extraCommands,
                                   serverDataDir = serverDataDir)
    elif type is 'binding':
        produce_foldx_pssm_jobs (mutations,
                                 pdbDir,
                                 outDir / 'data',
                                 parameters = foldxParam)
        produce_beluga_foldx_jobs (mutations,
                                   type,
                                   outDir / 'jobs',
                                   account = account,
                                   walltime = walltime,
                                   ntasks = ntasks,
                                   nodes = nodes,
                                   ntasks_per_node = ntasks_per_node,
                                   cpus_per_task = cpus_per_task,
                                   mem = mem,
                                   mem_per_cpu = mem_per_cpu,
                                   outputfile = outputfile,
                                   errorfile = errorfile,
                                   username = username,
                                   extraCommands = extraCommands,
                                   serverDataDir = serverDataDir)

def produce_foldx_buildmodel_jobs (mutations, pdbDir, outDir, parameters = None):
    """Produce jobs for the FoldX method for ∆∆G calculation using 'BuildModel' command.

    Args:
        mutations (dict): mutations associated with each structural model.
        pdbDir (Path): file directory containing PDB structures.
        outDir (Path): file directory to save FoldX jobs to.
        parameters (dict): parameters used for each FoldX job, otherwise default.

    """
    clear_structures()
    if not outDir.exists():
        os.makedirs(outDir)
    
    default_param = {'temp':298,
                     'ph':7,
                     'ionStrength':0.05,
                     'water':'-IGNORE',
                     'vdwDesign':2}
    if not parameters:
        parameters = {}
    mutations = mutations.items()
    n = len(mutations)
    print('Writing FoldX BuildModel files:')
    for i, (struc, mutList) in enumerate(mutations):
        sys.stdout.write('  Structure %d out of %d (%.2f%%) \r' % (i+1, n, 100*(i+1)/n))
        sys.stdout.flush()
        strucid = get_strucID (struc)
        strucDir = outDir / strucid
        if not strucDir.exists():
            os.makedirs(strucDir)
        mutListFile = strucDir / 'individual_list.txt'
        mutList = ['%s;' % mutList.pop(0)] + ['\n%s;' % mut for mut in mutList]
        with io.open(mutListFile, "w") as fout:
            for mut in mutList:
                fout.write(mut)
        pdbid, chainID = struc[:2]
        write_partial_structure (pdbid,
                                 [chainID],
                                 pdbDir,
                                 strucDir / (strucid + '.pdb'))
        write_foldx_config (strucDir / 'config_repairPDB.cfg',
                            'RepairPDB',
                            pdb_dir = '../data/%s' % strucid,
                            output_dir = '../data/%s' % strucid,
                            pdb_file = '%s.pdb' % strucid)
        mutParam = {k:v for k, v in default_param.items()}
        if struc in parameters:
            for k, v in parameters[struc].items():
                mutParam[k] = v
        write_foldx_config (strucDir / 'config_buildModel.cfg',
                            'BuildModel',
                            pdb_dir = '../data/%s' % strucid,
                            output_dir = '../data/%s' % strucid,
                            pdb_file = '%s_Repair.pdb' % strucid,
                            mutant_file = '../data/%s/individual_list.txt' % strucid,
                            temp = mutParam['temp'],
                            ph = mutParam['ph'],
                            ionStrength = mutParam['ionStrength'],
                            water = mutParam['water'],
                            vdwDesign = mutParam['vdwDesign'])
    print()

def produce_foldx_pssm_jobs (mutations, pdbDir, outDir, parameters = None):
    """Produce jobs for the FoldX method for ∆∆G calculation using PSSM command.

    Args:
        mutations (dict): mutations associated with each structural model.
        pdbDir (Path): file directory containing PDB structures.
        outDir (Path): file directory to save FoldX jobs to.
        parameters (dict): parameters used for each FoldX job, otherwise default.

    """
    clear_structures()
    if not outDir.exists():
        os.makedirs(outDir)
    
    default_param = {'temp':298,
                     'ph':7,
                     'ionStrength':0.05,
                     'water':'-IGNORE',
                     'vdwDesign':2}
    if not parameters:
        parameters = {}
    
    mutations = mutations.items()
    n = len(mutations)
    print('Writing FoldX PSSM files:')
    for i, (struc, mutList) in enumerate(mutations):
        sys.stdout.write('  Structure %d out of %d (%.2f%%) \r' % (i+1, n, 100*(i+1)/n))
        sys.stdout.flush()
        strucid = get_strucID (struc)
        pdbid, chainID1, chainID2 = struc[:3]
        for mut in mutList:
            mutID = '_'.join([strucid, mut])
            mutDir = outDir / mutID
            if not mutDir.exists():
                os.makedirs(mutDir)
            write_partial_structure (pdbid,
                                     [chainID1, chainID2],
                                     pdbDir,
                                     mutDir / (strucid + '.pdb'))
            write_foldx_config (mutDir / 'config_repairPDB.cfg',
                                'RepairPDB',
                                pdb_dir = '../data/%s' % mutID,
                                output_dir = '../data/%s' % mutID,
                                pdb_file = '%s.pdb' % strucid)
            mutParam = {k:v for k, v in default_param.items()}
            if (struc, mut) in parameters:
                for k, v in parameters[(struc, mut)].items():
                    mutParam[k] = v
            write_foldx_config (mutDir / 'config_pssm.cfg',
                                'Pssm',
                                other_cmd = ['analyseComplexChains=%s,%s' % (chainID1, chainID2),
                                             'aminoacids=%s' % mut[-1],
                                             'positions=%sa' % mut[:-1]],
                                pdb_dir = '../data/%s' % mutID,
                                output_dir = '../data/%s' % mutID,
                                pdb_file = '%s_Repair.pdb' % strucid,
                                temp = mutParam['temp'],
                                ph = mutParam['ph'],
                                ionStrength = mutParam['ionStrength'],
                                water = mutParam['water'],
                                vdwDesign = mutParam['vdwDesign'])
    print()

def write_foldx_config (outPath,
                        command,
                        other_cmd = None,
                        pdb_dir = None,
                        output_dir = None,
                        pdb_file = None,
                        mutant_file = None,
                        temp = 298,
                        ph = 7,
                        ionStrength = 0.05,
                        water = '-IGNORE',
                        vdwDesign = 2):
    """Produce FoldX configuration file.

    Args:
        outPath (Path): file path to save FoldX configurations.
        command (str): main command to be run by FoldX.
        other_cmd (list): additional commands in string format to be run by FoldX.
        pdb_dir (Path): file directory containing PDB structures.
        output_dir (Path): file directory where FoldX results are saved.
        pdb_file (Path): path to file containing PDB structure to be processed.
        mutant_file (Path): path to file containing list of mutations if applicable.
        temp (numeric): temperature in Kelvins.
        ph (numeric): PH value.
        ionStrength (numeric): ionic strength of the solution in Moles.
        water (str): how to handle water molecules, '-CRYSTAL', '-PREDICT', '-IGNORE' or '-COMPARE'.
        vdwDesign (numeric): VDW design of the experiment, 0 very soft, 1 medium soft, 2 strong.

    """
    with io.open(outPath, "w") as fout:
        fout.write('command=%s' % command)
        if other_cmd:
            fout.write('\n' + '\n'.join(other_cmd))
        if pdb_dir:
            fout.write('\n' + 'pdb-dir=%s' % pdb_dir)
        if output_dir:
            fout.write('\n' + 'output-dir=%s' % output_dir)
        if pdb_file:
            fout.write('\n' + 'pdb=%s' % pdb_file)
        if mutant_file:
            fout.write('\n' + 'mutant-file=%s' % mutant_file)
        fout.write('\n' + 'temperature=%.1f' % temp)
        fout.write('\n' + 'pH=%.1f' % ph)
        fout.write('\n' + 'ionStrength=%f' % ionStrength)
        fout.write('\n' + 'water=%s' % water)
        fout.write('\n' + 'vdwDesign=%d' % vdwDesign)

def produce_beluga_foldx_jobs (mutations,
                               type,
                               outDir,
                               account = 'ctb-yxia',
                               walltime = '1-00',
                               ntasks = 1,
                               nodes = 1,
                               ntasks_per_node = 1,
                               cpus_per_task = 1,
                               mem = None,
                               mem_per_cpu = '4G',
                               outputfile = '%x-%j.out',
                               errorfile = None,
                               username = '',
                               extraCommands = None,
                               serverDataDir = '../data'):
    """Produce Beluga server job files specific to FoldX ∆∆G calculations.
        See https://docs.computecanada.ca/wiki/Béluga/en

    Args:
        mutations (dict): mutations associated with each structural model.
        type (str): type of ∆∆G values; 'binding' for interface ∆∆G, 'folding' for protein folding ∆∆G.
        outDir (Path): file directory to save jobs to.
        account (str): project account name.
        walltime (str): maximum time allowed for job to run.
        ntasks (numeric): number of processes to be allocated.
        nodes (numeric): number of server nodes to be allocated.
        ntasks_per_node (numeric): number of processes to be allocated per node.
        cpus_per_task (numeric): number of nodes to be allocated per process.
        mem (str): memory per node.
        mem_per_cpu (str): memory per core.
        outputfile (Path): path to file where standard output is written.
        errorfile (Path): path to file where runtime error is written.
        username (str): user name to associate with Beluga server job.
        extraCommands (list): additional non-foldx commands to be written to job file.
        serverDataDir (Path): data directory used by Beluga server relative to job directory.

    """
    if not outDir.exists():
        os.makedirs(outDir)
    
    print('Writing Beluga job files')
    if type is 'folding':
        for struc, _ in mutations.items():
            strucid = get_strucID (struc)
            commands = ['../foldx -f %s/%s/config_repairPDB.cfg' % (serverDataDir, strucid),
                        '../foldx -f %s/%s/config_buildModel.cfg' % (serverDataDir, strucid)]
            if extraCommands:
                commands = extraCommands + commands
            write_beluga_job (outDir / (strucid + '_job.sh'),
                              account = account,
                              walltime = walltime,
                              ntasks = ntasks,
                              nodes = nodes,
                              ntasks_per_node = ntasks_per_node,
                              cpus_per_task = cpus_per_task,
                              mem = mem,
                              mem_per_cpu = mem_per_cpu,
                              outputfile = outputfile,
                              errorfile = errorfile,
                              jobname = '%s_foldx_buildmodel_%s' % (username, strucid),
                              commands = commands)
    elif type is 'binding':
        for struc, mutList in mutations.items():
            strucid = get_strucID (struc)
            for mut in mutList:
                mutID = '_'.join([strucid, mut])
                commands = ['../foldx -f %s/%s/config_repairPDB.cfg' % (serverDataDir, mutID),
                            '../foldx -f %s/%s/config_pssm.cfg' % (serverDataDir, mutID)]
                if extraCommands:
                    commands = extraCommands + commands
                write_beluga_job (outDir / (mutID + '_job.sh'),
                                  account = account,
                                  walltime = walltime,
                                  ntasks = ntasks,
                                  nodes = nodes,
                                  ntasks_per_node = ntasks_per_node,
                                  cpus_per_task = cpus_per_task,
                                  mem = mem,
                                  mem_per_cpu = mem_per_cpu,
                                  outputfile = outputfile,
                                  errorfile = errorfile,
                                  jobname = '%s_foldx_%s' % (username, mutID),
                                  commands = commands)

def read_foldx_results (inDir, type = 'binding'):
    """Read mutation ∆∆G results produced by FoldX commands.

    Args:
        inDir (Path): file directory containing FoldX results.
        type (str): type of ∆∆G values; 'binding' for interface ∆∆G, 'folding' for protein folding ∆∆G.

    Returns:
        dict, dict: processed and unprocessed mutations.

    """
    if type is 'folding':
        return read_foldx_buildmodel_results (inDir)
    elif type is 'binding':
        return read_foldx_pssm_results (inDir)

def read_foldx_buildmodel_results (inDir):
    """Read mutation ∆∆G results produced by FoldX BuildModel command.

    Args:
        inDir (Path): file directory containing FoldX results.

    Returns:
        dict, dict: processed and unprocessed mutations.

    """
    processed, unprocessed = {}, {}
    strucDir = os.listdir(inDir)
    strucDir = [dir for dir in strucDir if os.path.isdir(inDir / dir)]
    for strucID in strucDir:
        struc = tuple((solve_pdbfile_id(strucID)).split('_'))
        if re.match(r'\D\S\d+\D', struc[-1]):
            struc = struc[:-1]
        mutListFile = inDir / strucID / 'individual_list.txt'
        with io.open(mutListFile, "r") as f:
            mutList = list( map(str.strip, f.read().split(';')) )
        mutList.remove('')
        
        resultFile = None
        for filename in os.listdir(inDir / strucID):
            if filename.startswith('Average'):
                resultFile = filename
        
        results = []
        if resultFile:
            resultFile = inDir / strucID / resultFile
            headers = ["Pdb", "SD", "total energy", "Backbone Hbond", "Sidechain Hbond",
                       "Van der Waals", "Electrostatics", "Solvation Polar", "Solvation Hydrophobic", 
                       "Van der Waals clashes", "entropy sidechain", "entropy mainchain", 
                       "sloop_entropy", "mloop_entropy", "cis_bond", "torsional clash", "backbone clash", 
                       "helix dipole", "water bridge", "disulfide", "electrostatic kon", 
                       "partial covalent bonds", "energy Ionisation", "Entropy Complex"]
            with io.open(resultFile, "r") as f:
                for line in f:
                    linesplit = list(map(str.strip, line.split('\t')))
                    if linesplit == headers:
                        break
                for line in f:
                    linesplit = list(map(str.strip, line.split('\t')))
                    if len(linesplit) > 3:
                        results.append(linesplit[2])
                if len(results) == len(mutList):
                    for mut, ddg in zip(mutList, results):
                        processed[struc + (mut,)] = float(ddg)
        
        if not results:
            if len(mutList) == 1:
                processed[struc + (mutList.pop(),)] = 'X'
            elif len(mutList) > 1:
                for mut in mutList:
                    unprocessed[struc + (mut,)] = [mut]
    
    return processed, unprocessed

def read_foldx_pssm_results (inDir):
    """Read mutation ∆∆G results produced by FoldX PSSM command.

    Args:
        inDir (Path): file directory containing FoldX results.

    Returns:
        dict, dict: processed and unprocessed mutations.

    """
    processed, unprocessed = {}, {}
    strucDirs = os.listdir(inDir)
    strucDirs = [dir for dir in strucDirs if os.path.isdir(inDir / dir)]
    for strucID in strucDirs:
        strucDir = inDir / strucID
        struc = tuple((solve_pdbfile_id(strucID)).split('_'))
        if len(struc) > 3:
            struc = struc[:-1]
        
        mutListFile = resultFile = None
        for filename in os.listdir(strucDir):
            if filename.startswith('individual_list'):
                mutListFile = strucDir / filename
                with io.open(mutListFile, "r") as f:
                    mutList = list(map(str.strip, f.read().split(';')))
                    mutList.remove('')
            elif filename.startswith('_'.join(('PSSM',) + struc)):
                resultFile = strucDir / filename
        
        if resultFile:
            with io.open(resultFile, "r") as f:
                mutResidues = list(map(str.strip, f.readline().strip().split('\t')))
                for line in f:
                    ddgList = list(map(str.strip, line.split('\t')))
                    if len(ddgList) > 1:
                        wt = ddgList[0]
                        for mt, ddg in zip(mutResidues, ddgList[1:]):
                            processed[struc + (wt + mt,)] = float(ddg) if len(ddg) else 'X'
                            
            for mut in mutList:
                if struc + (mut,) not in processed:
                    processed[struc + (mut,)] = 'X'
        else:
            if len(mutList) == 1:
                processed[struc + (mutList.pop(),)] = 'X'
            elif len(mutList) > 1:
                for mut in mutList:
                    unprocessed[struc + (mut,)] = [mut]
    
    return processed, unprocessed

def produce_mCSM_jobs (mutations, pdbDir, outDir):
    """Produce jobs to be submitted to MCSM-PPI2 server for binding ∆∆G calculations.

    Args:
        mutations (dict): mutations associated with each structural model.
        pdbDir (Path): file directory containing PPI structures.
        outDir (Path): file directory to save MCSM-PPI2 jobs to.

    """
    clear_structures()
    if not outDir.exists():
        os.makedirs(outDir)
    
    mutations = mutations.items()
    n = len(mutations)
    print('Writing mCSM-PPI2 files:')
    for i, (struc, mutList) in enumerate(mutations):
        sys.stdout.write('  Structure %d out of %d (%.2f%%) \r' % (i+1, n, 100*(i+1)/n))
        sys.stdout.flush()
        strucid = get_strucID (struc)
        strucDir = outDir / strucid
        if not strucDir.exists():
            os.makedirs(strucDir)
        mutListFile = strucDir / 'mutList.txt'
        mutList = [mut[1] + ' ' + mut[0] + mut[2:] for mut in mutList]
        mutStr = '\n'.join(mutList)
        with io.open(mutListFile, "w") as fout:
            fout.write(mutStr)
        pdbid, chainID1, chainID2 = struc[:3]
        write_partial_structure (pdbid,
                                 [chainID1, chainID2],
                                 pdbDir,
                                 strucDir / (strucid + '.pdb'))

def read_mCSM_results (inPath):
    """Read mutation ∆∆G results produced by mCSM-PPI2.

    Args:
        inPath (Path): path to file containing mCSM-PPI2 results.

    Returns:
        dict: processed mutations.

    """
    processed = {}
    with io.open(inPath, "r", encoding="utf-8") as f:
        next(f)
        for line in f:
            linesplit = line.split(',')
            if len(linesplit) >= 8:
                strucID, chainID, wt, resNum, mt, dist, ddg, affinity = linesplit[:8]
                struc = tuple(strucID.split('_'))
                if len(struc) > 3:
                    struc = struc[:-1]
                wt, mt = toOneLetterAA(wt), toOneLetterAA(mt)
                processed[struc + (wt + chainID + resNum + mt,)] = -float(ddg)    
    return processed

def produce_dynamut_jobs (mutations, pdbDir, outDir):
    """Produce jobs to be submitted to DynaMut2 server for folding ∆∆G calculations.

    Args:
        mutations (dict): mutations associated with each structural model.
        pdbDir (Path): file directory containing PPI structures.
        outDir (Path): file directory to save DynaMut2 jobs to.

    """
    clear_structures()
    if not outDir.exists():
        os.makedirs(outDir)
    
    mutations = mutations.items()
    n = len(mutations)
    print('Writing DynaMut2 files:')
    for i, (struc, mutList) in enumerate(mutations):
        sys.stdout.write('  Structure %d out of %d (%.2f%%) \r' % (i+1, n, 100*(i+1)/n))
        sys.stdout.flush()
        strucid = get_strucID (struc)
        strucDir = outDir / strucid
        if not strucDir.exists():
            os.makedirs(strucDir)
        mutListFile = strucDir / 'mutList.txt'
        mutList = [mut[1] + ' ' + mut[0] + mut[2:] for mut in mutList]
        mutStr = '\n'.join(mutList)
        with io.open(mutListFile, "w") as fout:
            fout.write(mutStr)
        pdbid, chainID = struc[:2]
        write_partial_structure (pdbid,
                                 [chainID],
                                 pdbDir,
                                 strucDir / (strucid + '.pdb'))

def read_dynamut_results (inDir):
    """Read mutation ∆∆G results produced by DynaMut2.

    Args:
        inDir (Path): file directory containing DynaMut2 results.

    Returns:
        dict: processed mutations.

    """
    processed = {}
    strucDir = os.listdir(inDir)
    strucDir = [dir for dir in strucDir if os.path.isdir(inDir / dir)]
    for strucID in strucDir:
        struc = tuple((solve_pdbfile_id(strucID)).split('_'))
        if re.match(r'\D\S\d+\D', struc[-1]):
            struc = struc[:-1]
        
        for filename in os.listdir(inDir / strucID):
            if filename.startswith('results_'):
                with io.open(inDir / strucID / filename, "r") as f:
                    next(f)
                    for line in f:
                        linesplit = list(map(str.strip, line.split(',')))
                        if len(linesplit) > 2:
                            mut, chainID, ddg = linesplit[:3]
                            processed[struc + (mut[0] + chainID + mut[1:],)] = -float(ddg)   
    return processed

def read_protein_mutation_ddg (inPath, type = 'binding'):
    """Read protein mutation ∆∆G values from file.

    Args:
        inPath (Path): path to file containing mutations with ∆∆G values.
        type (str): type of ∆∆G values; 'binding' for interface ∆∆G, 'folding' for protein folding ∆∆G.

    Returns:
        dict

    """
    ddgDict = {}
    with io.open(inPath, "r", encoding="utf-8") as f:
        next(f)
        for line in f:
            linesplit = list( map ( str.strip, line.split('\t') ) )
            if len(linesplit) > 8:
                if linesplit[8] is not 'X':
                    protein, partner, pr_pos, pdbid, chainID, ch_pos, ch_mut, ch_partner, ddg = linesplit
                    if type is 'binding':
                        k = protein, partner, int(pr_pos), ch_mut[-1]
                        val = pdbid, chainID, ch_partner, ch_mut, float(ddg)
                    elif type is 'folding':
                        k = protein, int(pr_pos), ch_mut[-1]
                        val = pdbid, chainID, ch_mut, float(ddg) 
                    if k not in ddgDict:
                        ddgDict[k] = val
    return ddgDict

def read_chain_mutation_ddg (inPath, type = 'binding'):
    """Read PDB chain mutation ∆∆G values from file.

    Args:
        inPath (Path): path to file containing mutations with ∆∆G values.
        type (str): type of ∆∆G values; 'binding' for interface ∆∆G, 'folding' for protein folding ∆∆G.

    Returns:
        dict

    """
    ddgDict = {}
    with io.open(inPath, "r", encoding="utf-8") as f:
        next(f)
        for line in f:
            linesplit = list( map ( str.strip, line.split('\t') ) )
            if len(linesplit) > 8:
                protein, partner, pr_pos, pdbid, chainID, ch_pos, ch_mut, ch_partner, ddg = linesplit[:9]
                if type is 'binding':
                    k = (pdbid,) + tuple(sorted([chainID, ch_partner])) + (ch_mut,)
                elif type is 'folding':
                    k = pdbid, chainID, ch_mut
                ddgDict[k] = ddg
    return ddgDict

def copy_mutation_ddg (inPath1, inPath2, outPath, type = 'binding'):
    """Transfer mutation ∆∆G values from one file to another file based on similar mutation 
        structure mappings.

    Args:
        inPath1 (Path): path to file containing mutation ∆∆G values.
        inPath2 (Path): path to file containing mutations with missing ∆∆G values.
        outPath (Path): file path to save a copy of file in inPath2 with updated ∆∆G values.
        type (str): type of ∆∆G values; 'binding' for interface ∆∆G, 'folding' for protein folding ∆∆G.

    """
    ddg = read_chain_mutation_ddg (inPath1, type)
    write_mutation_ddg_tofile (ddg, inPath2, outPath, type)

def write_mutation_ddg_tofile (ddg, inPath, outPath, type = 'binding'):
    """Update file with mutation ∆∆G values.

    Args:
        ddg (dict): mutation ∆∆G values.
        inPath (Path): path to file whose mutations will be updated with ∆∆G values.
        outPath (Path): file path to save a copy of file in inPath with updated ∆∆G values.
        type (str): type of ∆∆G values; 'binding' for interface ∆∆G, 'folding' for protein folding ∆∆G.

    """
    with io.open(inPath, "r", encoding="utf-8") as f, io.open(outPath, "w") as fout:
        fout.write(f.readline().strip() + '\n')
        for line in f:
            strsplit = list(map(str.strip, line.split('\t')))
            if len(strsplit) == 8:
                protein, partner, pr_pos, pdbid, chainID, ch_pos, ch_mut, ch_partner = strsplit
                if type is 'binding':
                    k = (pdbid,) + tuple(sorted([chainID, ch_partner])) + (ch_mut,)
                elif type is 'folding':
                    k = pdbid, chainID, ch_mut
                if k in ddg:
                    strsplit.append(str(ddg[k]))
            fout.write('\t'.join(map(str, strsplit)) + '\n')

def append_mutation_ddg_files (inPath1, inPath2, outPath):
    """Append two ∆∆G files together and save to anothor file.

    Args:
        inPath1 (Path): path to first input file.
        inPath2 (Path): path to second input file.
        outPath (Path): file path to save appended output to.

    """
    with io.open(outPath, "w") as fout:
        with io.open(inPath1, "r", encoding="utf-8") as f1:
            fout.write(f1.readline().strip() + '\n')
            for line in f1:
                fout.write(line)
        with io.open(inPath2, "r", encoding="utf-8") as f2:
            next(f2)
            for line in f2:
                fout.write(line)

def get_strucID (struc):
    """Return structure file ID from structure, chains and mutation ID tuple.

    Args:
        strucid (tuple): structure, chains and mutation ID tuple.

    Returns:
        str: structure file ID.

    """
    if re.match(r'\D\S\d+\D', struc[-1]):
        return pdbfile_id ('_'.join(struc[:-1])) + '_' + struc[-1]
    else:
        return pdbfile_id ('_'.join(struc))
