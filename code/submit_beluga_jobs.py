from os import listdir
from pathlib import Path
from time import sleep
from subprocess import check_output

def main():

    # directory of server jobs
    jobDir = Path('./jobs')
    
    jobfiles = listdir(jobDir)
    for f in jobfiles:
        if f.endswith('.sh'):
            r = 'failed'
            while not r.startswith('Submitted batch job'):
                try:
                    r = check_output(['sbatch' , str(jobDir / f)])
                except:
                    pass
                sleep(3)

if __name__ == "__main__":
    main()