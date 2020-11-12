#----------------------------------------------------------------------------------------
# Run this script from inside a Beluga job folder, and it will submit all job scripts with
# extension .sh to the server.
#----------------------------------------------------------------------------------------

from os import listdir
from pathlib import Path
from time import sleep
from subprocess import check_output

def main():

    # directory of server jobs
    jobDir = Path('./')
    
    jobfiles = listdir(jobDir)
    jobfiles = [f for f in jobfiles if f.endswith('job.sh')]
    n = len(jobfiles)
    for i, f in enumerate(jobfiles):
        r = b'failed'
        while not r.startswith(b'Submitted batch job'):
            try:
                print('Submitting job %s, %d of %d (%.1f%%)' % (f, i+1, n, 100.*(i+1)/n))
                r = check_output(['sbatch' , str(jobDir / f)])
            except:
                print('Submission failed. Trying again ...')
                sleep(3)
        print('Submission succeeded\n')
        sleep(3)

if __name__ == "__main__":
    main()