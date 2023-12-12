import numpy as np
import pandas as pd
import subprocess
import os



ls_bam = os.listdir('/mnt/project2/tmp70T/gth/ICGC/coverage_10X')


for i,bam_name in enumerate(ls_bam):
    prefix = bam_name.replace('.bam','')
    command = "snakemake -s Snake_index.smk --cores 2 --config prefix="+prefix
    subprocess.run(command, shell=True)

    print('Indexation of file ',prefix, ' ended successfully')