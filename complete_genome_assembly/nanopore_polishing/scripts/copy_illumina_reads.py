import pandas as pd
import shutil
import os

df = pd.read_csv('../data/kit12/original_file_paths.csv')

if not os.path.exists('/n/holyscratch01/grad_lab/Users/qinqinyu/20240802_nanopore_polishing_kit12/illumina_reads/'):
    os.mkdir('/n/holyscratch01/grad_lab/Users/qinqinyu/20240802_nanopore_polishing_kit12/illumina_reads/')
    
for i, row in df.iterrows():
    shutil.copy(row['illumina_reads_1'], '/n/holyscratch01/grad_lab/Users/qinqinyu/20240802_nanopore_polishing_kit12/illumina_reads/' + row['strain'] + '_1.fastq.gz')
    shutil.copy(row['illumina_reads_2'], '/n/holyscratch01/grad_lab/Users/qinqinyu/20240802_nanopore_polishing_kit12/illumina_reads/' + row['strain'] + '_2.fastq.gz')