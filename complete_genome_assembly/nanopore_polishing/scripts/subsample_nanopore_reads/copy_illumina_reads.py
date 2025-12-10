import pandas as pd
import shutil
import os

df = pd.read_csv('../../data/kit12/original_file_paths.csv')

if not os.path.exists('/n/netscratch/grad_lab/Lab/qinqinyu/20241218_illumina_polishing_subsampled_reads/illumina_reads/'):
    os.mkdir('/n/netscratch/grad_lab/Lab/qinqinyu/20241218_illumina_polishing_subsampled_reads/illumina_reads/')
    
for i, row in df.iterrows():
    shutil.copy(row['illumina_reads_1'], '/n/netscratch/grad_lab/Lab/qinqinyu/20241218_illumina_polishing_subsampled_reads/illumina_reads/' + row['strain'] + '_1.fastq.gz')
    shutil.copy(row['illumina_reads_2'], '/n/netscratch/grad_lab/Lab/qinqinyu/20241218_illumina_polishing_subsampled_reads/illumina_reads/' + row['strain'] + '_2.fastq.gz')