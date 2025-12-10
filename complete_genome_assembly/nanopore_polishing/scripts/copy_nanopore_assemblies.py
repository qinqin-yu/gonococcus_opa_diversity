import pandas as pd
import shutil
import os

df = pd.read_csv('../data/kit12/original_file_paths.csv')

if not os.path.exists('../data/kit12/flye_assemblies'):
    os.mkdir('../data/kit12/flye_assemblies')
    
for i, row in df.iterrows():
    shutil.copy(row['flye_assembly'], '../data/kit12/flye_assemblies/' + row['strain'] + '.fa')
    shutil.copy(row['flye_assembly_medaka_polished'], '../data/kit12/flye_assemblies_medaka_polished/' + row['strain'] + '.medaka.fa')