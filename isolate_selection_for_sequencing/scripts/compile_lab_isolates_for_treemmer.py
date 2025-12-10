import pandas as pd
import shutil
import os
import glob
import numpy as np

# Get list of lab isolates and complete genomes

# Get strains that we have in the lab
gc_metadata = pd.read_csv('/n/holylfs05/LABS/grad_lab/Lab/repos/gc_genomics/metadata/Ng-Combined-Metadata.txt', sep = '\t')
lab_strains = gc_metadata[gc_metadata['isolate_in_lab']=='Y'][['wgs_id']]

# Get hybrid genome metadata and append to lab strains metadata table
hybrid_genome_metadata = pd.read_csv('../data/genome_metadata/hybrid_genomes_metadata.csv')

gc_metadata = pd.concat([gc_metadata, hybrid_genome_metadata], ignore_index = True)

# Add the UMASS-DGI strains
umass_filenames = glob.glob('/n/grad_lab2/Lab/gonococcus/datasets/umass_dgi/annotations/*.gff')
dgi_strains = []
for filename in umass_filenames:
    basename = filename.split('/')[-1]
    dgi_strains.append(basename[:basename.find('.gff')])
dgi_strains = np.sort(dgi_strains)
dgi_strains_df = pd.DataFrame({'wgs_id':dgi_strains})

# Add the lab strains
# Note that the only "lab strain" that is not represented elsewhere in this list is 28Bl
standard_lab_filenames = glob.glob('/n/grad_lab2/Lab/gonococcus/datasets/lab_strains/annotations/*.gff')
standard_lab_strains = []
for filename in standard_lab_filenames:
    basename = filename.split('/')[-1]
    standard_lab_strains.append(basename[:basename.find('.gff')])
standard_lab_strains = np.sort(standard_lab_strains)
standard_lab_strains_df = pd.DataFrame({'wgs_id':standard_lab_strains})

# Add the hybrid genomes
complete_genome_paths = pd.read_csv('../data/genome_metadata/hybrid_genome_paths.csv')
complete_genome_df = pd.DataFrame({'wgs_id':complete_genome_paths['sample_name'].values})

all_strains = pd.concat([lab_strains, dgi_strains_df, standard_lab_strains_df, complete_genome_df], ignore_index = True)
all_strains.drop_duplicates(inplace = True) # Deduplicate
all_strains.sort_values('wgs_id', inplace = True)
all_strains.reset_index(inplace = True, drop = True)

# Get the pseudogenome paths for the list of representative isolates + lab strains + hybrid genomes (will merge to only be lab strains and hybrid genomes later)

isolates = list(pd.read_csv('../data/gubbins/representative_isolates/representative_isolates.txt', header = None)[0].values)
pseudogenome_paths = list(pd.read_csv('../data/gubbins/representative_isolates/representative_isolates_pseudogenome_paths.txt', header = None)[0].values)

# Add additional pseudogenome paths that were not included in the representative isolates
additional_pseudogenome_paths = ['/n/grad_lab2/Lab/gonococcus/datasets/umass_dgi/',
                                '/n/grad_lab2/Lab/gonococcus/datasets/ethiopia_isolates_grad_lab/',
                                '/n/grad_lab2/Lab/gonococcus/datasets/lab_strains/',
                                '/n/grad_lab2/Lab/gonococcus/datasets/Liu_2022_NCIPpanel/']
for path in additional_pseudogenome_paths:
    for filename in glob.glob(path + 'pseudogenomes/*'):
        basename = filename.split('/')[-1]
        isolates.append(basename[:basename.find('_pseudogenome.fasta')])
        pseudogenome_paths.append(filename)
        
pseudogenome_paths_df = pd.DataFrame({'wgs_id':isolates, 'pseudogenome_path':pseudogenome_paths})
pseudogenome_paths_df.drop_duplicates(inplace = True, ignore_index = True)

# Merge all strains isolates list with pseudogenome paths and metadata and save to file
merged = all_strains.merge(gc_metadata, on = 'wgs_id', how = 'left')
merged = merged.merge(pseudogenome_paths_df, on = 'wgs_id', how = 'left')

# Get kit 12 strains (unsupported)
with open('../data/treemmer/kit12_unsupported.txt') as f:
    kit12 = f.read().splitlines()

# Drop strains in hybrid genome that were sequenced by kit12
complete_genome_df = complete_genome_df[~complete_genome_df.wgs_id.isin(kit12)]

# Get strains that were sequenced only by nanopore (two different basecallers, guppy and dorado)
guppy_assembly_paths = glob.glob('/n/grad_lab2/Lab/gonococcus/nanopore_genome_assemblies/guppy/*.fa*')
guppy_strains = []
for i in guppy_assembly_paths:
    filename = os.path.basename(i)
    strain_name = filename[:filename.find('.')]
    if strain_name.find('DGI')==0:
        strain_name = 'UMASS-DGI_' + str(int(strain_name[3:]))
    guppy_strains.append(strain_name)
    
dorado_assembly_paths = glob.glob('/n/grad_lab2/Lab/gonococcus/nanopore_genome_assemblies/dorado/*.fa*')
dorado_strains = []
for i in dorado_assembly_paths:
    filename = os.path.basename(i)
    strain_name = filename[:filename.find('.')]
    if strain_name.find('DGI')==0:
        strain_name = 'UMASS-DGI_' + str(int(strain_name[3:]))
    dorado_strains.append(strain_name)
    
nanopore_only_strains = guppy_strains + dorado_strains
nanopore_only_strains = list(set(nanopore_only_strains))
nanopore_only_strains = ['25818_3#234' if x=='NY0370' else x for x in nanopore_only_strains] # Replace name of NY strain to match that in the metadata table and database
nanopore_only_strains = ['Ng_188' if x=='Ng188' else x for x in nanopore_only_strains] # Replace name of strain to match that in the metadata table and database (SSS project)
nanopore_only_strains = ['Ng_118' if x=='Ng118' else x for x in nanopore_only_strains] # Replace name of strain to match that in the metadata table and database (SSS project)
nanopore_only_strains = ['FC428' if x=='H18208' else x for x in nanopore_only_strains] # Strain from David Eyre https://www.eurosurveillance.org/content/10.2807/1560-7917.ES.2019.24.10.1900147/?crawler=true (may be H18-209?). H18-208 and FC428 are compared in https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7545138/ they are very clonally related, obviously, but not identical

nanopore_only_strains = list(set(nanopore_only_strains) - set(kit12))
nanopore_only_strains.sort()

with open('../data/treemmer/nanopore_only_strains.txt', 'w') as f:
    for line in nanopore_only_strains:
        f.write(f"{line}\n")

# Add F62 strain - not in our metadata table (lab strain) and we don't have short read data
merged = pd.concat([merged, pd.DataFrame({'wgs_id':['F62']})], ignore_index = True)

# Drop duplicates (this keeps the entries from the Ng metadata table since they are first rather than the hybrid genomes metadata table)
merged.drop_duplicates(subset = ['wgs_id'], inplace = True, ignore_index = True)

merged['assembly'] = 'draft'
merged.loc[merged.wgs_id.isin(nanopore_only_strains), 'assembly'] = 'nanopore_only'
merged.loc[merged.wgs_id.isin(complete_genome_df['wgs_id']), 'assembly'] = 'hybrid'
# merged.loc[merged.wgs_id.isin(kit12), 'assembly'] = 'kit12'

if not os.path.exists('../data/treemmer/'):
    os.mkdir('../data/treemmer/')
merged.to_csv('../data/treemmer/lab_strains_hybrid_genomes_metadata.csv')

# Copy pseudogenome paths to run with Gubbins

paths_df = pd.read_csv('../data/treemmer/lab_strains_hybrid_genomes_metadata.csv')
paths_df.dropna(subset = ['pseudogenome_path'], inplace = True)
paths = paths_df['pseudogenome_path'].values

# for path in paths:
#     shutil.copy(path, '/n/holyscratch01/grad_lab/Users/qinqinyu/20240716_gubbins_lab_strains_hybrid_genomes/pseudogenomes')