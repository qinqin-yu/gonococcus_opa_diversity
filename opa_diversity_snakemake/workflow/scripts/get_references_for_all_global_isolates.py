import pandas as pd
import numpy as np

# Save table of number of isolates from each reference

metadata = pd.read_csv('../data/isolates_summary_and_qc/Ng-Combined-Metadata_20240620.txt', sep = '\t')

# Drop Ethiopian isolates
metadata = metadata[(metadata['reference']!='ethiopia_isolates')&(metadata['reference']!='ethiopia_isolates_grad_lab')]
metadata.reset_index(inplace = True, drop = True)

# Drop unpublished isolates
metadata = metadata[metadata['reference']!='unpublished']
print('Total number of global isolates: ', len(metadata))

# Basic filtering of QC
metadata_qc = metadata[
(metadata['reference_coverage']>40)&
(metadata['reference_percentage_mapped']>80)&
(metadata['assembly_length']>1.75*10**6)&
(metadata['assembly_length']<2.5*10**6)&
(metadata['percent_missing']<12)]
print('Total number of global isolates after QC filtering: ', len(metadata_qc))

# Get the number of isolates from each paper reference
references = []
num_isolates = []
for reference, df in metadata_qc.groupby('reference'):
    references.append(reference)
    num_isolates.append(len(df))
    
# Save results
global_isolates_summary = pd.DataFrame({'reference':references, 'num_isolates':num_isolates})
global_isolates_summary.to_csv('../../input_data/whole_genome_metadata/references_for_all_global_isolates.csv')

metadata_qc.to_csv('../../input_data/whole_genome_metadata/global_isolates_metadata_met_qc.csv', index=False)