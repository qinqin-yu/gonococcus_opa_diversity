import numpy as np
import pandas as pd

from Bio.Phylo.TreeConstruction import DistanceCalculator
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment

def parse_opa_name(name):
    split = name.split('_opa_')
    strain = split[0]
    opa_gene_num = split[1]
    return strain, opa_gene_num

def get_distance_matrix(aln, get_df = False):
    calculator = DistanceCalculator('identity')
    dm = calculator.get_distance(aln)

    names = np.array(dm.names)
    matrix = dm.matrix
    
    # Convert distance output list to array

    matrix_array = np.empty((len(matrix), len(matrix)))
    matrix_array.fill(np.nan)
    for i in range(len(matrix)):
        matrix_array[i, 0:i+1] = matrix[i]

    # Replace zeros on diagonal with nan
    np.fill_diagonal(matrix_array, np.nan)
    
    if get_df:
        matrix_array_t = np.transpose(matrix_array)

        # Get a df of the values
        id_A = []
        id_B = []
        opa_gene_num_A = []
        opa_gene_num_B = []
        strain_A = []
        strain_B = []
        distance = []
        for i in range(len(names)):
            for j in range(i+1, len(names)):
                id_A.append(names[i])
                id_B.append(names[j])
                strain_A.append(parse_opa_name(names[i])[0])
                strain_B.append(parse_opa_name(names[j])[0])
                opa_gene_num_A.append(parse_opa_name(names[i])[1])
                opa_gene_num_B.append(parse_opa_name(names[j])[1])
                distance.append(matrix_array_t[i,j])
        distance_df = pd.DataFrame({'id_A':id_A, 'id_B':id_B, 'strain_A':strain_A, 'strain_B':strain_B, 'opa_gene_num_A':opa_gene_num_A, 'opa_gene_num_B':opa_gene_num_B, 'distance':distance})
        return matrix_array, distance_df
    else:
        return matrix_array
    
def distance_between_different_alignments(aln1, aln2):
    calculator = DistanceCalculator('identity')

    # Combine alignments
    aln = []
    for i in range(len(aln1)):
        seq_id = aln1[i].id
        seq = aln1[i].seq
        aln.append(SeqRecord(seq,id=seq_id))
    for i in range(len(aln2)):
        seq_id = aln2[i].id
        seq = aln2[i].seq
        aln.append(SeqRecord(seq,id=seq_id))
    aln = MultipleSeqAlignment(aln)

    dm = calculator.get_distance(aln)

    names = np.array(dm.names)
    matrix = dm.matrix

    # Convert distance output list to array

    matrix_array = np.empty((len(matrix), len(matrix)))
    matrix_array.fill(np.nan)
    for i in range(len(matrix)):
        matrix_array[i, 0:i+1] = matrix[i]

    # Replace zeros on diagonal with nan
    np.fill_diagonal(matrix_array, np.nan)

    return matrix_array[-1*len(aln2):, :len(aln1)]

def opa_repertoire_distance(aln):
    # Get the distance between two opa repertoires by matching the closest opa genes in the two strains
    # (using a greedy algorithm) and then calculating total distance. If one strain has fewer opa genes, 
    # then add a distance of 1*(difference in number of opa genes between the two strains) to the total distance.
    # Save to a df.
    
    # First get the strain name from the opa gene name
    ids = []
    for i in range(len(aln)):
        ids.append(aln[i].id)

    strains = []
    opa_gene_nums = []
    for id_i in ids:
        strain, opa_gene_num = parse_opa_name(id_i)
        strains.append(strain)
        opa_gene_nums.append(opa_gene_num)

    ids_metadata = pd.DataFrame({'id':ids, 'strain':strains, 'opa_gene_nums':opa_gene_nums})
    
    strains = np.unique(ids_metadata['strain'])
    
    # Next, calculate the distances between different strains in the alignment
    
    id_strain1_list = []
    id_strain2_list = []
    distance_list = []
    strain1_list = []
    strain2_list = []

    # Loop through strain pairs
    for i in range(len(strains)):
        if i%10==0:
            print(str(i) + '/' + str(len(strains)) + ' strains processed')
        strain1 = strains[i]
        df1 = ids_metadata[ids_metadata['strain']==strain1]
        idx_strain1 = df1.index
        id_strain1 = df1['id'].values
        
        # Get the subset of the alignment for opa genes from this strain
        aln_strain1 = aln[np.min(idx_strain1):np.max(idx_strain1)+1]

        for j in range(i+1, len(strains)):
            strain2 = strains[j]
            df2 = ids_metadata[ids_metadata['strain']==strain2]
            idx_strain2 = df2.index
            id_strain2 = df2['id'].values
            aln_strain2 = aln[np.min(idx_strain2):np.max(idx_strain2)+1]

            # Calculate the distance between all opa genes in the two strains
            distances = distance_between_different_alignments(aln_strain1, aln_strain2)

            id_strain1_list_strain_pair = []
            id_strain2_list_strain_pair = []

            # Loop through and get the most closely related pair, mask those opa genes, and repeat
            for k in range(np.min(distances.shape)):
                strain1_list.append(strain1)
                strain2_list.append(strain2)
                flat_idx = np.nanargmin(distances)
                distance_idxs = np.unravel_index(flat_idx, distances.shape)
                id_strain1_list_strain_pair.append(id_strain1[distance_idxs[1]])
                id_strain2_list_strain_pair.append(id_strain2[distance_idxs[0]])
                distance_list.append(np.nanmin(distances))

                # Mask the rows and columns corresponding to the chosen opa genes
                distances[distance_idxs[0],:]=np.nan
                distances[:,distance_idxs[1]]=np.nan
            
            # Check if one strain has fewer opa genes than the other
            if distances.shape[0]!=distances.shape[1]:

                id_strain1_no_pair = list(set(id_strain1) - set(id_strain1_list_strain_pair))
                id_strain2_no_pair = list(set(id_strain2) - set(id_strain2_list_strain_pair))
                
                # Add a distance of 1 for each opa gene that is missing
                if len(id_strain1_no_pair)>0:
                    id_strain1_list_strain_pair = id_strain1_list_strain_pair + id_strain1_no_pair
                    id_strain2_list_strain_pair = id_strain2_list_strain_pair + [np.nan]*len(id_strain1_no_pair)
                    distance_list = distance_list + [1]*len(id_strain1_no_pair)
                    strain1_list = strain1_list + [strain1]*len(id_strain1_no_pair)
                    strain2_list = strain2_list + [strain2]*len(id_strain1_no_pair)
                elif len(id_strain2_no_pair)>0:
                    id_strain1_list_strain_pair = id_strain1_list_strain_pair + [np.nan]*len(id_strain2_no_pair)
                    id_strain2_list_strain_pair = id_strain2_list_strain_pair + id_strain2_no_pair
                    distance_list = distance_list + [1]*len(id_strain2_no_pair)
                    strain1_list = strain1_list + [strain1]*len(id_strain2_no_pair)
                    strain2_list = strain2_list + [strain2]*len(id_strain2_no_pair)

            id_strain1_list = id_strain1_list+id_strain1_list_strain_pair
            id_strain2_list = id_strain2_list+id_strain2_list_strain_pair
    # Save the best pairs and their distances into a df
    distance_df = pd.DataFrame({'id_A':id_strain1_list, 'id_B':id_strain2_list, 'strain_A':strain1_list, 'strain_B':strain2_list, 'distance':distance_list})
    
    # Calculate the total opa repertoire distance between two strains
    strains_A = []
    strains_B = []
    opa_distances = []
    num_opa_difference = []
    strain_A_num_opa = []
    strain_B_num_opa = []
    for strains, df in distance_df.groupby(['strain_A', 'strain_B']):
        strain_A = strains[0]
        strain_B = strains[1]
        strains_A.append(strain_A)
        strains_B.append(strain_B)
        opa_distances.append(np.sum(df['distance']))
        num_opa_difference.append(df[['id_A', 'id_B']].isna().values.sum())
        strain_A_num_opa.append(len(df[~df['id_A'].isna()]))
        strain_B_num_opa.append(len(df[~df['id_B'].isna()]))
    summary = pd.DataFrame({'strain_A':strains_A, 'strain_B':strains_B, 'opa_distance':opa_distances, 'strain_A_num_opa':strain_A_num_opa, 'strain_B_num_opa':strain_B_num_opa, 'num_opa_difference':num_opa_difference})

    return distance_df, summary