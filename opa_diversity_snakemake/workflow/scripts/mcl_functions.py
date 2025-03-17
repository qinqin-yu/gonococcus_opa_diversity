#!/usr/bin/env python

import pandas as pd

def parse_mcl_clusters(clusters_filename):
    # clusters filename is typically named "dump..."
    opa_ids = []
    mcl_clusters = []
    i = 1
    with open(clusters_filename, "r") as infile:
        for line in infile:
            opa_ids_cluster = line.strip('\n').split('\t')
            opa_ids_cluster = [s for s in opa_ids_cluster]
            opa_ids = opa_ids+opa_ids_cluster
            mcl_clusters = mcl_clusters + ([i]*len(opa_ids_cluster))
            i+=1
    annotation = 'cluster'
    df = pd.DataFrame({'id':opa_ids, annotation:mcl_clusters})
    return df
