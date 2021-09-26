import pandas as pd 
import numpy as np
import scipy.stats

import seaborn as sns
import matplotlib.pyplot as plt

from gseapy import enrichr
from sklearn.cluster import AgglomerativeClustering as hclust
from scipy.spatial import distance
from scipy.cluster import hierarchy

import GEOparse as geo 
import ensembl_loci

#basic methylation gse puller
#expects ID_REF and VALUE columns, but wont work on all GSEs
#requires more information to accurately split into control/case
def gse_df(gse):
    gse = geo.get_GEO(geo=gse,destdir=path)
    out = [x.table.set_index("ID_REF").VALUE for key,x in gse.gsms.items()]
    out = pd.DataFrame(out)
    out.index = list(gse.gsms.keys())
    return(out)

def gse_meta(gse):
    pass

def to_genes(gse):
    #load ensembl ids and probe mapping csvs
    #use those to collapse to mean across features
    pass

def load_data()
    pass

def collapse_data():
    pass
'
def get_de():
    pass
'
def correlate_de():
    pass

def create_clustermap(corr):
    correlations_array = np.asarray(corr)
    row_linkage = hierarchy.linkage(
        distance.pdist(correlations_array), method='average')

    col_linkage = hierarchy.linkage(
        distance.pdist(correlations_array.T), method='average')

    g = sns.clustermap(correlations, row_linkage=row_linkage, col_linkage=col_linkage, method="average", figsize=(10, 10),cmap=sns.color_palette("vlag", as_cmap=True),yticklabels=True,xticklabels=True)
    g.ax_row_dendrogram.set_visible(False)
    g.ax_cbar.set_position((0.1, .22, .03, .4))
    return(g)

def enrich_clusters(corr,n=3,gene_sets=["KEGG_2019_Human"]):
    clusters = hclust(n_clusters=3).fit_predict(corr)
    clusters = pd.DataFrame([clusters,list(toclust.index)])
    clusters = clusters.T
    clusters.columns = ["n","gene"]
    output = dict()
    for n in clusters.n.unique():
        output[n]=enrichr(gene_list=list(clusters[clusters.n==1].gene.unique()),background=b,gene_sets=gene_sets)
    return(output)