import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from scipy.cluster.hierarchy import dendrogram, linkage
from scipy.cluster.hierarchy import fcluster

clusters = pd.read_table('Kmedioid_larva.txt')

clusters = clusters.iloc[:,:clusters.shape[1] - 1] #Remove NaN from last column
nGene = clusters.shape[1]

#Normalize fpkm using hyperbolic sine transformation
def hypsine ( d ):
    return np.log(d + np.sqrt(d ** 2 + 1))

def quantileNormalize(df_input):
    df = df_input.copy()
    #compute rank
    dic = {}
    for col in df:
        dic.update({col : sorted(df[col])})
    sorted_df = pd.DataFrame(dic)
    rank = sorted_df.mean(axis = 1).tolist()
    #sort
    for col in df:
        t = np.searchsorted(np.sort(df[col]), df[col])
        df[col] = [rank[i] for i in t]
    return df

#Read in fpkm data
data = pd.read_csv("~/Data/Nurse_Larva/fpkm.csv",index_col = 0)
data = data.filter(regex='W_L',axis=1) #Keep only worker larvae samples
data = data[:50]
data = data.apply(hypsine)
data = quantileNormalize(data.transpose())
data = data.transpose()

#Goal is to find the configuration that maximizes inertia, by comparing profiles of centroids (a specific gene) to all genes in that cluster
#Calculate within-cluster sum of squares
dist = [np.sum([[(data.iloc[row,i] - data.iloc[clusters.iloc[b,row],i])**2 for i in range(data.shape[1])] for row in range(data.shape[0])]) for b in range(clusters.shape[0])]

optimal =  np.argmin(dist)
np.savetxt('finalClusters.txt',clusters.iloc[optimal])

