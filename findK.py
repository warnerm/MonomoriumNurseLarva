import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from joblib import Parallel, delayed
import multiprocessing

#from Bio.Cluster import kcluster
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_samples, silhouette_score
import random


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

#From https://github.com/salspaugh/machine_learning/blob/master/clustering/kmedoids.py
def cluster(distances, k=3):
    m = distances.shape[0]  # number of points

    # Pick k random medoids.
    curr_medoids = np.array([-1] * k)
    while not len(np.unique(curr_medoids)) == k:
        curr_medoids = np.array([random.randint(0, m - 1) for _ in range(k)])
    old_medoids = np.array([-1] * k)  # Doesn't matter what we initialize these to.
    new_medoids = np.array([-1] * k)

    # Until the medoids stop updating, do the following:
    while not ((old_medoids == curr_medoids).all()):
        # Assign each point to cluster with closest medoid.
        clusters = assign_points_to_clusters(curr_medoids, distances)

        # Update cluster medoids to be lowest cost point.
        for curr_medoid in curr_medoids:
            cluster = np.where(clusters == curr_medoid)[0]
            new_medoids[curr_medoids == curr_medoid] = compute_new_medoid(cluster, distances)

        old_medoids[:] = curr_medoids[:]
        curr_medoids[:] = new_medoids[:]

    return clusters, curr_medoids


def assign_points_to_clusters(medoids, distances):
    distances_to_medoids = distances[:, medoids]
    clusters = medoids[np.argmin(distances_to_medoids, axis=1)]
    clusters[medoids] = medoids
    return clusters


def compute_new_medoid(cluster, distances):
    mask = np.ones(distances.shape)
    mask[np.ix_(cluster, cluster)] = 0.
    cluster_distances = np.ma.masked_array(data=distances, mask=mask, fill_value=10e9)
    costs = cluster_distances.sum(axis=1)
    return costs.argmin(axis=0, fill_value=10e9)

#Calculate silhouette distribution
def calcSil(clusters, medioids, distances):
    distances_to_medioids = distances.iloc[:, medioids] #Distance each gene is from metioid
    orig = np.argmin(np.array(distances_to_medioids), axis=1)
    sil = []
    #For each gene, calculate silhouette
    for i in range(np.shape(clusters)[0]):
        distMed=[]
        for j in range(len(medioids)):
            x=np.mean(distances.iloc[i:,clusters[medioids[j]]])
            distMed.append(x)
        s = min(np.delete(distMed,orig[i]))
        #Calculate piecewise silhouette
        if s > distMed[orig[i]]:
            si = 1 - distMed[orig[i]]/s
        elif s == distMed[orig[i]]:
            si = 0
        else:
            si = s/distMed[orig[i]]
        sil.append(si)
    return sil

def clusterSil(K):
    clusters, medioids = cluster(np.array(dist),k = K)
    sil = calcSil(clusters, medioids, dist)
    sil = np.array(sil)
    return(clusters, medioids, sil)

def performCluster(b,kval):
    clusters, medioids, sil = clusterSil(kval)
    return np.mean(sil), clusters

pltN = 0
maxK = 4
minK = 2
boots = 5
nTest = maxK - minK + 1

#Read in fpkm data
data = pd.read_csv("~/Nurse_Larva/fpkm.csv",index_col = 0)
data = data.filter(regex='W_L',axis=1) #Keep only worker larvae samples
data = data.apply(hypsine) #hyperbolic sine, similar to log transform
data = quantileNormalize(data.transpose())
data = data.transpose()
pearson = data.transpose().corr(method = 'pearson')
dist = 1 - abs(pearson)

num_cores = multiprocessing.cpu_count()
f = open('findK_sil.txt','w')
f.write("K\tMeanSil\n")
f.close()

f = open('findK_cluster.txt','w')
for n in range(np.shape(data)[0]):
    f.write('gene' + str(n))
    if n < (np.shape(data)[0] - 1):
        f.write('\t')
f.write('\n')
f.close()

#Run algorithm for a range of K values. At each K value, return the configuration that minimizes sum of squares
for kval in range(minK,maxK+1):
    print kval
    result = Parallel(n_jobs=num_cores)(delayed(performCluster)(b,kval) for b in range(boots))
    # Calculate within-cluster sum of squares

    optimal = np.argmax(result[:][0])
    print optimal
    print max(result[:][0])
    f = open('findK_cluster.txt', 'a')
    for n in range(np.shape(data)[0]):
        f.write(result[optimal][n])
        if n < (np.shape(data)[0] - 1):
            f.write('\t')
    f.write('\n')
    f.close()
