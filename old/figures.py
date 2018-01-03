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

#Read in fpkm data
data = pd.read_csv("~/Data/Nurse_Larva/fpkm.csv",index_col = 0)
data = data.filter(regex='W_L',axis=1) #Keep only worker larvae samples
data = data[:50]

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
#data = data.apply(hypsine)
#data = quantileNormalize(data.transpose())
#data = data.transpose()
pearson = data.transpose().corr(method = 'pearson')
dist = 1 - abs(pearson)

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

def makeplt():
    ax = plt.subplot(nTest, 1, pltN)
    ax.set_xlim([min(sil), 1])

    # The (n_clusters+1)*10 is for inserting blank space between silhouette
    # plots of individual clusters, to demarcate them clearly.
    ax.set_ylim([0, len(clusters) + (n_clusters + 1) * 10])

    # The silhouette_score gives the average value for all the samples.
    # This gives a perspective into the density and separation of the formed
    # clusters
    # silhouette_avg = silhouette_score(X, cluster_labels)
    y_lower = 10
    for i in medioids:
        # Aggregate the silhouette scores for samples belonging to
        # cluster i, and sort them
        ith_cluster_silhouette_values = \
            sil[clusters == i]

        ith_cluster_silhouette_values.sort()

        size_cluster_i = ith_cluster_silhouette_values.shape[0]
        y_upper = y_lower + size_cluster_i

        color = cm.spectral(float(i) / n_clusters)
        ax.fill_betweenx(np.arange(y_lower, y_upper),
                         0, ith_cluster_silhouette_values,
                         facecolor=color, edgecolor=color, alpha=0.7)

        # Label the silhouette plots with their cluster numbers at the middle
        ax.text(-0.05, y_lower + 0.5 * size_cluster_i, str(i))

        # Compute the new y_lower for next plot
        y_lower = y_upper + 10  # 10 for the 0 samples

    ax.set_title("The silhouette plot for the various clusters.")
    ax.set_xlabel("The silhouette coefficient values")
    ax.set_ylabel("Cluster label")

    # The vertical line for average silhouette score of all the values
    # ax1.axvline(x=silhouette_avg, color="red", linestyle="--")

    # ax.set_yticks([])  # Clear the yaxis labels / ticks
    ax.set_xticks([-0.1, 0, 0.2, 0.4, 0.6, 0.8, 1])

def performCluster(b):
    clusters, medioids, sil = clusterSil(kval)
    return[kval,np.mean(sil)]

pltN = 0
maxK = 30
minK = 2
boots = 100
nTest = maxK - minK + 1
y_lower = 10
# plt.figure(figsize = (11,11))
#plt.subplot(nTest,1,1)
num_cores = multiprocessing.cpu_count()

for kval in range(minK,maxK+1):
    results = Parallel(n_jobs=num_cores)(delayed(performCluster)(b) for b in range(boots))

np.savetxt('findK.txt',results)

# plt.tight_layout()
# plt.subplots_adjust(hspace=0.5)
# plt.savefig('foo.png')
#plt.show()

#[clusterSil(k) for k in range(2,6)]






#Parallel processing for multiple numbers of clusters

#Implement Pham, Dimov, and Nguyen to find optimal number of clusters

#Return consensus clustering result

#X = np.array([(random.uniform(-1, 1), random.uniform(-1, 1)) for i in range(100)])
#clusterid, error, nfound = kcluster(X)
#print clusterid
