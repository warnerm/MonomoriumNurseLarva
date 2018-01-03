import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from joblib import Parallel, delayed
import multiprocessing as mp
import re


#Normalize fpkm using hyperbolic sine transformation
def hypsine ( d ):
    return np.log(d + np.sqrt(d ** 2 + 1))

#Perform quantile normalization
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


#Derive larva and nurse matrices
def getMat(nurse):
    dataH = data.filter(regex=nurse, axis=1)
    dataL = data.filter(regex='W_L', axis=1)
    if (nurse == 'RH') | (nurse == 'RG'):
        #RH is queen-right, but isn't labeled as such in the columns
        dataH.columns = [x[:3] + 'Q' for x in dataH.columns]
    else:
        dataH.columns = [x[:4] for x in dataH.columns]
    print dataH
    dataL.columns = [x[:4] for x in dataL.columns]
    dataL = dataL.iloc[meds, :] #Restrict to medoids
    l1 = np.array(dataL.columns)
    l2 = np.array(dataH.columns)
    #Get shared columns
    sharedCol = list(set(l1) & set(l2))
    dataL = dataL.loc[:, sharedCol]
    dataH = dataH.loc[:, sharedCol]
    return dataL, dataH

#Derive distance matrix, input distance matrix to callMod to derive module definitions
def sortNurse(dataL,nurseD,scramble=True):
    if (scramble):
        nurseD = nurseD.sample(n = nurseD.shape[1], axis = 1) #Scramble nurse stage labels
    pearson = [[np.corrcoef(dataL.iloc[i],nurseD.iloc[j])[1][0] for i in range(dataL.shape[0])] for j in range(nurseD.shape[0])]
    dnew = pd.DataFrame.from_records(pearson)
    dist = 1 - abs(dnew)
    modL = [np.argmin(dist.iloc[row, ]) for row in range(dist.shape[0])]
    modL = meds[modL]
    return modL

#Initalize data files with headers
def initialize(nurse,Dnurse):
    f = open('sortMods'+nurse+'.txt', 'w')
    for n in range(np.shape(Dnurse)[0]):
        f.write('gene'+str(n))
        if n < (np.shape(Dnurse)[0] - 1):
            f.write('\t')
    f.write('\n')
    f.close()

#Write a line of results
def writeRes(res,nurse,Dnurse):
    f = open('sortMods'+nurse+'.txt', 'a')
    for n in range(np.shape(Dnurse)[0]):
        f.write(str(res[n]))
        if n < (np.shape(Dnurse)[0] - 1):
            f.write('\t')
    f.write('\n')
    f.close()

def run(nurse,boots):
    dataL, Dnurse = getMat(nurse)
    initialize(nurse, Dnurse)  # write header for files
    res = sortNurse(dataL, Dnurse, scramble=False)
    writeRes(res, nurse, Dnurse)
    num_cores = mp.cpu_count()

    # Note: the backend="threading" is necessary so that the local variables defined in sortNurse aren't shared
    results = Parallel(n_jobs=num_cores, backend="threading")(delayed(sortNurse)(dataL, Dnurse) for i in range(boots))
    [writeRes(results[i], nurse, Dnurse) for i in range(len(results))]

if __name__ == '__main__':
    # Read in fpkm data
    data = pd.read_csv("~/Nurse_Larva/fpkm.csv", index_col=0)
    data = data.filter(regex='CG|CH|W_L|RH|RG', axis=1)  # Keep only worker larvae samples
    data = data.apply(hypsine) #hyperbolic sine, similar to log transform
    data = quantileNormalize(data.transpose())
    data = data.transpose() #quantile normalize operates on columns, we want to perform it on genes

    # Load larval module definitions
    mods = pd.read_table("~/Nurse_Larva/findK_cluster.txt")
    mods = mods.iloc[10, :]  # Based on SIL, K = 12, which is the 11th row, is the optimal number of medoids
    meds = pd.unique(mods)  # Get list of medoids
    #run('CH',100)
    #run('CG',100)
    run('RH',1)
    run('RG',100)
    run('QCH',100)
    run('QCG',100)