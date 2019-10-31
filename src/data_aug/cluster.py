#!/usr/bin/python
import numpy as np
from math import pi
from sklearn.cluster import DBSCAN
from sklearn.cluster import OPTICS
from sklearn.cluster import SpectralClustering
from sklearn.cluster import AffinityPropagation
from sklearn.cluster import AgglomerativeClustering
from sklearn.cluster import MeanShift


def test():
    print("helloworld")

def meanshift(params): 
    distance_path=''
    distance_path+=params["distance_path"]
    print(distance_path)
    distance=np.loadtxt(distance_path,dtype=np.float32)
    print(distance.shape)

    ms=MeanShift().fit(distance)

    #get labels
    labels = ms.labels_

    print(labels,labels.shape)
    #get number of clusters
    no_clusters = len(set(labels)) - (1 if -1 in labels else 0)
    print(no_clusters,"no_clusters")

    #for i in range(no_clusters):
        #print('Cluster  : ', np.nonzero(labels == i)[0])

    #print(type(labels))
    return_val=tuple(labels.tolist())
    #print(type(return_val))
    return return_val
    
def dbscan(params): 
    distance_path=''
    distance_path+=params["distance_path"]
    print(distance_path)
    distance=np.loadtxt(distance_path,dtype=np.float32)
    print(distance.shape)

    #using default values, set metric to 'precomputed'
    db = DBSCAN(eps=0.03, min_samples =5, metric='precomputed')
    #check db
    print(db)

    db.fit(distance)
    #get labels
    labels = db.labels_

    print(labels,labels.shape)
    #get number of clusters
    no_clusters = len(set(labels)) - (1 if -1 in labels else 0)
    print(no_clusters,"no_clusters")

    #for i in range(no_clusters):
        #print('Cluster  : ', np.nonzero(labels == i)[0])

    #print(type(labels))
    return_val=tuple(labels.tolist())
    #print(type(return_val))
    return return_val


def spectralclustering(params): 
    distance_path=''
    distance_path+=params["distance_path"]
    print(distance_path)
    distance=np.loadtxt(distance_path,dtype=np.float32)
    print(distance.shape)
    delta=2
    affinity=np.exp(-distance ** 2/ (2. * delta ** 2))

    #using default values, set metric to 'precomputed'
    sp=SpectralClustering(n_clusters=10,affinity='precomputed')
    print(sp)

    sp.fit(affinity)
    #get labels
    labels = sp.labels_

    print(labels,labels.shape)
    #get number of clusters
    no_clusters = len(set(labels)) - (1 if -1 in labels else 0)
    print(no_clusters,"no_clusters")

    #for i in range(no_clusters):
        #print('Cluster  : ', np.nonzero(labels == i)[0])

    #print(type(labels))
    return_val=tuple(labels.tolist())
    #print(type(return_val))
    return return_val

def agglomerativeclustering(params): 
    distance_path=''
    distance_path+=params["distance_path"]
    print(distance_path)
    distance=np.loadtxt(distance_path,dtype=np.float32)
    print(distance.shape)
    #delta=2
    #affinity=np.exp(-distance ** 2/ (2. * delta ** 2))

    #using default values, set metric to 'precomputed'
    # linkage: average complete single(result is bad)  result quality: average > complete > single
    agg=AgglomerativeClustering(n_clusters=8,affinity='precomputed',linkage='average').fit(distance)

    #get labels
    labels = agg.labels_

    print(labels,labels.shape)
    #get number of clusters
    no_clusters = len(set(labels)) - (1 if -1 in labels else 0)
    print(no_clusters,"no_clusters")

    #for i in range(no_clusters):
        #print('Cluster  : ', np.nonzero(labels == i)[0])

    #print(type(labels))
    return_val=tuple(labels.tolist())
    #print(type(return_val))
    return return_val


def affinitypropagation(params): 
    distance_path=''
    distance_path+=params["distance_path"]
    print(distance_path)
    distance=np.loadtxt(distance_path,dtype=np.float32)
    print(distance.shape)
    delta=2
    affinity=np.exp(-distance ** 2/ (2. * delta ** 2))

    #using default values, set metric to 'precomputed'
    aff=AffinityPropagation(affinity='precomputed')
    print(aff)

    aff.fit(affinity)
    #get labels
    labels = aff.labels_

    print(labels,labels.shape)
    #get number of clusters
    no_clusters = len(set(labels)) - (1 if -1 in labels else 0)
    print(no_clusters,"no_clusters")

    #for i in range(no_clusters):
        #print('Cluster  : ', np.nonzero(labels == i)[0])

    #print(type(labels))
    return_val=tuple(labels.tolist())
    #print(type(return_val))
    return return_val

def optics(params): 
    distance_path=''
    distance_path+=params["distance_path"]
    print(distance_path)
    distance=np.loadtxt(distance_path,dtype=np.float32)
    print(distance.shape)

    #using default values, set metric to 'precomputed'
    op = OPTICS(eps=0.03, min_samples =10, metric='precomputed')
    #check db
    print(op)

    op.fit(distance)
    #get labels
    labels = op.labels_

    print(labels,labels.shape)
    #get number of clusters
    no_clusters = len(set(labels)) - (1 if -1 in labels else 0)
    print(no_clusters,"no_clusters")

    #for i in range(no_clusters):
        #print('Cluster  : ', np.nonzero(labels == i)[0])

    #print(type(labels))
    return_val=tuple(labels.tolist())
    #print(type(return_val))
    return return_val

