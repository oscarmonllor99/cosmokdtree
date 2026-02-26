"""
Python wrapper for cosmokdtree. Interplay with f2py module
"""

#f2py-compiled module
import os
import numpy as np
import sys
sys.path.append('../')
import pycosmokdtree

#this file only containes calls for the F2PY-WRAPPED Fortran cosmokdtree routines

#WRAPPER FUNCTIONS
#build tree
def build_kdtree(points, leaf=None, boxsize=None):
    
    if leaf is None:
        leaf = 16

    assert leaf > 0, "Leaf size must be a positive integer"

    if boxsize is None:
        boxsize = -1. * np.ones(points.shape[1])

    assert len(boxsize) == points.shape[1], "Incorrect dimensions for boxsize var. Size must be [Lx, Ly, Lz, ...]"

    boxsize = np.array(boxsize, dtype=np.float64)

    tree = pycosmokdtree.pycosmokdtree.pybuild_kdtree(points, leaf, boxsize)
    return tree

#deallocation
def deallocate_kdtree(tree):
    pycosmokdtree.pycosmokdtree.pydeallocate_tree(tree)

#knn search
def knn_search(tree, target, k, sorted = True):
    assert type(k) in [int, np.int32, np.int64], "k must be an integer"
    target = np.array(target, dtype=np.float64)
    #convert target to 2D array (npoints, ndim) if only one point is given
    if len(target.shape) == 1:
        target = target.reshape(1, -1)

    ntar = target.shape[0]
    ndim = target.shape[1]

    dist, idx = pycosmokdtree.pycosmokdtree.pyknn_search(tree, target, ndim, ntar, k, sorted)
    idx -= 1 #convert to 0-based indexing
    return dist, idx

#ball search
def ball_search(tree, target, radius, sorted = False):
    assert radius > 0., "Radius must be a positive number"
    radius = float(radius)
    target = np.array(target, dtype=np.float64)
    #convert target to 2D array (npoints, ndim) if only one point is given
    if len(target.shape) == 1:
        target = target.reshape(1, -1)

    ntar = target.shape[0]
    ndim = target.shape[1]

    nball, nballmax, query = pycosmokdtree.pycosmokdtree.pyball_search_call_1(tree, target, ndim, ntar, 
                                                                    radius, sorted)
    dist, idx = pycosmokdtree.pycosmokdtree.pyball_search_call_2(ndim, ntar, query, nballmax, nball)
    idx -= 1 #convert to 0-based indexing

    #to list format
    dist_list = []
    idx_list = []
    for i in range(ntar):
        dist_list.append(dist[i, :nball[i]])
        idx_list.append(idx[i, :nball[i]])

    return dist_list, idx_list

#box search
def box_search(tree, box):
    box = np.array(box, dtype=np.float64)
    nbox, query = pycosmokdtree.pycosmokdtree.pybox_search_call_1(tree, box)
    idx = pycosmokdtree.pycosmokdtree.pybox_search_call_2(query, nbox)
    idx -= 1 #convert to 0-based indexing
    return idx