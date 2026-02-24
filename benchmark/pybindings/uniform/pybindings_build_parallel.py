import numpy as np
import sys
sys.path.append('/home/monllor/projects/gcc_kdtree/cosmokdtree-dev/src/')
from cosmokdtree import kdtree
import time
import timeit

#Set number of OMP threads
import os
nthreads = 16
os.environ["OMP_NUM_THREADS"] = str(nthreads)

#Building benchmark
L = 100.0
ndim = 3
narray = [1_000_000, 10_000_000, 100_000_000, 1_000_000_000]
nbuild = 10
means = []
stds = []
for n in narray:
    print(f"Building KDTree with {n} points...")
    points = np.random.rand(n, ndim) * L - L / 2
    times = np.zeros(nbuild)
    for i in range(nbuild):
        print(f"  Build {i+1}/{nbuild}...")
        t0 =  time.time()
        tree = kdtree.build_kdtree(points)
        t1 = time.time()
        times[i] = t1 - t0
        kdtree.deallocate_kdtree(tree)

    mean = np.mean(times)
    std = np.std(times)
    means.append(mean)
    stds.append(std)

means = np.array(means)
stds = np.array(stds)

np.savetxt("../../pybindings_16_build_times.txt", np.vstack((narray, means, stds)).T, header="npoints mean_time(s) stddev(s)")