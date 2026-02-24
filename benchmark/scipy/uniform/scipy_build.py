import numpy as np
from scipy.spatial import KDTree as KDTree
import time
import timeit

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
        tree = KDTree(points, balanced_tree=True)
        t1 = time.time()
        times[i] = t1 - t0

    mean = np.mean(times)
    std = np.std(times)
    means.append(mean)
    stds.append(std)

means = np.array(means)
stds = np.array(stds)

np.savetxt("../../scipy_build_times.txt", np.vstack((narray, means, stds)).T, header="npoints mean_time(s) stddev(s)")