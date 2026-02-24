import numpy as np
import sys
sys.path.append('/home/monllor/projects/gcc_kdtree/cosmokdtree-dev/src/')
from cosmokdtree import kdtree
import time
import timeit

#Tree construction
L = 100.
ndim = 3

n = 10_000_000
points = np.load(f"../../examples/points_7.npy")

tree = kdtree.build_kdtree(points)


#uniformly random between 0, L = 100
nquery = 100000
x = np.random.uniform(0, L, nquery)
y = np.random.uniform(0, L, nquery)
z = np.random.uniform(0, L, nquery)
query_pts = np.vstack((x, y, z)).T

# KNN
karray = [1, 10, 100, 1000, 10000]
means = []
stds = []
for k in karray:
    times = np.zeros(nquery)
    for i in range(nquery):
        t0 = time.time()
        kdtree.knn_search(tree, query_pts[i], k)
        t1 = time.time()
        times[i] = t1 - t0

    mean = np.mean(times)
    std = np.std(times)
    means.append(mean)
    stds.append(std)

mean = np.array(means)
std = np.array(stds)
np.savetxt("../../inhomo_pybindings_knn_times.txt", np.vstack((karray, means, stds)).T, header="k mean_time(s) stddev(s)")



# BALL
rad_array = [0.01, 0.1, 1., 5., 10.]
means = []
stds = []
for rad in rad_array:
    times = np.zeros(nquery)
    for i in range(nquery):
        t0 = time.time()
        kdtree.ball_search(tree, query_pts[i], rad)
        t1 = time.time()
        times[i] = t1 - t0

    mean = np.mean(times)
    std = np.std(times)
    means.append(mean)
    stds.append(std)

mean = np.array(means)
std = np.array(stds)
np.savetxt("../../inhomo_pybindings_ball_times.txt", np.vstack((rad_array, means, stds)).T, header="radius mean_time(s) stddev(s)")

kdtree.deallocate_kdtree(tree)