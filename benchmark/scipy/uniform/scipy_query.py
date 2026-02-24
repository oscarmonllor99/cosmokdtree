import numpy as np
from scipy.spatial import KDTree as KDTree
import time
import timeit

#Tree construction
L = 100.
ndim = 3

n = 10_000_000
points = np.random.rand(n, ndim) * L - L / 2

tree = KDTree(points)


#uniformly random between 0, L = 100
nquery = 10000
x = np.random.uniform(-L/2, L/2, nquery)
y = np.random.uniform(-L/2, L/2, nquery)
z = np.random.uniform(-L/2, L/2, nquery)
query_pts = np.vstack((x, y, z)).T

# SCIPY KNN
karray = [1, 10, 100, 1000, 10000]
means = []
stds = []
for k in karray:
    times = np.zeros(nquery)
    for i in range(nquery):
        t0 = time.time()
        tree.query(query_pts[i], k=k)
        t1 = time.time()
        times[i] = t1 - t0

    mean = np.mean(times)
    std = np.std(times)
    means.append(mean)
    stds.append(std)

mean = np.array(means)
std = np.array(stds)
np.savetxt("../../scipy_knn_times.txt", np.vstack((karray, means, stds)).T, header="k mean_time(s) stddev(s)")


# SCIPY BALL
rad_array = [0.01, 0.1, 1., 5., 10.]
means = []
stds = []
for rad in rad_array:
    times = np.zeros(nquery)
    for i in range(nquery):
        t0 = time.time()
        tree.query_ball_point(query_pts[i], r=rad)
        t1 = time.time()
        times[i] = t1 - t0

    mean = np.mean(times)
    std = np.std(times)
    means.append(mean)
    stds.append(std)

mean = np.array(means)
std = np.array(stds)
np.savetxt("../../scipy_ball_times.txt", np.vstack((rad_array, means, stds)).T, header="radius mean_time(s) stddev(s)")