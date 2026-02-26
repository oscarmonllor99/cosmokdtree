import torch
import torch_kdtree
import time
import numpy as np

#explicit 1 thread for cpu queries
torch.set_num_threads(1)

#Tree construction
L = 100.
ndim = 3
n = 2**23 #~10 million points
data = torch.rand([n, 3], device="cuda") * L - L/2
tree = torch_kdtree.torchBuildCUDAKDTree(data)

#queries run on cpu for this kdtree implementation
tree.cpu()
data = data.cpu()

#uniformly random for query points
nquery = 1000
query = torch.rand(nquery, 3) * L - L/2
ntrials = 5

# KNN
karray = [1, 10, 100, 1000, 10000]
means = []
stds = []
for k in karray:
    times = np.zeros(ntrials)
    for i in range(ntrials):
        t0 = time.perf_counter()
        tree.search_knn(query, k)
        t1 = time.perf_counter()
        times[i] = (t1 - t0) / nquery

    mean = np.mean(times)
    std = np.std(times)
    means.append(mean)
    stds.append(std)

mean = np.array(means)
std = np.array(stds)
np.savetxt("../../torch_knn_times.txt", np.vstack((karray, means, stds)).T, header="k mean_time(s) stddev(s)")

# BALL
rad_array = [0.01, 0.1, 1., 5., 10.]
means = []
stds = []
for rad in rad_array:
    times = np.zeros(ntrials)
    for i in range(ntrials):
        t0 = time.perf_counter()
        tree.search_radius(query, rad)
        t1 = time.perf_counter()
        times[i] = (t1 - t0) / nquery

    mean = np.mean(times)
    std = np.std(times)
    means.append(mean)
    stds.append(std)

mean = np.array(means)
std = np.array(stds)

del tree
torch.cuda.empty_cache()

np.savetxt("../../torch_ball_times.txt", np.vstack((rad_array, means, stds)).T, header="radius mean_time(s) stddev(s)")

