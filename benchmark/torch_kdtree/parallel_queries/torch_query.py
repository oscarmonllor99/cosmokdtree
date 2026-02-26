import torch
import torch_kdtree
import time
import numpy as np

#explicit 1 thread for cpu queries
torch.set_num_threads(20)

#Tree construction
L = 100.
ndim = 3
n = 2**23 #~10 million points
points = np.load(f"../../examples/points_2e23.npy")
data = torch.tensor(points, dtype=torch.float32, device="cuda")
data = data.T.contiguous()
tree = torch_kdtree.torchBuildCUDAKDTree(data)
#queries run on cpu for this kdtree implementation
tree.cpu()
data = data.cpu()

#uniformly random for query points
nquery = 100000
query = torch.rand(nquery, 3) * L
ntrials = 2

# KNN
k = 100
times = np.zeros(ntrials)
for i in range(ntrials):
    t0 = time.perf_counter()
    tree.search_knn(query, k)
    t1 = time.perf_counter()
    times[i] = (t1 - t0)

time_knn = np.mean(times)

# BALL
rad = 1.
times = np.zeros(ntrials)
for i in range(ntrials):
    t0 = time.perf_counter()
    tree.search_radius(query, rad)
    t1 = time.perf_counter()
    times[i] = (t1 - t0)

time_ball = np.mean(times)

del tree
torch.cuda.empty_cache()

print("KNN times (s):", time_knn)
print("Ball times (s):", time_ball)