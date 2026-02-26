import torch
import torch_kdtree
import time
import numpy as np
import gc
L = 100.0
ndim = 3
narray = [2**18, 2**22, 2**26]
nbuild = 10
means = []
stds = []
for n in narray:
    print(f"Building KDTree with {n} points...")
    points = torch.rand(n, ndim, device="cuda") * L - L / 2
    times = np.zeros(nbuild)
    for i in range(nbuild):
        print(f"  Build {i+1}/{nbuild}...")

        torch.cuda.synchronize()
        gc.disable()
        t0 = time.perf_counter()

        tree = torch_kdtree.torchBuildCUDAKDTree(points)
        torch.cuda.synchronize()
        t1 = time.perf_counter()
        gc.enable()

        times[i] = t1 - t0
        del tree
        torch.cuda.empty_cache()

    mean = np.mean(times)
    std = np.std(times)
    means.append(mean)
    stds.append(std)

means = np.array(means)
stds = np.array(stds)
np.savetxt("../../torch_build_times.txt", np.vstack((narray, means, stds)).T, header="npoints mean_time(s) stddev(s)")