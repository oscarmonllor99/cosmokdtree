import torch
import torch_kdtree
import time
import numpy as np
import gc

L = 100.0
ndim = 3
narray = [2**18, 2**22, 2**26]
exps = [18, 22, 26]
nbuild = 10
means = []
stds = []
for exp in exps:
    n = 2**exp
    print(f"Building KDTree with {n} points...")
    points_np = np.load(f"../../examples/points_2e{exp}.npy")
    points = torch.tensor(points_np, dtype=torch.float32, device="cuda")
    points = points.T.contiguous()
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

    del points
    torch.cuda.empty_cache()
    mean = np.mean(times)
    std = np.std(times)
    means.append(mean)
    stds.append(std)

means = np.array(means)
stds = np.array(stds)

np.savetxt("../../inhomo_torch_build_times.txt", np.vstack((narray, means, stds)).T, header="npoints mean_time(s) stddev(s)")