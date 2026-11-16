import json
import math
import os
import time
from itertools import product

import numpy as np


def new_rng(seed=None):
    """
        New numpy random generator.
        Without seed: seeded from the current time and process ID,
        so that parallel workers draw independent samples.
        With seed: reproducible.
    """
    if seed is None:
        seed = (time.time_ns() ^ os.getpid()) & 0xFFFFFFFF
    return np.random.default_rng(seed)


def default_num_cpus(max_jobs=None):
    """
        Use all but two CPUs (at least one), and no more CPUs than there are
        jobs.
    """
    num_cpus = max(1, os.cpu_count() - 2)
    if max_jobs is not None:
        num_cpus = min(num_cpus, max_jobs)
    return num_cpus


def split_samples(num_samples, num_cpus):
    """
        Split num_samples into num_cpus parts whose sizes differ by at most
        one.
    """
    base = num_samples // num_cpus
    remainder = num_samples % num_cpus
    return [base + 1] * remainder + [base] * (num_cpus - remainder)


def get_range(target, prec=100, step=2, steps_below=5, steps_above=2):
    """
        Values of al around target, in units of 1/prec:
        steps_below steps below target, steps_above steps above
        (negative values dropped).
    """
    digits = round(math.log10(prec))
    scaled = round(target * prec)
    R = [round(r / prec, digits)
         for r in range(scaled - (step * steps_below),
                        scaled + (step * steps_above) + 1, step)
         if r >= 0]
    return sorted(set(R))


def nonzero_coefficients(dim):
    """
        One non-zero vector c in F_3^dim of each pair {c, -c}, one per row.
        For a kernel basis B, the kernel elements are k = c * B. 

        The elements k and -k give the same collision pairs 
        (with x and x2 swapped), so only one of them is needed.
    """
    # All vectors (c_1, ..., c_dim) with entries in {0, 1, 2}, in counting
    # order: (0,...,0,0), (0,...,0,1), (0,...,0,2), ...
    all_coeffs = list(product(range(3), repeat=dim))

    # Keep the vectors whose first non-zero entry is 1: of each pair
    # {c, -c} exactly one (since -1 = 2 in F_3). This also drops (0, ..., 0).
    nonzero_coeffs = []
    for c in all_coeffs:
        nonzero_entries = [ci for ci in c if ci != 0]
        if nonzero_entries and nonzero_entries[0] == 1:
            nonzero_coeffs.append(c)

    # One coefficient vector per row: a (3^dim - 1) / 2 x dim matrix
    return np.array(nonzero_coeffs, dtype=np.int64)


def kernel_hamming_weights(G_check, B, C):
    """
        Hamming weights of G_check * k for the kernel elements k = c * B,
        c a row of C. Uses G_check * k = (G_check * B^T) * c, so only one
        matrix product per kernel.
    """
    # l x dim: images of the basis vectors
    M = G_check * B.transpose()
    M = np.array(M, dtype=np.int64)
    Z = (C @ M.T) % 3                               # row j: G_check * k_j
    return np.count_nonzero(Z, axis=1)


def dumpToFile(path, desc, dump_data):
    """
        Write dump_data to file.
    """
    with open(f"{path}/{desc}.data", 'w', encoding="utf-8") as file:
        json.dump(dump_data, file)


def readFromFile_distr(path, desc):
    """
        Read file to dict and convert keys (Hamming weights) to integers
        and values (average counts) to floats.
    """
    dump_data = readFromFile(path, desc)
    dump_data = {int(k): float(v) for k, v in dump_data.items()}
    return dump_data


def readFromFile(path, desc):
    """
        Read file to dict.
    """
    with open(f"{path}/{desc}", 'r', encoding="utf-8") as file:
        dump_data = json.load(file)
    return dump_data


def findXR(data_kernel, R):
    """
        Extract occurrences of Kernels for the R smallest Hamming weights
    """
    cnt_hw = 0
    data_kernel_R = {}
    x_ker = sorted(data_kernel.keys())
    for hw in x_ker:
        cnt_hw += data_kernel[hw]
        if cnt_hw < R:
            data_kernel_R[hw] = data_kernel[hw]
        else:  # elif cnt_hw >= R:
            data_kernel_R[hw] = R - (cnt_hw - data_kernel[hw])
            break
    return data_kernel_R
