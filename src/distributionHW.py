import multiprocessing as mp
import numpy as np
from collections import Counter

from sage.all import matrix, GF

import pcommon
import parameters
import plots

from paths import DIST_DIR, DATA_DIR

# File name (without extension) for saved data and plots
DESC = 'secpar-{secpar}_an-{an}_al-{al}_q-{q}_exps-{num_exps}'


def F3_HW(l, N):
    """
        Sample N vectors in F_3^l.
        Return a histogram: Hamming weight -> count.
    """
    rng = pcommon.new_rng()

    hw_to_count = Counter()
    for _ in range(N):
        hw = np.count_nonzero(rng.integers(0, 3, size=l))
        hw_to_count[int(hw)] += 1

    return hw_to_count


def ker_HW(n, l, N, q):
    """
        Draw random submatrices G1,G2,
        compute kernel of G1 and the Hamming weights of G2 * k for each
        kernel element k.
        Repeat until N elements have been processed.
        Return a histogram: Hamming weight -> count.

        Correspond to Line 1 - 6 of pseudocode
        (cf. Figure 1).
    """
    rng = pcommon.new_rng()

    d = n - l + q
    G = matrix(GF(3), rng.integers(0, 3, size=(n, d)))

    hw_to_count = Counter()

    cnt_T1 = 0
    while cnt_T1 < N:
        ker_rows = rng.choice(n, size=(n - l), replace=False).tolist()
        ker_set = set(ker_rows)
        check_rows = [i for i in range(n) if i not in ker_set]

        G_ker = G[ker_rows, :]
        G_check = G[check_rows, :]

        # Kernel basis, one basis vector per row.
        # Note: Sagemath leaks memory here, large numbers for n and N require
        # a lot of memory.
        B = G_ker.right_kernel_matrix()
        dim = B.nrows()

        # Trivial kernel: no non-zero elements
        if dim == 0:
            continue

        # Non-zero kernel elements +-k = c * B, at most N - cnt_T1 of them
        C = pcommon.nonzero_coefficients(dim)[:N - cnt_T1]
        hws = pcommon.kernel_hamming_weights(G_check, B, C)

        for hw in hws:
            hw_to_count[int(hw)] += 1

            cnt_T1 += 1
    return hw_to_count


def F3_histogram(l, samples_per_cpu, num_cpus):
    """
        F3: one run of sum(samples_per_cpu) samples, split over all CPUs.
    """
    F3_args = [(l, N_cpu) for N_cpu in samples_per_cpu]

    F3_hw_to_count = Counter()
    with mp.Pool(processes=num_cpus) as pool:
        for r in pool.starmap(F3_HW, F3_args):
            F3_hw_to_count.update(r)
    return F3_hw_to_count


def kernel_histogram(n, l, q, samples_per_cpu, num_cpus):
    """
        Kernel: one experiment of sum(samples_per_cpu) samples, split over all
        CPUs.
        New pool for every experiment, which frees memory leaked by SageMath.
    """
    ker_args = [(n, l, N_cpu, q) for N_cpu in samples_per_cpu]

    exp_hw_to_count = Counter()
    with mp.Pool(processes=num_cpus) as pool:
        for r in pool.starmap(ker_HW, ker_args):
            exp_hw_to_count.update(r)
    return exp_hw_to_count


def average(hw_to_count, num_exps):
    """
        Divide every count of a histogram by num_exps.
    """
    return {hw: count / num_exps for hw, count in hw_to_count.items()}


def stats(secpar, an, al, q, num_exps, T, num_cpus=None):
    """
        Manage experiments, multiprocessing, input and output.

        F3: one run of N = T samples.
        Kernel: num_exps experiments of N = T samples each.
        Each experiment keeps its own R lowest Hamming weights;
        the histograms and the largest of the R lowest Hamming weights (X_R)
        are averaged over the experiments.
    """
    if num_cpus is None:
        num_cpus = pcommon.default_num_cpus()

    n, l, d, m, t = parameters.sPRF_parameters(secpar, q, an=an, al=al)
    N, R_int, log3_ZR = parameters.attack_budget(T, T, d, l)

    # Each run of N samples is split over all CPUs
    samples_per_cpu = pcommon.split_samples(N, num_cpus)

    # F3: one run of N samples
    F3_hw_to_count = F3_histogram(l, samples_per_cpu, num_cpus)

    # Kernel: num_exps experiments of N samples, each keeps its own R lowest
    # Hamming weights
    sum_ker_hw_to_count = Counter()
    sum_ker_R_hw_to_count = Counter()
    exp_X_R = []
    for _ in range(num_exps):
        exp_hw_to_count = kernel_histogram(n, l, q, samples_per_cpu, num_cpus)
        exp_R_hw_to_count = pcommon.findXR(exp_hw_to_count, R_int)

        sum_ker_hw_to_count.update(exp_hw_to_count)
        sum_ker_R_hw_to_count.update(exp_R_hw_to_count)
        exp_X_R.append(max(exp_R_hw_to_count.keys()))

    # Average over experiments
    avg_ker_hw_to_count = average(sum_ker_hw_to_count, num_exps)
    avg_ker_R_hw_to_count = average(sum_ker_R_hw_to_count, num_exps)
    avg_X_R = sum(exp_X_R) / num_exps

    info = {
        'an': an,
        'al': al,
        'q': q,
        'n': n,
        'l': l,
        'd': d,
        'num_exps': num_exps,
        'T': T,
        'N': N,
        'R': parameters.get_R(T, d),
        'E[X_r]': f'{log3_ZR:.5f}',
        'Max(X_r)': f'{avg_X_R:.5f}',  # average over experiments
    }

    desc = DESC.format(secpar=secpar, an=an, al=al, q=q, num_exps=num_exps)

    pcommon.dumpToFile(DATA_DIR, 'F3_' + desc, F3_hw_to_count)
    pcommon.dumpToFile(DATA_DIR, 'ker_' + desc, avg_ker_hw_to_count)
    pcommon.dumpToFile(DATA_DIR, 'kerR_' + desc, avg_ker_R_hw_to_count)
    pcommon.dumpToFile(DATA_DIR, 'info_' + desc, info)

    plots.plot_dist(avg_ker_hw_to_count, avg_ker_R_hw_to_count,
                    F3_hw_to_count, log3_ZR, desc, DIST_DIR)


def plot_from_data(secpar, an, al, q, num_exps):
    """
        Load data saved with dumpToFile and regenerate plots.
    """
    desc = DESC.format(secpar=secpar, an=an, al=al, q=q, num_exps=num_exps)

    n, l, d, m, t = parameters.sPRF_parameters(secpar, q, an, al)
    F3_hw_to_count = pcommon.readFromFile_distr(
        DATA_DIR, 'F3_' + desc + '.data')
    ker_hw_to_count = pcommon.readFromFile_distr(
        DATA_DIR, 'ker_' + desc + '.data')
    ker_R_hw_to_count = pcommon.readFromFile_distr(
        DATA_DIR, 'kerR_' + desc + '.data')

    info = pcommon.readFromFile(DATA_DIR, 'info_' + desc + '.data')

    T1 = int(info['T'])
    T2 = T1
    N, R_int, log3_ZR = parameters.attack_budget(T1, T2, d, l)

    plots.plot_dist(ker_hw_to_count, ker_R_hw_to_count, F3_hw_to_count,
                    log3_ZR, desc, DIST_DIR)


if __name__ == '__main__':
    pass
