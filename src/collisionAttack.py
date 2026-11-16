import heapq
import math
import multiprocessing as mp
from itertools import product

import pcommon
import strongPRF
import parameters
import plots

from paths import ATK_DIR, DATA_DIR

from sage.all import vector, GF

F3 = GF(3)
F2 = GF(2)

# Step 1 keeps R_FACTOR * R kernel elements instead of R,
# in case the R lowest ones yield fewer collision pairs than expected.
R_FACTOR = 8

# File name (without extension) for saved data and plots
DESC = 'CollisionAttack_secpar-{secpar}_an-{an}_q-{q}_exps-{num_exps}'

def collision_pairs(k):
    """
        Generate all collision pairs (x, x2 = x + k) for a kernel element k,
        (cf. Definition 5)
            k_i = 1  ->  x_i = 0
            k_i = 2  ->  x_i = 1
            k_i = 0  ->  x_i in {0, 1}, all combinations
    """
    map_kernel_to_x = {1: 0, 2: 1}

    # x for fixed positions, None for free positions
    x_fixed = []
    zero_pos = []

    for i in range(len(k)):
        if k[i] == 0:
            zero_pos.append(i)
            x_fixed.append(None)
        else:
            x_fixed.append(map_kernel_to_x[k[i]])

    for bits in product([0, 1], repeat=len(zero_pos)):
        x = x_fixed[:]
        for i, b in zip(zero_pos, bits):
            x[i] = b
        x2 = [xi + ki for xi, ki in zip(x, k)]
        yield x, x2


def step_one(G_prime, n, l, N, R, seed=None):
    """
        Step (1) of the attack for the generalized regime.
        (cf. Fig. 1)

        Process N non-zero kernel elements and keep the R with the lowest
        Hamming weights. Return them as a list of (k, check_rows), sorted by
        increasing Hamming weight.

        seed: optional, makes the random choice of rows reproducible.
    """
    rng = pcommon.new_rng(seed)

    # Max-heap of the R lowest Hamming weights (heapq is a min-heap, so store
    # -hw). Entries: (-hw, cnt_T1, k, check_rows). cnt_T1 is unique, so k is
    # never compared.
    heap = []

    cnt_T1 = 0

    while cnt_T1 < N:
        ker_rows = rng.choice(n, size=(n - l), replace=False).tolist()
        ker_set = set(ker_rows)
        check_rows = [i for i in range(n) if i not in ker_set]

        G_ker = G_prime[ker_rows, :]
        G_check = G_prime[check_rows, :]

        # Kernel basis, one basis vector per row.
        # Note: Sagemath leaks memory here, large numbers for n and N require
        # a lot of memory.
        B = G_ker.right_kernel_matrix()
        dim = B.nrows()

        #print(f"cnt_T1 {cnt_T1}")

        # Trivial kernel: no non-zero elements
        if dim == 0:
            continue

        # Non-zero kernel elements +-k = c * B, at most N - cnt_T1 of them
        C = pcommon.nonzero_coefficients(dim)[:N - cnt_T1]
        hws = pcommon.kernel_hamming_weights(G_check, B, C)

        for c, hw in zip(C, hws):
            hw = int(hw)
            if len(heap) < R or hw < -heap[0][0]:
                # Build k only if it is kept
                k = vector(F3, c.tolist()) * B
                entry = (-hw, cnt_T1, k, check_rows)

                if len(heap) < R:
                    heapq.heappush(heap, entry)
                else:   # new hw < current largest kept hw: replace largest
                    heapq.heapreplace(heap, entry)
            # else: all R HWs are <= than current. Do nothing.

            cnt_T1 += 1

    # Sort by increasing Hamming weight (ties: in order of sampling)
    sorted_heap = sorted(heap, key=lambda entry: (-entry[0], entry[1]))
    return [(k, check_rows) for _, _, k, check_rows in sorted_heap]


def step_two(G, sparse, sorted_k_rows, T2):
    """
        Step (2) of the attack for the generalized regime.
        (cf. Fig. 1)

        Return the first collision found within T2 evaluations and its
        check_rows, or (None, None).
    """
    cnt_T2 = 0
    for k, check_rows in sorted_k_rows:
        G_check = G[check_rows, :]

        # z2 = G_check * (x2||1) = G_check * (x||1) + G_check * (k||0)
        #    = z1 + z_k.
        # Rows where z_k = 0 agree mod 3 already, so only the hw(z_k) rows in
        # its support need checking.
        z_k = G_check[:, :-1] * k
        support = z_k.support()

        G_supp = G_check[support, :]

        z_k_supp = vector(F3, [z_k[i] for i in support])

        # Search through the collision pairs of k until collision found
        for x, x2 in collision_pairs(k):
            cnt_T2 += 1

            z1 = strongPRF.G_x_F3(G_supp, x)
            z2 = z1 + z_k_supp

            if (strongPRF.vec_F3_to_F2(z1, sparse)
                    == strongPRF.vec_F3_to_F2(z2, sparse)):
                return [vector(F2, x), vector(F2, x2)], check_rows

            # Budget of T2 evaluations used up
            if cnt_T2 >= T2:
                return None, None

    return None, None


def collision_search(sPRF, l, N, R, T2, seed=None):
    """
        CollisionSearch algorithm for the generalized regime.
        (cf. Fig. 1)

        Return the collision (or None) and its check_rows.
    """
    G_prime = sPRF.G[:, :sPRF.d]

    sorted_k_rows = step_one(G_prime, sPRF.n, l, N, R, seed)

    collision, check_rows = step_two(sPRF.G, sPRF.sparse, sorted_k_rows, T2)



    # Sanity check
    if collision is not None:
        if sPRF.eval_sPRF(collision[0]) != sPRF.eval_sPRF(collision[1]):
            return None, None

    return collision, check_rows


def collision_search_kernel_regime(sPRF):
    """
        CollisionSearch algorithm for the kernel regime.
        (cf. Fig. 1)
    """
    assert sPRF.d > sPRF.n, "Parameter d > n required for kernel regime."

    G_prime = sPRF.G[:, :sPRF.d]
    K = G_prime.right_kernel()

    # Trivial kernel: no collision pair
    if K.dimension() == 0:
        return None, None

    # Any non-zero kernel element gives collision pairs, take the first one
    k = K.basis()[0]
    x, x2 = next(collision_pairs(k))
    return x, x2


def experiment(secpar, an, al, q, T, sparse=True):
    """
        One run of the attack on a fresh random sPRF. Return the collision or
        None.
    """
    n, l, d, m, t = parameters.sPRF_parameters(secpar, q, an, al)
    sPRF = strongPRF.simple_sPRF(secpar, d, n, m, t, sparse=sparse)

    N = parameters.get_N(T)
    R = parameters.get_R(T, d)
    T2 = T

    R_bound = R_FACTOR * math.ceil(R)

    collision, _ = collision_search(sPRF, l, N, R_bound, T2)

    return collision


def stats(secpar, an, max_al, q, T, num_exps, sparse=True, num_cpus=None):
    """
        Manage experiments, multiprocessing, input and output.
        Run num_exps experiments for each al in pcommon.get_range(max_al),
        save and plot the success rates.
    """

    data_dump = {
        'secpar': secpar,
        'T': T,
        'al': [],
        'max_al': max_al,
        'success': [],
        'num_exps': num_exps
    }

    range_al = pcommon.get_range(max_al, prec=100, step=2, steps_below=5, steps_above=2)

    # All experiments for all al in one pool
    pool_args = []
    for al in range_al:
        pool_args += [(secpar, an, al, q, T, sparse)] * num_exps

    if num_cpus is None:
        num_cpus = pcommon.default_num_cpus(max_jobs=len(pool_args))

    # https://stackoverflow.com/questions/53751050/
    #   multiprocessing-understanding-logic-behind-chunksize/54032744#54032744
    # maxtasksperchild=1, chunksize=1
    # starmap starts a new worker for maxtasksperchild*chunksize  
    # pool_args (ie, experiments), 
    # Each experiment leaks memory (by SageMath), 
    # which is freed after worker terminates. 
    with mp.Pool(num_cpus, maxtasksperchild=1) as pool:
        collisions = pool.starmap(experiment, pool_args, chunksize=1)

    # Results come back in order: num_exps per al
    for i, al in enumerate(range_al):
        collisions_al = collisions[i * num_exps:(i + 1) * num_exps]
        data_dump['success'].append(sum([1 for c in collisions_al if c is not None]))
        data_dump['al'].append(al)

    desc = DESC.format(secpar=secpar, an=an, q=q, num_exps=num_exps)

    pcommon.dumpToFile(DATA_DIR, desc, data_dump)
    plots.plot_attack(data_dump, desc, ATK_DIR)


def plot_from_data(secpar, an, q, num_exps):
    """
        Load data saved with dumpToFile and regenerate plots.
    """
    desc = DESC.format(secpar=secpar, an=an, q=q, num_exps=num_exps)
    print(f"read {DATA_DIR}/{desc}")
    data = pcommon.readFromFile(DATA_DIR, desc + '.data')
    plots.plot_attack(data, desc, ATK_DIR)


if __name__ == '__main__':
    pass
