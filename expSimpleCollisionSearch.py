"""
    Proof of concept: one collision in the kernel regime and one in the
    generalized regime.
"""
import argparse
import math
import sys

from paths import SRC_DIR
sys.path.append(str(SRC_DIR))

import parameters
import strongPRF
import collisionAttack

from sage.all import GF, vector

F2 = GF(2)
F3 = GF(3)


def print_vec(prefix, vec, width=1):
    print(prefix + " ".join(str(x).rjust(width) for x in vec))


def mark_vec(n, I, mark="x"):
    I = set(I)
    return [mark if i in I else "" for i in range(0, n)]


def print_collision(sPRF, x, x2, desc, mark_I=None):
    """
        Print the inputs x, x2 and the values G * (x||1) mod 3 and mod 2,
        as computed in the PRF evaluation.
    """
    n, d = sPRF.n, sPRF.d
    z1_col_F3 = strongPRF.G_x_F3(sPRF.G, x)
    z2_col_F3 = strongPRF.G_x_F3(sPRF.G, x2)

    z1_col_F2 = strongPRF.vec_F3_to_F2(z1_col_F3, sPRF.sparse)
    z2_col_F2 = strongPRF.vec_F3_to_F2(z2_col_F3, sPRF.sparse)

    print("-" * 50)
    print(f"Collision in {desc} regime for n={n}, d={d}")
    if mark_I is not None:
        print(f"x marks positions potentially != mod 3, but = mod 2")
    print(f"-" * 5)
    print_vec("x        : ", x)
    print_vec("x2       : ", x2)
    print(f"-" * 5)
    if mark_I is not None:
        print_vec(" " * 11, mark_vec(n, mark_I, mark="x"))
    print_vec("z1_col_F3: ", z1_col_F3)
    print_vec("z2_col_F3: ", z2_col_F3)
    print(f"-" * 5)
    if mark_I is not None:
        print_vec(" " * 11, mark_vec(n, mark_I, mark="x"))
    print_vec("z1_col_F2: ", z1_col_F2)
    print_vec("z2_col_F2: ", z2_col_F2)
    print("-" * 50)


def demonstrate_kernel_regime(secpar, an, seed):
    n, _, _, m, t = parameters.sPRF_parameters(secpar, q=0, an=an, al=1)
    d = n + 1

    sPRF = strongPRF.simple_sPRF(secpar, d, n, m, t, seed=seed)
    x, x2 = collisionAttack.collision_search_kernel_regime(sPRF)

    if x is None:
        print("Kernel is trivial, please try again.")
        return

    if sPRF.eval_sPRF(vector(F2, x)) == sPRF.eval_sPRF(vector(F2, x2)):
        print_collision(sPRF, x, x2, 'kernel')
    else:
        print("No collision in kernel regime, please try again.")


def demonstrate_generalized_regime(secpar, an, al, q, T, seed):
    T2 = T

    n, l, d, m, t = parameters.sPRF_parameters(secpar, q, an, al)
    sPRF = strongPRF.simple_sPRF(secpar, d, n, m, t, seed=seed)

    N = parameters.get_N(T)
    R = parameters.get_R(T, d)
    R_bound = collisionAttack.R_FACTOR * math.ceil(R)

    # Different seed than the sPRF, so the two random streams are independent
    collisions, check_rows = collisionAttack.collision_search(
        sPRF, l, N, R_bound, T2, seed=seed + 1)

    if collisions is None:
        print(f"No collision found within T2 = {T2} evaluations, "
              "please try again.")
        return

    if (sPRF.eval_sPRF(vector(F2, collisions[0]))
            == sPRF.eval_sPRF(vector(F2, collisions[1]))):
        print_collision(sPRF, collisions[0], collisions[1], 'generalized',
                        check_rows)
    else:
        print("No collision in generalized regime, please try again.")


def parse_args():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--an', type=float, default=4,
                        help='n = ceil(an * secpar), both regimes')
    parser.add_argument('--kernel_secpar', type=int, default=5,
                        help='security parameter, kernel regime (d = n + 1)')
    parser.add_argument('--secpar', type=int, default=10,
                        help='security parameter, generalized regime')
    parser.add_argument('--al', type=float, default=1.2,
                        help='l = ceil(al * secpar), generalized regime')
    parser.add_argument('--q', type=int, default=5,
                        help='d = n - l + q, generalized regime')
    parser.add_argument('--log2T', type=int, default=None,
                        help='T = 2**log2T, generalized regime; '
                             'None means secpar - 1')
    parser.add_argument('--seed', type=int, default=2026433,
                        help='random seed, both regimes')
    args = parser.parse_args()

    if args.log2T is None:
        args.log2T = args.secpar - 1
    return args


if __name__ == '__main__':
    args = parse_args()
    demonstrate_kernel_regime(args.kernel_secpar, args.an, args.seed)
    demonstrate_generalized_regime(args.secpar, args.an, args.al, args.q,
                                   2**(args.log2T), args.seed)
